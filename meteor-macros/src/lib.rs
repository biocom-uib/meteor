#![feature(iterator_try_collect)]

use darling::{ast, FromDeriveInput, FromField};
use inflector::cases::pascalcase::to_pascal_case;
use inflector::cases::screamingsnakecase::to_screaming_snake_case;
use proc_macro::TokenStream;
use proc_macro2::TokenStream as TokenStream2;
use proc_macro2::Ident as Ident2;
use proc_macro2::Span as Span2;
use quote::format_ident;
use quote::quote;
use quote::ToTokens;
use syn::{parse_macro_input, DeriveInput};

fn split_generics<F>(generics: &syn::Generics, add_bounds: F) -> (TokenStream2, TokenStream2)
where
    F: Fn(&syn::TypeParam) -> TokenStream2,
{
    let lifetimes: Vec<_> = generics.lifetimes().collect();

    let lt_params = quote! { #( #lifetimes , )* };

    let lt_args = lifetimes.iter().map(|lt| &lt.lifetime);
    let lt_args = quote! { #( #lt_args , )* };

    let type_params: Vec<_> = generics.type_params().collect();

    let generic_params = if lifetimes.is_empty() && type_params.is_empty() {
        TokenStream2::new()
    } else {
        let bounds = type_params.iter().copied().map(add_bounds);
        quote! { < #lt_params #( #bounds ),* > }
    };

    let generic_args = if lifetimes.is_empty() && type_params.is_empty() {
        TokenStream2::new()
    } else {
        let idents = type_params.iter().map(|param| &param.ident);
        quote! { < #lt_args #( #idents ),* > }
    };

    (generic_params, generic_args)
}

#[derive(Clone, Debug, FromField)]
#[darling(attributes(filter))]
struct FilterField {
    ident: Option<syn::Ident>,
    ty: syn::Type,

    #[darling(default)]
    as_string: bool,

    #[darling(default)]
    rename: Option<String>,

    #[darling(default)]
    skip: bool,
}

#[derive(Debug, FromDeriveInput)]
#[darling(supports(struct_named), attributes(filter))]
struct FilterableRecordOpts {
    vis: syn::Visibility,
    ident: syn::Ident,
    generics: syn::Generics,

    data: ast::Data<(), FilterField>,

    as_string: Option<syn::Path>,

    #[darling(default)]
    alias_field: bool,

    #[darling(default)]
    alias_filter: bool,

    #[darling(default)]
    polars_literal: bool,
}

enum VariantType {
    String,
    Other(syn::Type),
}

impl ToTokens for VariantType {
    fn to_tokens(&self, tokens: &mut TokenStream2) {
        tokens.extend(match self {
            VariantType::String => quote! { String },
            VariantType::Other(ty) => quote! { #ty },
        });
    }
}

struct Variant {
    field: Ident2,
    ident: Ident2,
    literal: String,
    ty: VariantType,
}

struct FieldsEnum {
    vis: syn::Visibility,
    ident: Ident2,
    variants: Vec<Variant>
}

impl ToTokens for FieldsEnum {
    fn to_tokens(&self, tokens: &mut TokenStream2) {
        let vis = &self.vis;
        let ident = &self.ident;

        let variant_toks = self
            .variants
            .iter()
            .map(|Variant { ident, ty, .. }| quote! { #ident(#ty) });

        tokens.extend(quote! {
            #[derive(Debug, Clone)]
            #vis enum #ident {
                #( #variant_toks ),*
            }
        });
    }
}

fn fields_enum_variants(opts: &FilterableRecordOpts, fields: &Vec<FilterField>) -> FieldsEnum {
    let vis = &opts.vis;

    let mut variants = Vec::new();

    for field in fields {
        if field.skip {
            continue;
        }

        let field_ident = field
            .ident
            .as_ref()
            .expect("FilterableRecord can only be derived for named structs");

        let field_name = field_ident.to_string();
        let variant_name = to_pascal_case(&field_name);
        let variant_ident = Ident2::new(&variant_name, Span2::call_site());

        let field_ty = &field.ty;

        let string_ty = opts.as_string.clone().map(|path| {
            syn::Type::from(syn::TypePath { qself: None, path })
        });

        let variant_ty = if field.as_string || string_ty.as_ref() == Some(field_ty) {
            VariantType::String
        } else {
            VariantType::Other(field_ty.clone())
        };

        variants.push(Variant {
            field: field_ident.clone(),
            ident: variant_ident,
            literal: field.rename.as_ref().unwrap_or(&field_name).to_owned(),
            ty: variant_ty,
        });
    }

    FieldsEnum {
        vis: vis.clone(),
        ident: Ident2::new("Fields", Span2::call_site()),
        variants,
    }
}

fn gen_fields_enum_from_str_field(fields_enum: &FieldsEnum) -> TokenStream2 {
    let field_names = fields_enum.variants.iter().map(|Variant { literal, .. }| {
        quote! { #literal }
    });

    let field_name_match_arms = fields_enum.variants.iter().map(|Variant { field: _, ident, ty: _, literal }| {
        quote! { Self::#ident(_) => #literal }
    });

    let try_from_match_arms =
        fields_enum
            .variants
            .iter()
            .map(|Variant { field: _, ident, literal, .. }| {
                let rhs = quote! { Ok(Self::#ident(value.parse()?)) };
                /*match ty {
                    VariantType::String => quote! { Ok(Self::#ident(value.to_owned())) },
                    VariantType::Other(_) => quote! { Ok(Self::#ident(value.parse()?)) }
                };*/

                quote! { #literal => #rhs }
            });

    let parse_apply_op_match_arms = fields_enum.variants.iter().map(|Variant { field: _, ident, ty, .. }| {
        let rhs = match ty {
            VariantType::String => quote! { Ok(op.apply(other, value.as_ref())) },
            VariantType::Other(_) => quote! { Ok(op.apply(&other.parse()?, value)) },
        };

        quote! { Self::#ident(value) => #rhs }
    });

    let known_fields = fields_enum.variants
        .iter()
        .map(|v| &v.literal )
        .map(|f| quote!{ #f });

    let fields_enum_ident = &fields_enum.ident;

    quote! {
        impl ::meteor::csv::filter::FromStrField for #fields_enum_ident {
            fn field_names() -> &'static [&'static str] {
                &[ #( #field_names ),* ]
            }

            fn field_name(&self) -> &str {
                match self {
                    #( #field_name_match_arms , )*
                }
            }

            fn parse_apply_op(
                other: &str,
                op: ::meteor::csv::filter::Op,
                this: &Self,
            ) -> Result<bool, Self::Err>
            {
                match this {
                    #( #parse_apply_op_match_arms , )*
                }
            }

            fn try_from_parts(key: &str, value: &str) -> Result<Self, Self::Err> {
                match key {
                    #( #try_from_match_arms , )*
                    _ => Err(Self::Err::UnknownField(key.to_owned(), &[ #( #known_fields ),* ])),
                }
            }
        }
    }
}

fn gen_fields_enum_polars_literal(fields_enum: &FieldsEnum) -> TokenStream2 {
    let fields_enum_ident = &fields_enum.ident;

    let match_arms = fields_enum.variants.iter().map(|Variant { field: _, ident, ty, .. }| {
        match ty {
            VariantType::String => quote! { Self::#ident(value) => value.to_string().lit() },
            VariantType::Other(_) => quote! { Self::#ident(value) => value.lit() },
        }
    });

    quote! {
        impl ::polars::prelude::Literal for #fields_enum_ident {
            fn lit(self) -> ::polars::prelude::Expr {
                match self {
                    #( #match_arms , )*
                }
            }
        }
    }
}

fn gen_filterable_record_impl(
    opts: &FilterableRecordOpts,
    fields_enum: &FieldsEnum,
) -> TokenStream2 {

    let struct_ident = &opts.ident;
    let fields_enum_ident = &fields_enum.ident;

    let (generic_params, generic_args) = split_generics(&opts.generics, |param| {
        let type_ident = &param.ident;
        quote! { #type_ident : ::std::convert::AsRef<str> }
    });

    let match_arms = fields_enum
        .variants
        .iter()
        .map(|Variant { field, ident, ty, .. }| {
            let lhs = quote! { #fields_enum_ident::#ident(value) };

            let rhs = match ty {
                VariantType::String => quote! { filter.op.apply(self.#field.as_ref(), value.as_str()) },
                VariantType::Other(_) => quote! { filter.op.apply(&self.#field, value) },
            };

            quote! { #lhs => #rhs }
        });

    quote! {
        impl #generic_params ::meteor::csv::filter::FilterableRecord for #struct_ident #generic_args {
            type Field = #fields_enum_ident;

            fn apply_filter(&self, filter: &::meteor::csv::filter::Filter<Self::Field>) -> bool {
                match &filter.field_value {
                    #( #match_arms , )*
                }
            }
        }
    }
}

fn gen_aliases(opts: &FilterableRecordOpts) -> TokenStream2 {
    let vis = &opts.vis;

    let lifetimes: Vec<_> = opts.generics.lifetimes().collect();
    let lt_statics = lifetimes.iter().map(|_| quote! { 'static });

    let type_params: Vec<_> = opts
        .generics
        .type_params()
        .filter(|param| param.default.is_none())
        .collect();

    let generic_params = if type_params.is_empty() {
        TokenStream2::new()
    } else {
        quote! { < #( #type_params ),* > }
    };

    let generic_args = type_params.iter().map(|param| &param.ident);
    let generic_args = if lifetimes.is_empty() && type_params.is_empty() {
        TokenStream2::new()
    } else {
        quote! { < #( #lt_statics , )* #( #generic_args ),* > }
    };

    let alias_field = if opts.alias_field {
        let struct_ident = &opts.ident;
        let ident = format_ident!("{}Field", &opts.ident);

        quote! {
            #vis type #ident #generic_params =
                <#struct_ident #generic_args as ::meteor::csv::filter::FilterableRecord>::Field;
        }
    } else {
        TokenStream2::new()
    };

    let alias_filter = if opts.alias_filter {
        let struct_ident = &opts.ident;
        let ident = format_ident!("{}Filter", &opts.ident);

        quote! {
            #vis type #ident #generic_params = ::meteor::csv::filter::Filter<
                <#struct_ident #generic_args as ::meteor::csv::filter::FilterableRecord>::Field
            >;
        }
    } else {
        TokenStream2::new()
    };

    quote! {
        #alias_field
        #alias_filter
    }
}

#[proc_macro_derive(FilterableRecord, attributes(filter))]
pub fn derive_filterable_record(input: TokenStream) -> TokenStream {
    let parsed = parse_macro_input!(input as DeriveInput);

    let (opts, fields) = match FilterableRecordOpts::from_derive_input(&parsed) {
        Ok(opts) => {
            let fields = match &opts.data {
                ast::Data::Enum(_) => {
                    panic!("FilterableRecord can only be derived for named structs")
                }
                ast::Data::Struct(fields) => fields.fields.clone(),
            };
            (opts, fields)
        }
        Err(e) => {
            return e.write_errors().into();
        }
    };

    let fields_enum = fields_enum_variants(&opts, &fields);
    let from_str_field_impl = gen_fields_enum_from_str_field(&fields_enum);

    let polars_literal_impl = if opts.polars_literal {
        gen_fields_enum_polars_literal(&fields_enum)
    } else {
        TokenStream2::new()
    };

    let filterable_record_impl = gen_filterable_record_impl(&opts, &fields_enum);

    let aliases = gen_aliases(&opts);

    let output = quote! {
        #[doc(hidden)]
        const _: () = {
            #fields_enum

            #from_str_field_impl

            #polars_literal_impl

            #filterable_record_impl
        };

        #aliases
    };

    output.into()
}

#[derive(Clone, Debug, FromField)]
#[darling(attributes(polars))]
struct PolarsSchemaField {
    ident: Option<syn::Ident>,
    ty: syn::Type,

    #[darling(default)]
    dtype: Option<syn::Path>,

    #[darling(default)]
    rename: Option<String>,

    #[darling(default)]
    skip: bool,
}

#[derive(Debug, FromDeriveInput)]
#[darling(supports(struct_named), attributes(polars))]
struct PolarsSchemaOpts {
    ident: syn::Ident,
    generics: syn::Generics,

    set_csv_options: Option<syn::Expr>,

    data: ast::Data<(), PolarsSchemaField>,
}

#[proc_macro_derive(PolarsSchema, attributes(polars))]
pub fn derive_polars_schema(input: TokenStream) -> TokenStream {
    let parsed = parse_macro_input!(input as DeriveInput);

    let (opts, fields) = match PolarsSchemaOpts::from_derive_input(&parsed) {
        Ok(opts) => {
            let fields = match &opts.data {
                ast::Data::Enum(_) => panic!("PolarsSchema can only be derived for named structs"),
                ast::Data::Struct(fields) => fields.fields.clone(),
            };
            (opts, fields)
        },
        Err(e) => {
            return e.write_errors().into();
        },
    };

    let record_ident = &opts.ident;

    let (generic_params, generic_args) = split_generics(&opts.generics, |type_param| {
        let type_ident = &type_param.ident;

        quote! { #type_ident : ::meteor::csv::polars::PolarsCompatible }
    });

    let (const_idents, const_decls): (Vec<_>, Vec<_>) = fields.iter().filter_map(|field| {
        if field.skip {
            return None
        }

        let field_name = field.ident.as_ref().unwrap().to_string();
        let literal = field.rename.as_ref().unwrap_or(&field_name).to_owned();
        let const_ident = Ident2::new(&to_screaming_snake_case(&field_name), Span2::call_site());

        let field_ty = &field.ty;

        let dtype = if let Some(chosen_dtype) = &field.dtype {
            quote! { &#chosen_dtype }
        } else {
            quote! { &<#field_ty as ::meteor::csv::polars::PolarsCompatible>::DTYPE }
        };

        let const_decl = quote! {
            pub const #const_ident: (&'static str, &'static ::polars::prelude::DataType) = (#literal, #dtype);
        };

        Some((const_ident, const_decl))
    })
    .unzip();

    let set_csv_options = if let Some(csv_options) = &opts.set_csv_options {
        quote! {
            fn set_csv_options(__options: ::polars::prelude::LazyCsvReader) -> ::polars::prelude::LazyCsvReader {
                let __setter = ( #csv_options );

                __setter(__options)
            }
        }
    } else {
        TokenStream2::new()
    };

    let output = quote! {
        impl #generic_params #record_ident #generic_args {
            #( #const_decls )*
        }

        impl #generic_params ::meteor::csv::polars::PolarsSchema for #record_ident #generic_args {
            const ALL_FIELDS: &'static [(&'static str, &'static ::polars::prelude::DataType)] = &[
                #( Self::#const_idents ),*
            ];

            fn polars_schema() -> ::polars::prelude::Schema {
                let mut schema = ::std::vec::Vec::<(::polars::datatypes::PlSmallStr, ::polars::prelude::DataType)>::new();

                for &(name, dtype) in Self::ALL_FIELDS {
                    schema.push((name.into(), dtype.clone()));
                }

                schema.into_iter().collect()
            }

            #set_csv_options
        }
    };

    output.into()
}

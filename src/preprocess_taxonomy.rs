use std::path::PathBuf;

use anyhow::{anyhow, Context};
use clap::{Args, ValueEnum};
use itertools::Itertools;

use crate::{
    preprocessed_taxonomy::{PreprocessedTaxonomy, PreprocessedTaxonomyFormat, SomeTaxonomy},
    taxonomy::{
        formats::{
            ncbi::{AllNames, NamesAssoc, NcbiTaxonomy},
            newick::NewickTaxonomy,
        },
        Taxonomy, TaxonomyMut,
    },
};

#[derive(ValueEnum, Debug, Copy, Clone)]
pub enum TaxonomyFormat {
    Ncbi,
    Newick,
}

impl Default for TaxonomyFormat {
    fn default() -> Self {
        Self::Ncbi
    }
}

/// Load, preprocess and serialize a taxonomy for further Meteor usage.
#[derive(Args)]
pub struct PreprocessTaxonomyArgs {
    /// Enables contraction. If provided, RANKS should be a comma-separated list of the ranks to keep after contraction.
    #[clap(short, long, value_name = "RANKS")]
    contract: Option<Option<String>>,

    /// Input taxonomy format
    #[clap(long, value_enum, default_value_t)]
    input_format: TaxonomyFormat,

    /// Path to the input taxonomy (in the case of NCBI, specify the extracted directory of taxdump).
    /// To download the NCBI taxonomy into "ncbi-taxonomy", run `meteor fetch ncbi-taxonomy -d ncbi-taxonomy`
    /// (see `meteor fetch --help` for more information)
    input_taxonomy: String,

    /// Comma-separated list of name classes to load (if INPUT_FORMAT is ncbi)
    #[clap(long, default_value = "scientific name,synonym,equivalent name,common name,acronym,blast name,genbank acronym,genbank common name,includes,in-part" )]
    ncbi_name_classes: String,

    /// Comma-separated list of rank names to associate to each level of the Newick taxonomy (if
    /// INPUT_FORMAT is newick)
    #[clap(long)]
    newick_ranks: Option<String>,

    /// Format in which the preprocessed taxonomy should be serialized.
    #[clap(long, value_enum, default_value_t)]
    output_format: PreprocessedTaxonomyFormat,

    /// Output taxonomy path
    output_taxonomy: String,
}

pub struct PreprocessNcbiTaxonomyArgs {
    taxdump_path: String,
    name_classes: Vec<String>,
    contraction_ranks: Option<Vec<String>>,
}

pub struct PreprocessNewickTaxonomyArgs {
    path: String,
    ranks: Vec<String>,
    contract: bool,
}

impl PreprocessTaxonomyArgs {
    pub fn get_contraction_ranks(&self) -> Option<Vec<String>> {
        const DEFAULT_CONTRACTION_RANKS: &str =
            "superkingdom,phylum,class,order,family,genus,species";

        self.contract.as_ref().map(|ranks| {
            ranks
                .as_deref()
                .unwrap_or(DEFAULT_CONTRACTION_RANKS)
                .split(',')
                .map(|rank| rank.trim().to_owned())
                .collect()
        })
    }

    pub fn into_ncbi_args(self) -> PreprocessNcbiTaxonomyArgs {
        let contraction_ranks = self.get_contraction_ranks();

        PreprocessNcbiTaxonomyArgs {
            taxdump_path: self.input_taxonomy,
            name_classes: self
                .ncbi_name_classes
                .split(',')
                .map(ToOwned::to_owned)
                .collect(),
            contraction_ranks,
        }
    }

    pub fn into_newick_args(self) -> anyhow::Result<PreprocessNewickTaxonomyArgs> {
        let ranks = self
            .newick_ranks
            .as_deref()
            .ok_or_else(|| anyhow!("Flag --newick-ranks is required for --input-format newick"))?
            .split(',')
            .map(ToOwned::to_owned)
            .collect();

        Ok(PreprocessNewickTaxonomyArgs {
            path: self.input_taxonomy,
            ranks,
            contract: self.contract.is_some(),
        })
    }
}

fn rank_syms_to_strings<Tax: Taxonomy>(
    tax: &Tax,
    syms: &[Tax::RankSym],
) -> anyhow::Result<Vec<String>> {
    syms.iter()
        .map(|&sym| tax.rank_sym_str(sym).map(ToOwned::to_owned))
        .collect::<Option<Vec<String>>>()
        .ok_or_else(|| anyhow!("Unable to translate ordered_ranks to Vec<String>. This is a bug"))
}

fn contract_and_order_ranks<Names: NamesAssoc + Send + 'static>(
    args: &PreprocessNcbiTaxonomyArgs,
    tax: &mut NcbiTaxonomy<Names>,
) -> anyhow::Result<Option<Vec<String>>> {
    let ordered_ranks = if let Some(ranks) = &args.contraction_ranks {
        let contraction_ranks_syms: Vec<_> = ranks
            .iter()
            .map(|rank| {
                tax.lookup_rank_sym(rank)
                    .ok_or_else(|| anyhow!("Unrecognized rank {rank}"))
            })
            .try_collect()
            .context("Error during contraction")?;

        tax.contract(contraction_ranks_syms.iter().copied().collect());

        let rank_ordering = tax.rank_ordering();

        if let Ok(rank_ordering) = &rank_ordering {
            assert_eq!(
                rank_syms_to_strings(tax, rank_ordering)?,
                rank_syms_to_strings(tax, &contraction_ranks_syms)?,
            );
        }

        assert!(tax.topology_health_check());

        rank_ordering
    } else {
        tax.rank_ordering()
    };

    let ordered_ranks = match ordered_ranks {
        Ok(ordered_ranks) => Some(
            rank_syms_to_strings(tax, &ordered_ranks)
                .context("Error mapping rank symbols back to strings. This is bug")?,
        ),

        Err(topsort_loop) => {
            eprintln!(
                "Warning: the taxonomy contains rank loops: {:?}",
                rank_syms_to_strings(tax, &topsort_loop.0),
            );
            None
        }
    };

    Ok(ordered_ranks)
}

fn preprocess_ncbi(args: &PreprocessNcbiTaxonomyArgs) -> anyhow::Result<PreprocessedTaxonomy> {
    let mut path = PathBuf::from(&args.taxdump_path);

    let mut tax = {
        path.push("nodes.dmp");
        let tax = NcbiTaxonomy::load_nodes(&path)?;
        path.pop();
        eprintln!("Loaded nodes.dmp");

        path.push("merged.dmp");
        let tax = tax.with_merged_taxids(&path)?;
        path.pop();
        eprintln!("Loaded merged.dmp");

        path.push("names.dmp");
        let names = AllNames::load_filtered_names_dmp(&args.name_classes, &path)?;
        path.pop();
        eprintln!("Loaded names.dmp with classes {:?}", &args.name_classes);

        tax.with_names(names)
    };

    let ordered_ranks = contract_and_order_ranks(args, &mut tax)?;

    Ok(PreprocessedTaxonomy::new(
        if args.name_classes.len() == 1 {
            match std::mem::take(&mut tax.names).only_of_class(&args.name_classes[0]) {
                Ok(names) => SomeTaxonomy::NcbiTaxonomyWithSingleClassNames(tax.with_names(names)),
                Err(names) => SomeTaxonomy::NcbiTaxonomyWithManyNames(tax.with_names(names)),
            }
        } else {
            SomeTaxonomy::NcbiTaxonomyWithManyNames(tax)
        },
        ordered_ranks,
        args.contraction_ranks.clone(),
    ))
}

fn preprocess_newick(args: &PreprocessNewickTaxonomyArgs) -> anyhow::Result<PreprocessedTaxonomy> {
    let tax = NewickTaxonomy::load_newick(&args.path, args.ranks.clone())?;

    if args.contract {
        anyhow::bail!("Contraction is not supported for Newick taxonomies");
    }

    let contraction_ranks = if args.contract {
        Some(args.ranks.clone())
    } else {
        None
    };

    Ok(PreprocessedTaxonomy::new(
        SomeTaxonomy::NewickTaxonomy(tax),
        Some(args.ranks.clone()),
        contraction_ranks,
    ))
}

pub fn preprocess_taxonomy(args: PreprocessTaxonomyArgs) -> anyhow::Result<()> {
    let output_taxonomy = args.output_taxonomy.clone();
    let output_format = args.output_format;

    let taxonomy = match args.input_format {
        TaxonomyFormat::Ncbi => preprocess_ncbi(&args.into_ncbi_args()),
        TaxonomyFormat::Newick => preprocess_newick(&args.into_newick_args()?),
    }?;

    taxonomy.serialize_with_format(&output_taxonomy, output_format)?;

    Ok(())
}

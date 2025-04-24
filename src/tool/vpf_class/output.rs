use std::{borrow::Cow, ops::Deref, path::Path};

use itertools::Itertools;
use lending_iterator::HKT;
use polars::{error::PolarsResult, prelude::LazyFrame, series::Series};
use serde::{ser::SerializeStruct, Deserialize, Serialize};

use crate::{
    csv::{
        filter::{Filter, FilterableRecord},
        flatten_fix::SerializeFlat,
        polars::{lazyframe_from_file, ExprExt, PolarsSchema},
    },
    taxonomy::{tree, LabeledTaxonomy, NodeId, Taxonomy},
};

#[derive(Debug, Serialize, Deserialize, FilterableRecord, PolarsSchema)]
#[filter(alias_field, alias_filter, polars_literal)]
pub struct VpfClassRecord<'a> {
    #[filter(as_string)]
    #[serde(borrow)]
    pub virus_name: Cow<'a, str>,
    #[filter(as_string)]
    #[serde(borrow)]
    pub class_name: Cow<'a, str>,
    pub membership_ratio: f32,
    pub virus_hit_score: f32,
    pub confidence_score: f32,
}

pub type VpfClassRecordHKT = HKT!(VpfClassRecord<'_>);

impl<'a> SerializeFlat for VpfClassRecord<'a> {
    const FIELD_COUNT: usize = 5;

    fn serialize_flat<Ser>(&self, row: &mut Ser) -> Result<(), Ser::Error>
    where
        Ser: SerializeStruct,
    {
        row.serialize_field("virus_name", &*self.virus_name)?;
        row.serialize_field("class_name", &*self.class_name)?;
        row.serialize_field("membership_ratio", &self.membership_ratio)?;
        row.serialize_field("virus_hit_score", &self.virus_hit_score)?;
        row.serialize_field("confidence_score", &self.confidence_score)?;
        Ok(())
    }
}

#[derive(Copy, Clone, Eq, PartialEq, PartialOrd, Ord, Hash, clap::ValueEnum)]
pub enum VpfClassVirusRank {
    Baltimore,
    Family,
    Genus
}

#[derive(Copy, Clone, Eq, PartialEq, PartialOrd, Ord, Hash, clap::ValueEnum)]
pub enum VpfClassHostRank {
    Phylum,
    Family,
    Genus
}

#[derive(Copy, Clone, Eq, PartialEq, PartialOrd, Ord, Hash)]
pub enum VpfClassRank {
    Virus(VpfClassVirusRank),
    Host(VpfClassHostRank),
}

impl From<VpfClassVirusRank> for VpfClassRank {
    fn from(value: VpfClassVirusRank) -> Self {
        Self::Virus(value)
    }
}

impl From<VpfClassHostRank> for VpfClassRank {
    fn from(value: VpfClassHostRank) -> Self {
        Self::Host(value)
    }
}

impl VpfClassRank {
    pub fn to_rank_str(self) -> &'static str {
        let self_str = self.to_str();
        self_str.strip_prefix("host_").unwrap_or(self_str)
    }

    pub fn to_str(self) -> &'static str {
        type V = VpfClassVirusRank;
        type H = VpfClassHostRank;

        match self {
            VpfClassRank::Virus(V::Baltimore) => "baltimore",
            VpfClassRank::Virus(V::Family) => "family",
            VpfClassRank::Virus(V::Genus) => "genus",
            VpfClassRank::Host(H::Phylum) => "host_phylum",
            VpfClassRank::Host(H::Family) => "host_family",
            VpfClassRank::Host(H::Genus) => "host_genus",
        }
    }

    pub fn lookup_sym<Tax: Taxonomy>(self, tax: &Tax) -> anyhow::Result<Tax::RankSym> {
        tax.lookup_rank_sym(self.to_rank_str()).ok_or_else(|| {
            anyhow::anyhow!(
                "Could not find VPF-Class rank {} in the taxonomy",
                self.to_rank_str()
            )
        })
    }
}

impl std::fmt::Display for VpfClassRank {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.write_str(self.to_str())
    }
}

impl clap::ValueEnum for VpfClassRank {
    fn value_variants<'a>() -> &'a [Self] {
        use VpfClassRank::{Virus, Host};

        &[
            Virus(VpfClassVirusRank::Baltimore),
            Virus(VpfClassVirusRank::Family),
            Virus(VpfClassVirusRank::Genus),
            Host(VpfClassHostRank::Phylum),
            Host(VpfClassHostRank::Family),
            Host(VpfClassHostRank::Genus),
        ]
    }

    fn to_possible_value(&self) -> Option<clap::builder::PossibleValue> {
        Some(clap::builder::PossibleValue::new(self.to_str()))
    }
}

pub fn lazyframe_from_rank(parent: &Path, rank: VpfClassRank) -> PolarsResult<LazyFrame> {
    let mut buf = parent.to_path_buf();
    let file = format!("{rank}.tsv");
    buf.push(&file);

    lazyframe_from_file::<VpfClassRecord>(&buf)
}

pub fn class_name_to_taxids<Tax>(tax: &Tax, rank_sym: Tax::RankSym, name: &str) -> Option<Vec<NodeId>>
where
    Tax: LabeledTaxonomy,
{
    tax.nodes_with_label_and_rank(rank_sym, name).map(Vec::from)
}

/**
 * Input columns:
 *   - `virus_name: str, {col_names}+: [t]+`
 *
 * Output columns:
 *   - `virus_name: str, {col_names}+: [t]+`
*/
pub fn join_lists_by_virus_name<IS, S>(df1: LazyFrame, df2: LazyFrame, col_names: IS) -> LazyFrame
where
    IS: IntoIterator<Item = S>,
    S: AsRef<str>,
{
    use polars::lazy::dsl::{col, cols};

    let col_names = col_names
        .into_iter()
        .map(|s| s.as_ref().to_owned())
        .collect_vec();

    let cols1 = cols(&col_names).name().suffix("1");
    let df1 = df1.select([col("virus_name"), cols1]);

    let cols2 = cols(&col_names).name().suffix("2");
    let df2 = df2.select([col("virus_name"), cols2]);

    let join = crate::csv::polars::full_join_builder(df1, df2)
        .on([col("virus_name")])
        .finish();

    let mut selection = col_names
        .iter()
        .map(|col_name| {
            let col1 = format!("{col_name}1");
            let col2 = format!("{col_name}2");
            crate::csv::polars::col_list_concat(&col1, &col2).alias(col_name)
        })
        .collect_vec();

    selection.insert(0, col("virus_name"));

    join.select(&selection)
}

/**
 * Output columns:
 *   - `<VpfClassRecord>, taxids: [u32]`
*/
#[expect(dead_code)]
pub fn lazyframe_with_taxids_from_rank<TaxPtr>(
    output_dir: &Path,
    filters: &[VpfClassRecordFilter],
    tax: TaxPtr,
    rank: VpfClassRank,
) -> anyhow::Result<LazyFrame>
where
    TaxPtr: Deref<Target: LabeledTaxonomy> + Send + Sync + 'static
{
    use polars::lazy::dsl::col;

    let mut df = lazyframe_from_rank(output_dir, rank)?;

    df = Filter::apply_polars(df, filters);

    let rank_sym = tax
        .lookup_rank_sym(rank.to_rank_str())
        .ok_or_else(|| anyhow::anyhow!("Could not find rank {rank} in the taxonomy"))?;

    df = df.with_column(
        col("class_name")
            .downcast_map_apply_to_vec(Series::str, move |class_name| {
                class_name_to_taxids(&*tax, rank_sym, class_name)
                    .map(tree::nodeid_to_u32_vec)
            })
            .alias("taxids"),
    );

    Ok(df)
}


#[cfg(test)]
pub(crate) mod tests {
    use std::sync::Arc;

    use polars::{prelude::LazyFrame, series::Series};

    use crate::{
        csv::polars::{lazyframe_from_string_strip, series_vec, ExprExt},
        taxonomy::{
            formats::newick::{self, NewickTaxonomy},
            tree, NodeId, Taxonomy,
        },
        tool::vpf_class::VpfClassRecord,
    };

    pub const FAMILY: &str = "\
        virus_name\tclass_name\tmembership_ratio\tvirus_hit_score\tconfidence_score
        v1        \t f5       \t 0.5            \t 0             \t 0.5
        v1        \t f13      \t 0.5            \t 0             \t 0.5
        v2        \t f13      \t 1              \t 0             \t 0.5
        v3        \t g6       \t 0.3            \t 0             \t 0.5
        v3        \t xxx      \t 0.3            \t 0             \t 0.5
        v3        \t f22      \t 0.3            \t 0             \t 0.5
        v3        \t f30      \t 0.4            \t 0             \t 0.5
        v4        \t f30      \t 1              \t 0             \t 0.5
        v4        \t f30      \t 1              \t 0             \t 0.5
    ";

    pub const GENUS: &str = "\
        virus_name\tclass_name\tmembership_ratio\tvirus_hit_score\tconfidence_score
        v1        \t f5       \t 0.5            \t 0             \t 0.5
        v1        \t g10      \t 0.5            \t 0             \t 0.5
        v3        \t g6       \t 0.3            \t 0             \t 0.5
        v3        \t g14      \t 0.3            \t 0             \t 0.5
        v3        \t g17      \t 0.3            \t 0             \t 0.5
        v3        \t g23      \t 0.1            \t 0             \t 0.5
        v4        \t g10      \t 0.5            \t 0             \t 0.5
        v4        \t g31      \t 0.5            \t 0             \t 0.5
        v5        \t g6       \t 0.5            \t 0             \t 0.5
        v5        \t s11      \t 0.5            \t 0             \t 0.5
        v6        \t yyy      \t 0.5            \t 0             \t 0.5
    ";

    #[test]
    fn test_class_name_to_taxid() {
        let tax = newick::tests::ambiguous_taxonomy();
        let genus_sym = tax.lookup_rank_sym("genus").unwrap();

        assert_eq!(super::class_name_to_taxids(&tax, genus_sym, "x"), Some(vec![NodeId(6)]));

        assert_eq!(super::class_name_to_taxids(&tax, genus_sym, "13"), Some(vec![NodeId(13)]));
    }

    fn test_agg_lazyframe(df: LazyFrame, tax: Arc<NewickTaxonomy>, rank_str: &str) -> LazyFrame {
        use polars::lazy::dsl::col;

        let rank_sym = tax.lookup_rank_sym(rank_str).unwrap();

        df.with_column(
            col("class_name")
                .downcast_map_apply_to_vec(Series::str, move |class_name| {
                    super::class_name_to_taxids(&*tax, rank_sym, class_name)
                        .map(tree::nodeid_to_u32_vec)
                })
                .alias("taxids"),
        )
        .group_by([col("virus_name")])
        .agg([col("taxids").explode()])
    }

    #[test]
    fn test_class_names_to_taxids() -> anyhow::Result<()> {
        let tax = Arc::new(newick::tests::sample_taxonomy());

        let family = lazyframe_from_string_strip::<VpfClassRecord>(FAMILY)?;

        let family = test_agg_lazyframe(family, Arc::clone(&tax), "family")
            .sort(["virus_name"], Default::default())
            .collect()?;

        let family_expected = polars::df!(
            "virus_name" => ["v1", "v2", "v3", "v4"],
            "taxids" => series_vec([
                vec![Some(5), Some(13)],
                vec![Some(13)],
                vec![Some(6), None, Some(22), Some(30)],
                vec![Some(30), Some(30)],
            ])
        )?;

        assert_eq!(family, family_expected);

        let genus = lazyframe_from_string_strip::<VpfClassRecord>(GENUS)?;

        let genus = test_agg_lazyframe(genus, Arc::clone(&tax), "genus")
            .sort(["virus_name"], Default::default())
            .collect()?;

        let genus_expected = polars::df!(
            "virus_name" => ["v1", "v3", "v4", "v5", "v6"],
            "taxids" => series_vec([
                vec![Some(5), Some(10)],
                vec![Some(6), Some(14), Some(17), Some(23)],
                vec![Some(10), Some(31)],
                vec![Some(6), Some(11)],
                vec![None],
            ])
        )?;

        assert_eq!(genus, genus_expected);

        Ok(())
    }

    #[test]
    fn test_join_ranks() -> anyhow::Result<()> {
        let tax = Arc::new(newick::tests::sample_taxonomy());

        let family = lazyframe_from_string_strip::<VpfClassRecord>(FAMILY)?;
        let family = test_agg_lazyframe(family, Arc::clone(&tax), "family");

        let genus = lazyframe_from_string_strip::<VpfClassRecord>(GENUS)?;
        let genus = test_agg_lazyframe(genus, Arc::clone(&tax), "genus");

        let joined = super::join_lists_by_virus_name(family, genus, ["taxids"]).collect()?;

        let expected = polars::df!(
            "virus_name" => ["v1", "v2", "v3", "v4", "v5", "v6"],
            "taxids"     => series_vec([
                vec![Some(5u32), Some(13), Some(5), Some(10)],
                vec![Some(13)],
                vec![Some(6), None, Some(22), Some(30), Some(6), Some(14), Some(17), Some(23)],
                vec![Some(30), Some(30), Some(10), Some(31)],
                vec![Some(6), Some(11)],
                vec![None],
            ])
        )?;

        assert_eq!(joined.sort(["virus_name"], Default::default())?, expected,);

        Ok(())
    }
}

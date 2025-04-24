use std::collections::HashSet;

use itertools::Itertools;
use serde::{ser::SerializeStruct, Serialize, Serializer};

use crate::{
    crispr_match::VirusHostMapping, csv::flatten_fix::{serialize_flat_struct, SerializeFlat}, taxonomy::tree::node_id::NodeIdSet, tool::vpf_class::VpfClassRecord
};

pub trait Enrichment<'a>: Default {
    type Context: Copy;
    type CsvFields: SerializeFlat + From<Self>;

    type SummaryClassData: Default + Clone;
    type SummaryClassStats: SerializeFlat + From<Self::SummaryClassData>;

    fn enrich<'r>(
        &mut self,
        context: Self::Context,
        vpf_class_record: &VpfClassRecord<'r>,
        assigned_taxids: &NodeIdSet,
        assigned_contigs: &HashSet<&'a str>,
    ) -> bool;

    fn add_class_data(&self, class_data: &mut Self::SummaryClassData);

    fn into_csv_fields(self) -> Self::CsvFields {
        self.into()
    }
}

#[derive(Copy, Clone)]
pub struct NoEnrichmentContext;

#[derive(Default, Clone)]
pub struct NoEnrichmentData;

impl SerializeFlat for NoEnrichmentData {
    const FIELD_COUNT: usize = 1;

    fn serialize_flat<Ser: SerializeStruct>(&self, _row: &mut Ser) -> Result<(), Ser::Error> {
        Ok(())
    }
}

#[derive(Default, Clone)]
pub struct NoEnrichment;

impl SerializeFlat for NoEnrichment {
    const FIELD_COUNT: usize = 1;

    fn serialize_flat<Ser: SerializeStruct>(&self, _row: &mut Ser) -> Result<(), Ser::Error> {
        Ok(())
    }
}

impl<'a> Enrichment<'a> for NoEnrichment {
    type Context = NoEnrichmentContext;
    type CsvFields = NoEnrichment;
    type SummaryClassData = NoEnrichmentData;
    type SummaryClassStats = NoEnrichmentData;

    #[must_use]
    fn enrich<'r>(
        &mut self,
        _context: Self::Context,
        _vpf_class_record: &VpfClassRecord<'r>,
        _assigned_taxids: &NodeIdSet,
        _assigned_contigs: &HashSet<&'a str>,
    ) -> bool {
        true
    }

    fn add_class_data(&self, _class_data: &mut Self::SummaryClassData) {
    }
}

#[derive(Default, Clone, Serialize)]
pub struct CrisprEnrichment<'a> {
    pub crispr_matches: HashSet<&'a str>,
}

pub struct CrisprMatchesCsvFields {
    pub crispr_matches: String,
}

impl<'a> From<CrisprEnrichment<'a>> for CrisprMatchesCsvFields {
    fn from(value: CrisprEnrichment) -> Self {
        Self {
            crispr_matches: value.crispr_matches.iter().join(";"),
        }
    }
}

impl SerializeFlat for CrisprMatchesCsvFields {
    const FIELD_COUNT: usize = 1;

    fn serialize_flat<Ser: SerializeStruct>(&self, row: &mut Ser) -> Result<(), Ser::Error> {
        row.serialize_field("crispr_matches", &self.crispr_matches)
    }
}

#[derive(Default, Clone)]
pub struct CrisprEnrichmentSummaryClassData<'a> {
    pub crispr_match_contigs: HashSet<&'a str>,
}

#[derive(Default, Copy, Clone)]
pub struct CrisprEnrichmentSummaryClassStats {
    pub num_crispr_matched_contigs: usize,
}

impl<'a> From<CrisprEnrichmentSummaryClassData<'a>> for CrisprEnrichmentSummaryClassStats {
    fn from(value: CrisprEnrichmentSummaryClassData<'a>) -> Self {
        Self {
            num_crispr_matched_contigs: value.crispr_match_contigs.len(),
        }
    }
}

impl SerializeFlat for CrisprEnrichmentSummaryClassStats {
    const FIELD_COUNT: usize = 1;

    fn serialize_flat<Ser: SerializeStruct>(&self, row: &mut Ser) -> Result<(), Ser::Error> {
        row.serialize_field("num_crispr_matched_contigs", &self.num_crispr_matched_contigs)
    }
}

impl<'a> Enrichment<'a> for CrisprEnrichment<'a> {
    type Context = &'a VirusHostMapping;
    type CsvFields = CrisprMatchesCsvFields;

    fn enrich<'r>(
        &mut self,
        context: &'a VirusHostMapping,
        vpf_class_record: &VpfClassRecord<'r>,
        _assigned_taxids: &NodeIdSet,
        assigned_contigs: &HashSet<&'a str>,
    ) -> bool {
        self.crispr_matches.extend(
            context
                .0
                .lookup(vpf_class_record.virus_name.as_ref())
                .filter(|contig| assigned_contigs.contains(contig))
        );

        !self.crispr_matches.is_empty()
    }

    type SummaryClassData = CrisprEnrichmentSummaryClassData<'a>;
    type SummaryClassStats = CrisprEnrichmentSummaryClassStats;

    fn add_class_data(&self, class_data: &mut Self::SummaryClassData) {
        class_data.crispr_match_contigs.extend(&self.crispr_matches);
    }
}


pub struct EnrichedVpfClassRecord<'a, 'r, CE> {
    pub vpf_class_record: VpfClassRecord<'r>,
    pub assigned_taxids: NodeIdSet,
    pub assigned_contigs: HashSet<&'a str>,
    pub crispr_enrichment: CE,
}

impl<'a, 'r, CE: Default> EnrichedVpfClassRecord<'a, 'r, CE> {
    pub fn new(vpf_class_record: VpfClassRecord<'r>) -> Self {
        EnrichedVpfClassRecord {
            vpf_class_record,
            assigned_taxids: NodeIdSet::new(),
            assigned_contigs: HashSet::new(),
            crispr_enrichment: CE::default(),
        }
    }

    pub fn has_assignments(&self) -> bool {
        !self.assigned_taxids.is_empty() || !self.assigned_contigs.is_empty()
    }
}

#[derive(Debug)]
pub struct CsvEnrichedVpfClassRecord<'a, 'r, CE: Enrichment<'a>> {
    pub vpf_class_record: VpfClassRecord<'r>,

    pub assigned_taxids: String,
    pub num_assigned_contigs: usize,

    pub crispr_match_fields: CE::CsvFields,
}

impl<'a, 'r, CE> Serialize for CsvEnrichedVpfClassRecord<'a, 'r, CE>
where
    CE: Enrichment<'a>,
{
    fn serialize<Ser: Serializer>(&self, serializer: Ser) -> Result<Ser::Ok, Ser::Error> {
        serialize_flat_struct(serializer, "CsvEnrichedVpfClassRecord", self)
    }
}

impl<'a, 'r, CE> SerializeFlat for CsvEnrichedVpfClassRecord<'a, 'r, CE>
where
    CE: Enrichment<'a>,
{
    const FIELD_COUNT: usize = VpfClassRecord::FIELD_COUNT + 2 + CE::CsvFields::FIELD_COUNT;

    fn serialize_flat<Ser: SerializeStruct>(&self, row: &mut Ser) -> Result<(), Ser::Error> {
        self.vpf_class_record.serialize_flat(row)?;

        row.serialize_field("assigned_taxids", &self.assigned_taxids)?;
        row.serialize_field("num_assigned_contigs", &self.num_assigned_contigs)?;

        self.crispr_match_fields.serialize_flat(row)?;

        Ok(())
    }
}

impl<'a, 'r, CE> From<EnrichedVpfClassRecord<'a, 'r, CE>> for CsvEnrichedVpfClassRecord<'a, 'r, CE>
where
    CE: Enrichment<'a>,
{
    fn from(other: EnrichedVpfClassRecord<'a, 'r, CE>) -> Self {
        CsvEnrichedVpfClassRecord {
            vpf_class_record: other.vpf_class_record,
            assigned_taxids: other.assigned_taxids.iter().join(";"),
            num_assigned_contigs: other.assigned_contigs.len(),
            crispr_match_fields: other.crispr_enrichment.into_csv_fields(),
        }
    }
}

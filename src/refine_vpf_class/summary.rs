use std::{
    collections::{HashMap, HashSet},
    io,
};

use clap::ValueEnum;
use itertools::Itertools;
use serde::{ser::SerializeStruct, Serialize, Serializer};

use crate::{
    csv::flatten_fix::{serialize_flat_struct, SerializeFlat},
    taxonomy::tree::node_id::NodeIdSet,
};

use super::enrichment::{CrisprEnrichment, EnrichedVpfClassRecord, Enrichment, NoEnrichment};


#[derive(Clone, Copy, Eq, PartialEq, Ord, PartialOrd, Hash)]
pub enum RecordDropReason {
    Filtered,
    UnknownHostTaxId,
    NotInAssignment,
    NoCrisprInfo,
}

#[derive(Default)]
pub struct RecordDropStats {
    num_filtered: u32,
    num_unknown_host_taxid: u32,
    num_not_in_assignment: u32,
    num_no_crispr_info: u32,
}

struct ClassData<'a, CE: Enrichment<'a>> {
    virus_count: u32,
    assigned_taxids: NodeIdSet,
    assigned_contigs: HashSet<&'a str>,
    crispr_enrichment_data: CE::SummaryClassData,
}

struct ClassStats<'a, CE: Enrichment<'a>> {
    class_name: String,
    virus_count: i32,
    num_assigned_taxids: i32,
    num_assigned_contigs: i32,
    crispr_enrichment_stats: CE::SummaryClassStats,
}

impl RecordDropStats {
    pub fn total_dropped(&self) -> u32 {
        self.num_filtered + self.num_unknown_host_taxid + self.num_not_in_assignment + self.num_no_crispr_info
    }
}

impl<'a, CE: Enrichment<'a>> Default for ClassData<'a, CE> {
    fn default() -> Self {
        Self {
            virus_count: Default::default(),
            assigned_taxids: Default::default(),
            assigned_contigs: Default::default(),
            crispr_enrichment_data: Default::default(),
        }
    }
}

impl<'a, CE: Enrichment<'a>> ClassData<'a, CE> {
    fn add<'r>(&mut self, record: &EnrichedVpfClassRecord<'a, 'r, CE>) {
        self.virus_count += 1;
        self.assigned_taxids.extend(&record.assigned_taxids);
        self.assigned_contigs.extend(&record.assigned_contigs);

        record.crispr_enrichment.add_class_data(&mut self.crispr_enrichment_data);
    }
}

impl<'a, CE: Enrichment<'a>> SerializeFlat for ClassStats<'a, CE> {
    const FIELD_COUNT: usize = 4 + CE::SummaryClassStats::FIELD_COUNT;

    fn serialize_flat<Ser: SerializeStruct>(&self, row: &mut Ser) -> Result<(), Ser::Error> {
        row.serialize_field("class_name", &self.class_name)?;
        row.serialize_field("virus_count", &self.virus_count)?;
        row.serialize_field("num_assigned_taxids", &self.num_assigned_taxids)?;
        row.serialize_field("num_assigned_contigs", &self.num_assigned_contigs)?;
        self.crispr_enrichment_stats.serialize_flat(row)?;
        Ok(())
    }
}

impl<'a, CE: Enrichment<'a>> Serialize for ClassStats<'a, CE> {
    fn serialize<S: Serializer>(&self, serializer: S) -> Result<S::Ok, S::Error> {
        serialize_flat_struct(serializer, "ClassStats", self)
    }
}


#[derive(ValueEnum, Debug, Copy, Clone, PartialEq, Eq, PartialOrd, Ord)]
pub enum SummarySortBy {
    ClassName,
    VirusCount,
    AssignedTaxIds,
    AssignedContigs,
}

pub struct EnrichmentSummary<'a, CE: Enrichment<'a>> {
    class_data: HashMap<String, ClassData<'a, CE>>,
    dropped: RecordDropStats,
    num_records: u32,
}

impl<'a, CE: Enrichment<'a>> Default for EnrichmentSummary<'a, CE> {
    fn default() -> Self {
        Self {
            class_data: Default::default(),
            dropped: Default::default(),
            num_records: 0,
        }
    }
}

impl<'a, CE: Enrichment<'a>> EnrichmentSummary<'a, CE> {
    pub fn num_records(&self) -> u32 {
        self.num_records
    }

    pub fn account_kept(&mut self, _verbose: bool) {
        self.num_records += 1;
    }

    pub fn account_dropped(&mut self, reason: RecordDropReason, verbose: bool) {
        self.num_records += 1;

        use RecordDropReason::*;

        macro_rules! explain {
            ($fmt_text:literal) => {{
                if verbose {
                    eprintln!($fmt_text, self.num_records);
                }
            }}
        }

        match reason {
            Filtered => {
                self.dropped.num_filtered += 1;
                explain!("Record #{} dropped by --filter");
            }
            UnknownHostTaxId => {
                self.dropped.num_unknown_host_taxid += 1;
                explain!("Record #{} dropped because the host prediction had no matching taxids");
            }
            NotInAssignment => {
                self.dropped.num_not_in_assignment += 1;
                explain!("Record #{} dropped because it had no relatives in the assignment");
            }
            NoCrisprInfo => {
                self.dropped.num_no_crispr_info += 1;
                explain!("Record #{} dropped because it had no relatives in the assignment");
            }
        }
    }

    pub fn account_classes<'r>(&mut self, record: &EnrichedVpfClassRecord<'a, 'r, CE>) {
        let class_name = &*record.vpf_class_record.class_name;

        self
            .class_data
            .raw_entry_mut()
            .from_key(class_name)
            .or_insert_with(|| (class_name.to_owned(), ClassData::default()))
            .1
            .add(record);
    }

    pub fn write_classes<W: io::Write>(&self, writer: W, sort: SummarySortBy) -> csv::Result<()> {
        let mut class_stats = self
            .class_data
            .iter()
            .map(|(class_name, data)| ClassStats::<CE> {
                class_name: class_name.clone(),
                virus_count: data.virus_count as i32,
                num_assigned_taxids: data.assigned_taxids.len() as i32,
                num_assigned_contigs: data.assigned_contigs.len() as i32,
                crispr_enrichment_stats: data.crispr_enrichment_data.clone().into(),
            })
            .collect_vec();

        match sort {
            SummarySortBy::ClassName => class_stats.sort_by(|r1, r2| r1.class_name.cmp(&r2.class_name)),
            SummarySortBy::VirusCount => class_stats.sort_by_key(|r| -r.virus_count),
            SummarySortBy::AssignedTaxIds => class_stats.sort_by_key(|r| -r.num_assigned_taxids),
            SummarySortBy::AssignedContigs => class_stats.sort_by_key(|r| -r.num_assigned_contigs),
        }

        let mut csv_writer = csv::WriterBuilder::new()
            .has_headers(true)
            .delimiter(b'\t')
            .from_writer(writer);

        for row in class_stats {
            csv_writer.serialize(&row)?;
        }

        csv_writer.flush()?;
        Ok(())
    }
}

fn write_basic_drop_stats<W: io::Write>(mut w: W, dropped: &RecordDropStats) -> io::Result<()> {
    writeln!(
        &mut w,
        "{} records dropped due to --filter",
        dropped.num_filtered
    )?;

    writeln!(
        &mut w,
        "{} records dropped because the host prediction was not found in the reference taxonomy",
        dropped.num_unknown_host_taxid
    )?;

    writeln!(
        &mut w,
        "{} records unrelated to the metagenomic assignment dropped",
        dropped.num_not_in_assignment
    )?;

    Ok(())
}

impl<'a> EnrichmentSummary<'a, NoEnrichment> {
    pub fn write_drop_stats<W: io::Write>(&self, mut w: W) -> io::Result<()> {
        write_basic_drop_stats(&mut w, &self.dropped)?;

        writeln!(
            &mut w,
            "{} records kept",
            self.num_records - self.dropped.total_dropped()
        )?;

        Ok(())
    }
}

impl<'a> EnrichmentSummary<'a, CrisprEnrichment<'a>> {
    pub fn write_drop_stats<W: io::Write>(&self, mut w: W) -> io::Result<()> {
        write_basic_drop_stats(&mut w, &self.dropped)?;

        writeln!(
            &mut w,
            "{} records with no CRISPR matches dropped",
            self.dropped.num_no_crispr_info
        )?;

        writeln!(
            &mut w,
            "{} records kept",
            self.num_records - self.dropped.total_dropped()
        )?;

        Ok(())
    }
}

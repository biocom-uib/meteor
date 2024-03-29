use std::{io::BufRead, iter};

use anyhow::Context;
use clap::Args;
use itertools::Itertools;

use crate::{
    preprocessed_taxonomy::{with_some_taxonomy, PreprocessedTaxonomyArgs},
    taxonomy::LabelledTaxonomy, util,
};

/// Obtain lineage information from a list of taxids. The input is provided as a taxid per line on
/// STDIN.
#[derive(Args)]
pub struct GetLineageArgs {
    #[clap(flatten)]
    taxonomy: PreprocessedTaxonomyArgs,
}

pub fn get_lineage_with_taxonomy(taxo: &impl LabelledTaxonomy) -> anyhow::Result<()> {
    let mut output = csv::WriterBuilder::new()
        .delimiter(b'\t')
        .from_writer(std::io::stdout());

    output.write_record(["taxid", "lineage_taxids", "lineage_names", "lineage_ranks"])?;

    let stdin_lines = std::io::stdin().lock().lines();

    itertools::process_results(stdin_lines, |lines| -> anyhow::Result<()> {
        for line in lines {
            let line = line.trim();

            if line.is_empty() {
                continue;
            }

            let taxid = line
                .parse()
                .with_context(|| format!("Unable to parse input taxid: {line:?}"))?;

            let taxid = if let Some(taxid) = taxo.fixup_node(taxid) {
                taxid
            } else {
                eprintln!("Warning: Unknown input taxid {taxid}, skipping");
                continue;
            };

            let root = taxo.get_root();

            let ancestors = iter::once(taxid)
                .chain(taxo.strict_ancestors(taxid).filter(|ancestor| *ancestor != root))
                .collect_vec()
                .into_iter()
                .rev();

            let (lineage_taxids, lineage_labels, lineage_ranks): (Vec<_>, Vec<_>, Vec<_>) =
                ancestors
                    .map(|ancestor_id| {
                        let label = taxo.some_label_of(ancestor_id).unwrap_or("");
                        let rank = taxo.find_rank_str(ancestor_id).unwrap_or("");
                        (ancestor_id.to_string(), label, rank)
                    })
                    .multiunzip();

            let lineage_taxids = lineage_taxids.into_iter().join(";");
            let lineage_names = lineage_labels.into_iter().join(";");
            let lineage_ranks = lineage_ranks.into_iter().join(";");

            output.write_record([
                &taxid.to_string(),
                &lineage_taxids,
                &lineage_names,
                &lineage_ranks,
            ])?;
        }

        Ok(())
    })?
}

pub fn get_lineage(args: GetLineageArgs) -> anyhow::Result<()> {
    let now = std::time::Instant::now();

    let taxonomy = args.taxonomy.deserialize()?;

    eprintln!("Loaded taxonomy in {} seconds", now.elapsed().as_secs());

    with_some_taxonomy!(&taxonomy.tree, tax => {
        util::ignore_broken_pipe_anyhow(get_lineage_with_taxonomy(tax))
    })
}

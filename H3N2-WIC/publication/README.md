# Public sequence identifiers

`gisaid_isolate_identifiers.tsv` contains one `isolate_id` column with the 2,501 unique GISAID identifiers recorded for the GISAID-sourced rows of the finalized WIC sequence-source table. It includes the four records whose accession-database annotation remains IRD but whose selected sequence payload came from GISAID. This is the finalized source inventory before passage filtering, not a claim that every identifier enters every model variant.

The exported identifiers come from `gisaid_selected_isolate_id` in the local `H3N2-WIC/data/final/H3N2_WIC_seq_data_final.tsv`, restricted to rows where `sequence_source` is `GISAID`, then deduplicated. The public column is named `isolate_id`; it is not copied from the internal column with that name.

This file is an identifier list for locating records through GISAID EpiFlu. It contains no sequences, raw FASTA headers, collection dates, locations, host fields, laboratory fields, or internal selection diagnostics. The reported GISAID searches were manual, as confirmed by the author on 19 September 2026.

The complete internal mapping, canonical-isolate table, filled manifests, final sequence-source table, and download-helper CSV remain local and are ignored by Git. They have not been replaced with reduced-schema files at their pipeline-input paths. The manuscript repository therefore does not provide an unrestricted copy of every input needed to rerun the analysis. Authorized users must obtain restricted records through GISAID under its [access terms](https://gisaid.org/terms-of-use/).

The [GISAID publication guide](https://gisaid.org/publish/) permits accession identifiers and describes contributor acknowledgment and EPI_SET records. This identifier list is not an official GISAID acknowledgment table or an EPI_SET record. A study-specific acknowledgment artifact has not yet been supplied for this package; it remains a submission item. Any official acknowledgment table supplied later must be retained in its original form.

These exclusions apply to the current repository tree. Previously committed metadata remains in Git history; no history rewrite was performed.

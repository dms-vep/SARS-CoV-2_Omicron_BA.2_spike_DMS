"""Top-level ``snakemake`` file that runs analysis."""


import os


configfile: "config.yaml"


# include `dms-vep-pipeline` pipeline Snakemake file
include: os.path.join(config["pipeline_path"], "pipeline.smk")


rule all:
    input:
        variant_count_files,
        rules.check_adequate_variant_counts.output.passed,
        # antibody escape output
        prob_escape_files,
        expand(
            rules.fit_polyclonal.output.pickle,
            antibody_selection_group=antibody_selections["selection_group"].unique(),
        ),
        # functional effects of mutations output
        func_score_files,
        expand(
            rules.fit_globalepistasis.output.muteffects_latent,
            func_selection=func_selections["selection_name"],
        ),
        expand(
            rules.fit_globalepistasis.output.muteffects_observed,
            func_selection=func_selections["selection_name"],
        ),
        (
            [config["muteffects_observed"], config["muteffects_latent"]]
            if len(func_selections)
            else []
        ),
        # HTML documentation output
        config["docs"],


# Arbitrary other rules should be added here
rule site_numbering_map:
    """Map sequential numbering of protein in experiments to standard reference."""
    input:
        prot=config["gene_sequence_protein"],
    output:
        reference="results/site_numbering/numbering_reference.fa",
        alignment="results/site_numbering/alignment.fa",
        to_align="results/site_numbering/to_align.fa",
        site_numbering_map=config["site_numbering_map"],
    params:
        numbering_reference_accession=config["numbering_reference_accession"],
    log:
        os.path.join(config["logdir"], "site_numbering_map.txt"),
    conda:
        "dms-vep-pipeline/environment.yml"
    script:
        "scripts/site_numbering_map.py"


# Add any extra data/results files for docs with name: file
extra_data_files = {
    "sequential to reference site numbering": config["site_numbering_map"],
}


# include `dms-vep-pipeline` docs building Snakemake file
include: os.path.join(config["pipeline_path"], "docs.smk")

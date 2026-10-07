# Environment notes

This workflow was developed and validated using a mixed software stack rather than a single monolithic Conda environment.

## Workflow engine

The workflow was orchestrated with:

- Snakemake 7.32.4

## External software

The analysis relies on the following command-line tools:

- fastp
- Bowtie2
- SAMtools
- MEGAHIT
- CD-HIT-DUP
- CD-HIT-EST
- NCBI BLAST+
- Docker
- geNomad

Exact software versions used for the validated analysis are listed in:

`docs/software_versions.md`

## geNomad

geNomad was executed using the pinned Docker image:

```text
community.wave.seqera.io/library/genomad:1.12.0--27836e6e665e84b5
```

The geNomad database is not distributed with this repository and must be provided locally.

## Plotting environment

The plotting scripts require Python, matplotlib, and NumPy.

A reference Conda environment is provided at:

```text
workflow/envs/plotting.yaml
```

It can be created with:

```bash
conda env create -f workflow/envs/plotting.yaml
conda activate wastewater-capture-plotting
```

Alternatively, an existing compatible Python interpreter can be supplied through the `PLOT_PYTHON` environment variable:

```bash
export PLOT_PYTHON=/path/to/python
snakemake --cores 8
```

The validated plotting environment used:

- Python 3.10.20
- matplotlib 3.10.9
- NumPy 2.2.6

## Reference databases

Large reference databases and indexes are not distributed with this repository.

The workflow requires local access to:

- a Bowtie2 GRCh38 human-reference index;
- the geNomad database;
- the NCBI `nt_viruses` BLAST database.

Database locations should be configured locally and should not be committed to the repository.

## Reproducibility

The end-to-end Snakemake DAG was validated from preprocessing through the final portfolio heatmap.

A reproducibility test was performed by deleting the final portfolio matrix and heatmap outputs and regenerating them through Snakemake.

The regenerated TSV and PNG outputs were byte-identical to the original outputs. The regenerated PDF had different binary metadata but identical rendered content.

Because several external tools are provided by system installations, existing Conda environments, or containers, the current workflow does not require every Snakemake rule to use an individual Conda environment.

Validated software versions are documented in:

`docs/software_versions.md`

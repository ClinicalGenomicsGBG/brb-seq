# ClinicalGenomicsGBG/brb-seq

## Introduction

**ClinicalGenomicsGBG/brb-seq** is a bioinformatics pipeline that preprocesses raw sequencing data from BRB-seq and computes count matrices.

## Usage

> [!NOTE]
> If you are new to Nextflow and nf-core, please refer to [this page](https://nf-co.re/docs/usage/installation) on how to set-up Nextflow.

This pipeline demultiplexes a BRB-seq sequencing run with BCL Convert, then performs read QC, alignment and count-matrix generation. Prepare a CSV samplesheet with one row for each biological sample.

`samplesheet.csv`:

```csv
sample_id,udi,barcode
SAMPLE_A,MQ-UDI-1,AAACCCAAGAAACACT
SAMPLE_B,MQ-UDI-1,AAACCCAAGAAACCAT
```

The three columns are mandatory and must be provided in this order:

- `sample_id`: biological sample name used in the count matrices and demultiplexed FASTQ filenames
- `udi`: unique dual index name matching the sample sheet used for the sequencing run
- `barcode`: BRB-seq barcode for the biological sample

The samplesheet must be a CSV file. The sequencing run directory (or `.tar.gz` archive) is supplied separately with `--rundir`; the BCL Convert sample sheet defaults to `assets/bclconvert/samplesheet_forward.csv`.

Run the pipeline with:

```bash
nextflow run ClinicalGenomicsGBG/brb-seq \
   -profile <docker/singularity/.../institute> \
   --input samplesheet.csv \
   --rundir /path/to/sequencing/run \
   --fasta reference_genome.fa \
   --gtf reference_genome.gtf \
   --outdir <OUTDIR>
```

If you have a pre-computed STAR index for your genome, supply it using `--star_index` and omit `--fasta` and `--gtf`.

> [!WARNING]
> Please provide pipeline parameters via the CLI or Nextflow `-params-file` option. Custom config files including those provided by the `-c` Nextflow option can be used to provide any configuration _**except for parameters**_; see [docs](https://nf-co.re/docs/usage/getting_started/configuration#custom-configuration-files).

## Credits

ClinicalGenomicsGBG/brb-seq was originally written by Daniel Schmitz.

## Contributions and Support

If you would like to contribute to this pipeline, please see the [contributing guidelines](.github/CONTRIBUTING.md).

## Citations

An extensive list of references for the tools used by the pipeline can be found in the [`CITATIONS.md`](CITATIONS.md) file.

This pipeline uses code and infrastructure developed and maintained by the [nf-core](https://nf-co.re) community, reused here under the [MIT license](https://github.com/nf-core/tools/blob/main/LICENSE).

> **The nf-core framework for community-curated bioinformatics pipelines.**
>
> Philip Ewels, Alexander Peltzer, Sven Fillinger, Harshil Patel, Johannes Alneberg, Andreas Wilm, Maxime Ulysse Garcia, Paolo Di Tommaso & Sven Nahnsen.
>
> _Nat Biotechnol._ 2020 Feb 13. doi: [10.1038/s41587-020-0439-x](https://dx.doi.org/10.1038/s41587-020-0439-x).

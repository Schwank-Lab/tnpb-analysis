# tnpb-analysis

Analysis code for the directed evolution of IsDra2 TnpB: PacBio long-read genotyping of
evolving phage populations, and Illumina-based pooled self-targeting screens in mammalian
cells.

## Layout

```
genotyping/               alignment, mutation calling, PacBio demultiplexing
sequencing/               FASTQ anchoring and read-part extraction
screening/                editing-rate, indel and normalisation routines
pipelines/                per-sample read-processing task

pacbio_data_analysis/     PacBio genotyping
  A8-03/                    PACE campaign 1
  A8-04/                    PACE campaign 2
cell_screening_analysis/  pooled self-targeting screens
  A9_screen/ A9_07_screen/  5N TAM libraries
  A13_screen/               canonical-TAM library
  A14_02/ A14_03/           benchmarking library (A14_03 also covers indels)
```

The four top-level packages are the shared analysis library. Each is an extract containing
only the functions the notebooks here call, plus their helpers.

## Environment

All dependencies, including the `bwa` and `fastqc` binaries, are declared in
`environment.yml`:

```bash
conda env create -f environment.yml
conda activate tnpb-analysis
```

Versions that affect reported numbers are pinned; `bwa` in particular, because alignment of
long, noisy reads is sensitive to the release.

Nothing needs installing: the scripts and notebooks add the repository root to `sys.path`
themselves. In VS Code, also set `"jupyter.notebookFileRoot": "${fileDirname}"` so the
working directory follows the notebook.

## Data

Raw FASTQ files and the bulk intermediates derived from them are not in this repository;
see the data availability statement in the accompanying publication. Each experiment's
config points at its data directory at the top of the file — edit those paths for your
machine.

The PacBio per-sample summary tables are included, so the mutation analysis in
`pacbio_data_analysis/` runs without the raw data. The alignment, demultiplexing and screen
read-processing steps do not.

The cell-screening figure notebooks in `cell_screening_analysis/Figures/` are the exception:
each reads one per-sample editing-rate table that is too large to include here, named at the
top of the notebook. Those tables come from the data availability statement as well. The
collected per-replicate tables the notebooks write are included.

## Running

Each experiment directory is self-contained. Start Jupyter with that directory as the
working directory, or run the pipeline script directly:

```bash
cd cell_screening_analysis/A14_02
python Priscilla_A14_02.py
```

The PacBio analyses read their reference and metadata from the experiment directory and
write to `out/` within it.

## License

MIT. See [LICENSE](LICENSE).

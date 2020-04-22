# Analysis of OMM sequence files

## Start pipeline

```bash
snakemake --jobs 40
```

This will create:
- per sample *fastp* reports written to docs/reports/fastp
- per sample indexed alignments written to processed/$sample_id


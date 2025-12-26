---
name: mirna-seq
description: "miRNA sequencing analysis pipeline. Process raw FASTQ.gz files, trim adapters, align to miRBase, quantify expression, perform differential expression analysis, and visualize results."
---

# miRNA-seq Analysis

## Overview

Complete miRNA sequencing analysis pipeline for processing raw FASTQ.gz files through to differential expression results. Perform quality control, adapter trimming, alignment to miRBase, expression quantification, and differential expression analysis with publication-quality visualizations.

## When to Use This Skill

This skill should be used when:
- Processing raw miRNA-seq data in FASTQ or FASTQ.gz format
- Analyzing small RNA sequencing experiments
- Quantifying mature miRNA expression levels
- Performing differential expression analysis of miRNAs
- Comparing miRNA profiles between conditions (e.g., disease vs healthy)
- Users mention "miRNA-seq", "microRNA sequencing", "small RNA-seq", or "miRNA expression"

## Quick Start Workflow

For users who want to perform a standard miRNA-seq analysis:

```python
from mirna_analysis import MiRNASeqPipeline

# Initialize pipeline
pipeline = MiRNASeqPipeline(
    output_dir="mirna_results",
    species="human",  # or "mouse", "rat"
    adapter="TGGAATTCTCGGGTGCCAAGG"  # Illumina small RNA adapter
)

# Run complete pipeline
results = pipeline.run(
    fastq_files=["sample1.fastq.gz", "sample2.fastq.gz"],
    sample_names=["control_1", "treated_1"],
    metadata="metadata.csv"
)

# Access results
print(f"Detected {len(results.mirnas)} miRNAs")
print(results.count_matrix.head())
```

## Pipeline Stages

### Stage 1: Quality Control

Assess raw sequencing data quality before processing.

```python
from mirna_analysis import quality_control

# Run FastQC-like quality assessment
qc_report = quality_control.analyze_fastq(
    "sample.fastq.gz",
    output_dir="qc_reports"
)

print(f"Total reads: {qc_report.total_reads:,}")
print(f"Mean quality: {qc_report.mean_quality:.1f}")
print(f"Read length range: {qc_report.min_length}-{qc_report.max_length} bp")
print(f"GC content: {qc_report.gc_content:.1f}%")
```

**Quality metrics:**
- Per-base sequence quality (Phred scores)
- Per-sequence quality distribution
- Read length distribution (miRNAs are typically 18-25 bp)
- GC content analysis
- Adapter contamination detection
- Overrepresented sequences

### Stage 2: Adapter Trimming

Remove sequencing adapters from reads. Critical for miRNA-seq due to short insert sizes.

```python
from mirna_analysis import adapter_trimming

# Trim adapters with common presets
trimmed_fastq = adapter_trimming.trim_adapters(
    input_fastq="sample.fastq.gz",
    output_fastq="sample_trimmed.fastq.gz",
    adapter="illumina_smallrna",  # Preset or custom sequence
    min_length=18,                # Minimum read length after trimming
    max_length=30,                # Maximum read length
    quality_cutoff=20,            # Trim low-quality bases
    trim_n=True                   # Remove reads with N bases
)

print(f"Reads before trimming: {trimmed_fastq.reads_input:,}")
print(f"Reads after trimming: {trimmed_fastq.reads_output:,}")
print(f"Adapter detected in: {trimmed_fastq.adapter_rate:.1f}% of reads")
```

**Common adapter sequences:**
| Platform | Adapter Name | Sequence |
|----------|--------------|----------|
| Illumina TruSeq Small RNA | 3' adapter | `TGGAATTCTCGGGTGCCAAGG` |
| Illumina NEBNext | 3' adapter | `AGATCGGAAGAGCACACGTCT` |
| Qiagen QIAseq | 3' adapter | `AACTGTAGGCACCATCAAT` |
| Lexogen Small RNA | 3' adapter | `TGGAATTCTCGGGTGCCAAGGAACTCCAGTCAC` |

### Stage 3: Alignment to miRBase

Map trimmed reads to mature miRNA sequences from miRBase.

```python
from mirna_analysis import alignment

# Build or load miRBase index
mirbase_index = alignment.prepare_mirbase_index(
    species="hsa",           # Human (hsa), Mouse (mmu), Rat (rno)
    mirbase_version="22.1",  # miRBase version
    include_hairpins=True    # Also index precursor sequences
)

# Align reads
alignment_results = alignment.align_to_mirbase(
    fastq="sample_trimmed.fastq.gz",
    index=mirbase_index,
    output_bam="sample_aligned.bam",
    allow_mismatches=1,      # Allow 1 mismatch for isomiRs
    seed_length=18,          # Minimum seed length
    report_all=False         # Report best alignment only
)

print(f"Mapped reads: {alignment_results.mapped_reads:,}")
print(f"Mapping rate: {alignment_results.mapping_rate:.1f}%")
print(f"Unique miRNAs detected: {alignment_results.unique_mirnas}")
```

**Alignment considerations:**
- miRNAs are short (18-25 bp), requiring sensitive alignment parameters
- Allow 1-2 mismatches to capture isomiRs (miRNA variants)
- Consider both 5p and 3p arms of precursor hairpins
- Multi-mapping reads are common due to miRNA families

### Stage 4: Expression Quantification

Count reads mapping to each miRNA to generate expression matrix.

```python
from mirna_analysis import quantification

# Quantify miRNA expression
count_matrix = quantification.count_mirnas(
    bam_files=["sample1.bam", "sample2.bam", "sample3.bam"],
    sample_names=["ctrl_1", "ctrl_2", "treat_1"],
    count_mode="unique",     # Count only uniquely mapped reads
    min_count=5,             # Filter miRNAs with < 5 total counts
    normalize=True           # Return both raw and normalized counts
)

# Access results
print(f"miRNAs quantified: {count_matrix.shape[0]}")
print(f"Samples: {count_matrix.shape[1]}")

# Save count matrix
count_matrix.to_csv("mirna_counts.csv")

# Get normalized counts (RPM/CPM)
rpm_matrix = quantification.normalize_counts(count_matrix, method="rpm")
```

**Normalization methods:**
- **RPM/CPM**: Reads/Counts per million mapped reads
- **TMM**: Trimmed Mean of M-values (recommended for DE analysis)
- **RLE**: Relative Log Expression (used by DESeq2)
- **Upper Quartile**: Robust to highly expressed outliers

### Stage 5: Differential Expression Analysis

Identify differentially expressed miRNAs between conditions.

```python
from mirna_analysis import differential_expression

# Load count matrix and metadata
import pandas as pd

counts = pd.read_csv("mirna_counts.csv", index_col=0)
metadata = pd.read_csv("metadata.csv", index_col=0)

# Run differential expression analysis
de_results = differential_expression.run_deseq2(
    counts=counts,
    metadata=metadata,
    design="~condition",
    contrast=["condition", "treated", "control"],
    alpha=0.05
)

# Get significant miRNAs
significant = de_results[de_results.padj < 0.05]
print(f"Significant miRNAs: {len(significant)}")
print(f"Upregulated: {len(significant[significant.log2FoldChange > 0])}")
print(f"Downregulated: {len(significant[significant.log2FoldChange < 0])}")

# Save results
de_results.to_csv("differential_expression_results.csv")
```

### Stage 6: Visualization

Generate publication-quality plots for miRNA-seq results.

```python
from mirna_analysis import visualization
import matplotlib.pyplot as plt

# Volcano plot
fig, ax = visualization.volcano_plot(
    de_results,
    padj_threshold=0.05,
    lfc_threshold=1.0,
    label_top_n=10,
    title="miRNA Differential Expression"
)
plt.savefig("volcano_plot.png", dpi=300)

# Expression heatmap
fig, ax = visualization.expression_heatmap(
    counts=rpm_matrix,
    mirnas=significant.index[:50],  # Top 50 significant miRNAs
    metadata=metadata,
    cluster_samples=True,
    cluster_mirnas=True,
    scale="row"  # Z-score normalization per miRNA
)
plt.savefig("expression_heatmap.png", dpi=300)

# Read length distribution
fig, ax = visualization.length_distribution(
    fastq="sample_trimmed.fastq.gz",
    title="Read Length Distribution"
)
plt.savefig("length_distribution.png", dpi=300)

# miRNA family analysis
fig, ax = visualization.family_barplot(
    de_results,
    top_families=20,
    show_members=True
)
plt.savefig("mirna_families.png", dpi=300)
```

## Complete Pipeline Example

### Single-Command Analysis

```python
from mirna_analysis import MiRNASeqPipeline

# Define samples
samples = {
    "control_1": "data/ctrl_1.fastq.gz",
    "control_2": "data/ctrl_2.fastq.gz",
    "control_3": "data/ctrl_3.fastq.gz",
    "treated_1": "data/treat_1.fastq.gz",
    "treated_2": "data/treat_2.fastq.gz",
    "treated_3": "data/treat_3.fastq.gz",
}

# Create metadata
metadata = pd.DataFrame({
    "condition": ["control", "control", "control", "treated", "treated", "treated"]
}, index=samples.keys())

# Initialize and run pipeline
pipeline = MiRNASeqPipeline(
    output_dir="mirna_analysis_results",
    species="human",
    adapter="TGGAATTCTCGGGTGCCAAGG",
    threads=4
)

results = pipeline.run(
    fastq_files=list(samples.values()),
    sample_names=list(samples.keys()),
    metadata=metadata,
    design="~condition",
    contrast=["condition", "treated", "control"]
)

# Generate report
pipeline.generate_report(results, format="html")
```

### Using the Command-Line Script

```bash
# Basic usage
python scripts/run_mirna_analysis.py \
  --fastq data/*.fastq.gz \
  --metadata metadata.csv \
  --species human \
  --adapter illumina_smallrna \
  --output results/

# With all options
python scripts/run_mirna_analysis.py \
  --fastq data/*.fastq.gz \
  --metadata metadata.csv \
  --species human \
  --adapter TGGAATTCTCGGGTGCCAAGG \
  --design "~condition" \
  --contrast condition treated control \
  --output results/ \
  --threads 8 \
  --min-length 18 \
  --max-length 30 \
  --min-count 10 \
  --plots
```

## Working with Raw FASTQ.gz Files

### Reading Compressed FASTQ Files

```python
import gzip
from Bio import SeqIO

# Read FASTQ.gz file
def read_fastq_gz(filepath):
    """Read sequences from gzipped FASTQ file."""
    with gzip.open(filepath, "rt") as handle:
        for record in SeqIO.parse(handle, "fastq"):
            yield record

# Process reads
for record in read_fastq_gz("sample.fastq.gz"):
    print(f"ID: {record.id}")
    print(f"Sequence: {record.seq}")
    print(f"Quality: {record.letter_annotations['phred_quality']}")
    break
```

### Using Pysam for FASTQ Processing

```python
import pysam

# Read FASTQ with pysam (handles gzip automatically)
with pysam.FastxFile("sample.fastq.gz") as fastq:
    for entry in fastq:
        name = entry.name
        sequence = entry.sequence
        quality = entry.quality

        # Filter by length (typical miRNA size)
        if 18 <= len(sequence) <= 30:
            print(f"{name}: {sequence}")
```

### Streaming Large Files

```python
from mirna_analysis import utils

# Process large FASTQ.gz files in chunks
for chunk in utils.stream_fastq("large_file.fastq.gz", chunk_size=100000):
    # Process each chunk of 100,000 reads
    for read in chunk:
        # Your processing logic
        pass
```

## IsomiR Detection

Identify miRNA variants (isomiRs) with 5' or 3' modifications.

```python
from mirna_analysis import isomirs

# Detect isomiRs
isomir_results = isomirs.detect_isomirs(
    bam="sample_aligned.bam",
    mirbase_mature="mature.fa",
    min_reads=10
)

# Summarize isomiR types
for mirna, variants in isomir_results.items():
    print(f"\n{mirna}:")
    for variant in variants:
        print(f"  {variant.type}: {variant.sequence} ({variant.count} reads)")
```

**IsomiR types detected:**
- **5' trimming/addition**: Variants at the 5' end
- **3' trimming/addition**: Variants at the 3' end
- **3' non-templated additions**: A/U tailing
- **Internal modifications**: SNPs or editing events

## Target Prediction Integration

Predict miRNA target genes using external databases.

```python
from mirna_analysis import targets

# Get validated targets from miRTarBase
validated_targets = targets.get_mirtarbase_targets(
    mirnas=["hsa-miR-21-5p", "hsa-miR-155-5p"],
    evidence="strong"  # Strong, weak, or all
)

# Predict targets using TargetScan
predicted_targets = targets.predict_targetscan(
    mirnas=["hsa-miR-21-5p"],
    species="human",
    context_score_threshold=-0.2
)

# Perform pathway enrichment on targets
enrichment = targets.pathway_enrichment(
    genes=validated_targets["hsa-miR-21-5p"],
    database="KEGG"
)
```

## Integration with Other Skills

### With PyDESeq2

```python
from pydeseq2.dds import DeseqDataSet
from pydeseq2.ds import DeseqStats
import pandas as pd

# Load miRNA count matrix
counts = pd.read_csv("mirna_counts.csv", index_col=0).T
metadata = pd.read_csv("metadata.csv", index_col=0)

# Run DESeq2 analysis
dds = DeseqDataSet(
    counts=counts,
    metadata=metadata,
    design="~condition"
)
dds.deseq2()

ds = DeseqStats(dds, contrast=["condition", "treated", "control"])
ds.summary()

results = ds.results_df
significant_mirnas = results[results.padj < 0.05]
```

### With Scanpy (Single-Cell miRNA)

```python
import scanpy as sc
import anndata as ad

# Create AnnData object from miRNA counts
adata = ad.AnnData(X=count_matrix.values.T)
adata.obs_names = count_matrix.columns
adata.var_names = count_matrix.index

# Standard scanpy workflow
sc.pp.filter_genes(adata, min_cells=3)
sc.pp.normalize_total(adata, target_sum=1e6)
sc.pp.log1p(adata)
sc.pp.highly_variable_genes(adata)
sc.tl.pca(adata)
sc.pp.neighbors(adata)
sc.tl.umap(adata)
sc.pl.umap(adata, color="condition")
```

## Troubleshooting

### Low Mapping Rate

**Problem:** Less than 50% of reads map to miRBase

**Possible causes:**
1. Incorrect adapter sequence
2. Low quality data
3. High rRNA/tRNA contamination
4. Wrong species selected

**Solutions:**
```python
# Check adapter content
from mirna_analysis import quality_control
adapter_check = quality_control.detect_adapters("sample.fastq.gz")
print(f"Detected adapter: {adapter_check.best_match}")

# Check RNA composition
composition = quality_control.rna_composition("sample_trimmed.fastq.gz")
print(f"miRNA: {composition.mirna_pct:.1f}%")
print(f"rRNA: {composition.rrna_pct:.1f}%")
print(f"tRNA: {composition.trna_pct:.1f}%")
```

### No Significant miRNAs

**Problem:** Differential expression finds no significant miRNAs

**Diagnostics:**
```python
# Check sample clustering
from mirna_analysis import visualization
visualization.pca_plot(count_matrix, metadata, color_by="condition")

# Check dispersion
visualization.dispersion_plot(de_results)

# Check for batch effects
visualization.batch_effect_plot(count_matrix, metadata, batch_col="batch")
```

### Memory Issues with Large Files

**Solution:** Use streaming mode
```python
pipeline = MiRNASeqPipeline(
    output_dir="results",
    streaming=True,      # Enable streaming mode
    chunk_size=500000    # Process 500K reads at a time
)
```

## Reference Documentation

For comprehensive details beyond this workflow-oriented guide:

- **API Reference** (`references/api_reference.md`): Complete documentation of miRNA-seq classes, methods, and parameters

- **Workflow Guide** (`references/workflow_guide.md`): In-depth guide covering analysis workflows, quality control, adapter trimming, alignment strategies, and best practices

- **Database Reference** (`references/database_reference.md`): Information about miRBase, isomiR databases, and target prediction resources

## Installation and Requirements

```bash
# Core dependencies
uv pip install biopython pysam pandas numpy matplotlib seaborn

# For differential expression
uv pip install pydeseq2

# For target prediction (optional)
uv pip install requests  # For API access to miRTarBase/TargetScan
```

**System requirements:**
- Python 3.9+
- 8GB RAM minimum (16GB recommended for large datasets)
- 50GB disk space for miRBase indices

## Key Reminders

1. **Adapter trimming is critical:** miRNA inserts are shorter than read length, so nearly all reads contain adapter sequence.

2. **Quality filtering matters:** Use Phred quality cutoff of 20+ to remove low-quality reads.

3. **Length filtering:** Mature miRNAs are 18-25 bp. Filter reads outside 18-30 bp range after trimming.

4. **Multi-mapping:** Many miRNAs belong to families with similar sequences. Decide how to handle multi-mappers before quantification.

5. **Normalization for DE:** Use DESeq2's RLE normalization or TMM, not simple RPM, for differential expression.

6. **Biological replicates:** Aim for at least 3 biological replicates per condition for reliable DE analysis.

7. **IsomiRs matter:** Consider isomiR variants in your analysis—they can have different biological functions.

8. **Validate findings:** Use qRT-PCR to validate key miRNA findings from sequencing.

## Additional Resources

- **miRBase:** https://www.mirbase.org/
- **miRTarBase:** https://mirtarbase.cuhk.edu.cn/
- **TargetScan:** https://www.targetscan.org/
- **isomiRage:** https://cru.genomics.iit.it/isomirage/
- **sRNAtoolbox:** https://bioinfo5.ugr.es/srnatoolbox/

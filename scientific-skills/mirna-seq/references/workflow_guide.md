# miRNA-seq Workflow Guide

This guide provides comprehensive workflows for miRNA sequencing analysis, from raw data to biological interpretation.

## Table of Contents

1. [Complete Analysis Workflow](#complete-analysis-workflow)
2. [Quality Control](#quality-control)
3. [Adapter Trimming](#adapter-trimming)
4. [Alignment Strategies](#alignment-strategies)
5. [Quantification Methods](#quantification-methods)
6. [Differential Expression](#differential-expression)
7. [Advanced Topics](#advanced-topics)
8. [Best Practices](#best-practices)

---

## Complete Analysis Workflow

### Standard Pipeline Overview

```
Raw FASTQ.gz
    │
    ▼
┌─────────────────┐
│ Quality Control │  → QC reports, read statistics
└────────┬────────┘
         │
         ▼
┌─────────────────┐
│ Adapter Trimming│  → Trimmed FASTQ, adapter stats
└────────┬────────┘
         │
         ▼
┌─────────────────┐
│ Length Filtering│  → Size-selected reads (18-30 bp)
└────────┬────────┘
         │
         ▼
┌─────────────────┐
│ Alignment       │  → BAM files, mapping stats
└────────┬────────┘
         │
         ▼
┌─────────────────┐
│ Quantification  │  → Count matrix
└────────┬────────┘
         │
         ▼
┌─────────────────┐
│ Normalization   │  → Normalized expression values
└────────┬────────┘
         │
         ▼
┌─────────────────┐
│ Differential    │  → DE results, significant miRNAs
│ Expression      │
└────────┬────────┘
         │
         ▼
┌─────────────────┐
│ Visualization & │  → Plots, reports, target predictions
│ Interpretation  │
└─────────────────┘
```

### Implementation

```python
import os
import pandas as pd
import gzip
from pathlib import Path
from Bio import SeqIO
from collections import defaultdict

class MiRNASeqPipeline:
    """Complete miRNA-seq analysis pipeline."""

    def __init__(self, output_dir, species="human", adapter=None, threads=4):
        self.output_dir = Path(output_dir)
        self.output_dir.mkdir(parents=True, exist_ok=True)
        self.species = species
        self.adapter = adapter or self._get_default_adapter()
        self.threads = threads

        # Create subdirectories
        for subdir in ["qc", "trimmed", "aligned", "counts", "results", "plots"]:
            (self.output_dir / subdir).mkdir(exist_ok=True)

    def _get_default_adapter(self):
        """Return default Illumina small RNA adapter."""
        return "TGGAATTCTCGGGTGCCAAGG"

    def run(self, fastq_files, sample_names, metadata, design="~condition",
            contrast=None):
        """Execute complete pipeline."""
        results = {}

        # Stage 1: Quality Control
        print("Stage 1: Quality Control")
        qc_results = self._run_qc(fastq_files, sample_names)
        results["qc"] = qc_results

        # Stage 2: Adapter Trimming
        print("\nStage 2: Adapter Trimming")
        trimmed_files = self._run_trimming(fastq_files, sample_names)
        results["trimming"] = trimmed_files

        # Stage 3: Alignment
        print("\nStage 3: Alignment to miRBase")
        aligned_files = self._run_alignment(trimmed_files, sample_names)
        results["alignment"] = aligned_files

        # Stage 4: Quantification
        print("\nStage 4: Expression Quantification")
        count_matrix = self._run_quantification(aligned_files, sample_names)
        results["counts"] = count_matrix

        # Stage 5: Differential Expression
        if contrast:
            print("\nStage 5: Differential Expression Analysis")
            de_results = self._run_de_analysis(count_matrix, metadata, design, contrast)
            results["de"] = de_results

        return results

    def _run_qc(self, fastq_files, sample_names):
        """Run quality control on input files."""
        qc_results = {}

        for fastq, name in zip(fastq_files, sample_names):
            print(f"  Processing {name}...")
            stats = self._calculate_fastq_stats(fastq)
            qc_results[name] = stats

            # Save individual QC report
            report_path = self.output_dir / "qc" / f"{name}_qc.txt"
            self._write_qc_report(stats, report_path)

        return qc_results

    def _calculate_fastq_stats(self, fastq_path):
        """Calculate quality statistics for FASTQ file."""
        stats = {
            "total_reads": 0,
            "total_bases": 0,
            "mean_length": 0,
            "min_length": float("inf"),
            "max_length": 0,
            "mean_quality": 0,
            "gc_content": 0,
            "length_distribution": defaultdict(int)
        }

        quality_sum = 0
        gc_count = 0

        opener = gzip.open if str(fastq_path).endswith(".gz") else open
        mode = "rt" if str(fastq_path).endswith(".gz") else "r"

        with opener(fastq_path, mode) as handle:
            for record in SeqIO.parse(handle, "fastq"):
                seq = str(record.seq)
                qual = record.letter_annotations["phred_quality"]

                stats["total_reads"] += 1
                stats["total_bases"] += len(seq)
                stats["min_length"] = min(stats["min_length"], len(seq))
                stats["max_length"] = max(stats["max_length"], len(seq))
                stats["length_distribution"][len(seq)] += 1

                quality_sum += sum(qual)
                gc_count += seq.count("G") + seq.count("C")

        if stats["total_reads"] > 0:
            stats["mean_length"] = stats["total_bases"] / stats["total_reads"]
            stats["mean_quality"] = quality_sum / stats["total_bases"]
            stats["gc_content"] = (gc_count / stats["total_bases"]) * 100

        return stats

    def _write_qc_report(self, stats, output_path):
        """Write QC statistics to file."""
        with open(output_path, "w") as f:
            f.write("miRNA-seq Quality Control Report\n")
            f.write("=" * 40 + "\n\n")
            f.write(f"Total reads: {stats['total_reads']:,}\n")
            f.write(f"Total bases: {stats['total_bases']:,}\n")
            f.write(f"Mean read length: {stats['mean_length']:.1f} bp\n")
            f.write(f"Read length range: {stats['min_length']}-{stats['max_length']} bp\n")
            f.write(f"Mean quality score: {stats['mean_quality']:.1f}\n")
            f.write(f"GC content: {stats['gc_content']:.1f}%\n")

    def _run_trimming(self, fastq_files, sample_names):
        """Trim adapters from reads."""
        trimmed_files = {}

        for fastq, name in zip(fastq_files, sample_names):
            print(f"  Trimming {name}...")
            output_path = self.output_dir / "trimmed" / f"{name}_trimmed.fastq.gz"
            stats = self._trim_adapters(fastq, output_path)
            trimmed_files[name] = {
                "path": output_path,
                "stats": stats
            }

        return trimmed_files

    def _trim_adapters(self, input_path, output_path, min_length=18, max_length=30):
        """Trim adapter sequences from reads."""
        stats = {
            "reads_input": 0,
            "reads_output": 0,
            "reads_too_short": 0,
            "reads_too_long": 0,
            "adapter_found": 0
        }

        opener = gzip.open if str(input_path).endswith(".gz") else open
        mode = "rt" if str(input_path).endswith(".gz") else "r"

        with opener(input_path, mode) as infile, \
             gzip.open(output_path, "wt") as outfile:

            for record in SeqIO.parse(infile, "fastq"):
                stats["reads_input"] += 1
                seq = str(record.seq)

                # Find adapter position
                adapter_pos = seq.find(self.adapter[:8])  # Use first 8 bp of adapter

                if adapter_pos != -1:
                    stats["adapter_found"] += 1
                    seq = seq[:adapter_pos]
                    record = record[:adapter_pos]

                # Length filtering
                if len(seq) < min_length:
                    stats["reads_too_short"] += 1
                    continue
                elif len(seq) > max_length:
                    stats["reads_too_long"] += 1
                    continue

                stats["reads_output"] += 1
                SeqIO.write(record, outfile, "fastq")

        return stats

    def _run_alignment(self, trimmed_files, sample_names):
        """Align reads to miRBase."""
        # In practice, this would use bowtie/bowtie2 or similar
        # This is a simplified implementation
        aligned_files = {}

        for name in sample_names:
            print(f"  Aligning {name}...")
            input_path = trimmed_files[name]["path"]
            output_path = self.output_dir / "aligned" / f"{name}.bam"

            # Placeholder for actual alignment
            # In real implementation, would call bowtie2 or similar
            aligned_files[name] = {
                "path": output_path,
                "stats": {"mapped_reads": 0, "mapping_rate": 0}
            }

        return aligned_files

    def _run_quantification(self, aligned_files, sample_names):
        """Quantify miRNA expression."""
        # In practice, this would parse BAM files and count reads per miRNA
        # This is a simplified placeholder
        count_matrix = pd.DataFrame()
        return count_matrix

    def _run_de_analysis(self, count_matrix, metadata, design, contrast):
        """Run differential expression analysis."""
        from pydeseq2.dds import DeseqDataSet
        from pydeseq2.ds import DeseqStats

        # Ensure proper format
        counts_df = count_matrix.T  # Samples × miRNAs

        dds = DeseqDataSet(
            counts=counts_df,
            metadata=metadata,
            design=design
        )
        dds.deseq2()

        ds = DeseqStats(dds, contrast=contrast)
        ds.summary()

        return ds.results_df
```

---

## Quality Control

### Pre-Alignment QC Metrics

Key metrics to assess before processing:

| Metric | Good | Warning | Poor |
|--------|------|---------|------|
| Total reads | >5M | 1-5M | <1M |
| Mean quality | >30 | 25-30 | <25 |
| Adapter content | >80% | 50-80% | <50% |
| Read length mode | 21-23 bp | 18-25 bp | Outside range |
| GC content | 45-55% | 40-60% | Outside range |

### Quality Score Analysis

```python
import numpy as np
import matplotlib.pyplot as plt

def plot_quality_profile(fastq_path, output_path=None):
    """Plot per-position quality scores."""
    position_qualities = defaultdict(list)

    opener = gzip.open if str(fastq_path).endswith(".gz") else open
    mode = "rt" if str(fastq_path).endswith(".gz") else "r"

    # Sample first 100,000 reads
    with opener(fastq_path, mode) as handle:
        for i, record in enumerate(SeqIO.parse(handle, "fastq")):
            if i >= 100000:
                break
            for pos, qual in enumerate(record.letter_annotations["phred_quality"]):
                position_qualities[pos].append(qual)

    # Calculate statistics per position
    positions = sorted(position_qualities.keys())
    means = [np.mean(position_qualities[p]) for p in positions]
    q25 = [np.percentile(position_qualities[p], 25) for p in positions]
    q75 = [np.percentile(position_qualities[p], 75) for p in positions]

    # Plot
    fig, ax = plt.subplots(figsize=(12, 6))
    ax.fill_between(positions, q25, q75, alpha=0.3, color="blue", label="IQR")
    ax.plot(positions, means, color="blue", linewidth=2, label="Mean")
    ax.axhline(y=30, color="green", linestyle="--", label="Q30 threshold")
    ax.axhline(y=20, color="orange", linestyle="--", label="Q20 threshold")

    ax.set_xlabel("Position in Read (bp)")
    ax.set_ylabel("Quality Score (Phred)")
    ax.set_title("Per-Position Quality Score Distribution")
    ax.legend()
    ax.set_ylim(0, 42)

    if output_path:
        plt.savefig(output_path, dpi=300, bbox_inches="tight")

    return fig, ax
```

### Adapter Detection

```python
def detect_adapter_sequence(fastq_path, known_adapters=None):
    """Detect adapter sequence in FASTQ file."""
    if known_adapters is None:
        known_adapters = {
            "illumina_truseq": "TGGAATTCTCGGGTGCCAAGG",
            "illumina_nebnext": "AGATCGGAAGAGCACACGTCT",
            "qiagen_qiaseq": "AACTGTAGGCACCATCAAT",
        }

    adapter_counts = defaultdict(int)
    total_reads = 0

    opener = gzip.open if str(fastq_path).endswith(".gz") else open
    mode = "rt" if str(fastq_path).endswith(".gz") else "r"

    with opener(fastq_path, mode) as handle:
        for i, record in enumerate(SeqIO.parse(handle, "fastq")):
            if i >= 100000:  # Sample first 100K reads
                break
            total_reads += 1
            seq = str(record.seq)

            for name, adapter in known_adapters.items():
                if adapter[:12] in seq:  # Check first 12 bp
                    adapter_counts[name] += 1

    results = {
        name: count / total_reads * 100
        for name, count in adapter_counts.items()
    }

    best_match = max(results, key=results.get) if results else None

    return {
        "total_reads_sampled": total_reads,
        "adapter_percentages": results,
        "best_match": best_match,
        "best_match_sequence": known_adapters.get(best_match)
    }
```

---

## Adapter Trimming

### Trimming Parameters

| Parameter | Recommended | Description |
|-----------|-------------|-------------|
| min_length | 18 | Minimum length after trimming |
| max_length | 30 | Maximum length after trimming |
| quality_cutoff | 20 | Trim low-quality 3' bases |
| error_rate | 0.1 | Max error rate for adapter matching |
| min_overlap | 3 | Minimum adapter overlap |

### Cutadapt Integration

```python
import subprocess

def trim_with_cutadapt(input_fastq, output_fastq, adapter,
                       min_length=18, max_length=30, quality_cutoff=20):
    """Trim adapters using cutadapt."""
    cmd = [
        "cutadapt",
        "-a", adapter,           # 3' adapter sequence
        "-m", str(min_length),   # Minimum length
        "-M", str(max_length),   # Maximum length
        "-q", str(quality_cutoff),  # Quality cutoff
        "--trim-n",              # Trim N bases
        "-o", output_fastq,      # Output file
        input_fastq              # Input file
    ]

    result = subprocess.run(cmd, capture_output=True, text=True)

    # Parse cutadapt statistics from stderr
    stats = parse_cutadapt_output(result.stderr)

    return stats

def parse_cutadapt_output(output_text):
    """Parse cutadapt statistics."""
    stats = {}

    for line in output_text.split("\n"):
        if "Total reads processed:" in line:
            stats["total_reads"] = int(line.split(":")[1].strip().replace(",", ""))
        elif "Reads with adapters:" in line:
            stats["reads_with_adapter"] = int(line.split(":")[1].split("(")[0].strip().replace(",", ""))
        elif "Reads written" in line:
            stats["reads_written"] = int(line.split(":")[1].split("(")[0].strip().replace(",", ""))
        elif "Total basepairs processed:" in line:
            stats["bp_processed"] = int(line.split(":")[1].strip().replace(",", "").replace(" bp", ""))

    return stats
```

### Quality-Based Trimming

```python
def quality_trim(sequence, qualities, min_quality=20, window_size=4):
    """Trim low-quality bases from 3' end using sliding window."""
    if len(sequence) != len(qualities):
        raise ValueError("Sequence and quality lengths must match")

    # Sliding window from 3' end
    for i in range(len(sequence) - window_size, -1, -1):
        window_qual = sum(qualities[i:i+window_size]) / window_size
        if window_qual >= min_quality:
            return sequence[:i+window_size], qualities[:i+window_size]

    # If no window passes, return empty
    return "", []
```

---

## Alignment Strategies

### Bowtie Alignment

```python
def build_bowtie_index(fasta_path, index_prefix):
    """Build bowtie index for miRNA sequences."""
    cmd = ["bowtie-build", fasta_path, index_prefix]
    subprocess.run(cmd, check=True)

def align_with_bowtie(fastq_path, index_prefix, output_sam,
                      mismatches=1, seed_length=18):
    """Align reads using bowtie."""
    cmd = [
        "bowtie",
        "-x", index_prefix,      # Index prefix
        "-q", fastq_path,        # Input FASTQ
        "-S", output_sam,        # Output SAM
        "-v", str(mismatches),   # Mismatches allowed
        "-l", str(seed_length),  # Seed length
        "-n", "0",               # Mismatches in seed
        "--best",                # Report best alignment
        "--strata",              # Only report best stratum
        "-p", "4",               # Threads
        "-k", "1"                # Report up to k alignments
    ]

    result = subprocess.run(cmd, capture_output=True, text=True)
    return parse_bowtie_output(result.stderr)
```

### Multi-Mapping Strategy

miRNAs often have family members with similar sequences, leading to multi-mapping reads.

**Handling strategies:**

1. **Unique only**: Count only uniquely mapped reads
   - Conservative, may miss family-wide effects
   - Best for specific miRNA quantification

2. **Random assignment**: Randomly assign multi-mappers to one location
   - Simple, but adds noise
   - Acceptable for exploratory analysis

3. **Fractional assignment**: Distribute counts equally among mappings
   - Preserves total read count
   - Good for family-level analysis

4. **EM algorithm**: Iteratively estimate expression and reassign reads
   - Most accurate, computationally intensive
   - Best for publication-quality results

```python
def handle_multimappers(alignments, strategy="fractional"):
    """Handle multi-mapping reads."""
    if strategy == "unique":
        return [a for a in alignments if a.mapping_quality >= 10]

    elif strategy == "random":
        import random
        grouped = group_by_read(alignments)
        return [random.choice(group) for group in grouped.values()]

    elif strategy == "fractional":
        counts = defaultdict(float)
        grouped = group_by_read(alignments)
        for read_id, aligns in grouped.items():
            weight = 1.0 / len(aligns)
            for align in aligns:
                counts[align.reference_name] += weight
        return counts
```

---

## Quantification Methods

### Count Matrix Generation

```python
import pysam

def count_mirnas(bam_files, sample_names, min_mapq=10):
    """Generate miRNA count matrix from BAM files."""
    all_mirnas = set()
    sample_counts = {}

    for bam_path, sample_name in zip(bam_files, sample_names):
        counts = defaultdict(int)

        with pysam.AlignmentFile(bam_path, "rb") as bam:
            for read in bam:
                if read.is_unmapped or read.mapping_quality < min_mapq:
                    continue

                mirna_name = read.reference_name
                counts[mirna_name] += 1
                all_mirnas.add(mirna_name)

        sample_counts[sample_name] = counts

    # Build count matrix
    count_matrix = pd.DataFrame(
        index=sorted(all_mirnas),
        columns=sample_names
    ).fillna(0).astype(int)

    for sample_name, counts in sample_counts.items():
        for mirna, count in counts.items():
            count_matrix.loc[mirna, sample_name] = count

    return count_matrix
```

### Normalization Methods

```python
import numpy as np
from scipy import stats

def normalize_counts(count_matrix, method="rpm"):
    """Normalize count matrix."""
    if method == "rpm":
        # Reads Per Million
        library_sizes = count_matrix.sum(axis=0)
        return count_matrix / library_sizes * 1e6

    elif method == "tmm":
        # Trimmed Mean of M-values
        return tmm_normalize(count_matrix)

    elif method == "rle":
        # Relative Log Expression (DESeq2 method)
        return rle_normalize(count_matrix)

    elif method == "uq":
        # Upper Quartile
        upper_quartiles = count_matrix[count_matrix > 0].quantile(0.75)
        return count_matrix / upper_quartiles * np.median(upper_quartiles)

def rle_normalize(count_matrix):
    """DESeq2-style RLE normalization."""
    # Compute geometric mean per gene
    log_counts = np.log(count_matrix.replace(0, np.nan))
    geo_means = log_counts.mean(axis=1, skipna=True)

    # Compute size factors
    size_factors = []
    for sample in count_matrix.columns:
        sample_log = np.log(count_matrix[sample].replace(0, np.nan))
        ratios = sample_log - geo_means
        size_factor = np.exp(np.nanmedian(ratios))
        size_factors.append(size_factor)

    # Normalize
    return count_matrix / size_factors
```

---

## Differential Expression

### DESeq2 Analysis for miRNAs

```python
from pydeseq2.dds import DeseqDataSet
from pydeseq2.ds import DeseqStats

def run_mirna_deseq2(count_matrix, metadata, design="~condition",
                     contrast=None, alpha=0.05, min_count=10):
    """Run DESeq2 differential expression analysis for miRNAs."""
    # Filter low-count miRNAs
    counts_filtered = count_matrix[count_matrix.sum(axis=1) >= min_count]
    print(f"miRNAs after filtering: {len(counts_filtered)}")

    # Transpose to samples × miRNAs
    counts_df = counts_filtered.T

    # Ensure sample order matches
    common_samples = counts_df.index.intersection(metadata.index)
    counts_df = counts_df.loc[common_samples]
    metadata = metadata.loc[common_samples]

    # Initialize DESeq2
    dds = DeseqDataSet(
        counts=counts_df,
        metadata=metadata,
        design=design,
        refit_cooks=True
    )

    # Run DESeq2 pipeline
    dds.deseq2()

    # Statistical testing
    if contrast:
        ds = DeseqStats(dds, contrast=contrast, alpha=alpha)
    else:
        ds = DeseqStats(dds, alpha=alpha)

    ds.summary()

    # Apply LFC shrinkage for visualization
    ds.lfc_shrink()

    results = ds.results_df.copy()

    # Add significance classification
    results["significant"] = results["padj"] < alpha
    results["regulation"] = "NS"
    results.loc[(results["significant"]) & (results["log2FoldChange"] > 0), "regulation"] = "Up"
    results.loc[(results["significant"]) & (results["log2FoldChange"] < 0), "regulation"] = "Down"

    return results, dds
```

### Multi-Group Comparisons

```python
def pairwise_comparisons(count_matrix, metadata, condition_col="condition",
                         reference="control"):
    """Perform all pairwise comparisons against reference."""
    conditions = metadata[condition_col].unique()
    test_conditions = [c for c in conditions if c != reference]

    all_results = {}

    for test_cond in test_conditions:
        # Subset to relevant samples
        subset_samples = metadata[
            metadata[condition_col].isin([reference, test_cond])
        ].index

        subset_counts = count_matrix[subset_samples]
        subset_metadata = metadata.loc[subset_samples]

        # Run DE analysis
        results, _ = run_mirna_deseq2(
            subset_counts,
            subset_metadata,
            design=f"~{condition_col}",
            contrast=[condition_col, test_cond, reference]
        )

        all_results[f"{test_cond}_vs_{reference}"] = results

    return all_results
```

---

## Advanced Topics

### IsomiR Analysis

```python
def detect_isomirs(aligned_reads, mature_sequences, min_reads=10):
    """Detect isomiR variants from aligned reads."""
    isomirs = defaultdict(lambda: defaultdict(int))

    for read in aligned_reads:
        mirna = read.reference_name
        read_seq = read.query_sequence
        ref_seq = mature_sequences.get(mirna, "")

        if not ref_seq:
            continue

        # Classify isomiR type
        isomir_type = classify_isomir(read_seq, ref_seq, read)
        isomirs[mirna][isomir_type] += 1

    # Filter by minimum reads
    filtered = {}
    for mirna, variants in isomirs.items():
        filtered_variants = {k: v for k, v in variants.items() if v >= min_reads}
        if filtered_variants:
            filtered[mirna] = filtered_variants

    return filtered

def classify_isomir(read_seq, ref_seq, alignment):
    """Classify isomiR variant type."""
    ref_start = alignment.reference_start
    ref_end = alignment.reference_end

    # 5' variation
    if ref_start > 0:
        return f"5p-{ref_start}"
    elif alignment.query_alignment_start > 0:
        return f"5p+{alignment.query_alignment_start}"

    # 3' variation
    len_diff = len(read_seq) - len(ref_seq)
    if len_diff != 0:
        # Check for non-templated additions
        if len_diff > 0:
            tail = read_seq[-len_diff:]
            if all(b in "AU" for b in tail):
                return f"3p-NTA-{tail}"
            else:
                return f"3p+{len_diff}"
        else:
            return f"3p{len_diff}"

    return "canonical"
```

### Arm Switching Analysis

```python
def analyze_arm_switching(count_matrix, precursor_arms):
    """Analyze 5p/3p arm expression ratios."""
    arm_ratios = {}

    for precursor, arms in precursor_arms.items():
        if len(arms) != 2:
            continue

        mirna_5p = arms.get("5p")
        mirna_3p = arms.get("3p")

        if mirna_5p in count_matrix.index and mirna_3p in count_matrix.index:
            counts_5p = count_matrix.loc[mirna_5p]
            counts_3p = count_matrix.loc[mirna_3p]

            # Calculate log2 ratio (5p/3p)
            # Add pseudocount to avoid division by zero
            ratio = np.log2((counts_5p + 1) / (counts_3p + 1))

            arm_ratios[precursor] = {
                "5p_mirna": mirna_5p,
                "3p_mirna": mirna_3p,
                "5p_counts": counts_5p.values,
                "3p_counts": counts_3p.values,
                "log2_ratio": ratio.values,
                "dominant_arm": "5p" if ratio.mean() > 0 else "3p"
            }

    return pd.DataFrame(arm_ratios).T
```

---

## Best Practices

### Experimental Design

1. **Biological replicates**: Minimum 3 per condition, ideally 5-6
2. **Sequencing depth**: 5-10 million reads per sample
3. **Library prep**: Use validated small RNA library prep kits
4. **Size selection**: Ensure proper gel-based or bead-based size selection

### Data Processing

1. **Always trim adapters**: miRNA inserts are shorter than read length
2. **Quality filtering**: Use Phred ≥20 cutoff
3. **Length filtering**: 18-30 bp for mature miRNAs
4. **Document all parameters**: Ensure reproducibility

### Analysis

1. **Use appropriate normalization**: RLE or TMM for DE analysis
2. **Correct for multiple testing**: Use FDR (Benjamini-Hochberg)
3. **Consider batch effects**: Include batch in design formula if applicable
4. **Validate findings**: qRT-PCR for key miRNAs

### Reporting

1. **Report filtering criteria**: How many miRNAs/reads removed at each step
2. **Show QC metrics**: Include in supplementary materials
3. **Provide count matrices**: Enable reproducibility
4. **Deposit data**: Submit to GEO/SRA with metadata

---

## Troubleshooting

### Common Issues

| Issue | Possible Cause | Solution |
|-------|----------------|----------|
| Low mapping rate | Wrong adapter | Detect adapter sequence |
| Few miRNAs detected | Contamination | Check RNA composition |
| High duplication | Low input | Increase RNA input |
| Batch effects | Technical variation | Include batch in design |
| No DE miRNAs | Small effects | Increase sample size |

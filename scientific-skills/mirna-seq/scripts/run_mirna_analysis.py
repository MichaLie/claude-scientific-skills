#!/usr/bin/env python3
"""
miRNA-seq Analysis Pipeline

Complete pipeline for analyzing miRNA sequencing data from raw FASTQ files
through differential expression analysis.

Usage:
    python run_mirna_analysis.py --fastq data/*.fastq.gz --metadata metadata.csv \
           --species human --output results/

Requirements:
    - biopython
    - pandas
    - numpy
    - matplotlib (optional, for plots)
    - pydeseq2 (optional, for differential expression)
"""

import argparse
import gzip
import os
import sys
from collections import defaultdict
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Iterator, List, Optional, Tuple

import numpy as np
import pandas as pd

try:
    from Bio import SeqIO
    from Bio.SeqRecord import SeqRecord
except ImportError:
    print("Error: biopython not installed. Install with: pip install biopython")
    sys.exit(1)


# ============================================================================
# Data Classes
# ============================================================================

@dataclass
class QCStats:
    """Quality control statistics for a FASTQ file."""
    total_reads: int
    total_bases: int
    mean_length: float
    min_length: int
    max_length: int
    mean_quality: float
    gc_content: float
    length_distribution: Dict[int, int]


@dataclass
class TrimStats:
    """Statistics from adapter trimming."""
    reads_input: int
    reads_output: int
    reads_with_adapter: int
    reads_too_short: int
    reads_too_long: int
    adapter_rate: float


# ============================================================================
# Adapter Sequences
# ============================================================================

ADAPTER_PRESETS = {
    "illumina_smallrna": "TGGAATTCTCGGGTGCCAAGG",
    "illumina_truseq": "TGGAATTCTCGGGTGCCAAGG",
    "illumina_nebnext": "AGATCGGAAGAGCACACGTCT",
    "qiagen_qiaseq": "AACTGTAGGCACCATCAAT",
    "lexogen_smallrna": "TGGAATTCTCGGGTGCCAAGGAACTCCAGTCAC",
    "nextflex": "TGGAATTCTCGGGTGCCAAGG",
}

SPECIES_CODES = {
    "human": "hsa",
    "mouse": "mmu",
    "rat": "rno",
    "zebrafish": "dre",
    "chicken": "gga",
}


# ============================================================================
# FASTQ Processing
# ============================================================================

def read_fastq(filepath: str, limit: Optional[int] = None) -> Iterator[SeqRecord]:
    """Read FASTQ file, handling gzip compression automatically."""
    opener = gzip.open if str(filepath).endswith(".gz") else open
    mode = "rt" if str(filepath).endswith(".gz") else "r"

    with opener(filepath, mode) as handle:
        for i, record in enumerate(SeqIO.parse(handle, "fastq")):
            if limit and i >= limit:
                break
            yield record


def calculate_qc_stats(fastq_path: str, sample_limit: int = 500000) -> QCStats:
    """Calculate quality control statistics for a FASTQ file."""
    total_reads = 0
    total_bases = 0
    quality_sum = 0
    gc_count = 0
    min_length = float("inf")
    max_length = 0
    length_dist = defaultdict(int)

    for record in read_fastq(fastq_path, limit=sample_limit):
        seq = str(record.seq)
        quals = record.letter_annotations.get("phred_quality", [])

        total_reads += 1
        total_bases += len(seq)
        min_length = min(min_length, len(seq))
        max_length = max(max_length, len(seq))
        length_dist[len(seq)] += 1

        if quals:
            quality_sum += sum(quals)
        gc_count += seq.count("G") + seq.count("C")

    mean_length = total_bases / total_reads if total_reads > 0 else 0
    mean_quality = quality_sum / total_bases if total_bases > 0 else 0
    gc_content = (gc_count / total_bases * 100) if total_bases > 0 else 0

    return QCStats(
        total_reads=total_reads,
        total_bases=total_bases,
        mean_length=mean_length,
        min_length=min_length if min_length != float("inf") else 0,
        max_length=max_length,
        mean_quality=mean_quality,
        gc_content=gc_content,
        length_distribution=dict(length_dist)
    )


def detect_adapter(fastq_path: str, sample_limit: int = 50000) -> Tuple[str, float]:
    """Detect adapter sequence in FASTQ file."""
    adapter_counts = defaultdict(int)
    total_reads = 0

    for record in read_fastq(fastq_path, limit=sample_limit):
        total_reads += 1
        seq = str(record.seq)

        for name, adapter in ADAPTER_PRESETS.items():
            # Check for first 10bp of adapter
            if adapter[:10] in seq:
                adapter_counts[name] += 1

    if not adapter_counts:
        return "illumina_smallrna", 0.0

    best_adapter = max(adapter_counts, key=adapter_counts.get)
    detection_rate = adapter_counts[best_adapter] / total_reads * 100

    return best_adapter, detection_rate


def trim_adapters(
    input_path: str,
    output_path: str,
    adapter: str,
    min_length: int = 18,
    max_length: int = 30,
    quality_cutoff: int = 20
) -> TrimStats:
    """Trim adapter sequences from FASTQ file."""
    # Get adapter sequence from preset or use directly
    adapter_seq = ADAPTER_PRESETS.get(adapter, adapter)

    stats = TrimStats(
        reads_input=0,
        reads_output=0,
        reads_with_adapter=0,
        reads_too_short=0,
        reads_too_long=0,
        adapter_rate=0.0
    )

    # Open output file
    with gzip.open(output_path, "wt") as out_handle:
        for record in read_fastq(input_path):
            stats.reads_input += 1
            seq = str(record.seq)

            # Find adapter (using first 8bp as seed)
            adapter_pos = seq.find(adapter_seq[:8])

            if adapter_pos != -1:
                stats.reads_with_adapter += 1
                # Trim at adapter position
                record = record[:adapter_pos]
                seq = str(record.seq)

            # Quality trim from 3' end
            if quality_cutoff > 0 and record.letter_annotations.get("phred_quality"):
                quals = record.letter_annotations["phred_quality"]
                # Find position where quality drops below threshold
                trim_pos = len(seq)
                for i in range(len(quals) - 1, -1, -1):
                    if quals[i] >= quality_cutoff:
                        trim_pos = i + 1
                        break
                record = record[:trim_pos]
                seq = str(record.seq)

            # Length filtering
            if len(seq) < min_length:
                stats.reads_too_short += 1
                continue
            if len(seq) > max_length:
                stats.reads_too_long += 1
                continue

            # Write passing read
            stats.reads_output += 1
            SeqIO.write(record, out_handle, "fastq")

    stats.adapter_rate = (stats.reads_with_adapter / stats.reads_input * 100
                          if stats.reads_input > 0 else 0)

    return stats


# ============================================================================
# Quantification (Simplified - Sequence Counting)
# ============================================================================

def count_sequences(
    fastq_files: List[str],
    sample_names: List[str],
    min_count: int = 5
) -> pd.DataFrame:
    """
    Count unique sequences across samples.

    This is a simplified quantification method that counts exact sequence matches.
    For production use, align to miRBase reference.
    """
    all_sequences = defaultdict(lambda: defaultdict(int))

    for fastq, sample in zip(fastq_files, sample_names):
        print(f"  Counting sequences in {sample}...")
        for record in read_fastq(fastq):
            seq = str(record.seq)
            all_sequences[seq][sample] += 1

    # Build count matrix
    sequences = list(all_sequences.keys())
    count_matrix = pd.DataFrame(
        index=sequences,
        columns=sample_names,
        dtype=int
    ).fillna(0)

    for seq, counts in all_sequences.items():
        for sample, count in counts.items():
            count_matrix.loc[seq, sample] = count

    # Filter low-count sequences
    total_counts = count_matrix.sum(axis=1)
    count_matrix = count_matrix[total_counts >= min_count]

    # Add sequence length
    count_matrix["length"] = [len(seq) for seq in count_matrix.index]

    print(f"  Unique sequences: {len(count_matrix)}")

    return count_matrix


def normalize_counts(
    count_matrix: pd.DataFrame,
    method: str = "rpm"
) -> pd.DataFrame:
    """Normalize count matrix."""
    # Exclude non-count columns
    count_cols = [c for c in count_matrix.columns if c != "length"]
    counts = count_matrix[count_cols]

    if method == "rpm":
        # Reads per million
        library_sizes = counts.sum(axis=0)
        normalized = counts / library_sizes * 1e6
    elif method == "rle":
        # DESeq2-style normalization
        log_counts = np.log(counts.replace(0, np.nan))
        geo_means = log_counts.mean(axis=1, skipna=True)
        size_factors = []
        for col in counts.columns:
            sample_log = np.log(counts[col].replace(0, np.nan))
            ratios = sample_log - geo_means
            sf = np.exp(np.nanmedian(ratios))
            size_factors.append(sf)
        normalized = counts / size_factors
    else:
        normalized = counts

    # Add back length column if present
    if "length" in count_matrix.columns:
        normalized["length"] = count_matrix["length"]

    return normalized


# ============================================================================
# Differential Expression
# ============================================================================

def run_differential_expression(
    count_matrix: pd.DataFrame,
    metadata: pd.DataFrame,
    design: str,
    contrast: List[str],
    alpha: float = 0.05
) -> Optional[pd.DataFrame]:
    """Run differential expression analysis using PyDESeq2."""
    try:
        from pydeseq2.dds import DeseqDataSet
        from pydeseq2.ds import DeseqStats
    except ImportError:
        print("Warning: pydeseq2 not installed. Skipping differential expression.")
        print("Install with: pip install pydeseq2")
        return None

    # Prepare counts (exclude non-count columns)
    count_cols = [c for c in count_matrix.columns if c not in ["length"]]
    counts = count_matrix[count_cols].T  # Transpose to samples × features

    # Ensure sample order matches
    common_samples = counts.index.intersection(metadata.index)
    counts = counts.loc[common_samples]
    metadata = metadata.loc[common_samples]

    print(f"  Running DESeq2 with {len(counts)} samples, {counts.shape[1]} features")

    # Initialize and run DESeq2
    dds = DeseqDataSet(
        counts=counts.astype(int),
        metadata=metadata,
        design=design,
        refit_cooks=True
    )
    dds.deseq2()

    # Statistical testing
    ds = DeseqStats(dds, contrast=contrast, alpha=alpha)
    ds.summary()

    # Get results
    results = ds.results_df.copy()

    # Add classification
    results["significant"] = results["padj"] < alpha
    results["regulation"] = "NS"
    results.loc[(results["significant"]) & (results["log2FoldChange"] > 0), "regulation"] = "Up"
    results.loc[(results["significant"]) & (results["log2FoldChange"] < 0), "regulation"] = "Down"

    return results


# ============================================================================
# Visualization
# ============================================================================

def create_plots(
    count_matrix: pd.DataFrame,
    de_results: Optional[pd.DataFrame],
    output_dir: Path,
    metadata: pd.DataFrame = None
):
    """Create visualization plots."""
    try:
        import matplotlib.pyplot as plt
        import seaborn as sns
    except ImportError:
        print("Warning: matplotlib/seaborn not installed. Skipping plots.")
        return

    plots_dir = output_dir / "plots"
    plots_dir.mkdir(exist_ok=True)

    # 1. Read length distribution
    if "length" in count_matrix.columns:
        plt.figure(figsize=(10, 6))
        lengths = count_matrix["length"]
        count_cols = [c for c in count_matrix.columns if c != "length"]
        weights = count_matrix[count_cols].sum(axis=1)

        plt.hist(lengths, bins=range(15, 35), weights=weights, edgecolor="black", alpha=0.7)
        plt.xlabel("Sequence Length (bp)")
        plt.ylabel("Read Count")
        plt.title("Sequence Length Distribution")
        plt.axvline(x=22, color="red", linestyle="--", label="Typical miRNA (22bp)")
        plt.legend()
        plt.savefig(plots_dir / "length_distribution.png", dpi=300, bbox_inches="tight")
        plt.close()
        print(f"  Saved: {plots_dir / 'length_distribution.png'}")

    # 2. Sample correlation heatmap
    count_cols = [c for c in count_matrix.columns if c != "length"]
    if len(count_cols) >= 2:
        plt.figure(figsize=(10, 8))
        corr = count_matrix[count_cols].corr()
        sns.heatmap(corr, annot=True, cmap="RdYlBu_r", vmin=0.5, vmax=1,
                    square=True, linewidths=0.5)
        plt.title("Sample Correlation (Pearson)")
        plt.tight_layout()
        plt.savefig(plots_dir / "sample_correlation.png", dpi=300, bbox_inches="tight")
        plt.close()
        print(f"  Saved: {plots_dir / 'sample_correlation.png'}")

    # 3. Volcano plot (if DE results available)
    if de_results is not None and len(de_results) > 0:
        plt.figure(figsize=(10, 8))

        # Add -log10(padj)
        de_results["-log10(padj)"] = -np.log10(de_results["padj"].fillna(1))

        # Color by significance
        colors = []
        for _, row in de_results.iterrows():
            if row["regulation"] == "Up":
                colors.append("red")
            elif row["regulation"] == "Down":
                colors.append("blue")
            else:
                colors.append("gray")

        plt.scatter(
            de_results["log2FoldChange"],
            de_results["-log10(padj)"],
            c=colors,
            alpha=0.6,
            s=20
        )

        plt.axhline(-np.log10(0.05), color="black", linestyle="--", alpha=0.5)
        plt.axvline(-1, color="gray", linestyle="--", alpha=0.3)
        plt.axvline(1, color="gray", linestyle="--", alpha=0.3)

        plt.xlabel("Log2 Fold Change")
        plt.ylabel("-Log10(Adjusted P-value)")
        plt.title("Differential Expression Volcano Plot")

        # Add legend
        from matplotlib.lines import Line2D
        legend_elements = [
            Line2D([0], [0], marker='o', color='w', markerfacecolor='red',
                   markersize=10, label='Upregulated'),
            Line2D([0], [0], marker='o', color='w', markerfacecolor='blue',
                   markersize=10, label='Downregulated'),
            Line2D([0], [0], marker='o', color='w', markerfacecolor='gray',
                   markersize=10, label='Not significant'),
        ]
        plt.legend(handles=legend_elements, loc='upper right')

        plt.savefig(plots_dir / "volcano_plot.png", dpi=300, bbox_inches="tight")
        plt.close()
        print(f"  Saved: {plots_dir / 'volcano_plot.png'}")

        # 4. MA plot
        plt.figure(figsize=(10, 6))
        de_results["log_baseMean"] = np.log10(de_results["baseMean"] + 1)

        plt.scatter(
            de_results["log_baseMean"],
            de_results["log2FoldChange"],
            c=colors,
            alpha=0.6,
            s=20
        )

        plt.axhline(0, color="black", linestyle="-", alpha=0.5)
        plt.xlabel("Log10(Base Mean + 1)")
        plt.ylabel("Log2 Fold Change")
        plt.title("MA Plot")
        plt.savefig(plots_dir / "ma_plot.png", dpi=300, bbox_inches="tight")
        plt.close()
        print(f"  Saved: {plots_dir / 'ma_plot.png'}")


# ============================================================================
# Report Generation
# ============================================================================

def generate_report(
    qc_results: Dict[str, QCStats],
    trim_results: Dict[str, TrimStats],
    count_matrix: pd.DataFrame,
    de_results: Optional[pd.DataFrame],
    output_dir: Path
):
    """Generate analysis report."""
    report_path = output_dir / "analysis_report.txt"

    with open(report_path, "w") as f:
        f.write("=" * 70 + "\n")
        f.write("miRNA-seq Analysis Report\n")
        f.write("=" * 70 + "\n\n")

        # QC Summary
        f.write("QUALITY CONTROL SUMMARY\n")
        f.write("-" * 40 + "\n")
        for sample, stats in qc_results.items():
            f.write(f"\n{sample}:\n")
            f.write(f"  Total reads: {stats.total_reads:,}\n")
            f.write(f"  Mean length: {stats.mean_length:.1f} bp\n")
            f.write(f"  Mean quality: {stats.mean_quality:.1f}\n")
            f.write(f"  GC content: {stats.gc_content:.1f}%\n")

        # Trimming Summary
        f.write("\n\nADAPTER TRIMMING SUMMARY\n")
        f.write("-" * 40 + "\n")
        for sample, stats in trim_results.items():
            f.write(f"\n{sample}:\n")
            f.write(f"  Input reads: {stats.reads_input:,}\n")
            f.write(f"  Output reads: {stats.reads_output:,}\n")
            f.write(f"  Adapter detection rate: {stats.adapter_rate:.1f}%\n")
            f.write(f"  Reads too short: {stats.reads_too_short:,}\n")
            f.write(f"  Reads too long: {stats.reads_too_long:,}\n")
            retention = stats.reads_output / stats.reads_input * 100
            f.write(f"  Retention rate: {retention:.1f}%\n")

        # Quantification Summary
        f.write("\n\nQUANTIFICATION SUMMARY\n")
        f.write("-" * 40 + "\n")
        count_cols = [c for c in count_matrix.columns if c != "length"]
        f.write(f"Unique sequences: {len(count_matrix)}\n")
        f.write(f"Samples: {len(count_cols)}\n")
        for col in count_cols:
            total = count_matrix[col].sum()
            f.write(f"  {col}: {total:,} reads\n")

        # Length distribution
        if "length" in count_matrix.columns:
            f.write("\nLength distribution (by read count):\n")
            length_counts = count_matrix.groupby("length")[count_cols[0]].sum()
            for length, count in sorted(length_counts.items()):
                if 18 <= length <= 30:
                    f.write(f"  {length} bp: {count:,}\n")

        # DE Summary
        if de_results is not None:
            f.write("\n\nDIFFERENTIAL EXPRESSION SUMMARY\n")
            f.write("-" * 40 + "\n")
            total = len(de_results)
            sig = de_results["significant"].sum()
            up = (de_results["regulation"] == "Up").sum()
            down = (de_results["regulation"] == "Down").sum()

            f.write(f"Total features tested: {total}\n")
            f.write(f"Significant (padj < 0.05): {sig}\n")
            f.write(f"  Upregulated: {up}\n")
            f.write(f"  Downregulated: {down}\n")

            if sig > 0:
                f.write("\nTop 20 significant features:\n")
                top = de_results[de_results["significant"]].nsmallest(20, "padj")
                for seq, row in top.iterrows():
                    f.write(f"  {seq[:30]}... LFC={row['log2FoldChange']:.2f} "
                            f"padj={row['padj']:.2e}\n")

        f.write("\n" + "=" * 70 + "\n")
        f.write("Analysis complete.\n")

    print(f"Report saved: {report_path}")


# ============================================================================
# Main Pipeline
# ============================================================================

def main():
    parser = argparse.ArgumentParser(
        description="miRNA-seq analysis pipeline",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Basic analysis
  python run_mirna_analysis.py \\
    --fastq data/*.fastq.gz \\
    --metadata metadata.csv \\
    --species human \\
    --output results/

  # With differential expression
  python run_mirna_analysis.py \\
    --fastq data/*.fastq.gz \\
    --metadata metadata.csv \\
    --species human \\
    --design "~condition" \\
    --contrast condition treated control \\
    --output results/ \\
    --plots
        """
    )

    parser.add_argument("--fastq", nargs="+", required=True,
                        help="Input FASTQ or FASTQ.gz files")
    parser.add_argument("--metadata", required=True,
                        help="Sample metadata CSV file")
    parser.add_argument("--species", default="human",
                        choices=["human", "mouse", "rat", "zebrafish"],
                        help="Species (default: human)")
    parser.add_argument("--adapter", default=None,
                        help="Adapter sequence or preset name (auto-detected if not specified)")
    parser.add_argument("--output", default="mirna_results",
                        help="Output directory (default: mirna_results)")
    parser.add_argument("--design", default=None,
                        help="Design formula for DE analysis (e.g., '~condition')")
    parser.add_argument("--contrast", nargs=3, default=None,
                        metavar=("VARIABLE", "TEST", "REFERENCE"),
                        help="Contrast specification for DE analysis")
    parser.add_argument("--min-length", type=int, default=18,
                        help="Minimum read length after trimming (default: 18)")
    parser.add_argument("--max-length", type=int, default=30,
                        help="Maximum read length after trimming (default: 30)")
    parser.add_argument("--min-count", type=int, default=10,
                        help="Minimum total count for feature filtering (default: 10)")
    parser.add_argument("--plots", action="store_true",
                        help="Generate visualization plots")
    parser.add_argument("--threads", type=int, default=4,
                        help="Number of threads (default: 4)")

    args = parser.parse_args()

    # Setup
    output_dir = Path(args.output)
    output_dir.mkdir(parents=True, exist_ok=True)
    (output_dir / "trimmed").mkdir(exist_ok=True)
    (output_dir / "counts").mkdir(exist_ok=True)

    # Load metadata
    print("Loading metadata...")
    metadata = pd.read_csv(args.metadata, index_col=0)

    # Get sample names from filenames
    sample_names = []
    for fastq in args.fastq:
        name = Path(fastq).stem
        if name.endswith(".fastq"):
            name = name[:-6]
        sample_names.append(name)

    print(f"Samples: {len(sample_names)}")

    # Step 1: Quality Control
    print("\n" + "=" * 60)
    print("STEP 1: Quality Control")
    print("=" * 60)

    qc_results = {}
    for fastq, name in zip(args.fastq, sample_names):
        print(f"  Analyzing {name}...")
        qc_results[name] = calculate_qc_stats(fastq)
        print(f"    Reads: {qc_results[name].total_reads:,}")
        print(f"    Mean quality: {qc_results[name].mean_quality:.1f}")

    # Step 2: Detect/Set Adapter
    print("\n" + "=" * 60)
    print("STEP 2: Adapter Detection")
    print("=" * 60)

    if args.adapter:
        adapter = args.adapter
        print(f"  Using specified adapter: {adapter}")
    else:
        adapter, rate = detect_adapter(args.fastq[0])
        print(f"  Auto-detected adapter: {adapter} ({rate:.1f}% detection rate)")

    adapter_seq = ADAPTER_PRESETS.get(adapter, adapter)
    print(f"  Adapter sequence: {adapter_seq}")

    # Step 3: Adapter Trimming
    print("\n" + "=" * 60)
    print("STEP 3: Adapter Trimming")
    print("=" * 60)

    trim_results = {}
    trimmed_files = []

    for fastq, name in zip(args.fastq, sample_names):
        print(f"  Trimming {name}...")
        output_path = output_dir / "trimmed" / f"{name}_trimmed.fastq.gz"
        trim_results[name] = trim_adapters(
            fastq,
            str(output_path),
            adapter_seq,
            min_length=args.min_length,
            max_length=args.max_length
        )
        trimmed_files.append(str(output_path))
        print(f"    Input: {trim_results[name].reads_input:,} reads")
        print(f"    Output: {trim_results[name].reads_output:,} reads")
        print(f"    Adapter rate: {trim_results[name].adapter_rate:.1f}%")

    # Step 4: Quantification
    print("\n" + "=" * 60)
    print("STEP 4: Sequence Quantification")
    print("=" * 60)

    count_matrix = count_sequences(trimmed_files, sample_names, min_count=args.min_count)

    # Save count matrix
    count_path = output_dir / "counts" / "raw_counts.csv"
    count_matrix.to_csv(count_path)
    print(f"  Saved raw counts: {count_path}")

    # Normalize
    normalized = normalize_counts(count_matrix, method="rpm")
    norm_path = output_dir / "counts" / "normalized_counts_rpm.csv"
    normalized.to_csv(norm_path)
    print(f"  Saved normalized counts: {norm_path}")

    # Step 5: Differential Expression (optional)
    de_results = None
    if args.design and args.contrast:
        print("\n" + "=" * 60)
        print("STEP 5: Differential Expression Analysis")
        print("=" * 60)

        de_results = run_differential_expression(
            count_matrix,
            metadata,
            args.design,
            args.contrast
        )

        if de_results is not None:
            de_path = output_dir / "de_results.csv"
            de_results.to_csv(de_path)
            print(f"  Saved DE results: {de_path}")

            sig_count = de_results["significant"].sum()
            print(f"  Significant features: {sig_count}")

    # Step 6: Visualization (optional)
    if args.plots:
        print("\n" + "=" * 60)
        print("STEP 6: Generating Plots")
        print("=" * 60)

        create_plots(count_matrix, de_results, output_dir, metadata)

    # Step 7: Generate Report
    print("\n" + "=" * 60)
    print("STEP 7: Generating Report")
    print("=" * 60)

    generate_report(qc_results, trim_results, count_matrix, de_results, output_dir)

    print("\n" + "=" * 60)
    print("ANALYSIS COMPLETE")
    print("=" * 60)
    print(f"Results saved to: {output_dir}")


if __name__ == "__main__":
    main()

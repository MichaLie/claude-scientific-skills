# miRNA-seq API Reference

Complete API documentation for the miRNA-seq analysis module.

## Table of Contents

1. [Core Classes](#core-classes)
2. [Quality Control Module](#quality-control-module)
3. [Adapter Trimming Module](#adapter-trimming-module)
4. [Alignment Module](#alignment-module)
5. [Quantification Module](#quantification-module)
6. [Differential Expression Module](#differential-expression-module)
7. [Visualization Module](#visualization-module)
8. [Utility Functions](#utility-functions)

---

## Core Classes

### MiRNASeqPipeline

Main pipeline class for complete miRNA-seq analysis.

```python
class MiRNASeqPipeline:
    """Complete miRNA sequencing analysis pipeline."""

    def __init__(
        self,
        output_dir: str,
        species: str = "human",
        adapter: str = None,
        threads: int = 4,
        mirbase_version: str = "22.1"
    ):
        """
        Initialize the miRNA-seq analysis pipeline.

        Parameters
        ----------
        output_dir : str
            Directory for output files. Created if doesn't exist.
        species : str, default "human"
            Species for miRBase reference. Options: "human", "mouse", "rat",
            or 3-letter code (hsa, mmu, rno).
        adapter : str, optional
            3' adapter sequence. If None, uses Illumina TruSeq Small RNA adapter.
        threads : int, default 4
            Number of threads for parallel processing.
        mirbase_version : str, default "22.1"
            miRBase version for reference sequences.
        """

    def run(
        self,
        fastq_files: List[str],
        sample_names: List[str],
        metadata: Union[str, pd.DataFrame],
        design: str = "~condition",
        contrast: List[str] = None,
        skip_steps: List[str] = None
    ) -> Dict:
        """
        Execute the complete analysis pipeline.

        Parameters
        ----------
        fastq_files : List[str]
            Paths to input FASTQ or FASTQ.gz files.
        sample_names : List[str]
            Sample names corresponding to each FASTQ file.
        metadata : str or DataFrame
            Path to metadata CSV or pandas DataFrame with sample information.
        design : str, default "~condition"
            Design formula for differential expression analysis.
        contrast : List[str], optional
            Contrast specification [variable, test, reference].
            Required for differential expression.
        skip_steps : List[str], optional
            Pipeline steps to skip. Options: "qc", "trim", "align", "count", "de".

        Returns
        -------
        Dict
            Dictionary containing results from each pipeline step:
            - "qc": Quality control results
            - "trimming": Trimming statistics
            - "alignment": Alignment results
            - "counts": Count matrix (pandas DataFrame)
            - "de": Differential expression results (if contrast provided)
        """

    def generate_report(
        self,
        results: Dict,
        format: str = "html",
        output_path: str = None
    ) -> str:
        """
        Generate analysis report.

        Parameters
        ----------
        results : Dict
            Results dictionary from run() method.
        format : str, default "html"
            Report format. Options: "html", "pdf", "markdown".
        output_path : str, optional
            Custom output path. If None, saves to output_dir/report.{format}.

        Returns
        -------
        str
            Path to generated report.
        """
```

### QCResult

Container for quality control results.

```python
@dataclass
class QCResult:
    """Quality control result container."""

    total_reads: int
    """Total number of reads in file."""

    total_bases: int
    """Total number of bases."""

    mean_length: float
    """Mean read length in bp."""

    min_length: int
    """Minimum read length."""

    max_length: int
    """Maximum read length."""

    mean_quality: float
    """Mean Phred quality score."""

    gc_content: float
    """GC content percentage (0-100)."""

    length_distribution: Dict[int, int]
    """Dictionary mapping read length to count."""

    quality_by_position: List[float]
    """Mean quality score at each position."""

    adapter_content: float
    """Percentage of reads containing adapter (0-100)."""

    overrepresented_sequences: List[Tuple[str, int]]
    """List of (sequence, count) for overrepresented sequences."""
```

### TrimmingResult

Container for adapter trimming results.

```python
@dataclass
class TrimmingResult:
    """Adapter trimming result container."""

    reads_input: int
    """Number of input reads."""

    reads_output: int
    """Number of output reads after filtering."""

    reads_with_adapter: int
    """Number of reads where adapter was found."""

    reads_too_short: int
    """Number of reads discarded (too short after trimming)."""

    reads_too_long: int
    """Number of reads discarded (too long)."""

    bases_trimmed: int
    """Total number of bases trimmed."""

    adapter_rate: float
    """Percentage of reads with adapter detected."""

    output_path: str
    """Path to trimmed FASTQ file."""
```

### AlignmentResult

Container for alignment results.

```python
@dataclass
class AlignmentResult:
    """Alignment result container."""

    total_reads: int
    """Total reads in input."""

    mapped_reads: int
    """Number of mapped reads."""

    uniquely_mapped: int
    """Number of uniquely mapped reads."""

    multi_mapped: int
    """Number of multi-mapped reads."""

    unmapped_reads: int
    """Number of unmapped reads."""

    mapping_rate: float
    """Mapping rate percentage (0-100)."""

    unique_mirnas: int
    """Number of unique miRNAs with at least 1 read."""

    output_bam: str
    """Path to output BAM file."""
```

---

## Quality Control Module

### Functions

```python
def analyze_fastq(
    fastq_path: str,
    output_dir: str = None,
    sample_limit: int = None
) -> QCResult:
    """
    Perform quality control analysis on FASTQ file.

    Parameters
    ----------
    fastq_path : str
        Path to FASTQ or FASTQ.gz file.
    output_dir : str, optional
        Directory for QC output files. If None, only returns results.
    sample_limit : int, optional
        Limit analysis to first N reads. If None, analyzes all reads.

    Returns
    -------
    QCResult
        Container with all QC metrics.

    Examples
    --------
    >>> qc = analyze_fastq("sample.fastq.gz")
    >>> print(f"Total reads: {qc.total_reads:,}")
    >>> print(f"Mean quality: {qc.mean_quality:.1f}")
    """


def detect_adapters(
    fastq_path: str,
    known_adapters: Dict[str, str] = None,
    sample_limit: int = 100000
) -> Dict:
    """
    Detect adapter sequences in FASTQ file.

    Parameters
    ----------
    fastq_path : str
        Path to FASTQ or FASTQ.gz file.
    known_adapters : Dict[str, str], optional
        Dictionary mapping adapter names to sequences.
        If None, uses common small RNA adapters.
    sample_limit : int, default 100000
        Number of reads to sample for detection.

    Returns
    -------
    Dict
        Dictionary containing:
        - "best_match": Name of most likely adapter
        - "best_match_sequence": Sequence of best match
        - "adapter_percentages": Dict mapping each adapter to detection rate
        - "total_reads_sampled": Number of reads analyzed

    Examples
    --------
    >>> adapters = detect_adapters("sample.fastq.gz")
    >>> print(f"Detected: {adapters['best_match']}")
    >>> print(f"Sequence: {adapters['best_match_sequence']}")
    """


def rna_composition(
    fastq_path: str,
    annotation_db: str = None
) -> Dict[str, float]:
    """
    Estimate RNA composition (miRNA, rRNA, tRNA, etc.).

    Parameters
    ----------
    fastq_path : str
        Path to trimmed FASTQ file.
    annotation_db : str, optional
        Path to annotation database. If None, uses built-in.

    Returns
    -------
    Dict[str, float]
        Dictionary mapping RNA type to percentage.

    Examples
    --------
    >>> comp = rna_composition("sample_trimmed.fastq.gz")
    >>> print(f"miRNA: {comp['mirna']:.1f}%")
    >>> print(f"rRNA: {comp['rrna']:.1f}%")
    """
```

---

## Adapter Trimming Module

### Functions

```python
def trim_adapters(
    input_fastq: str,
    output_fastq: str,
    adapter: str = "illumina_smallrna",
    min_length: int = 18,
    max_length: int = 30,
    quality_cutoff: int = 20,
    error_rate: float = 0.1,
    min_overlap: int = 3,
    trim_n: bool = True,
    discard_untrimmed: bool = False
) -> TrimmingResult:
    """
    Trim adapter sequences from reads.

    Parameters
    ----------
    input_fastq : str
        Path to input FASTQ or FASTQ.gz file.
    output_fastq : str
        Path for output trimmed FASTQ file.
    adapter : str, default "illumina_smallrna"
        Adapter sequence or preset name. Presets:
        - "illumina_smallrna": TGGAATTCTCGGGTGCCAAGG
        - "illumina_nebnext": AGATCGGAAGAGCACACGTCT
        - "qiagen_qiaseq": AACTGTAGGCACCATCAAT
    min_length : int, default 18
        Minimum read length after trimming.
    max_length : int, default 30
        Maximum read length after trimming.
    quality_cutoff : int, default 20
        Trim low-quality bases from 3' end below this threshold.
    error_rate : float, default 0.1
        Maximum allowed error rate for adapter matching.
    min_overlap : int, default 3
        Minimum overlap between read and adapter for trimming.
    trim_n : bool, default True
        Discard reads containing N bases.
    discard_untrimmed : bool, default False
        Discard reads where no adapter was found.

    Returns
    -------
    TrimmingResult
        Container with trimming statistics.

    Examples
    --------
    >>> result = trim_adapters(
    ...     "sample.fastq.gz",
    ...     "sample_trimmed.fastq.gz",
    ...     adapter="illumina_smallrna"
    ... )
    >>> print(f"Retained: {result.reads_output:,} reads")
    """


def batch_trim(
    input_files: List[str],
    output_dir: str,
    adapter: str = "illumina_smallrna",
    threads: int = 4,
    **kwargs
) -> Dict[str, TrimmingResult]:
    """
    Trim adapters from multiple files in parallel.

    Parameters
    ----------
    input_files : List[str]
        List of input FASTQ file paths.
    output_dir : str
        Directory for trimmed output files.
    adapter : str, default "illumina_smallrna"
        Adapter sequence or preset name.
    threads : int, default 4
        Number of parallel processes.
    **kwargs
        Additional arguments passed to trim_adapters().

    Returns
    -------
    Dict[str, TrimmingResult]
        Dictionary mapping input filename to trimming results.
    """
```

### Adapter Presets

```python
ADAPTER_PRESETS = {
    "illumina_smallrna": "TGGAATTCTCGGGTGCCAAGG",
    "illumina_truseq": "TGGAATTCTCGGGTGCCAAGG",
    "illumina_nebnext": "AGATCGGAAGAGCACACGTCT",
    "qiagen_qiaseq": "AACTGTAGGCACCATCAAT",
    "lexogen_smallrna": "TGGAATTCTCGGGTGCCAAGGAACTCCAGTCAC",
    "takara_smarter": "AAAAAAAAAA",  # Poly-A
    "nextflex": "TGGAATTCTCGGGTGCCAAGG",
}
```

---

## Alignment Module

### Functions

```python
def prepare_mirbase_index(
    species: str,
    mirbase_version: str = "22.1",
    include_hairpins: bool = True,
    output_dir: str = None
) -> str:
    """
    Download and prepare miRBase reference index.

    Parameters
    ----------
    species : str
        Species identifier. Options:
        - Full name: "human", "mouse", "rat"
        - 3-letter code: "hsa", "mmu", "rno"
    mirbase_version : str, default "22.1"
        miRBase version to download.
    include_hairpins : bool, default True
        Include precursor (hairpin) sequences in addition to mature.
    output_dir : str, optional
        Directory for index files. If None, uses default cache.

    Returns
    -------
    str
        Path to index prefix.
    """


def align_to_mirbase(
    fastq: str,
    index: str,
    output_bam: str,
    allow_mismatches: int = 1,
    seed_length: int = 18,
    report_all: bool = False,
    threads: int = 4
) -> AlignmentResult:
    """
    Align reads to miRBase reference.

    Parameters
    ----------
    fastq : str
        Path to trimmed FASTQ file.
    index : str
        Path to bowtie index prefix.
    output_bam : str
        Path for output BAM file.
    allow_mismatches : int, default 1
        Number of mismatches allowed in alignment.
    seed_length : int, default 18
        Minimum seed length for alignment.
    report_all : bool, default False
        If True, report all valid alignments (for multi-mapper analysis).
        If False, report only best alignment.
    threads : int, default 4
        Number of alignment threads.

    Returns
    -------
    AlignmentResult
        Container with alignment statistics.
    """


def align_to_genome(
    fastq: str,
    genome_index: str,
    output_bam: str,
    annotation_gtf: str = None,
    threads: int = 4
) -> AlignmentResult:
    """
    Align reads to full genome (alternative to miRBase-only alignment).

    Parameters
    ----------
    fastq : str
        Path to trimmed FASTQ file.
    genome_index : str
        Path to genome index (STAR or bowtie2).
    output_bam : str
        Path for output BAM file.
    annotation_gtf : str, optional
        Path to annotation GTF for feature assignment.
    threads : int, default 4
        Number of alignment threads.

    Returns
    -------
    AlignmentResult
        Container with alignment statistics.
    """
```

---

## Quantification Module

### Functions

```python
def count_mirnas(
    bam_files: List[str],
    sample_names: List[str],
    count_mode: str = "unique",
    min_mapq: int = 10,
    min_count: int = 5,
    normalize: bool = False
) -> pd.DataFrame:
    """
    Generate miRNA count matrix from BAM files.

    Parameters
    ----------
    bam_files : List[str]
        Paths to aligned BAM files.
    sample_names : List[str]
        Sample names for count matrix columns.
    count_mode : str, default "unique"
        How to count reads. Options:
        - "unique": Count only uniquely mapped reads
        - "fractional": Distribute multi-mappers equally
        - "all": Count all mapped reads (multi-mappers counted multiple times)
    min_mapq : int, default 10
        Minimum mapping quality for counting.
    min_count : int, default 5
        Filter miRNAs with less than this total count across samples.
    normalize : bool, default False
        If True, return normalized counts (RPM) instead of raw.

    Returns
    -------
    pd.DataFrame
        Count matrix with miRNAs as rows and samples as columns.
    """


def normalize_counts(
    count_matrix: pd.DataFrame,
    method: str = "rpm",
    size_factors: List[float] = None
) -> pd.DataFrame:
    """
    Normalize count matrix.

    Parameters
    ----------
    count_matrix : pd.DataFrame
        Raw count matrix (miRNAs × samples).
    method : str, default "rpm"
        Normalization method. Options:
        - "rpm": Reads Per Million
        - "cpm": Counts Per Million (same as RPM)
        - "tmm": Trimmed Mean of M-values
        - "rle": Relative Log Expression (DESeq2 method)
        - "uq": Upper Quartile
        - "custom": Use provided size_factors
    size_factors : List[float], optional
        Custom size factors for normalization. Required if method="custom".

    Returns
    -------
    pd.DataFrame
        Normalized count matrix.
    """


def filter_low_counts(
    count_matrix: pd.DataFrame,
    min_count: int = 10,
    min_samples: int = None,
    min_cpm: float = None
) -> pd.DataFrame:
    """
    Filter low-count miRNAs.

    Parameters
    ----------
    count_matrix : pd.DataFrame
        Raw count matrix.
    min_count : int, default 10
        Minimum total count across all samples.
    min_samples : int, optional
        Minimum number of samples with non-zero count.
    min_cpm : float, optional
        Minimum CPM in at least one sample.

    Returns
    -------
    pd.DataFrame
        Filtered count matrix.
    """
```

---

## Differential Expression Module

### Functions

```python
def run_deseq2(
    counts: pd.DataFrame,
    metadata: pd.DataFrame,
    design: str = "~condition",
    contrast: List[str] = None,
    alpha: float = 0.05,
    shrink_lfc: bool = True,
    min_count: int = 10
) -> pd.DataFrame:
    """
    Run DESeq2 differential expression analysis.

    Parameters
    ----------
    counts : pd.DataFrame
        Raw count matrix (miRNAs × samples).
    metadata : pd.DataFrame
        Sample metadata with design variables.
    design : str, default "~condition"
        Design formula for modeling.
    contrast : List[str], optional
        Contrast specification [variable, test_level, reference_level].
    alpha : float, default 0.05
        Significance threshold for adjusted p-values.
    shrink_lfc : bool, default True
        Apply apeGLM LFC shrinkage for visualization.
    min_count : int, default 10
        Filter miRNAs with less than this total count.

    Returns
    -------
    pd.DataFrame
        Results table with columns:
        - baseMean: Mean normalized count
        - log2FoldChange: Log2 fold change
        - lfcSE: Standard error of LFC
        - stat: Wald statistic
        - pvalue: Raw p-value
        - padj: Adjusted p-value (BH)
        - significant: Boolean significance flag
        - regulation: "Up", "Down", or "NS"
    """


def run_edger(
    counts: pd.DataFrame,
    metadata: pd.DataFrame,
    design: str = "~condition",
    contrast: List[str] = None,
    alpha: float = 0.05
) -> pd.DataFrame:
    """
    Run edgeR differential expression analysis (alternative to DESeq2).

    Parameters
    ----------
    counts : pd.DataFrame
        Raw count matrix.
    metadata : pd.DataFrame
        Sample metadata.
    design : str, default "~condition"
        Design formula.
    contrast : List[str], optional
        Contrast specification.
    alpha : float, default 0.05
        Significance threshold.

    Returns
    -------
    pd.DataFrame
        Results table with similar columns to run_deseq2().
    """


def multi_group_comparison(
    counts: pd.DataFrame,
    metadata: pd.DataFrame,
    condition_col: str = "condition",
    reference: str = None,
    method: str = "deseq2"
) -> Dict[str, pd.DataFrame]:
    """
    Perform all pairwise comparisons.

    Parameters
    ----------
    counts : pd.DataFrame
        Raw count matrix.
    metadata : pd.DataFrame
        Sample metadata.
    condition_col : str, default "condition"
        Column name containing condition labels.
    reference : str, optional
        Reference condition for comparisons. If None, uses first level.
    method : str, default "deseq2"
        DE method. Options: "deseq2", "edger".

    Returns
    -------
    Dict[str, pd.DataFrame]
        Dictionary mapping comparison name to results DataFrame.
    """
```

---

## Visualization Module

### Functions

```python
def volcano_plot(
    de_results: pd.DataFrame,
    padj_threshold: float = 0.05,
    lfc_threshold: float = 1.0,
    label_top_n: int = 10,
    point_size: int = 20,
    alpha: float = 0.7,
    title: str = "Volcano Plot",
    figsize: Tuple[int, int] = (10, 8)
) -> Tuple[plt.Figure, plt.Axes]:
    """
    Create volcano plot of differential expression results.

    Parameters
    ----------
    de_results : pd.DataFrame
        Differential expression results with padj and log2FoldChange columns.
    padj_threshold : float, default 0.05
        Significance threshold for coloring.
    lfc_threshold : float, default 1.0
        Log2 fold change threshold for labeling.
    label_top_n : int, default 10
        Number of top significant miRNAs to label.
    point_size : int, default 20
        Size of points.
    alpha : float, default 0.7
        Point transparency.
    title : str, default "Volcano Plot"
        Plot title.
    figsize : Tuple[int, int], default (10, 8)
        Figure size in inches.

    Returns
    -------
    Tuple[plt.Figure, plt.Axes]
        Matplotlib figure and axes objects.
    """


def expression_heatmap(
    counts: pd.DataFrame,
    mirnas: List[str] = None,
    metadata: pd.DataFrame = None,
    cluster_samples: bool = True,
    cluster_mirnas: bool = True,
    scale: str = "row",
    cmap: str = "RdBu_r",
    figsize: Tuple[int, int] = (12, 10)
) -> Tuple[plt.Figure, plt.Axes]:
    """
    Create expression heatmap.

    Parameters
    ----------
    counts : pd.DataFrame
        Normalized count matrix.
    mirnas : List[str], optional
        Subset of miRNAs to plot. If None, uses all.
    metadata : pd.DataFrame, optional
        Sample metadata for annotation bar.
    cluster_samples : bool, default True
        Hierarchically cluster samples.
    cluster_mirnas : bool, default True
        Hierarchically cluster miRNAs.
    scale : str, default "row"
        Scaling method. Options: "row", "column", None.
    cmap : str, default "RdBu_r"
        Colormap name.
    figsize : Tuple[int, int], default (12, 10)
        Figure size in inches.

    Returns
    -------
    Tuple[plt.Figure, plt.Axes]
        Matplotlib figure and axes objects.
    """


def length_distribution(
    fastq: str,
    sample_limit: int = 100000,
    title: str = "Read Length Distribution",
    figsize: Tuple[int, int] = (10, 6)
) -> Tuple[plt.Figure, plt.Axes]:
    """
    Plot read length distribution.

    Parameters
    ----------
    fastq : str
        Path to FASTQ file.
    sample_limit : int, default 100000
        Number of reads to sample.
    title : str, default "Read Length Distribution"
        Plot title.
    figsize : Tuple[int, int], default (10, 6)
        Figure size.

    Returns
    -------
    Tuple[plt.Figure, plt.Axes]
        Matplotlib figure and axes objects.
    """


def pca_plot(
    counts: pd.DataFrame,
    metadata: pd.DataFrame,
    color_by: str = "condition",
    shape_by: str = None,
    n_components: int = 2,
    figsize: Tuple[int, int] = (10, 8)
) -> Tuple[plt.Figure, plt.Axes]:
    """
    Create PCA plot of samples.

    Parameters
    ----------
    counts : pd.DataFrame
        Normalized count matrix.
    metadata : pd.DataFrame
        Sample metadata.
    color_by : str, default "condition"
        Metadata column for point colors.
    shape_by : str, optional
        Metadata column for point shapes.
    n_components : int, default 2
        Number of PCA components.
    figsize : Tuple[int, int], default (10, 8)
        Figure size.

    Returns
    -------
    Tuple[plt.Figure, plt.Axes]
        Matplotlib figure and axes objects.
    """


def ma_plot(
    de_results: pd.DataFrame,
    padj_threshold: float = 0.05,
    figsize: Tuple[int, int] = (10, 6)
) -> Tuple[plt.Figure, plt.Axes]:
    """
    Create MA plot (log fold change vs mean expression).

    Parameters
    ----------
    de_results : pd.DataFrame
        Differential expression results.
    padj_threshold : float, default 0.05
        Significance threshold.
    figsize : Tuple[int, int], default (10, 6)
        Figure size.

    Returns
    -------
    Tuple[plt.Figure, plt.Axes]
        Matplotlib figure and axes objects.
    """
```

---

## Utility Functions

### File I/O

```python
def read_fastq(
    filepath: str,
    limit: int = None
) -> Iterator[SeqRecord]:
    """
    Read FASTQ file (handles gzip automatically).

    Parameters
    ----------
    filepath : str
        Path to FASTQ or FASTQ.gz file.
    limit : int, optional
        Maximum number of records to yield.

    Yields
    ------
    SeqRecord
        Biopython sequence record.
    """


def write_fastq(
    records: Iterator[SeqRecord],
    filepath: str,
    compress: bool = True
):
    """
    Write sequences to FASTQ file.

    Parameters
    ----------
    records : Iterator[SeqRecord]
        Sequence records to write.
    filepath : str
        Output file path.
    compress : bool, default True
        Gzip compress output.
    """


def stream_fastq(
    filepath: str,
    chunk_size: int = 100000
) -> Iterator[List[SeqRecord]]:
    """
    Stream FASTQ file in chunks for memory-efficient processing.

    Parameters
    ----------
    filepath : str
        Path to FASTQ file.
    chunk_size : int, default 100000
        Number of records per chunk.

    Yields
    ------
    List[SeqRecord]
        Chunk of sequence records.
    """
```

### miRBase Utilities

```python
def download_mirbase(
    species: str,
    version: str = "22.1",
    output_dir: str = None
) -> Dict[str, str]:
    """
    Download miRBase sequences for a species.

    Parameters
    ----------
    species : str
        Species identifier (e.g., "hsa", "mmu").
    version : str, default "22.1"
        miRBase version.
    output_dir : str, optional
        Output directory for files.

    Returns
    -------
    Dict[str, str]
        Dictionary with paths to downloaded files:
        - "mature": Path to mature miRNA FASTA
        - "hairpin": Path to precursor FASTA
    """


def parse_mirbase_id(mirna_id: str) -> Dict[str, str]:
    """
    Parse miRBase identifier into components.

    Parameters
    ----------
    mirna_id : str
        miRNA identifier (e.g., "hsa-miR-21-5p").

    Returns
    -------
    Dict[str, str]
        Dictionary with:
        - "species": Species code (e.g., "hsa")
        - "name": miRNA name (e.g., "miR-21")
        - "arm": Arm designation (e.g., "5p")
        - "full": Original identifier
    """


def get_mirna_family(mirna_id: str) -> str:
    """
    Get miRNA family from identifier.

    Parameters
    ----------
    mirna_id : str
        miRNA identifier.

    Returns
    -------
    str
        Family name (e.g., "miR-21" -> "miR-21").
    """
```

---

## Constants

```python
# Species codes
SPECIES_CODES = {
    "human": "hsa",
    "mouse": "mmu",
    "rat": "rno",
    "zebrafish": "dre",
    "drosophila": "dme",
    "c_elegans": "cel",
}

# Default quality thresholds
DEFAULT_MIN_QUALITY = 20
DEFAULT_MIN_LENGTH = 18
DEFAULT_MAX_LENGTH = 30

# miRBase current version
MIRBASE_CURRENT_VERSION = "22.1"
```

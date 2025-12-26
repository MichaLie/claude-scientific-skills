# miRNA Database Reference

This guide covers miRNA databases and resources used in miRNA-seq analysis.

## Table of Contents

1. [miRBase](#mirbase)
2. [Target Prediction Databases](#target-prediction-databases)
3. [Validated Target Databases](#validated-target-databases)
4. [Disease Association Databases](#disease-association-databases)
5. [IsomiR Databases](#isomir-databases)
6. [API Access](#api-access)

---

## miRBase

The primary database for miRNA sequences and annotations.

### Overview

- **URL**: https://www.mirbase.org/
- **Current Version**: 22.1 (March 2018)
- **Content**: 38,589 hairpin precursors, 48,860 mature miRNAs across 271 species

### File Downloads

```python
import requests
from pathlib import Path

def download_mirbase_files(species="hsa", version="22.1", output_dir="."):
    """Download miRBase files for a species."""
    base_url = "https://www.mirbase.org/download"
    output_dir = Path(output_dir)
    output_dir.mkdir(exist_ok=True)

    files = {
        "mature": f"{base_url}/mature.fa.gz",
        "hairpin": f"{base_url}/hairpin.fa.gz",
        "aliases": f"{base_url}/aliases.txt.gz",
        "mirna_dat": f"{base_url}/miRNA.dat.gz"
    }

    downloaded = {}
    for name, url in files.items():
        output_path = output_dir / f"{name}.fa.gz"
        response = requests.get(url, stream=True)
        with open(output_path, "wb") as f:
            for chunk in response.iter_content(chunk_size=8192):
                f.write(chunk)
        downloaded[name] = output_path

    return downloaded
```

### Species-Specific Extraction

```python
import gzip
from Bio import SeqIO

def extract_species_mirnas(fasta_path, species="hsa"):
    """Extract miRNAs for a specific species from miRBase FASTA."""
    mirnas = {}

    opener = gzip.open if str(fasta_path).endswith(".gz") else open
    mode = "rt" if str(fasta_path).endswith(".gz") else "r"

    with opener(fasta_path, mode) as handle:
        for record in SeqIO.parse(handle, "fasta"):
            if record.id.startswith(species):
                mirnas[record.id] = str(record.seq)

    return mirnas

# Example usage
human_mature = extract_species_mirnas("mature.fa.gz", "hsa")
print(f"Human mature miRNAs: {len(human_mature)}")
```

### miRNA Naming Convention

| Component | Example | Description |
|-----------|---------|-------------|
| Species | hsa | 3-letter species code |
| Type | miR | Mature miRNA (vs mir for precursor) |
| Number | 21 | Unique identifier |
| Letter suffix | a/b/c | Paralog variants |
| Arm | 5p/3p | Precursor arm of origin |

**Examples:**
- `hsa-miR-21-5p`: Human, mature, miR-21, 5' arm
- `hsa-mir-21`: Human, precursor hairpin
- `mmu-miR-155-3p`: Mouse, mature, miR-155, 3' arm

### Species Codes

| Species | Code | Common Name |
|---------|------|-------------|
| Homo sapiens | hsa | Human |
| Mus musculus | mmu | Mouse |
| Rattus norvegicus | rno | Rat |
| Danio rerio | dre | Zebrafish |
| Drosophila melanogaster | dme | Fruit fly |
| Caenorhabditis elegans | cel | Nematode |
| Arabidopsis thaliana | ath | Thale cress |
| Gallus gallus | gga | Chicken |

---

## Target Prediction Databases

### TargetScan

Predicts biological targets of miRNAs by searching for conserved sites matching the seed region.

- **URL**: https://www.targetscan.org/
- **Species**: Human, mouse, fly, worm, fish
- **Method**: Seed matching + context scores

```python
import requests
import pandas as pd

def query_targetscan(mirna, species="human"):
    """Query TargetScan for miRNA targets."""
    # Note: TargetScan requires downloading data files
    # This is a conceptual example

    species_map = {"human": "Human", "mouse": "Mouse"}
    base_url = f"https://www.targetscan.org/cgi-bin/targetscan/vert_80/targetscan.cgi"

    # In practice, download and parse their data files
    # Available at: https://www.targetscan.org/vert_80/vert_80_data_download/

    return pd.DataFrame()  # Placeholder
```

**TargetScan File Downloads:**
```bash
# Download conserved family info
wget https://www.targetscan.org/vert_80/vert_80_data_download/Conserved_Family_Info.txt.zip

# Download predicted targets
wget https://www.targetscan.org/vert_80/vert_80_data_download/Predicted_Targets_Info.default_predictions.txt.zip
```

### miRDB

- **URL**: http://mirdb.org/
- **Content**: Computationally predicted miRNA targets
- **Method**: Machine learning (MirTarget)

```python
def download_mirdb_predictions(species="human", output_path="mirdb_predictions.txt"):
    """Download miRDB prediction file."""
    urls = {
        "human": "http://mirdb.org/download/miRDB_v6.0_prediction_result_human.txt.gz",
        "mouse": "http://mirdb.org/download/miRDB_v6.0_prediction_result_mouse.txt.gz",
        "rat": "http://mirdb.org/download/miRDB_v6.0_prediction_result_rat.txt.gz",
    }

    if species not in urls:
        raise ValueError(f"Species must be one of: {list(urls.keys())}")

    response = requests.get(urls[species])
    with open(output_path, "wb") as f:
        f.write(response.content)

    return output_path
```

### DIANA-microT

- **URL**: http://diana.imis.athena-innovation.gr/DianaTools/
- **Method**: MicroT-CDS algorithm
- **Features**: 3'UTR and CDS binding sites

---

## Validated Target Databases

### miRTarBase

Experimentally validated miRNA-target interactions.

- **URL**: https://mirtarbase.cuhk.edu.cn/
- **Content**: >2.2 million validated interactions
- **Evidence Types**: Reporter assay, Western blot, qPCR, etc.

```python
import pandas as pd

def download_mirtarbase(version="9.0", output_path="mirtarbase.xlsx"):
    """Download miRTarBase validated interactions."""
    url = f"https://mirtarbase.cuhk.edu.cn/~miRTarBase/miRTarBase_2022/cache/download/9.0/miRTarBase_MTI.xlsx"

    response = requests.get(url)
    with open(output_path, "wb") as f:
        f.write(response.content)

    return pd.read_excel(output_path)

def query_mirtarbase(mirna, mirtarbase_df, evidence="strong"):
    """Query miRTarBase for validated targets."""
    # Filter by miRNA
    targets = mirtarbase_df[mirtarbase_df["miRNA"] == mirna]

    # Filter by evidence strength
    if evidence == "strong":
        strong_methods = ["Luciferase reporter assay", "Western blot", "qRT-PCR"]
        targets = targets[targets["Experiments"].str.contains("|".join(strong_methods), na=False)]

    return targets
```

**Evidence Categories:**
| Category | Methods |
|----------|---------|
| Strong | Reporter assay, Western blot, qRT-PCR |
| Functional | Knockdown, Overexpression |
| Computational | CLIP-seq, Degradome-seq |

### TarBase

- **URL**: https://dianalab.e-ce.uth.gr/html/diana/web/index.php?r=tarbasev8
- **Content**: Experimentally supported interactions
- **Integration**: Links to DIANA tools

---

## Disease Association Databases

### HMDD (Human microRNA Disease Database)

- **URL**: https://www.cuilab.cn/hmdd
- **Content**: miRNA-disease associations
- **Evidence**: Literature-curated

```python
def query_hmdd(mirna=None, disease=None):
    """Query HMDD for miRNA-disease associations."""
    # HMDD provides downloadable files
    # Manual download required from website

    # Example data structure
    associations = pd.DataFrame({
        "miRNA": ["hsa-miR-21-5p", "hsa-miR-155-5p"],
        "Disease": ["Breast Cancer", "Inflammation"],
        "PMID": ["12345678", "23456789"],
        "Description": ["Upregulated in tumor", "Immune response"]
    })

    if mirna:
        associations = associations[associations["miRNA"] == mirna]
    if disease:
        associations = associations[associations["Disease"].str.contains(disease, case=False)]

    return associations
```

### miR2Disease

- **URL**: http://www.mir2disease.org/
- **Focus**: Human miRNA-disease relationships
- **Content**: Manually curated from literature

### DBDEMC (Database of Differentially Expressed MiRNAs in Human Cancers)

- **URL**: https://www.biosino.org/dbDEMC/
- **Focus**: Cancer-specific miRNA expression
- **Content**: DE miRNAs from cancer studies

---

## IsomiR Databases

### isomiRage

- **URL**: https://cru.genomics.iit.it/isomirage/
- **Content**: IsomiR profiles from sRNA-seq data
- **Features**: 5'/3' variations, non-templated additions

### YM500

- **URL**: http://ngs.ym.edu.tw/ym500/
- **Content**: Small RNA-seq analysis platform
- **Features**: IsomiR detection, arm switching

### IsomiR Classification

| Type | Description | Example |
|------|-------------|---------|
| 5' isomiR | Different 5' start position | +1, -1 nucleotide |
| 3' isomiR | Different 3' end position | +1, -2 nucleotides |
| NTA | Non-templated addition | 3' uridylation |
| SNP | Internal modification | A>G editing |

```python
def classify_isomir(read_seq, canonical_seq, alignment):
    """Classify isomiR variant type."""
    classification = {
        "type": "canonical",
        "5p_offset": 0,
        "3p_offset": 0,
        "nta": None,
        "modifications": []
    }

    # Check 5' offset
    if alignment.reference_start != 0:
        classification["type"] = "5p_isomir"
        classification["5p_offset"] = alignment.reference_start

    # Check 3' offset
    len_diff = len(read_seq) - len(canonical_seq)
    if len_diff != 0:
        classification["type"] = "3p_isomir"
        classification["3p_offset"] = len_diff

        # Check for non-templated additions
        if len_diff > 0:
            tail = read_seq[-len_diff:]
            if all(b in "AU" for b in tail):
                classification["nta"] = tail

    return classification
```

---

## API Access

### miRBase REST API

```python
import requests

def mirbase_lookup(mirna_id):
    """Look up miRNA information from miRBase."""
    # Note: miRBase doesn't have a public REST API
    # This shows the conceptual approach using web scraping
    # or downloaded data files

    url = f"https://www.mirbase.org/cgi-bin/mirna_entry.pl?acc={mirna_id}"
    # Parse response...
    return {}
```

### Ensembl REST API for miRNA Annotations

```python
def ensembl_mirna_lookup(mirna_id, species="human"):
    """Look up miRNA from Ensembl."""
    server = "https://rest.ensembl.org"
    ext = f"/lookup/symbol/{species}/{mirna_id}"

    response = requests.get(
        server + ext,
        headers={"Content-Type": "application/json"}
    )

    if response.ok:
        return response.json()
    return None
```

### NCBI Gene for miRNA Information

```python
from Bio import Entrez

def ncbi_mirna_search(mirna_name, email="your@email.com"):
    """Search NCBI Gene for miRNA information."""
    Entrez.email = email

    # Search for miRNA
    handle = Entrez.esearch(db="gene", term=f"{mirna_name}[gene] AND human[organism]")
    record = Entrez.read(handle)
    handle.close()

    if record["IdList"]:
        gene_id = record["IdList"][0]

        # Fetch gene details
        handle = Entrez.efetch(db="gene", id=gene_id, rettype="xml")
        details = Entrez.read(handle)
        handle.close()

        return details

    return None
```

---

## Local Database Setup

### SQLite Database for miRNA Data

```python
import sqlite3
import pandas as pd

def create_mirna_database(db_path="mirna_data.db"):
    """Create local SQLite database for miRNA data."""
    conn = sqlite3.connect(db_path)

    # Create tables
    conn.execute("""
        CREATE TABLE IF NOT EXISTS mirnas (
            id TEXT PRIMARY KEY,
            name TEXT,
            species TEXT,
            sequence TEXT,
            arm TEXT,
            precursor_id TEXT
        )
    """)

    conn.execute("""
        CREATE TABLE IF NOT EXISTS targets (
            mirna_id TEXT,
            gene_symbol TEXT,
            gene_id TEXT,
            evidence TEXT,
            score REAL,
            source TEXT,
            FOREIGN KEY (mirna_id) REFERENCES mirnas(id)
        )
    """)

    conn.execute("""
        CREATE TABLE IF NOT EXISTS diseases (
            mirna_id TEXT,
            disease TEXT,
            association TEXT,
            pmid TEXT,
            FOREIGN KEY (mirna_id) REFERENCES mirnas(id)
        )
    """)

    conn.commit()
    return conn

def populate_from_mirbase(conn, mirbase_fasta, species="hsa"):
    """Populate database from miRBase FASTA."""
    mirnas = extract_species_mirnas(mirbase_fasta, species)

    data = []
    for mirna_id, sequence in mirnas.items():
        parts = mirna_id.split("-")
        species_code = parts[0]
        name = "-".join(parts[1:-1]) if len(parts) > 2 else parts[1]
        arm = parts[-1] if parts[-1] in ["5p", "3p"] else None

        data.append({
            "id": mirna_id,
            "name": name,
            "species": species_code,
            "sequence": sequence,
            "arm": arm,
            "precursor_id": f"{species_code}-{name}"
        })

    df = pd.DataFrame(data)
    df.to_sql("mirnas", conn, if_exists="replace", index=False)

    return len(data)
```

---

## Data Integration Example

```python
class MiRNADatabase:
    """Unified interface for miRNA database queries."""

    def __init__(self, cache_dir=".mirna_cache"):
        self.cache_dir = Path(cache_dir)
        self.cache_dir.mkdir(exist_ok=True)

    def get_mirna_info(self, mirna_id):
        """Get comprehensive information about a miRNA."""
        info = {
            "id": mirna_id,
            "sequence": self._get_sequence(mirna_id),
            "targets": {
                "predicted": self._get_predicted_targets(mirna_id),
                "validated": self._get_validated_targets(mirna_id)
            },
            "diseases": self._get_disease_associations(mirna_id),
            "family": self._get_family_members(mirna_id)
        }
        return info

    def _get_sequence(self, mirna_id):
        """Get miRNA sequence from miRBase."""
        # Implementation
        pass

    def _get_predicted_targets(self, mirna_id):
        """Get predicted targets from TargetScan/miRDB."""
        # Implementation
        pass

    def _get_validated_targets(self, mirna_id):
        """Get validated targets from miRTarBase."""
        # Implementation
        pass

    def _get_disease_associations(self, mirna_id):
        """Get disease associations from HMDD."""
        # Implementation
        pass

    def _get_family_members(self, mirna_id):
        """Get miRNA family members."""
        # Implementation
        pass
```

---

## Best Practices

1. **Version Control**: Always record database versions used in analysis
2. **Regular Updates**: miRBase updates periodically; maintain latest reference
3. **Cross-Reference**: Validate findings across multiple databases
4. **Evidence Filtering**: For targets, prefer experimentally validated over predicted
5. **Species Matching**: Ensure database species matches your samples
6. **Cache Data**: Download and cache database files for reproducibility

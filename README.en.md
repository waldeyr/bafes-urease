# BAFES-urease-mining

[![Language: Portuguese](https://img.shields.io/badge/Language-Portugu%C3%AAs--BR-blue)](README.md)
[![Language: English](https://img.shields.io/badge/Language-English--US-green)](README.en.md)

> 🌐 **Language Switch / Alternar Idioma**: Read this documentation in [Português-BR](README.md).

---

## 1. General Project Description

**BAFES-urease-mining** is an automated bioinformatics pipeline built in Nextflow DSL2 and Python 3.11, fully containerized via Docker. It is designed for genomic mining, functional identification, and evolutionary analysis of genes involved in urea catabolism across bacterial genomes.

The pipeline analyses whatever list of genomes you provide — GenBank/RefSeq accessions and, optionally, local FASTA assemblies (see [Input genomes](#input-genomes)) — and investigates three metabolic avenues in each one:

1. **Canonical Urease Pathway**: Catalytic subunits ($\gamma$ - UreA, $\beta$ - UreB, $\alpha$ - UreC) and accessory maturation proteins (UreD, UreE, UreF, UreG).
2. **Active Urea Transport System**: ABC transporter complex $urtABCDE$ (UrtA, UrtB, UrtC, UrtD, UrtE).
3. **Nickel-Independent Alternative Pathway**: Urea carboxylase (`uc`, EC 6.3.4.6) coupled with allophanate hydrolase (`ah`, EC 3.5.1.54).

**Example dataset.** The repository ships configured with the dataset that motivated the project:
10 aerobic endospore-forming bacteria (**BAFES**, from the Brazilian Portuguese *Bactérias
Aeróbicas Formadoras de Endósporos*) from the **CBAFES/UnB collection**, isolated from soils of
Distrito Federal, Brazil (accessions `VKHW00000000.1` to `VKIC00000000.1`, in
[data/accessions.tsv](data/accessions.tsv)). Throughout this README, counts and strain names (S1,
S2–S10) refer to that dataset and are for illustration only.

### Input genomes

The pipeline runs on the genomes listed in up to two tab-separated files with a header row:

| File | Option | Columns | Purpose |
| :--- | :--- | :--- | :--- |
| `accessions.tsv` | `--accessions` (default: `data/accessions.tsv`) | `strain`, `species`, `genbank_accession` | Public genomes, downloaded from NCBI. Takes an assembly (`GCA_`/`GCF_`) or a WGS master (`VKHW00000000.1`), which is resolved to its assembly before download. |
| Local genomes TSV | `--local-genomes` (optional) | `strain`, `species`, `fasta` | Assemblies that are not on NCBI. The FASTA path is relative to the launch directory; the genome skips the download and goes straight to QC and Bakta. Example: [data/local_genomes.tsv](data/local_genomes.tsv). |

Rules:

- `strain` is the short identifier that names every output (`results/00_genomes/{strain}.fna`,
  `{strain}_candidates.tsv`, …). It must not repeat, within a file or across the two: a repeated
  ID aborts the run before any task starts.
- `species` is informative (it shows in the report); it may be `Unknown`.
- Every genome, downloaded or local, goes through the same validation: contigs, >= 500 kb and the
  nucleotide alphabet. An invalid genome, or one that fails to download, stops the whole run
  (Nextflow's default error strategy); fix or drop the row and continue with `--resume`.
- The pre-check requires at least one accession in `accessions.tsv`.

```bash
# Other public genomes
./run.sh --exec --accessions my_list.tsv --outdir results_my_list --verbose

# Public genomes + local assemblies
./run.sh --exec --accessions my_list.tsv --local-genomes my_fastas.tsv --verbose
```

Two things depend on the taxonomic group of the input genomes:

- **BLASTp references.** The default set is Firmicutes, the phylum of the example dataset. For
  genomes from other groups, change the slice with `fetch_references.py --query` (see the reference
  set paragraph in [Data Integrity](#7-data-integrity)).
- **Resources.** Run time grows with the number of genomes. The CPU and memory limits in
  [nextflow.config](nextflow.config) are sized for a 96-core host, so that up to 10 strains clear
  hmmscan in a single wave.
---

## 2. Quick-Start

### Prerequisites

The whole pipeline runs **inside the `bafes-urease` container**. The host only needs:

- Docker Engine >= 24.0 (with your user in the `docker` group)
- Git, Bash, and `curl`
- Disk space: **~160 GB free** for the full Bakta database (use `--bakta-db-type light` for ~10 GB,
  or `--skip-bakta-db` if you already have one)

Nextflow, Bakta, BLAST+, HMMER, MAFFT, trimAl, IQ-TREE, QUAST, and clinker all ship inside the
image — none of them need to be installed on the host. To run without Docker (with a preconfigured
Conda/Mamba environment), use `--runtime local`.

### Three-Step Command Sequence

```bash
# Clone the repository and enable the wrapper
git clone https://github.com/waldeyr/bafes-urease.git
cd bafes-urease
chmod +x run.sh

# Step A — Bootstrap: builds the image, downloads Pfam and the full Bakta database.
#    First run takes hours (image build ~20 min + Bakta full DB ~75 GB). Idempotent:
#    every completed step is skipped on subsequent runs.
#    Upgrading from an install that used the old `bafes_urease` tag? Retag first to
#    avoid a ~20-minute rebuild:
#        docker tag bafes_urease bafes-urease && docker rmi bafes_urease
./run.sh --bootstrap --runtime docker --container-image bafes-urease

# Step B — Build & Pre-flight Check (validates environment & resource reachability)
./run.sh --build --runtime docker --container-image bafes-urease

# Step C — Pipeline Execution with real-time tqdm progress tracking & detailed logging
./run.sh --exec --runtime docker --container-image bafes-urease --verbose
```

Every `run.sh` invocation returns the prompt immediately and keeps running inside a detached
`screen` session — see [Automatic `screen` Session](#87-automatic-screen-session) for following the
log and reading the real exit code. Further parameter combinations (reduced database, resuming an
interrupted run, alternative paths, running without a container) are in
[Installation & Usage](#8-installation--usage).

---

## 3. Pipeline Overview

The pipeline executes a modular, multi-tier analysis workflow:

```mermaid
flowchart TD
    Preflight["00_preflight: Resource Pre-checks & Network Validation"] --> Download["00_genomes: NCBI Genome Download"]
    Download --> QC["01_qc: QUAST & CheckM2 Quality Control"]
    QC --> Bakta["02_bakta: Functional Annotation (Bakta)"]
    Bakta --> Candidates["03_candidates: 3-Layer Screening (Bakta + BLASTp + HMMER3)"]
    Candidates --> Summary["04_summary: Presence/Absence Matrix & Unified Summary"]
    Candidates --> Synteny["05_synteny: Gene Neighborhood Analysis (clinker ±10 kb)"]
    Candidates --> Phylogeny["06_phylogeny: UreC Phylogeny (MAFFT + trimAl + IQ-TREE 2)"]
```

---

## 4. Phase-by-Phase Detailed View, Rules & Filters

| Phase ID & Name | Process / Script | Description | Rules & Applied Filters | Output Artifacts |
| :--- | :--- | :--- | :--- | :--- |
| **00_preflight** | `resource_checker.py` / `PREFLIGHT_CHECK` | Pre-flight validation of external resources & network endpoints | Validates HTTP Status 200/JSON for every input NCBI accession, Pfam FTP URL, and Bakta Zenodo DB. Aborts with Exit Code 4 if critical resources are unreachable. | `results/00_preflight/preflight_status.txt` |
| **00_genomes** | `DOWNLOAD_GENOME` | Automated acquisition of bacterial assembly FASTA files | Queries NCBI Datasets CLI v2 by accession. Automatically falls back to NCBI Entrez `efetch` if Datasets API fails. Local genomes (`--local-genomes`) skip the download and go through the same validation in `STAGE_LOCAL_GENOME`. | `results/00_genomes/{strain}.fna` |
| **01_qc** | `QC_QUAST` & `QC_CHECKM2` | Genomic assembly quality control & contamination check | Evaluates N50, GC content, total length (QUAST), completeness >= 95%, and contamination <= 5% (CheckM2). | `results/01_qc/quast/`, `results/01_qc/checkm2/` |
| **02_bakta** | `BAKTA_ANNOTATE` | Structural & functional genome annotation | Uses Bakta v1.9+ with full database. Assigns RefSeq/UniRef IDs, EC numbers (3.5.1.5, 6.3.4.6, 3.5.1.54), and gene symbols. | `results/02_bakta/{strain}/` (.json, .gff3, .faa, .gbff) |
| **03_candidates** | `extract_urease_candidates.py` | 3-Layer candidate screening & evidence consolidation | **Layer A (Bakta)**: EC & gene symbol filtering.<br>**Layer B (BLASTp)**: E-value <= 1e-5, Identity >= 30%, Coverage >= 50%.<br>**Layer C (HMMER3)**: Domain score >= 20.0 against Pfam profiles.<br>**Score**: see [Candidate screening score](#candidate-screening-score) — High (>= 2 points), Medium (1), Low (TrEMBL BLASTp only). | `results/03_candidates/{strain}/{strain}_candidates.tsv` |
| **04_summary** | `merge_strain_results.py` | Cross-strain matrix aggregation | Merges candidate TSVs across all strains. Constructs binary (1/0) matrix for 11 target genes ($ureA$–$ureH$, $urtA$, $uc$, $ah$). | `results/04_summary/urease_presence_absence.tsv`, `summary_unified.tsv` |
| **05_synteny** | `synteny_analysis.py` | Gene neighborhood & synteny visualization | Extracts $\pm 10$ kb genomic window centered on $ureC$. Runs `clinker` to produce an interactive SVG/HTML map. | `results/05_synteny/synteny.html` |
| **06_phylogeny** | `build_phylogeny.py` | UreC evolutionary tree reconstruction | Requires >= 3 UreC protein sequences. Adds the reference UreC and the outgroup (`data/urec_outgroup.fasta`). Runs MAFFT `--auto` alignment $\rightarrow$ trimAl `-automated1` $\rightarrow$ IQ-TREE 2 (1000 ultrafast bootstraps, rooted on the outgroup with `-o`). | `results/06_phylogeny/ureC.treefile`, `ureC.aln` |

### Candidate screening score

Each layer is worth 1 point, with one exception: a BLASTp whose **best hit** is a TrEMBL
(unreviewed) reference is recorded but scores nothing. TrEMBL annotation is automatic, like
Bakta's, so agreeing with it is not independent evidence. The status of the best hit goes in the
`blast_ref_status` column, and the points in `evidence_points`.

| Evidence | Points |
| :--- | :--- |
| Layer A — Bakta: target gene, EC number or product | 1 |
| Layer B — BLASTp with its best hit on a reviewed (Swiss-Prot) reference | 1 |
| Layer B — BLASTp with its best hit on an unreviewed (TrEMBL) reference | 0 (recorded) |
| Layer C — HMMER: Pfam profile from the list, domain score >= 20 | 1 |

| Points | Confidence |
| :--- | :--- |
| >= 2 | High |
| 1 | Medium |
| 0, with a TrEMBL BLASTp hit | Low |

| Example | Points | Confidence |
| :--- | :--- | :--- |
| Bakta + BLASTp (Swiss-Prot) + HMMER | 3 | High |
| Bakta + BLASTp (Swiss-Prot) | 2 | High |
| Bakta + HMMER + BLASTp (TrEMBL) | 2 | High |
| Bakta + BLASTp (TrEMBL) | 1 | Medium |
| HMMER + BLASTp (TrEMBL) | 1 | Medium |
| Bakta only, or HMMER only | 1 | Medium |
| BLASTp (TrEMBL) only | 0 | Low |

Today only UrtA–E, `uca` and `atzF` come from TrEMBL. UrtA–E and `atzF` have no profile in the
Layer C list, so in practice they reach Medium at most. Urea carboxylase (`uca`) carries the CT_C_D
domain (PF02682), which is on the list, so it can still reach High with Bakta + HMMER. A custom reference FASTA (`--references`) whose
headers do not follow the `ref|gene|UniProt:accession|organism|status|Naa` format gets status
`unknown` and scores as reviewed.

---

## 5. Data Provenance & Lineage Diagram

The diagram below maps data flow from primary external databases to downstream analytical artifacts:

```mermaid
flowchart LR
    subgraph Primary_Sources ["Primary Data Sources"]
        NCBI["NCBI GenBank (nuccore / WGS)"]
        PfamDB["EMBL-EBI Pfam-A Database"]
        BaktaDB["Bakta Database (Zenodo DB)"]
        RefDB["UniProt Reference FASTA (data/urease_references.fasta)"]
    end

    subgraph Internal_Processing ["Pipeline Processing Engine"]
        Accessions["data/accessions.tsv"] -->|Accessions| Preflight["Pre-flight Check"]
        NCBI -->|FASTA Download| Genomes["00_genomes/"]
        Genomes -->|Assembly FASTA| QC["01_qc (QUAST & CheckM2)"]
        Genomes -->|Assembly FASTA| Bakta["02_bakta (Annotation)"]
        BaktaDB -->|RefSeq/UniRef| Bakta
        
        Bakta -->|FAA & JSON| Screening["03_candidates (BLASTp + HMMER3)"]
        PfamDB -->|PF00547, etc.| Screening
        RefDB -->|Reference Proteins| Screening
    end

    subgraph Derived_Outputs ["Derived Analytical Products"]
        Screening -->|Consolidated Candidates| Summary["04_summary (urease_presence_absence.tsv)"]
        Screening -->|GBFF Loci (±10 kb)| Synteny["05_synteny (synteny.html)"]
        Screening -->|UreC FASTA| Phylogeny["06_phylogeny (ureC.treefile)"]
    end
```

---

## 6. Technology & Software Versions

| Technology / Software | Version | Role in Pipeline | Reference / Source |
| :--- | :--- | :--- | :--- |
| **Linux Ubuntu (Jammy)** | 22.04 LTS | Container OS base | Canonical |
| **Python** | 3.11.x | Core scripting, consolidation, tqdm tracking | Python Software Foundation |
| **OpenJDK JRE** | 17.0+ | Runtime environment for Nextflow | Eclipse Temurin / OpenJDK |
| **Nextflow** | >= 24.04.0 | DSL2 Workflow Orchestrator | Seqera / Nextflow |
| **Docker Engine** | >= 24.0+ | Containerization & runtime isolation | Docker Inc. |
| **NCBI Datasets CLI / Entrez** | v2 / E-utilities | Automated genome assembly retrieval | NCBI / NLM / NIH |
| **Bakta** | >= 1.9.0 | Functional genome annotation engine | Schwengers et al. (2021) |
| **QUAST** | >= 5.2.0 | Genome assembly quality evaluation | Gurevich et al. (2013) |
| **CheckM2** | >= 1.0.2 | Genome completeness & contamination estimation | Chklovski et al. (2023) |
| **BLAST+ (blastp)** | >= 2.14.0 | Protein sequence alignment | Camacho et al. (2009) |
| **HMMER3 (hmmscan)** | >= 3.3.2 | Profile Hidden Markov Model search | Eddy et al. (2011) |
| **clinker** | >= 0.0.28 | Gene cluster synteny comparison & visualization | Gilchrist & Chooi (2021) |
| **MAFFT** | >= 7.520 | Multiple protein sequence alignment | Katoh & Standley (2013) |
| **trimAl** | >= 1.4.1 | Alignment trimming & gap filtering | Capella-Gutiérrez et al. (2009) |
| **IQ-TREE 2** | >= 2.2.0 | Maximum Likelihood phylogenetic tree inference | Minh et al. (2020) |
| **pandas / tqdm / Biopython** | Latest stable | Dataframes, progress bars, and FASTA handling | PyPI Community |

---

## 7. Data Integrity

The pipeline **never fabricates data**. Any step that cannot produce a real result from the
input genomes fails explicitly (exit != 0) instead of substituting a plausible-looking value:

- **Genomes.** WGS master accessions (like the example dataset's, `VKHW00000000.1`) carry no
  sequence — an `efetch` on them returns zero bytes.
  [download_genome.py](bin/download_genome.py) resolves the accession to its assembly
  (`GCA_`/`GCF_`) before downloading, then validates contigs, base count (>= 500 kb), and the
  nucleotide alphabet. Without a real genome, the strain does not continue.
- **QC.** CheckM2 actually runs against its database. There are no hardcoded
  completeness/contamination values.
- **Annotation.** Without the Bakta database, `BAKTA_ANNOTATE` fails. There is no substitute
  annotation with invented genes.
- **References.** The sequences in [urease_references.fasta](data/urease_references.fasta) and
  the outgroup in [urec_outgroup.fasta](data/urec_outgroup.fasta) are downloaded from UniProt via
  REST queries (described below). Nothing is synthesised locally — a network failure, a
  truncated response or an empty result set aborts without writing either file. Regenerate with
  `python3 bin/fetch_references.py`.
- **Synteny and phylogeny.** No placeholder HTML and no partial trees. With fewer than 3 UreC
  sequences, that is recorded in `phylogeny_status.txt` as a legitimate outcome — no tree is
  possible — and no tree file is written.

Reference set — query `(taxonomy_id:1239) AND (((protein_name:urease) AND (reviewed:true)) OR
(gene_exact:urtA OR … OR gene_exact:urtE OR gene_exact:uca OR gene_exact:atzF))` against
`https://rest.uniprot.org/uniprotkb/stream`. All entries are Firmicutes: the urease subunits and
accessory proteins (UreA–I) are reviewed (Swiss-Prot) only, while the transporter (UrtA–E), urea
carboxylase (`uca`) and allophanate hydrolase (`atzF`) may be unreviewed (TrEMBL). The set keeps
one sequence per gene and species — reviewed first, then the lowest accession; `--keep-redundant`
turns this off. A sanity floor of 100 sequences applies to the raw response (below that it is
treated as truncated and aborts). The UreC tree outgroup — the reviewed UreC of *Klebsiella
aerogenes* (P18314) and *Proteus mirabilis* (P17086) — goes to `data/urec_outgroup.fasta`, used
only by the tree and kept out of BLASTp. The set is **live**: it tracks UniProt and changes with
every release (on release 2026_03, 960 sequences after deduplication — 109 reviewed and 851
unreviewed — out of 1,074 in the response). Reproducibility does not rest on a fixed accession
list but on the provenance records written next to each FASTA
(`data/urease_references.provenance.txt`, `data/urec_outgroup.provenance.txt` — release, date,
counts and `sha256`). To pin a different slice, use `--query` or `--url`.

---

## 8. Installation & Usage

### Prerequisites

Same as the [Quick-Start](#2-quick-start): Docker Engine >= 24.0 with your user in the `docker`
group, Git, Bash, `curl`, and ~160 GB free for the full Bakta database. Every bioinformatics tool
ships inside the image.

### 8.1 Clone the repository
```bash
git clone https://github.com/waldeyr/bafes-urease.git
cd bafes-urease
```

### 8.2 Make the wrapper script executable
```bash
chmod +x run.sh
```

### 8.3 Quick-Start sequence

The recommended sequence, repeated here in full:

```bash
# Step A — Bootstrap: builds the image, downloads Pfam and the full Bakta database.
#    First run takes hours (image build ~20 min + Bakta full DB ~75 GB). Idempotent:
#    every completed step is skipped on subsequent runs.
#    Upgrading from an install that used the old `bafes_urease` tag? Retag first to
#    avoid a ~20-minute rebuild:
#        docker tag bafes_urease bafes-urease && docker rmi bafes_urease
./run.sh --bootstrap --runtime docker --container-image bafes-urease

# Step B — Build & Pre-flight Check (validates environment & resource reachability)
./run.sh --build --runtime docker --container-image bafes-urease

# Step C — Pipeline Execution with real-time tqdm progress tracking & detailed logging
./run.sh --exec --runtime docker --container-image bafes-urease --verbose
```

### 8.4 Further parameter combinations

All three modes (`--bootstrap`, `--build`, `--exec`) accept the same options, and the defaults
already match the Quick-Start (`--runtime docker`, `--container-image bafes-urease`,
`--bakta-db-type full`). The examples below cover the most common deviations.

**Bootstrap with the reduced database** — when ~160 GB is not available; the light Bakta DB takes
~10 GB and annotates fewer RefSeq/UniRef identifiers:

```bash
./run.sh --bootstrap --bakta-db-type light
```

**Bootstrap reusing a Bakta database already on the host** — avoids re-downloading 75 GB by
pointing at the copy you already have:

```bash
./run.sh --bootstrap --skip-bakta-db --bakta-db /data/bakta/db
./run.sh --exec      --bakta-db /data/bakta/db --verbose
```

**Bootstrap forcing an image rebuild** — after editing `Dockerfile` or `environment.yml`:

```bash
./run.sh --bootstrap --force-rebuild
```

**Check and update drifted databases** — without `--update-db`, bootstrap only warns about the
version drift and keeps the local copy:

```bash
./run.sh --bootstrap --update-db
```

**Chain build and execution propagating the exit code** — `--wait` blocks until each phase ends,
so the `&&` only fires if the previous one exited 0:

```bash
./run.sh --build --wait && ./run.sh --exec --wait --verbose
```

**Debug in the foreground** — no `screen` session, output straight to the terminal:

```bash
./run.sh --exec --no-screen --verbose
```

**Resume an interrupted run** — reuses the Nextflow cache instead of redoing completed phases:

```bash
./run.sh --exec --resume --verbose
```

**Run a subset of strains into a separate directory** — useful to smoke-test the pipeline with 1
or 2 genomes before the full round, without overwriting `results/`:

```bash
head -3 data/accessions.tsv > data/accessions_test.tsv
./run.sh --exec --accessions data/accessions_test.tsv --outdir results_test --verbose
```

**Point every input at alternative paths** — databases and references outside the project
directory (for example, on a shared volume):

```bash
./run.sh --exec \
  --accessions /data/bafes/accessions.tsv \
  --references /data/bafes/urease_references.fasta \
  --pfam       /data/databases/Pfam-A.hmm \
  --bakta-db   /data/databases/bakta/db \
  --outdir     /data/bafes/results \
  --verbose
```

**Run without Docker** — requires a preconfigured Conda/Mamba environment with the tools on the
host:

```bash
./run.sh --exec --runtime local --verbose
```

**Run with Singularity/Apptainer** — on HPC clusters where Docker is unavailable:

```bash
./run.sh --exec --runtime singularity --container-image bafes-urease --verbose
```

**Full unattended round** — reduced database, everything in the foreground and chained, suitable
for launching inside a queue job:

```bash
./run.sh --bootstrap --bakta-db-type light --no-screen \
  && ./run.sh --build --no-screen \
  && ./run.sh --exec --no-screen --verbose
```

### 8.5 Available options

| Option | Effect |
| :--- | :--- |
| `--bootstrap` | Prepares inputs, image, and databases (idempotent) |
| `--build` | Validates environment and resources (pre-flight check) |
| `--exec` | Runs the Nextflow pipeline |
| `--verbose` | Enables detailed reports and logs |
| `--resume` | Resumes an interrupted run from the Nextflow cache |
| `--force-rebuild` | Rebuilds the Docker image even if it already exists |
| `--bakta-db-type full\|light` | Full database (~75 GB, default) or reduced one (~10 GB) |
| `--skip-bakta-db` | Skips the Bakta database download (use if you already have one at `--bakta-db`) |
| `--update-db` | Re-downloads databases whose version drifted from the available one (without it, only warns) |
| `--runtime docker\|singularity\|local` | Execution engine (default: `docker`) |
| `--container-image IMAGE` | Container image tag (default: `bafes-urease`) |
| `--no-screen` | Runs in the foreground, without creating a session |
| `--wait` | Blocks until the session finishes and propagates the real exit code |
| `--accessions PATH` | Path to `accessions.tsv` (default: `data/accessions.tsv`) |
| `--local-genomes PATH` | TSV of local FASTA assemblies (`strain`, `species`, `fasta`), in addition to the accessions |
| `--references PATH` | Path to the references FASTA (default: `data/urease_references.fasta`) |
| `--pfam PATH` | Path to `Pfam-A.hmm` (default: `Pfam-A.hmm`) |
| `--bakta-db PATH` | Path to the Bakta database (default: `db/db`) |
| `--outdir PATH` | Output directory (default: `results`) |
| `-h`, `--help` | Displays the help message |

### 8.6 Database Version Checking

Databases are not only checked for existence — the installed version is checked too. This matters
because the source URLs are **rolling**: `Pfam/current_release` serves different releases over
time, so "the file exists" does not say *which* release sits on disk.

`run.sh` records the downloaded version under `db/.versions/` and compares it against what is
available:

| Database | Local version | Available version |
| :--- | :--- | :--- |
| Pfam | `db/.versions/pfam` (release recorded at download) | `Pfam.version.gz` at EMBL-EBI |
| Bakta | `major.minor` from `version.json` + type in `db/.versions/bakta_type` | `bakta_db list` |
| CheckM2 | DIAMOND db filename (e.g. `uniref100.KO.1.dmnd`) | no queryable index |
| AMRFinderPlus | target of the `amrfinderplus-db/latest` symlink | queried by `amrfinder_update` |
| UniProt references | `uniprot_release` in `data/urease_references.provenance.txt` | `x-uniprot-release` header from the REST API |

**Nothing is re-downloaded while the database exists**, not even on version drift: `--bootstrap`
only warns and keeps the local copy. Re-fetching the full Bakta DB means ~75 GB and hours of
network, so updating is always explicit, via `--update-db`. With no network to query the remote
version, the local copy is kept without error.

The UniProt references come from a **live** query (`fetch_references.py`) whose result set changes
with every release. What identifies the local FASTA is therefore the release stamped into
`data/urease_references.provenance.txt` — alongside the record count and `sha256`. A FASTA
generated before that record has no known release: `--bootstrap` warns and keeps the file.

When updating (Bakta and CheckM2), the current database is moved to `.old` and only discarded once
the new one arrives intact; if the download fails, the working database is restored, not lost.

```bash
# Check versions and update only what drifted
./run.sh --bootstrap --update-db
```

### 8.7 Automatic `screen` Session

Every `run.sh` invocation creates its own **detached** `screen` session and returns the prompt
immediately, so an SSH drop does not kill the run — which matters because bootstrap takes hours
and `--exec` takes longer still.

```
>>> Sessão / Session : bafes-bootstrap (screen)
>>> Reanexar / Attach: screen -r bafes-bootstrap        (detach with Ctrl-A D)
>>> Log              : .logs/bootstrap-20260727-1430.log
>>> Código de saída  : .logs/bootstrap-20260727-1430.exitcode
```

**Mind the exit code.** Because the run continues in the background, the outer invocation returns
`0` even when the pipeline fails. The real code is written to the `.exitcode` file. To chain
commands or automate, use `--wait`, which blocks until completion and propagates the code:

```bash
./run.sh --build --runtime docker --wait && ./run.sh --exec --runtime docker --wait
```

If `screen` is unavailable, `run.sh` falls back to `tmux`; with neither, it runs in the foreground
with a warning. A second invocation of the same phase is refused (exit 8) rather than duplicated.

### 8.8 Common Problem: `permission denied` on the Docker socket

If the command works in a fresh terminal but fails inside an older `screen` session, Docker is not
the cause. A process's supplementary groups are fixed at login by `setgroups()` and are never
re-read from `/etc/group`; a session created **before** `usermod -aG docker` carries the old set,
and every shell inside it inherits that staleness. `run.sh` detects this and re-executes through
`sg docker`, which re-reads `/etc/group` — there is no need to kill the session.

If it still fails, the user genuinely is not in the group. The admin must run
`sudo usermod -aG docker $USER`, and you must open a new SSH session.

### 8.9 Configuration via `.env`

Copy the template and fill it in. `run.sh` loads `.env` automatically and forwards the
variables into the container. `.env` is in `.gitignore` — **never commit your key**.

```bash
cp .env.example .env
```

| Variable | Effect |
| :--- | :--- |
| `NCBI_API_KEY` | Raises the E-utilities rate limit from 3 to 10 req/s. **Recommended**: without it the accession pre-check and the genome downloads are throttled to 3 req/s and are more prone to HTTP 429. Generate one at [account.ncbi.nlm.nih.gov](https://account.ncbi.nlm.nih.gov/) → Account Settings → API Key Management |
| `NCBI_EMAIL` | Contact e-mail required by the NCBI E-utilities policy |
| `BAFES_INSECURE_SSL=1` | Disables TLS verification (only for networks behind a TLS-intercepting proxy) |

### 8.10 Direct Nextflow Execution (Optional)

Requires Nextflow and the tools on the host, or the `docker` profile (one container per process):

```bash
nextflow run main.nf \
   --accessions data/accessions.tsv \
   --bakta_db ./db/db \
   --pfam_hmm ./Pfam-A.hmm \
   --references data/urease_references.fasta \
   --outdir results \
   -profile docker
```

The reports (`nf_report.html`, `nf_trace.txt`, `nf_timeline.html`, `nf_dag.html`) are already
enabled in `nextflow.config` — there is no need to pass the `-with-*` flags.

### 8.11 Generated Output Structure
```
results/
├── 00_preflight/    # Pre-flight validation status
├── 00_genomes/      # Assembly FASTA files downloaded from NCBI
├── 01_qc/           # QC reports (QUAST + CheckM2)
├── 02_bakta/        # Annotations in GFF3, FAA, JSON, GBFF per strain
├── 03_candidates/   # Consolidated candidates per strain (Bakta + BLASTp + HMMER)
├── 04_summary/      # urease_presence_absence.tsv matrix and unified summary
├── 05_synteny/      # Interactive synteny viewer (synteny.html)
└── 06_phylogeny/    # UreC phylogenetic tree (ureC.treefile and alignments)
```

---

## 9. Database & Data Dictionary

### Database Diagram

The pipeline persists tab-delimited tabular files (`.tsv`), structured JSON documents, FASTA
alignments, and HTML viewers on a directory-structured filesystem:

```mermaid
erDiagram
    ACCESSIONS_TSV {
        string strain PK
        string species
        string genbank_accession
    }
    GENOME_FASTA {
        string strain FK
        string fasta_header
        string sequence
    }
    BAKTA_JSON {
        string strain FK
        string feature_id PK
        string gene_symbol
        string product_name
        string ec_number
    }
    CANDIDATES_TSV {
        string strain FK
        string feature_id PK
        string gene_symbol
        boolean layer_bakta
        boolean layer_blast
        boolean layer_hmmer
        string confidence_level
    }
    PRESENCE_ABSENCE_MATRIX {
        string strain PK
        int ureC
        int ureA
        int ureB
        int ureD
        int ureE
        int ureF
        int ureG
        int urtA
        int uc
        int ah
    }

    ACCESSIONS_TSV ||--|| GENOME_FASTA : "specifies download"
    GENOME_FASTA ||--|| BAKTA_JSON : "annotated into"
    BAKTA_JSON ||--|{ CANDIDATES_TSV : "screened into"
    CANDIDATES_TSV }|--|| PRESENCE_ABSENCE_MATRIX : "consolidated into"
```

### Data Dictionary

#### Table: `data/accessions.tsv`
| Field | Type | Description |
| :--- | :--- | :--- |
| `strain` | String (PK) | Short unique identifier of the bacterial strain (e.g. S1, S2) |
| `species` | String | Bacterial species binomial (e.g. Lysinibacillus fusiformis) |
| `genbank_accession` | String | Assembly (`GCA_`/`GCF_`) or WGS master GenBank accession (e.g. VKHW00000000.1) |

#### Table: `--local-genomes` TSV (e.g. `data/local_genomes.tsv`)
| Field | Type | Description |
| :--- | :--- | :--- |
| `strain` | String (PK) | Short unique identifier; must not repeat one from `accessions.tsv` |
| `species` | String | Species binomial, or `Unknown` |
| `fasta` | String | Path to the FASTA assembly, relative to the launch directory |

#### Table: `results/03_candidates/{strain}/{strain}_candidates.tsv`
| Field | Type | Description |
| :--- | :--- | :--- |
| `strain` | String (FK) | Bacterial strain identifier |
| `feature_id` | String (PK) | Identifier of the gene/protein annotated by Bakta |
| `gene` | String | Gene symbol assigned (e.g. ureC, urtA, uc) |
| `product` | String | Functional description of the protein |
| `layer_bakta` | Boolean | Detection flag for Layer A (Bakta) |
| `layer_blast` | Boolean | Detection flag for Layer B (BLASTp) |
| `layer_hmmer` | Boolean | Detection flag for Layer C (Pfam HMMER) |
| `blast_ref_status` | String | Status of the reference behind the best BLASTp hit: `reviewed` (Swiss-Prot), `unreviewed` (TrEMBL) or `unknown` |
| `evidence_points` | Integer | Evidence points (see [Candidate screening score](#candidate-screening-score)) |
| `confidence` | String | Confidence derived from the points: `Alta/High`, `Media/Medium` or `Baixa/Low` |

#### Table: `results/04_summary/urease_presence_absence.tsv`
| Field | Type | Description |
| :--- | :--- | :--- |
| `strain` | String (PK) | Unique strain identifier |
| `ureC` | Integer (0/1) | Presence (1) or absence (0) of the urease $\alpha$ subunit |
| `ureA` | Integer (0/1) | Presence (1) or absence (0) of the urease $\gamma$ subunit |
| `ureB` | Integer (0/1) | Presence (1) or absence (0) of the urease $\beta$ subunit |
| `ureD` | Integer (0/1) | Presence (1) or absence (0) of the accessory protein UreD |
| `ureE` | Integer (0/1) | Presence (1) or absence (0) of the accessory protein UreE |
| `ureF` | Integer (0/1) | Presence (1) or absence (0) of the accessory protein UreF |
| `ureG` | Integer (0/1) | Presence (1) or absence (0) of the accessory GTPase UreG |
| `ureH` | Integer (0/1) | Presence (1) or absence (0) of the accessory protein UreH |
| `urtA` | Integer (0/1) | Presence (1) or absence (0) of the urea transporter UrtA |
| `uc` | Integer (0/1) | Presence (1) or absence (0) of urea carboxylase |
| `ah` | Integer (0/1) | Presence (1) or absence (0) of allophanate hydrolase |

---

## 10. Pipeline Decision Guide (Q&A)

### Which databases, and why?

Each database answers a different question, and candidate screening deliberately uses three
independent sources: no single one (annotation, similarity or domain) is enough to claim a gene is
present. A candidate is **High** confidence only when at least two layers agree.

| Database | Step | Why |
| :--- | :--- | :--- |
| **NCBI GenBank / Assembly** | `DOWNLOAD_GENOME` | Public source of the genomes listed in `accessions.tsv`. WGS master accessions (like the example dataset's, `VKHW00000000.1`) hold no sequence; [download_genome.py](bin/download_genome.py) resolves them to the `GCA_`/`GCF_` assembly before downloading. Local genomes (`--local-genomes`) do not go through here. |
| **CheckM2** (`uniref100.KO.1.dmnd`) | `QC_CHECKM2` | Estimates completeness and contamination with a machine-learning model that does not depend on lineage-specific markers — a good fit for strains with no species assigned, like the example dataset's *Bacillus* sp. |
| **Bakta DB full** | `BAKTA_ANNOTATE` | Standardised bacterial annotation. Its JSON carries gene, product, EC numbers and `db_xrefs`, which make up Layer A. The *full* DB (~75 GB) assigns more RefSeq/UniRef identifiers than *light*. It bundles AMRFinderPlus, an internal Bakta dependency — which is why it shows up in the version check. |
| **UniProtKB** (Swiss-Prot + TrEMBL) | `EXTRACT_CANDIDATES`, Layer B (BLASTp) | Restricted by default to **Firmicutes** (`taxonomy_id:1239`), the example dataset's phylum, which keeps identity high against the 30 % identity and 50 % coverage thresholds. For urease, **reviewed** entries only (Swiss-Prot), with curator-verified function. For UrtA–E, `uca` and `atzF` Swiss-Prot has no Firmicutes entry (and no UrtA–E in any organism), so TrEMBL entries come in, searched by exact gene name. The query is live, with release and `sha256` recorded in [urease_references.provenance.txt](data/urease_references.provenance.txt). |
| **UniProtKB/Swiss-Prot** (outgroup) | `BUILD_PHYLOGENY` | Two Gammaproteobacteria UreC that root the UreC tree (see the species question). They live in [urec_outgroup.fasta](data/urec_outgroup.fasta), apart from the BLASTp references. |
| **Pfam-A** | `EXTRACT_CANDIDATES`, Layer C (hmmscan) | HMM profiles catch remote homology that BLASTp misses and do not depend on Bakta's annotation. Only 10 profiles count as evidence (`UREASE_PFAM_ACCESSIONS` in [extract_urease_candidates.py](bin/extract_urease_candidates.py)). PF01979 (Amidohydro_1) and PF00005 (ABC_tran) are left out because they span whole superfamilies. |

**Known limitation.** Pfam and Bakta come from rolling URLs: "the file exists" does not say which
release is on disk. Reproducibility relies on the versions recorded in `db/.versions/` (see
[Database Version Checking](#86-database-version-checking)).

### What is the execution order, and why?

Nextflow runs as a dataflow: the "order" is the chain of dependencies, and each strain moves on
independently. The linear diagram in the [Pipeline Overview](#3-pipeline-overview) simplifies this.

1. **`PREFLIGHT_CHECK`** — Tests NCBI, Pfam, Bakta and the references first, so a run fails in
   seconds rather than after hours of annotation and hmmscan.
2. **`DOWNLOAD_GENOME` / `STAGE_LOCAL_GENOME`** — Fetches and validates the FASTA (contigs,
   >= 500 kb, alphabet). Without a real genome, the strain stops here.
3. **`QC_QUAST` and `QC_CHECKM2`**, in parallel with **`BAKTA_ANNOTATE`** — QC is informative and
   does **not** filter strains: the genome goes straight to Bakta. Running in parallel keeps QC
   from delaying the longest step.
4. **`BAKTA_ANNOTATE`** before screening — BLASTp and hmmscan run on the `.faa` that Bakta writes.
   That way all three layers refer to the same locus tag, and consensus is computed by identifier.
5. **`EXTRACT_CANDIDATES`** — BLASTp, hmmscan and the three-layer consolidation, per strain. This
   is the pipeline's bottleneck (hmmscan against the whole of Pfam-A).
6. **`MERGE_ALL_STRAINS`**, **`EXTRACT_LOCUS`** and **`BUILD_PHYLOGENY`** after screening — All
   three consume the candidates: the matrix consolidates them, the locus slice anchors on them and
   the phylogeny takes UreC from them, never from the whole proteome. They are independent of each
   other.
7. **`EXTRACT_LOCUS`** before **`ANALYZE_SYNTENY`** — clinker aligns every protein against every
   other. In the example dataset, with the 10 whole genomes (~46,500 CDS), that was ~10⁹ alignments
   and the step blew its 8 h limit (exit 143). With the locus sliced out (±5 flanking genes), the
   step takes seconds — and its cost now grows with the size of the locus, not of the genome.

### Which reference species are used at each step, and why?

| Step | Reference | Why |
| :--- | :--- | :--- |
| **Input genomes** | Those listed in `accessions.tsv` and in the `--local-genomes` TSV. In the example dataset: S1 *Lysinibacillus fusiformis* and S2–S10 *Bacillus* sp. (CBAFES/UnB collection), plus the local genomes SDF0008, SDF0018 and SDF0030 from [local_genomes.tsv](data/local_genomes.tsv). | They are the subject of the analysis, not references. |
| **QC** | None. | QUAST runs without `-r`: the pipeline does not assume the species of the input genomes, so there is no fixed reference genome. CheckM2 uses a universal model. |
| **Annotation (Bakta)** | No species fixed (no `--genus`/`--species`). | Bakta's databases are universal; fixing a species would only change the name in the output. |
| **Layer B (BLASTp), urease** | 109 Swiss-Prot proteins from 21 Firmicutes species, one per gene and species (release 2026_03): UreA 18, UreB 19, UreC 19, UreD 12, UreE 16, UreF 11, UreG 12, UreH 1, UreI 1. Among them *Sporosarcina pasteurii*, *Lysinibacillus sphaericus*, *Bacillus subtilis*, *B. cereus*, *B. velezensis*, *Geobacillus kaustophilus*, *Halalkalibacterium halodurans*, plus *Staphylococcus*, *Streptococcus* and clostridia. | Same phylum as the example dataset, with verified function. *S. pasteurii* is the model Firmicutes urease; in the example, *L. sphaericus* shares S1's genus, and *Bacillus* and *Geobacillus* are close relatives of S2–S10. For genomes from another phylum, change the slice with `--query`. Deduplication took *S. aureus* from 94 entries down to 7 (one per gene). |
| **Layer B (BLASTp), transport and alternative pathway** | 851 TrEMBL Firmicutes proteins, one per gene and species: UrtA 130, UrtB 139, UrtC 136, UrtD 137, UrtE 137 (141 species), `uca` 84 and `atzF` 88. Most represented genera: *Paenibacillus*, *Cohnella*, *Bacillus*, *Metabacillus*, *Geobacillus*. | They are the only Firmicutes carrying these genes in UniProt. `uca` and `atzF` are searched by gene name rather than EC number, because EC 3.5.1.54 also brings the KipI/KipA "allophanate hydrolases", which in *Bacillus* are the 5-oxoprolinase PxpB/PxpC. |
| **Layer C (Pfam)** | No species. | Family profiles built from alignments spanning many taxa. |
| **Synteny** | Input strains only. | The goal is to compare operon organisation across the isolates. A strain with no urease is left off the map, and no locus is invented for it. |
| **UreC phylogeny** | The strains' UreC + 19 reference UreC, one per Firmicutes species (including *S. pasteurii*, *B. subtilis*, *B. cereus*, *B. velezensis*, *G. kaustophilus* and *H. halodurans*) + outgroup: the UreC of *Klebsiella aerogenes* (P18314) and *Proteus mirabilis* (P17086). | IQ-TREE needs >= 4 sequences for bootstrap, and the references place the strains within the phylum. The outgroup are Gammaproteobacteria model ureases with solved structures, clearly outside Firmicutes; IQ-TREE roots the tree on them (`-o`). Can be turned off with `--phylo_include_refs false` and `--phylo_outgroup false`. |

**Known limitations.**

- **Unreviewed references for transport and the alternative pathway.** UrtA–E, `uca` and `atzF`
  come from TrEMBL, which is annotated automatically, like Bakta. For these genes a BLASTp hit
  confirms similarity to proteins annotated the same way, not verified function. That is why such
  a hit is flagged (`blast_ref_status = unreviewed`, outlined "L" badge in the report) and scores
  nothing — see [Candidate screening score](#candidate-screening-score).
- **Composition bias in transport.** *Paenibacillus* accounts for almost half of the UrtA–E
  references, because it is the Firmicutes genus with the most genomes carrying the operon.
- **Few species with reviewed urease.** Deduplication removed the *S. aureus* redundancy but adds
  no species: Swiss-Prot holds reviewed urease from only 21 Firmicutes species.
- **The presence/absence matrix reads Bakta only.** It uses the annotated gene and product; a
  protein found only by BLASTp or Pfam changes the confidence in the candidate table but does not
  mark presence in the matrix.

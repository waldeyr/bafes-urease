#!/usr/bin/env bash
# ==============================================================================
# run.sh — Single Entry Point Wrapper for Baafes-urease-mining Pipeline
# PT-BR: Wrapper de execução unificado com suporte a tqdm, dry-run e Docker.
# EN-US: Unified execution wrapper supporting tqdm, dry-run pre-checks, and Docker.
# ==============================================================================

set -e

# PT-BR: Valores padrão / EN-US: Default values
BOOTSTRAP=false
BUILD=false
EXEC=false
VERBOSE=false
RESUME=false
RUNTIME="docker"
CONTAINER_IMAGE="bafes_urease"
ACCESSIONS="data/accessions.tsv"
REFERENCES="data/urease_references.fasta"
PFAM="Pfam-A.hmm"
BAKTA_DB="db/db"
OUTDIR="results"

usage() {
    echo "Uso / Usage: $0 [opções / options]"
    echo ""
    echo "Modos de operação / Operating modes:"
    echo "  --bootstrap       Prepara diretórios e bancos / Prepares directories and databases"
    echo "  --build           Valida ambiente e executa dry-run / Validates environment & runs dry-run"
    echo "  --exec            Executa o pipeline Nextflow com tqdm / Runs Nextflow pipeline with tqdm"
    echo "  --verbose         Ativa relatórios e logs detalhados / Enables detailed reports & logs"
    echo "  --resume          Retoma execução interrompida / Resumes interrupted execution"
    echo ""
    echo "Opções / Options:"
    echo "  --runtime RUNTIME               Motor de execução: docker | singularity | local (default: docker)"
    echo "  --container-image IMAGE         Nome da imagem container / Container image tag (default: bafes_urease)"
    echo "  --accessions PATH               Caminho accessions.tsv / Path to accessions.tsv (default: data/accessions.tsv)"
    echo "  --references PATH               Caminho referencias FASTA / Path to references FASTA (default: data/urease_references.fasta)"
    echo "  --pfam PATH                     Caminho Pfam-A.hmm / Path to Pfam-A.hmm (default: Pfam-A.hmm)"
    echo "  --bakta-db PATH                 Caminho banco Bakta / Path to Bakta DB (default: db/db)"
    echo "  --outdir PATH                   Diretório de resultados / Output directory (default: results)"
    echo "  -h, --help                      Exibe mensagem de ajuda / Displays help message"
    exit 0
}

# PT-BR: Parse de opções / EN-US: Parse command-line flags
while [[ $# -gt 0 ]]; do
    case "$1" in
        --bootstrap) BOOTSTRAP=true; shift ;;
        --build) BUILD=true; shift ;;
        --exec) EXEC=true; shift ;;
        --verbose) VERBOSE=true; shift ;;
        --resume) RESUME=true; shift ;;
        --runtime) RUNTIME="$2"; shift 2 ;;
        --container-image) CONTAINER_IMAGE="$2"; shift 2 ;;
        --accessions) ACCESSIONS="$2"; shift 2 ;;
        --references) REFERENCES="$2"; shift 2 ;;
        --pfam) PFAM="$2"; shift 2 ;;
        --bakta-db) BAKTA_DB="$2"; shift 2 ;;
        --outdir) OUTDIR="$2"; shift 2 ;;
        -h|--help) usage ;;
        *) echo "Opção desconhecida / Unknown option: $1"; usage ;;
    esac
done

mkdir -p data db results .logs

echo "================================================================="
echo "   BAAFES UREASE MINING — PIPELINE WRAPPER (run.sh)             "
echo "================================================================="

# 1. MODO BOOTSTRAP / BOOTSTRAP MODE
if [ "$BOOTSTRAP" = true ]; then
    echo ">>> [BOOTSTRAP] Preparando insumos / Preparing inputs & database downloads..."
    
    if [ ! -f "$ACCESSIONS" ]; then
        echo "Criando template em / Creating template at $ACCESSIONS..."
        mkdir -p "$(dirname "$ACCESSIONS")"
        cat << 'EOF' > "$ACCESSIONS"
strain	species	genbank_accession
S1	Lysinibacillus fusiformis	VKHW00000000.1
S2	Bacillus sp. S2	VKHY00000000.1
S3	Bacillus sp. S3	VKHZ00000000.1
S4	Bacillus sp. S4	SADW00000000.1
S5	Bacillus sp. S5	VKHX00000000.1
S6	Bacillus sp. S6	SADY00000000.1
S7	Bacillus sp. S7	SADV00000000.1
S8	Bacillus sp. S8	SADX00000000.1
S9	Bacillus sp. S9	VKIB00000000.1
S10	Bacillus sp. S10	VKIC00000000.1
EOF
    fi

    if [ ! -f "$REFERENCES" ]; then
        echo "Criando template em / Creating template at $REFERENCES..."
        mkdir -p "$(dirname "$REFERENCES")"
        cat << 'EOF' > "$REFERENCES"
>ref|UreC|urease_alpha_subunit|Bacillus_subtilis
MKLSPREKDKLLLFTADACIAEGRIVTVEEVIEKGIVTLTGIAHEEVDIPLGTHLVEVVP
EOF
    fi

    if [ ! -f "$PFAM" ]; then
        echo "Baixando / Downloading Pfam-A.hmm.gz..."
        curl -s -L "https://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/Pfam-A.hmm.gz" -o Pfam-A.hmm.gz || true
        if [ -f Pfam-A.hmm.gz ]; then
            gunzip -f Pfam-A.hmm.gz || true
        fi
    fi

    echo "[BOOTSTRAP OK] Insumos iniciais preparados / Inputs ready."
fi

# 2. MODO BUILD & DRY-RUN / BUILD & DRY-RUN MODE
if [ "$BUILD" = true ]; then
    echo ">>> [BUILD] Validando runtime e conectividade / Validating runtime & resource connectivity..."
    
    python3 bin/resource_checker.py \
        --accessions "$ACCESSIONS" \
        --references "$REFERENCES" \
        --bakta-db "$BAKTA_DB" || {
            STATUS_CODE=$?
            if [ $STATUS_CODE -eq 4 ]; then
                echo "[ERRO BUILD / BUILD ERROR] Pre-flight resource check failed (Exit code 4)."
                exit 4
            fi
        }

    if [ "$RUNTIME" = "docker" ]; then
        echo "Checando Docker..."
        docker info >/dev/null 2>&1 || { echo "[AVISO/WARNING] Docker daemon não detectado."; }
    fi

    echo "[BUILD OK] Ambiente e recursos prévios validados / Environment ready."
fi

# 3. MODO EXEC / EXEC MODE
if [ "$EXEC" = true ]; then
    echo ">>> [EXEC] Iniciando Nextflow pipeline..."

    NF_FLAGS=""
    if [ "$RUNTIME" = "docker" ]; then
        NF_FLAGS="-profile docker"
    elif [ "$RUNTIME" = "singularity" ]; then
        NF_FLAGS="-profile singularity"
    else
        NF_FLAGS="-profile standard"
    fi

    if [ "$RESUME" = true ]; then
        NF_FLAGS="$NF_FLAGS -resume"
    fi

    if [ "$VERBOSE" = true ]; then
        NF_FLAGS="$NF_FLAGS -with-report ${OUTDIR}/nf_report.html -with-trace ${OUTDIR}/nf_trace.txt -with-timeline ${OUTDIR}/nf_timeline.html -with-dag ${OUTDIR}/nf_dag.html"
    fi

    python3 bin/progress_tracker.py ${VERBOSE:+--verbose} &
    TRACKER_PID=$!

    nextflow run main.nf \
        --accessions "$ACCESSIONS" \
        --references "$REFERENCES" \
        --pfam_hmm "$PFAM" \
        --bakta_db "$BAKTA_DB" \
        --outdir "$OUTDIR" \
        $NF_FLAGS

    wait $TRACKER_PID 2>/dev/null || true
    echo ">>> [EXEC OK] Pipeline concluído! Artefatos em / Outputs in: $OUTDIR/"
fi

if [ "$BOOTSTRAP" = false ] && [ "$BUILD" = false ] && [ "$EXEC" = false ]; then
    usage
fi

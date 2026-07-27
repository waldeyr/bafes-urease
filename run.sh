#!/usr/bin/env bash
# ==============================================================================
# run.sh — Single Entry Point Wrapper for Baafes-urease-mining Pipeline
# PT-BR: Wrapper de execução unificado. Todo o pipeline roda DENTRO do container.
# EN-US: Unified execution wrapper. The whole pipeline runs INSIDE the container.
# ==============================================================================

set -e

# PT-BR: Valores padrão / EN-US: Default values
BOOTSTRAP=false
BUILD=false
EXEC=false
VERBOSE=false
RESUME=false
FORCE_REBUILD=false
SKIP_BAKTA_DB=false
NO_SCREEN=false
WAIT_FOR_SESSION=false
RUNTIME="docker"
CONTAINER_IMAGE="bafes-urease"
BAKTA_DB_TYPE="full"
ACCESSIONS="data/accessions.tsv"
REFERENCES="data/urease_references.fasta"
REFERENCE_ACCESSIONS="data/urease_reference_accessions.tsv"
PFAM="Pfam-A.hmm"
BAKTA_DB="db/db"
CHECKM2_DB_DIR="db/checkm2"
OUTDIR="results"

# PT-BR: Carrega o .env, se existir. Ver .env.example / EN-US: Loads .env if present. See .env.example
if [ -f .env ]; then
    set -a
    # shellcheck disable=SC1091
    . ./.env
    set +a
fi

# PT-BR: Espaço mínimo para o banco Bakta (KiB) / EN-US: Minimum free disk for Bakta DB (KiB)
BAKTA_FULL_MIN_KB=$((160 * 1024 * 1024))   # ~160 GB
BAKTA_LIGHT_MIN_KB=$((10 * 1024 * 1024))   # ~10 GB

usage() {
    echo "Uso / Usage: $0 [opções / options]"
    echo ""
    echo "Modos de operação / Operating modes:"
    echo "  --bootstrap       Prepara insumos, imagem e bancos / Prepares inputs, image and databases"
    echo "  --build           Valida ambiente e recursos / Validates environment & resources"
    echo "  --exec            Executa o pipeline Nextflow / Runs the Nextflow pipeline"
    echo "  --verbose         Ativa relatórios e logs detalhados / Enables detailed reports & logs"
    echo "  --resume          Retoma execução interrompida / Resumes interrupted execution"
    echo ""
    echo "Opções / Options:"
    echo "  --runtime RUNTIME               Motor: docker | singularity | local (default: docker)"
    echo "  --container-image IMAGE         Imagem container / Container image tag (default: bafes-urease)"
    echo "  --no-screen                     Executa em primeiro plano, sem criar sessão screen"
    echo "                                  Runs in the foreground, without creating a screen session"
    echo "  --wait                          Aguarda a sessão terminar e propaga o código de saída real"
    echo "                                  Waits for the session to finish and propagates the real exit code"
    echo "  --force-rebuild                 Reconstrói a imagem mesmo se existir / Rebuilds image even if present"
    echo "  --bakta-db-type TYPE            Banco Bakta: full | light (default: full)"
    echo "  --skip-bakta-db                 Não baixa o banco Bakta / Skips Bakta DB download"
    echo "  --accessions PATH               Caminho accessions.tsv / Path to accessions.tsv (default: data/accessions.tsv)"
    echo "  --references PATH               Caminho referencias FASTA / Path to references FASTA (default: data/urease_references.fasta)"
    echo "  --pfam PATH                     Caminho Pfam-A.hmm / Path to Pfam-A.hmm (default: Pfam-A.hmm)"
    echo "  --bakta-db PATH                 Caminho banco Bakta / Path to Bakta DB (default: db/db)"
    echo "  --outdir PATH                   Diretório de resultados / Output directory (default: results)"
    echo "  -h, --help                      Exibe mensagem de ajuda / Displays help message"
    echo ""
    echo "Variáveis de ambiente / Environment variables:"
    echo "  NCBI_API_KEY      Chave da API do NCBI; eleva o limite de 3 para 10 req/s"
    echo "                    NCBI API key; raises rate limit from 3 to 10 req/s"
    echo "  NCBI_EMAIL        E-mail de contato exigido pela política do NCBI E-utilities"
    echo "                    Contact e-mail required by NCBI E-utilities policy"
    exit 0
}

# PT-BR: Guarda os argumentos originais ANTES do parsing — o laço abaixo os consome
#        com 'shift', e tanto a reexecução via 'sg docker' quanto a criação da sessão
#        screen precisam repassar a linha de comando intacta ao processo interno.
# EN-US: Save the original arguments BEFORE parsing — the loop below consumes them with
#        'shift', and both the 'sg docker' re-execution and the screen session need to
#        forward the untouched command line to the inner process.
ORIGINAL_ARGS=("$@")

# PT-BR: Parse de opções / EN-US: Parse command-line flags
while [[ $# -gt 0 ]]; do
    case "$1" in
        --bootstrap) BOOTSTRAP=true; shift ;;
        --build) BUILD=true; shift ;;
        --exec) EXEC=true; shift ;;
        --verbose) VERBOSE=true; shift ;;
        --resume) RESUME=true; shift ;;
        --force-rebuild) FORCE_REBUILD=true; shift ;;
        --skip-bakta-db) SKIP_BAKTA_DB=true; shift ;;
        --no-screen) NO_SCREEN=true; shift ;;
        --wait) WAIT_FOR_SESSION=true; shift ;;
        --runtime) RUNTIME="$2"; shift 2 ;;
        --container-image) CONTAINER_IMAGE="$2"; shift 2 ;;
        --bakta-db-type) BAKTA_DB_TYPE="$2"; shift 2 ;;
        --accessions) ACCESSIONS="$2"; shift 2 ;;
        --references) REFERENCES="$2"; shift 2 ;;
        --pfam) PFAM="$2"; shift 2 ;;
        --bakta-db) BAKTA_DB="$2"; shift 2 ;;
        --outdir) OUTDIR="$2"; shift 2 ;;
        -h|--help) usage ;;
        *) echo "Opção desconhecida / Unknown option: $1"; usage ;;
    esac
done

if [ "$BOOTSTRAP" = false ] && [ "$BUILD" = false ] && [ "$EXEC" = false ]; then
    usage
fi

case "$RUNTIME" in
    docker|singularity|local) ;;
    *) echo "[ERRO/ERROR] --runtime inválido / invalid: $RUNTIME (docker | singularity | local)"; exit 2 ;;
esac

case "$BAKTA_DB_TYPE" in
    full|light) ;;
    *) echo "[ERRO/ERROR] --bakta-db-type inválido / invalid: $BAKTA_DB_TYPE (full | light)"; exit 2 ;;
esac

# ==============================================================================
# PT-BR: Credencial de grupo defasada.
#        Os grupos suplementares de um processo são fixados no login por setgroups()
#        e nunca são relidos de /etc/group. Uma sessão screen/tmux iniciada ANTES do
#        'usermod -aG docker' carrega o conjunto antigo, e todo shell criado dentro
#        dela herda a defasagem — inclusive um screen novo criado a partir dela.
#        O 'sg' relê /etc/group, então reexecutar por ele corrige a sessão velha.
# EN-US: Stale group credential.
#        A process's supplementary groups are fixed at login by setgroups() and are
#        never re-read from /etc/group. A screen/tmux session started BEFORE
#        'usermod -aG docker' carries the old set, and every shell created inside it
#        inherits that staleness — including a new screen spawned from it.
#        'sg' re-reads /etc/group, so re-executing through it repairs the old session.
# ==============================================================================
if [ "$RUNTIME" = "docker" ] && [ -z "${BAFES_SG_RETRIED:-}" ] \
   && ! id -nG | tr ' ' '\n' | grep -qx docker \
   && getent group docker 2>/dev/null | cut -d: -f4 | tr ',' '\n' | grep -qx "$USER"; then
    echo ">>> Credencial de grupo defasada nesta sessão / Stale group credential in this session."
    echo "    Esta sessão é anterior ao 'usermod -aG docker'; reexecutando via 'sg docker'."
    echo "    This session predates 'usermod -aG docker'; re-executing through 'sg docker'."
    export BAFES_SG_RETRIED=1
    # PT-BR: Sob 'sg docker' o GID primário vira o do grupo docker. Preservamos o GID
    #        original para que os artefatos não saiam com grupo 'docker'.
    # EN-US: Under 'sg docker' the primary GID becomes docker's. We preserve the original
    #        GID so artifacts do not end up group-owned by 'docker'.
    export BAFES_HOST_GID="$(id -g)"
    exec sg docker -c "$(printf '%q ' "$0" "${ORIGINAL_ARGS[@]}")"
fi

# ==============================================================================
# PT-BR: Sessão screen — protege execuções longas contra queda do SSH.
#        O --bootstrap leva horas (imagem + Bakta full + CheckM2) e o --exec mais ainda.
# EN-US: Screen session — protects long runs against SSH disconnection.
#        --bootstrap takes hours (image + full Bakta + CheckM2) and --exec even longer.
# ==============================================================================
session_alive() {
    case "${1:-}" in
        tmux)   tmux has-session -t "$SESSION" 2>/dev/null ;;
        *)      screen -ls 2>/dev/null | grep -qE "[0-9]+\.${SESSION}[[:space:]]" ;;
    esac
}

if [ -z "${BAFES_IN_SCREEN:-}" ] && [ "$NO_SCREEN" = false ]; then
    PHASES=""
    [ "$BOOTSTRAP" = true ] && PHASES="${PHASES}bootstrap-"
    [ "$BUILD" = true ]     && PHASES="${PHASES}build-"
    [ "$EXEC" = true ]      && PHASES="${PHASES}exec-"
    PHASES="${PHASES%-}"

    SESSION="bafes-${PHASES}"
    mkdir -p .logs
    STAMP="$(date +%Y%m%d-%H%M%S)"
    LOGFILE=".logs/${PHASES}-${STAMP}.log"
    EXITFILE=".logs/${PHASES}-${STAMP}.exitcode"

    BACKEND=""
    if command -v screen >/dev/null 2>&1; then
        BACKEND="screen"
    elif command -v tmux >/dev/null 2>&1; then
        BACKEND="tmux"
    fi

    if [ -n "$BACKEND" ]; then
        if session_alive "$BACKEND"; then
            echo "[ERRO/ERROR] Já existe uma sessão ativa / A session is already running: $SESSION"
            echo "             Reanexar / Attach:  $BACKEND $([ "$BACKEND" = tmux ] && echo "attach -t" || echo "-r") $SESSION"
            echo "             Encerrar / Kill:    $([ "$BACKEND" = tmux ] && echo "tmux kill-session -t $SESSION" || echo "screen -S $SESSION -X quit")"
            exit 8
        fi

        # PT-BR: tee mantém a saída visível ao reanexar E grava o log; PIPESTATUS
        #        recupera o código real do run.sh, não o do tee.
        # EN-US: tee keeps output visible on attach AND writes the log; PIPESTATUS
        #        recovers run.sh's real exit code, not tee's.
        WRAPPER='BAFES_IN_SCREEN=1 "$@" 2>&1 | tee '"$(printf '%q' "$LOGFILE")"'
code=${PIPESTATUS[0]}
echo "$code" > '"$(printf '%q' "$EXITFILE")"'
exit "$code"'

        if [ "$BACKEND" = "tmux" ]; then
            tmux new-session -d -s "$SESSION" bash -c "$WRAPPER" bash "$0" "${ORIGINAL_ARGS[@]}"
            ATTACH="tmux attach -t $SESSION"
        else
            screen -dmS "$SESSION" bash -c "$WRAPPER" bash "$0" "${ORIGINAL_ARGS[@]}"
            ATTACH="screen -r $SESSION"
        fi

        echo "================================================================="
        echo "   BAAFES UREASE MINING — SESSÃO INICIADA / SESSION STARTED      "
        echo "================================================================="
        echo ">>> Sessão / Session : $SESSION ($BACKEND)"
        echo ">>> Reanexar / Attach: $ATTACH        (solte com / detach with Ctrl-A D)"
        echo ">>> Log              : $LOGFILE"
        echo ">>> Código de saída  : $EXITFILE"

        if [ "$WAIT_FOR_SESSION" = true ]; then
            echo ">>> Aguardando a sessão terminar / Waiting for the session to finish..."
            while session_alive "$BACKEND"; do sleep 5; done
            REAL_CODE="$(cat "$EXITFILE" 2>/dev/null || echo 1)"
            echo ">>> Sessão concluída / Session finished (exit=$REAL_CODE)."
            exit "$REAL_CODE"
        fi

        echo ""
        echo "[ATENÇÃO/NOTE] Este comando retorna 0 porque a execução seguiu em segundo plano."
        echo "               O código de saída REAL fica em $EXITFILE ao final."
        echo "               This command returns 0 because the run continues in the background."
        echo "               The REAL exit code lands in $EXITFILE when it completes."
        echo "               Use --wait para bloquear e propagar o código / to block and propagate it."
        exit 0
    fi

    echo "[AVISO/WARNING] Nem 'screen' nem 'tmux' encontrados / Neither 'screen' nor 'tmux' found."
    echo "                Executando em primeiro plano; a queda do SSH interrompe a execução."
    echo "                Running in the foreground; an SSH drop will kill the run."
fi

mkdir -p data db results .logs .nfhome

echo "================================================================="
echo "   BAAFES UREASE MINING — PIPELINE WRAPPER (run.sh)             "
echo "================================================================="

# ==============================================================================
# PT-BR: Camada de container — tudo roda dentro da imagem quando RUNTIME=docker.
# EN-US: Container layer — everything runs inside the image when RUNTIME=docker.
# ==============================================================================

check_docker() {
    if ! command -v docker >/dev/null 2>&1; then
        echo "[ERRO/ERROR] Docker não encontrado no PATH / Docker not found on PATH."
        echo "             Instale o Docker ou use --runtime local / Install Docker or use --runtime local."
        exit 3
    fi

    local err
    if err=$(docker info 2>&1 >/dev/null); then
        return 0
    fi

    echo "[ERRO/ERROR] Não foi possível falar com o daemon Docker / Cannot reach the Docker daemon."
    if echo "$err" | grep -qi "permission denied"; then
        # PT-BR: A defasagem de credencial já foi tratada antes; se chegamos aqui, o
        #        usuário realmente não pertence ao grupo docker.
        # EN-US: The stale-credential case was already handled earlier; reaching here
        #        means the user genuinely does not belong to the docker group.
        echo "             Causa provável / Likely cause: usuário sem acesso ao socket / user lacks socket access."
        echo "             O admin precisa executar / The admin must run:"
        echo "                 sudo usermod -aG docker $USER"
        echo "             Depois, abra uma sessão SSH nova / Then open a new SSH session."
    else
        echo "             Causa provável / Likely cause: daemon parado / daemon not running."
        echo "             Solução / Fix: sudo systemctl start docker"
    fi
    echo "             Detalhe / Detail: $err"
    exit 3
}

# PT-BR: Executa um comando dentro do container (ou no host se RUNTIME=local).
# EN-US: Runs a command inside the container (or on the host when RUNTIME=local).
in_container() {
    if [ "$RUNTIME" = "local" ]; then
        "$@"
        return $?
    fi

    # PT-BR: BAFES_HOST_GID preserva o GID original quando reexecutamos via 'sg docker',
    #        que troca o GID primário do processo. Sem isso os artefatos em results/,
    #        db/ e work/ sairiam com grupo 'docker'.
    # EN-US: BAFES_HOST_GID preserves the original GID when we re-execute through
    #        'sg docker', which swaps the process's primary GID. Without it, artifacts in
    #        results/, db/ and work/ would end up group-owned by 'docker'.
    docker run --rm \
        -v "$PWD":/workspace \
        -w /workspace \
        -u "$(id -u):${BAFES_HOST_GID:-$(id -g)}" \
        -e HOME=/workspace/.nfhome \
        -e NXF_HOME=/workspace/.nfhome \
        -e NCBI_API_KEY="${NCBI_API_KEY:-}" \
        -e NCBI_EMAIL="${NCBI_EMAIL:-}" \
        -e BAFES_INSECURE_SSL="${BAFES_INSECURE_SSL:-}" \
        "$CONTAINER_IMAGE" "$@"
}

# PT-BR: Espaço livre em KiB no diretório atual / EN-US: Free space in KiB on the current directory
free_space_kb() {
    df -Pk . | awk 'NR==2 {print $4}'
}

human_gb() {
    awk -v kb="$1" 'BEGIN { printf "%.1f", kb / 1024 / 1024 }'
}

if [ "$RUNTIME" = "docker" ]; then
    check_docker
elif [ "$RUNTIME" = "singularity" ]; then
    echo "[AVISO/WARNING] RUNTIME=singularity ainda não implementa a camada de container."
    echo "                RUNTIME=singularity does not implement the container layer yet."
    echo "                Comandos rodarão no host / Commands will run on the host."
    RUNTIME="local"
fi

# ==============================================================================
# 1. MODO BOOTSTRAP / BOOTSTRAP MODE
# ==============================================================================
if [ "$BOOTSTRAP" = true ]; then
    echo ">>> [BOOTSTRAP] Preparando insumos / Preparing inputs & database downloads..."

    # --- 1.1 Templates de entrada / Input templates -------------------------
    if [ ! -f "$ACCESSIONS" ]; then
        echo "  [1/5] Criando template em / Creating template at $ACCESSIONS..."
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
    else
        echo "  [1/5] Accessions já presente / already present: $ACCESSIONS"
    fi

    # PT-BR: As referências são BAIXADAS do UniProt a partir de accessions fixadas.
    #        Nunca são escritas localmente: um stub inventado geraria evidência falsa
    #        de BLASTp e contaminaria toda a triagem de candidatos.
    # EN-US: References are DOWNLOADED from UniProt using pinned accessions. They are
    #        never written locally: a made-up stub would generate false BLASTp evidence
    #        and contaminate the whole candidate screening.
    if [ ! -s "$REFERENCES" ]; then
        echo "        Baixando referências curadas do UniProt / Downloading curated references from UniProt..."
        mkdir -p "$(dirname "$REFERENCES")"
        if ! in_container python3 bin/fetch_references.py \
            --accessions "$REFERENCE_ACCESSIONS" \
            --out "$REFERENCES"; then
            echo "  [ERRO/ERROR] Falha ao obter as referências / Failed to obtain references."
            echo "               Nenhuma sequência será inventada / No sequence will be invented."
            exit 7
        fi
    else
        echo "        Referências já presentes / already present: $REFERENCES"
    fi

    # --- 1.2 Imagem Docker / Docker image -----------------------------------
    if [ "$RUNTIME" = "docker" ]; then
        if [ "$FORCE_REBUILD" = true ]; then
            echo "  [2/5] Reconstruindo imagem / Rebuilding image '$CONTAINER_IMAGE' (--force-rebuild)..."
            docker build --pull -t "$CONTAINER_IMAGE" .
        elif docker image inspect "$CONTAINER_IMAGE" >/dev/null 2>&1; then
            echo "  [2/5] Imagem já existe / Image already present: $CONTAINER_IMAGE"
            echo "        Use --force-rebuild para reconstruir / to rebuild."
        else
            echo "  [2/5] Construindo imagem / Building image '$CONTAINER_IMAGE' (pode levar ~20 min)..."
            docker build -t "$CONTAINER_IMAGE" .
        fi
    else
        echo "  [2/5] RUNTIME=local: build de imagem ignorado / image build skipped."
    fi

    # --- 1.3 Pfam-A.hmm ------------------------------------------------------
    if [ -s "$PFAM" ]; then
        echo "  [3/5] Pfam já presente / already present: $PFAM"
    else
        echo "  [3/5] Baixando / Downloading Pfam-A.hmm.gz -> $PFAM ..."
        PFAM_DIR="$(dirname "$PFAM")"
        [ "$PFAM_DIR" != "." ] && mkdir -p "$PFAM_DIR"
        PFAM_GZ="${PFAM}.gz"

        if ! curl -fL -# --retry 3 --retry-delay 5 \
            "https://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/Pfam-A.hmm.gz" \
            -o "$PFAM_GZ"; then
            echo "  [ERRO/ERROR] Falha ao baixar o Pfam / Failed to download Pfam."
            rm -f "$PFAM_GZ"
            exit 5
        fi

        if ! gunzip -t "$PFAM_GZ" 2>/dev/null; then
            echo "  [ERRO/ERROR] Arquivo Pfam corrompido / Corrupted Pfam archive: $PFAM_GZ"
            rm -f "$PFAM_GZ"
            exit 5
        fi

        gunzip -f "$PFAM_GZ"
        echo "        [OK] Pfam pronto / ready: $PFAM"
    fi

    # PT-BR: O hmmscan exige o Pfam indexado pelo hmmpress. Sem os .h3* ele falha
    #        em toda estirpe — indexar uma vez aqui evita repetir por processo.
    # EN-US: hmmscan requires Pfam pressed by hmmpress. Without the .h3* files it fails
    #        on every strain — pressing once here avoids repeating it per process.
    if [ -f "${PFAM}.h3i" ]; then
        echo "        Pfam já indexado / already pressed: ${PFAM}.h3i"
    else
        echo "        Indexando Pfam com hmmpress / Pressing Pfam with hmmpress..."
        if ! in_container hmmpress -f "$PFAM"; then
            echo "  [ERRO/ERROR] hmmpress falhou / failed."
            exit 5
        fi
    fi

    # --- 1.4 Banco do Bakta / Bakta database --------------------------------
    if [ "$SKIP_BAKTA_DB" = true ]; then
        echo "  [4/5] Download do banco Bakta ignorado / Bakta DB download skipped (--skip-bakta-db)."
    elif [ -f "$BAKTA_DB/version.json" ]; then
        echo "  [4/5] Banco Bakta já presente / already present: $BAKTA_DB"
    else
        if [ "$BAKTA_DB_TYPE" = "full" ]; then
            MIN_KB=$BAKTA_FULL_MIN_KB
        else
            MIN_KB=$BAKTA_LIGHT_MIN_KB
        fi
        AVAIL_KB=$(free_space_kb)

        echo "  [4/5] Banco Bakta '$BAKTA_DB_TYPE': necessário / required ~$(human_gb "$MIN_KB") GB; disponível / available $(human_gb "$AVAIL_KB") GB."
        if [ "$AVAIL_KB" -lt "$MIN_KB" ]; then
            echo "  [ERRO/ERROR] Espaço em disco insuficiente / Insufficient disk space."
            echo "               Libere espaço, use --bakta-db-type light ou --skip-bakta-db."
            echo "               Free up space, use --bakta-db-type light or --skip-bakta-db."
            exit 6
        fi

        echo "        [AVISO/WARNING] O download do banco '$BAKTA_DB_TYPE' leva horas / takes hours."
        BAKTA_OUT="$(dirname "$BAKTA_DB")"
        mkdir -p "$BAKTA_OUT"

        if ! in_container bakta_db download --output "$BAKTA_OUT" --type "$BAKTA_DB_TYPE"; then
            echo "  [ERRO/ERROR] Falha ao baixar o banco do Bakta / Failed to download Bakta DB."
            exit 6
        fi

        echo "        Atualizando AMRFinderPlus DB / Updating AMRFinderPlus DB..."
        if ! in_container amrfinder_update --force_update --database "$BAKTA_DB/amrfinderplus-db"; then
            echo "        [AVISO/WARNING] amrfinder_update falhou / failed; a anotação segue sem AMRFinderPlus."
        fi
        echo "        [OK] Banco Bakta pronto / ready: $BAKTA_DB"
    fi

    # --- 1.4b Banco do CheckM2 / CheckM2 database ---------------------------
    # PT-BR: O QC_CHECKM2 executa o CheckM2 de verdade e exige este banco (~3 GB).
    #        Antes o processo gravava completude/contaminação fixas no código.
    # EN-US: QC_CHECKM2 actually runs CheckM2 and requires this database (~3 GB).
    #        The process used to write hardcoded completeness/contamination values.
    if [ -n "$(find "$CHECKM2_DB_DIR" -name '*.dmnd' 2>/dev/null | head -1)" ]; then
        echo "        Banco CheckM2 já presente / already present: $CHECKM2_DB_DIR"
    else
        echo "        Baixando banco do CheckM2 (~3 GB) / Downloading CheckM2 database (~3 GB)..."
        mkdir -p "$CHECKM2_DB_DIR"
        if ! in_container checkm2 database --download --path "$CHECKM2_DB_DIR"; then
            echo "  [ERRO/ERROR] Falha ao baixar o banco do CheckM2 / Failed to download the CheckM2 database."
            exit 6
        fi
    fi

    # --- 1.5 Nextflow no host (best-effort) / Host Nextflow (best-effort) ----
    if command -v nextflow >/dev/null 2>&1; then
        echo "  [5/5] Nextflow já disponível no host / already available on host."
    elif [ -x "./bin/nextflow" ]; then
        echo "  [5/5] Nextflow já instalado em / already installed at ./bin/nextflow"
    elif command -v java >/dev/null 2>&1; then
        echo "  [5/5] Instalando Nextflow no host / Installing Nextflow on host (./bin/nextflow)..."
        # PT-BR: O instalador pode retornar != 0 mesmo criando o launcher; o critério
        #        de sucesso é o arquivo existir e ser executável, não o exit code.
        # EN-US: The installer may return != 0 even after creating the launcher; the
        #        success criterion is the file existing and being executable, not the exit code.
        (cd bin && curl -s https://get.nextflow.io | bash) >/dev/null 2>&1 || true
        if [ -x "./bin/nextflow" ]; then
            echo "        [OK] Nextflow instalado / installed: ./bin/nextflow"
        else
            echo "        [AVISO/WARNING] Instalação do Nextflow no host falhou / host install failed."
            echo "        Não é bloqueante: o pipeline usa o Nextflow do container."
            echo "        Not blocking: the pipeline uses the container's Nextflow."
        fi
    else
        echo "  [5/5] [AVISO/WARNING] Java não encontrado; Nextflow não será instalado no host."
        echo "        Java not found; skipping host Nextflow install."
        echo "        Não é bloqueante: o pipeline usa o Nextflow do container."
        echo "        Not blocking: the pipeline uses the container's Nextflow."
    fi

    echo "[BOOTSTRAP OK] Insumos iniciais preparados / Inputs ready."
fi

# ==============================================================================
# 2. MODO BUILD / BUILD MODE
# ==============================================================================
if [ "$BUILD" = true ]; then
    echo ">>> [BUILD] Validando runtime e conectividade / Validating runtime & resource connectivity..."

    if [ "$RUNTIME" = "docker" ] && ! docker image inspect "$CONTAINER_IMAGE" >/dev/null 2>&1; then
        echo "[ERRO BUILD / BUILD ERROR] Imagem não encontrada / Image not found: $CONTAINER_IMAGE"
        echo "                           Rode primeiro / Run first: $0 --bootstrap"
        exit 3
    fi

    CHECKER_FLAGS=""
    [ "$RUNTIME" != "local" ] && CHECKER_FLAGS="--in-container"

    set +e
    in_container python3 bin/resource_checker.py \
        --accessions "$ACCESSIONS" \
        --references "$REFERENCES" \
        --bakta-db "$BAKTA_DB" \
        --pfam "$PFAM" \
        $CHECKER_FLAGS
    STATUS_CODE=$?
    set -e

    if [ $STATUS_CODE -ne 0 ]; then
        echo "[ERRO BUILD / BUILD ERROR] Pre-flight resource check failed (Exit code $STATUS_CODE)."
        exit $STATUS_CODE
    fi

    echo "[BUILD OK] Ambiente e recursos prévios validados / Environment ready."
fi

# ==============================================================================
# 3. MODO EXEC / EXEC MODE
# ==============================================================================
if [ "$EXEC" = true ]; then
    echo ">>> [EXEC] Iniciando Nextflow pipeline..."

    # PT-BR: Dentro do container usamos o executor local; o container JÁ é o ambiente.
    # EN-US: Inside the container we use the local executor; the container IS the environment.
    NF_FLAGS="-profile standard"

    if [ "$RESUME" = true ]; then
        NF_FLAGS="$NF_FLAGS -resume"
    fi

    if [ "$VERBOSE" = true ]; then
        # PT-BR: report/trace/timeline/dag já vêm habilitados via nextflow.config.
        # EN-US: report/trace/timeline/dag are already enabled through nextflow.config.
        NF_FLAGS="$NF_FLAGS -ansi-log false"
    fi

    # PT-BR: Barra de progresso roda no host; ausência de python3 não bloqueia.
    # EN-US: Progress bar runs on the host; a missing python3 is not blocking.
    TRACKER_PID=""
    if command -v python3 >/dev/null 2>&1; then
        TRACKER_ARGS=""
        [ "$VERBOSE" = true ] && TRACKER_ARGS="--verbose"
        python3 bin/progress_tracker.py $TRACKER_ARGS &
        TRACKER_PID=$!
    fi

    # PT-BR: O layout do banco do CheckM2 varia por versão; localiza o .dmnd real.
    # EN-US: The CheckM2 database layout varies by version; locate the actual .dmnd.
    CHECKM2_DB="$(find "$CHECKM2_DB_DIR" -name '*.dmnd' 2>/dev/null | head -1)"
    if [ -z "$CHECKM2_DB" ]; then
        echo "[ERRO/ERROR] Banco do CheckM2 não encontrado em / CheckM2 database not found in: $CHECKM2_DB_DIR"
        echo "             Rode primeiro / Run first: $0 --bootstrap"
        exit 6
    fi

    set +e
    in_container nextflow run main.nf \
        --accessions "$ACCESSIONS" \
        --references "$REFERENCES" \
        --pfam_hmm "$PFAM" \
        --bakta_db "$BAKTA_DB" \
        --checkm2_db "$CHECKM2_DB" \
        --outdir "$OUTDIR" \
        $NF_FLAGS
    NF_STATUS=$?
    set -e

    [ -n "$TRACKER_PID" ] && { wait "$TRACKER_PID" 2>/dev/null || true; }

    if [ $NF_STATUS -ne 0 ]; then
        echo ">>> [EXEC FALHOU / EXEC FAILED] Nextflow retornou / returned $NF_STATUS."
        exit $NF_STATUS
    fi

    echo ">>> [EXEC OK] Pipeline concluído! Artefatos em / Outputs in: $OUTDIR/"
fi

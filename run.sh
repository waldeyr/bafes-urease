#!/usr/bin/env bash
# ==============================================================================
# run.sh — Single Entry Point Wrapper for BAFES-urease-mining Pipeline
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
RESUME_SESSION=""
FORCE_REBUILD=false
SKIP_BAKTA_DB=false
UPDATE_DB=false
NO_SCREEN=false
WAIT_FOR_SESSION=false
RUNTIME="docker"
CONTAINER_IMAGE="bafes-urease"
BAKTA_DB_TYPE="full"
ACCESSIONS="data/accessions.tsv"
REFERENCES="data/urease_references.fasta"
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
    echo "  --resume [ID]     Retoma execução interrompida; sem ID retoma a ÚLTIMA sessão do"
    echo "                    .nextflow/history (um 'nextflow -preview' também conta e traz"
    echo "                    cache vazio) / Resumes interrupted execution; with no ID it"
    echo "                    resumes the LAST session in .nextflow/history (a"
    echo "                    'nextflow -preview' counts too, and brings an empty cache)"
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
    echo "  --update-db                     Rebaixa bancos cuja versão divergiu da disponível"
    echo "                                  Re-downloads databases whose version drifted from the available one"
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
        # PT-BR: --resume aceita um ID de sessão opcional. Sem ele o Nextflow retoma a
        #        ÚLTIMA sessão do .nextflow/history — e um `nextflow run -preview` no mesmo
        #        diretório também entra nesse histórico, com cache vazio. Quando isso
        #        acontece, o `-resume` seguinte não reaproveita nada e o pipeline inteiro
        #        (Bakta incluso) refaz do zero. Passe o ID explícito para retomar a sessão
        #        que você realmente quer: ./run.sh --exec --resume <uuid>
        # EN-US: --resume takes an optional session ID. Without it Nextflow resumes the LAST
        #        session in .nextflow/history — and a `nextflow run -preview` in the same
        #        directory also lands in that history, with an empty cache. When that
        #        happens the next `-resume` reuses nothing and the whole pipeline (Bakta
        #        included) reruns from scratch. Pass the explicit ID to resume the session
        #        you actually want: ./run.sh --exec --resume <uuid>
        --resume)
            RESUME=true
            shift
            if [[ $# -gt 0 && "$1" != --* ]]; then
                RESUME_SESSION="$1"
                shift
            fi
            ;;
        --force-rebuild) FORCE_REBUILD=true; shift ;;
        --skip-bakta-db) SKIP_BAKTA_DB=true; shift ;;
        --update-db) UPDATE_DB=true; shift ;;
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
        echo "   BAFES UREASE MINING — SESSÃO INICIADA / SESSION STARTED      "
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
echo "   BAFES UREASE MINING — PIPELINE WRAPPER (run.sh)             "
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

# PT-BR: Diretório do projeto no HOST, canonizado. O container monta este mesmo caminho
#        NO MESMO lugar (mount identidade), então todo caminho absoluto — o work dir do
#        Nextflow, o ${projectDir}, os bancos — significa a mesma coisa dentro e fora.
#        'pwd -P' resolve symlinks por dois motivos: o 'docker run -v' precisa de um
#        caminho real do host, e o Nextflow deriva o launchDir do getcwd() da JVM, que já
#        é o caminho físico — canonizar aqui faz os dois lados coincidirem exatamente.
# EN-US: The project directory on the HOST, canonicalized. The container mounts this very
#        path at the SAME location (identity mount), so every absolute path — Nextflow's
#        work dir, ${projectDir}, the databases — means the same thing inside and outside.
#        'pwd -P' resolves symlinks for two reasons: 'docker run -v' needs a real host
#        path, and Nextflow derives launchDir from the JVM's getcwd(), which is already the
#        physical path — canonicalizing here makes both sides agree exactly.
BAFES_WORKDIR="$(pwd -P)"

# PT-BR: O mount identidade impõe duas condições ao caminho do projeto. Falhamos aqui,
#        antes de qualquer execução de horas, e não no meio de uma tarefa do Nextflow.
# EN-US: The identity mount imposes two conditions on the project path. We fail here,
#        before any hours-long run, instead of in the middle of a Nextflow task.
check_workdir_mount() {
    # PT-BR: ':' é o separador de campos do '-v'; espaço/tab quebram os scripts gerados
    #        pelo Nextflow, onde as interpolações de ${projectDir} e absPath() não são
    #        citadas. Hoje o /workspace esconde o caminho do host do container.
    # EN-US: ':' is '-v''s field separator; whitespace breaks the scripts Nextflow
    #        generates, where the ${projectDir} and absPath() interpolations are unquoted.
    #        Today /workspace hides the host path from the container.
    case "$BAFES_WORKDIR" in
        *:*|*[[:space:]]*)
            echo "[ERRO/ERROR] Caminho do projeto contém ':' ou espaço / path has ':' or whitespace:" >&2
            echo "             $BAFES_WORKDIR" >&2
            echo "             Mova o projeto para um caminho sem esses caracteres." >&2
            echo "             Move the project to a path without those characters." >&2
            exit 9
            ;;
    esac

    # PT-BR: Montar sobre um diretório que existe na imagem esconderia o conteúdo dela —
    #        /opt/conda guarda TODOS os ambientes das ferramentas.
    # EN-US: Mounting over a directory that exists in the image would shadow its content —
    #        /opt/conda holds ALL the tool environments.
    case "$BAFES_WORKDIR" in
        /|/bin|/bin/*|/boot|/boot/*|/dev|/dev/*|/etc|/etc/*|/lib|/lib/*|/lib64|/lib64/*|\
        /opt|/opt/conda|/opt/conda/*|/proc|/proc/*|/root|/root/*|/sbin|/sbin/*|\
        /sys|/sys/*|/usr|/usr/*|/var|/var/*)
            echo "[ERRO/ERROR] Projeto em caminho que existe dentro da imagem:" >&2
            echo "             Project sits at a path that exists inside the image:" >&2
            echo "             $BAFES_WORKDIR" >&2
            echo "             O mount esconderia o conteúdo da imagem / it would shadow the image." >&2
            exit 9
            ;;
    esac
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
    # PT-BR: Mount identidade — o caminho DENTRO do container é idêntico ao do host.
    #        NXF_HOME é o que realmente importa: é lá que o Nextflow grava plugins e
    #        cache, e precisa cair dentro do mount para persistir. O HOME abaixo é inerte
    #        na prática — o entrypoint do micromamba o reescreve para /home/mambauser
    #        sempre que passamos '-u' (verificado na lisina) — mas é derivado da MESMA
    #        variável para não divergir se a imagem base mudar esse comportamento.
    # EN-US: Identity mount — the path INSIDE the container is identical to the host's.
    #        NXF_HOME is the one that matters: Nextflow writes its plugins and cache there
    #        and it must land inside the mount to persist. The HOME below is inert in
    #        practice — micromamba's entrypoint rewrites it to /home/mambauser whenever we
    #        pass '-u' (verified on lisina) — but it derives from the SAME variable so the
    #        two cannot drift apart if the base image ever stops doing that.
    docker run --rm \
        -v "$BAFES_WORKDIR":"$BAFES_WORKDIR" \
        -w "$BAFES_WORKDIR" \
        -u "$(id -u):${BAFES_HOST_GID:-$(id -g)}" \
        -e HOME="$BAFES_WORKDIR/.nfhome" \
        -e NXF_HOME="$BAFES_WORKDIR/.nfhome" \
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

# ==============================================================================
# PT-BR: Registro de versões dos bancos.
#        Os bancos são grandes (Pfam ~1,5 GB, Bakta full ~75 GB, CheckM2 ~3 GB) e as
#        URLs de origem são ROLANTES: o mesmo endereço serve releases diferentes ao
#        longo do tempo. Só "o arquivo existe" não diz QUAL versão está no disco.
#        Guardamos a versão baixada aqui para comparar sem rebaixar nada.
# EN-US: Database version registry.
#        The databases are large (Pfam ~1.5 GB, Bakta full ~75 GB, CheckM2 ~3 GB) and the
#        source URLs are ROLLING: the same address serves different releases over time.
#        "The file exists" alone does not say WHICH version sits on disk. We record the
#        downloaded version here so it can be compared without re-downloading anything.
# ==============================================================================
DB_VERSION_DIR="db/.versions"

# PT-BR: Lê a versão registrada de um banco / EN-US: Reads a database's recorded version
db_version_get() {
    local file="$DB_VERSION_DIR/$1"
    [ -s "$file" ] || return 1
    head -n1 "$file"
}

# PT-BR: Registra a versão instalada / EN-US: Records the installed version
db_version_set() {
    mkdir -p "$DB_VERSION_DIR"
    printf '%s\n' "$2" > "$DB_VERSION_DIR/$1"
}

# PT-BR: Compara versão local x remota e define DB_DECISION=skip|update.
#        Padrão conservador: se o banco já existe, NADA é rebaixado — nem quando a
#        versão divergiu. Uma divergência só vira download com --update-db explícito,
#        porque rebaixar o Bakta full são ~75 GB e horas de rede.
#        $1 rótulo, $2 versão local (vazio = desconhecida), $3 versão remota (vazio = indisponível)
# EN-US: Compares local vs remote version and sets DB_DECISION=skip|update.
#        Conservative default: if the database is already there, NOTHING is re-downloaded —
#        not even on version drift. A drift only becomes a download with an explicit
#        --update-db, because re-fetching the full Bakta DB means ~75 GB and hours of network.
#        $1 label, $2 local version (empty = unknown), $3 remote version (empty = unavailable)
db_decide() {
    local label="$1" local_ver="$2" remote_ver="$3"
    DB_DECISION="skip"

    if [ -n "$local_ver" ] && [ -n "$remote_ver" ] && [ "$local_ver" = "$remote_ver" ]; then
        echo "        Versão local igual à disponível / Local version matches the available one: $local_ver"
        echo "        Nada a baixar / Nothing to download."
        return 0
    fi

    if [ -z "$local_ver" ]; then
        echo "        Versão local desconhecida / Unknown local version (baixado antes do registro de versões)."
        [ -n "$remote_ver" ] && echo "        Versão disponível / Available version: $remote_ver"
    elif [ -z "$remote_ver" ]; then
        echo "        Versão remota indisponível / Remote version unavailable; mantendo a local / keeping local: $local_ver"
        return 0
    else
        echo "        [AVISO/WARNING] Versão divergente / Version drift: local=$local_ver, disponível/available=$remote_ver"
    fi

    if [ "$UPDATE_DB" = true ]; then
        echo "        --update-db: rebaixando / re-downloading $label."
        DB_DECISION="update"
    else
        echo "        Mantendo a cópia local / Keeping the local copy. Use --update-db para atualizar / to update."
    fi
}

# PT-BR: Release do Pfam publicada no EMBL-EBI. Best-effort: sem rede, ecoa vazio.
#        O Pfam-A.hmm não carrega a release no cabeçalho — ela só existe neste arquivo.
# EN-US: Pfam release published at EMBL-EBI. Best-effort: with no network it echoes empty.
#        Pfam-A.hmm does not carry the release in its header — it only lives in this file.
pfam_remote_version() {
    curl -fsSL --connect-timeout 10 --max-time 30 \
        "https://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/Pfam.version.gz" 2>/dev/null \
        | gunzip -c 2>/dev/null \
        | awk -F':[[:space:]]*' '/^Pfam release/ { gsub(/[[:space:]]/, "", $2); print $2; exit }'
}

# PT-BR: Versão do banco Bakta local, lida do version.json que o próprio bakta_db escreve.
# EN-US: Local Bakta DB version, read from the version.json that bakta_db itself writes.
bakta_local_version() {
    local vj="$1/version.json" major minor
    [ -s "$vj" ] || return 1
    major="$(tr -d ' \n' < "$vj" | sed -n 's/.*"major":\([0-9][0-9]*\).*/\1/p')"
    minor="$(tr -d ' \n' < "$vj" | sed -n 's/.*"minor":\([0-9][0-9]*\).*/\1/p')"
    [ -n "$major" ] && [ -n "$minor" ] || return 1
    printf '%s.%s\n' "$major" "$minor"
}

# PT-BR: Maior versão do banco compatível com o bakta instalado, via 'bakta_db list'.
#        O parsing é tolerante: só considera colunas no formato N.N.
# EN-US: Highest DB version compatible with the installed bakta, via 'bakta_db list'.
#        Parsing is lenient: it only considers columns shaped like N.N.
bakta_remote_version() {
    in_container bakta_db list 2>/dev/null \
        | awk '$1 ~ /^[0-9]+\.[0-9]+$/ { print $1 }' \
        | sort -t. -k1,1n -k2,2n \
        | tail -n1
}

# PT-BR: Release do UniProt registrada na procedência do FASTA local. O conjunto de
#        referências vem de uma consulta VIVA — muda a cada release —, então a release
#        gravada pelo fetch_references.py é a única "versão" que o FASTA carrega.
# EN-US: UniProt release recorded in the local FASTA's provenance. The reference set comes
#        from a LIVE query — it changes with every release — so the release stamped by
#        fetch_references.py is the only "version" the FASTA carries.
references_local_release() {
    local prov="${1%.*}.provenance.txt" release
    [ -s "$prov" ] || return 1
    release="$(awk -F'\t' '$1 == "uniprot_release" { print $2; exit }' "$prov")"
    [ -n "$release" ] && [ "$release" != "unknown" ] || return 1
    printf '%s\n' "$release"
}

# PT-BR: Release corrente do UniProt, do cabeçalho x-uniprot-release da API REST.
#        Best-effort: sem rede, ecoa vazio e a cópia local é mantida.
# EN-US: Current UniProt release, from the REST API's x-uniprot-release header.
#        Best-effort: with no network it echoes empty and the local copy is kept.
uniprot_remote_release() {
    curl -fsSI --connect-timeout 10 --max-time 30 \
        "https://rest.uniprot.org/uniprotkb/P41022.fasta" 2>/dev/null \
        | awk 'tolower($1) == "x-uniprot-release:" { gsub(/\r/, "", $2); print $2; exit }'
}

if [ "$RUNTIME" = "docker" ]; then
    check_docker
    check_workdir_mount
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

    # PT-BR: As referências são BAIXADAS do UniProt por consulta REST (ureases de Firmicutes
    #        revisadas). Nunca são escritas localmente: um stub inventado geraria evidência
    #        falsa de BLASTp e contaminaria toda a triagem de candidatos.
    # EN-US: References are DOWNLOADED from UniProt via a REST query (reviewed Firmicutes
    #        ureases). They are never written locally: a made-up stub would generate false
    #        BLASTp evidence and contaminate the whole candidate screening.
    # PT-BR: A "versão" das referências é a release do UniProt gravada na procedência
    #        pelo fetch_references.py. Sem procedência (FASTA de uma versão anterior do
    #        script) a release é desconhecida e a cópia local é mantida.
    # EN-US: The references' "version" is the UniProt release stamped into the provenance
    #        by fetch_references.py. With no provenance (a FASTA from an earlier version of
    #        the script) the release is unknown and the local copy is kept.
    REF_DOWNLOAD=true
    if [ -s "$REFERENCES" ]; then
        REF_LOCAL_REL="$(references_local_release "$REFERENCES" || true)"
        echo "        Referências já presentes / already present: $REFERENCES${REF_LOCAL_REL:+ (UniProt $REF_LOCAL_REL)}"
        db_decide "referências UniProt / UniProt references" \
            "$REF_LOCAL_REL" "$(uniprot_remote_release || true)"
        [ "$DB_DECISION" = "update" ] || REF_DOWNLOAD=false
    fi

    if [ "$REF_DOWNLOAD" = true ]; then
        echo "        Baixando referências de urease do UniProt / Downloading urease references from UniProt..."
        mkdir -p "$(dirname "$REFERENCES")"
        if ! in_container python3 bin/fetch_references.py --out "$REFERENCES"; then
            echo "  [ERRO/ERROR] Falha ao obter as referências / Failed to obtain references."
            echo "               Nenhuma sequência será inventada / No sequence will be invented."
            exit 7
        fi
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
    PFAM_DOWNLOAD=true
    PFAM_REMOTE_VER="$(pfam_remote_version || true)"

    if [ -s "$PFAM" ]; then
        PFAM_LOCAL_VER="$(db_version_get pfam || true)"
        echo "  [3/5] Pfam já presente / already present: $PFAM${PFAM_LOCAL_VER:+ (release $PFAM_LOCAL_VER)}"
        db_decide "Pfam" "$PFAM_LOCAL_VER" "$PFAM_REMOTE_VER"
        [ "$DB_DECISION" = "update" ] || PFAM_DOWNLOAD=false
    fi

    if [ "$PFAM_DOWNLOAD" = true ]; then
        echo "  [3/5] Baixando / Downloading Pfam-A.hmm.gz -> $PFAM ${PFAM_REMOTE_VER:+(release $PFAM_REMOTE_VER)}..."
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
        # PT-BR: Só registramos a versão DEPOIS do arquivo estar íntegro no disco.
        # EN-US: We only record the version AFTER the file is intact on disk.
        db_version_set pfam "${PFAM_REMOTE_VER:-desconhecida/unknown}"
        echo "        [OK] Pfam pronto / ready: $PFAM${PFAM_REMOTE_VER:+ (release $PFAM_REMOTE_VER)}"
    fi

    # PT-BR: O hmmscan exige o Pfam indexado pelo hmmpress. Sem os .h3* ele falha
    #        em toda estirpe — indexar uma vez aqui evita repetir por processo.
    #        Um Pfam recém-baixado invalida os índices antigos: reindexa sempre.
    # EN-US: hmmscan requires Pfam pressed by hmmpress. Without the .h3* files it fails
    #        on every strain — pressing once here avoids repeating it per process.
    #        A freshly downloaded Pfam invalidates the old indices: always re-press.
    if [ "$PFAM_DOWNLOAD" = false ] && [ -f "${PFAM}.h3i" ]; then
        echo "        Pfam já indexado / already pressed: ${PFAM}.h3i"
    else
        echo "        Indexando Pfam com hmmpress / Pressing Pfam with hmmpress..."
        if ! in_container hmmpress -f "$PFAM"; then
            echo "  [ERRO/ERROR] hmmpress falhou / failed."
            exit 5
        fi
    fi

    # --- 1.4 Banco do Bakta / Bakta database --------------------------------
    # PT-BR: A versão local vem do version.json escrito pelo próprio bakta_db; o tipo
    #        (full/light) não está lá, então guardamos à parte. Trocar de light para full
    #        conta como divergência de versão — os bancos não são intercambiáveis.
    # EN-US: The local version comes from the version.json bakta_db itself writes; the type
    #        (full/light) is not in it, so we record it separately. Switching light to full
    #        counts as version drift — the databases are not interchangeable.
    BAKTA_DOWNLOAD=true
    BAKTA_LOCAL_VER=""

    if [ "$SKIP_BAKTA_DB" = true ]; then
        echo "  [4/5] Download do banco Bakta ignorado / Bakta DB download skipped (--skip-bakta-db)."
        BAKTA_DOWNLOAD=false
    elif [ -f "$BAKTA_DB/version.json" ]; then
        BAKTA_LOCAL_VER="$(bakta_local_version "$BAKTA_DB" || true)"
        BAKTA_LOCAL_TYPE="$(db_version_get bakta_type || true)"
        echo "  [4/5] Banco Bakta já presente / already present: $BAKTA_DB${BAKTA_LOCAL_VER:+ (v$BAKTA_LOCAL_VER)}"

        BAKTA_LOCAL_ID=""
        [ -n "$BAKTA_LOCAL_VER" ] && BAKTA_LOCAL_ID="${BAKTA_LOCAL_VER}/${BAKTA_LOCAL_TYPE:-$BAKTA_DB_TYPE}"
        BAKTA_REMOTE_VER="$(bakta_remote_version || true)"
        BAKTA_REMOTE_ID=""
        [ -n "$BAKTA_REMOTE_VER" ] && BAKTA_REMOTE_ID="${BAKTA_REMOTE_VER}/${BAKTA_DB_TYPE}"

        db_decide "Bakta DB" "$BAKTA_LOCAL_ID" "$BAKTA_REMOTE_ID"
        [ "$DB_DECISION" = "update" ] || BAKTA_DOWNLOAD=false
    fi

    if [ "$BAKTA_DOWNLOAD" = true ]; then
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

        # PT-BR: O bakta_db extrai por cima; um banco antigo deixaria arquivos órfãos de
        #        outra release misturados aos novos. Movemos para .old em vez de apagar:
        #        se o download de horas falhar, o banco que já funcionava é restaurado em
        #        vez de perdido. A checagem de espaço acima já garante folga para os dois.
        # EN-US: bakta_db extracts in place; an old database would leave orphan files from
        #        another release mixed with the new ones. We move it to .old instead of
        #        deleting: if the hours-long download fails, the working database is restored
        #        rather than lost. The disk check above already guarantees room for both.
        BAKTA_OLD=""
        if [ -n "$BAKTA_LOCAL_VER" ] && [ -d "$BAKTA_DB" ]; then
            BAKTA_OLD="${BAKTA_DB}.old"
            rm -rf "${BAKTA_OLD:?}"
            echo "        Guardando o banco antigo / Setting the old database aside: $BAKTA_DB (v$BAKTA_LOCAL_VER) -> $BAKTA_OLD"
            mv "$BAKTA_DB" "$BAKTA_OLD"
        fi

        echo "        [AVISO/WARNING] O download do banco '$BAKTA_DB_TYPE' leva horas / takes hours."
        BAKTA_OUT="$(dirname "$BAKTA_DB")"
        mkdir -p "$BAKTA_OUT"

        if ! in_container bakta_db download --output "$BAKTA_OUT" --type "$BAKTA_DB_TYPE"; then
            echo "  [ERRO/ERROR] Falha ao baixar o banco do Bakta / Failed to download Bakta DB."
            if [ -n "$BAKTA_OLD" ] && [ -d "$BAKTA_OLD" ]; then
                echo "        Restaurando o banco anterior / Restoring the previous database: $BAKTA_OLD -> $BAKTA_DB"
                rm -rf "${BAKTA_DB:?}"
                mv "$BAKTA_OLD" "$BAKTA_DB"
            fi
            exit 6
        fi

        [ -n "$BAKTA_OLD" ] && rm -rf "${BAKTA_OLD:?}"
        db_version_set bakta_type "$BAKTA_DB_TYPE"
        echo "        [OK] Banco Bakta pronto / ready: $BAKTA_DB ($(bakta_local_version "$BAKTA_DB" || echo '?')/$BAKTA_DB_TYPE)"
    fi

    # PT-BR: O amrfinder_update baixa a release MAIS RECENTE sempre que a chamamos, o que
    #        mudaria a anotação entre execuções sem o usuário pedir. Então só chamamos
    #        quando o banco ainda não existe, quando o Bakta acabou de ser baixado, ou
    #        com --update-db. A versão instalada é o alvo do symlink 'latest'.
    # EN-US: amrfinder_update downloads the LATEST release every time we call it, which
    #        would change annotation between runs without the user asking. So we only call
    #        it when the database is missing, when Bakta was just downloaded, or with
    #        --update-db. The installed version is the target of the 'latest' symlink.
    if [ "$SKIP_BAKTA_DB" = false ] && [ -d "$BAKTA_DB" ]; then
        AMR_DB="$BAKTA_DB/amrfinderplus-db"
        AMR_LOCAL_VER=""
        if [ -e "$AMR_DB/latest" ]; then
            AMR_LINK="$(readlink "$AMR_DB/latest" 2>/dev/null || true)"
            AMR_LOCAL_VER="$(basename "${AMR_LINK:-latest}")"
        fi

        if [ -n "$AMR_LOCAL_VER" ] && [ "$UPDATE_DB" = false ] && [ "$BAKTA_DOWNLOAD" = false ]; then
            echo "        AMRFinderPlus DB já presente / already present: $AMR_LOCAL_VER — nada a baixar / nothing to download."
        else
            AMRFINDER_FLAGS=""
            [ "$UPDATE_DB" = true ] && AMRFINDER_FLAGS="--force_update"
            echo "        Atualizando AMRFinderPlus DB / Updating AMRFinderPlus DB..."
            if ! in_container amrfinder_update $AMRFINDER_FLAGS --database "$AMR_DB"; then
                echo "        [AVISO/WARNING] amrfinder_update falhou / failed; a anotação segue sem AMRFinderPlus."
            fi
        fi
    fi

    # --- 1.4b Banco do CheckM2 / CheckM2 database ---------------------------
    # PT-BR: O QC_CHECKM2 executa o CheckM2 de verdade e exige este banco (~3 GB).
    #        Antes o processo gravava completude/contaminação fixas no código.
    #        O CheckM2 não publica índice de versões consultável: a identidade do banco é
    #        o nome do DIAMOND db (ex.: uniref100.KO.1.dmnd). Presente => não baixa.
    # EN-US: QC_CHECKM2 actually runs CheckM2 and requires this database (~3 GB).
    #        The process used to write hardcoded completeness/contamination values.
    #        CheckM2 publishes no queryable version index: the database's identity is the
    #        DIAMOND db filename (e.g. uniref100.KO.1.dmnd). Present => no download.
    CHECKM2_DMND="$(find "$CHECKM2_DB_DIR" -name '*.dmnd' 2>/dev/null | head -1)"
    if [ -n "$CHECKM2_DMND" ] && [ "$UPDATE_DB" = false ]; then
        echo "        Banco CheckM2 já presente / already present: $CHECKM2_DB_DIR"
        echo "        Versão / Version: $(db_version_get checkm2 || basename "$CHECKM2_DMND")"
    else
        # PT-BR: Mesma proteção do Bakta: o banco atual só é descartado depois que o
        #        novo chega inteiro. Ver o bloco 1.4.
        # EN-US: Same protection as Bakta: the current database is only discarded after
        #        the new one arrives intact. See block 1.4.
        CHECKM2_OLD=""
        if [ -n "$CHECKM2_DMND" ]; then
            CHECKM2_OLD="${CHECKM2_DB_DIR}.old"
            rm -rf "${CHECKM2_OLD:?}"
            echo "        --update-db: rebaixando o banco do CheckM2 / re-downloading the CheckM2 database."
            mv "$CHECKM2_DB_DIR" "$CHECKM2_OLD"
        fi
        echo "        Baixando banco do CheckM2 (~3 GB) / Downloading CheckM2 database (~3 GB)..."
        mkdir -p "$CHECKM2_DB_DIR"
        if ! in_container checkm2 database --download --path "$CHECKM2_DB_DIR" --no_write_json_db; then
            echo "  [ERRO/ERROR] Falha ao baixar o banco do CheckM2 / Failed to download the CheckM2 database."
            if [ -n "$CHECKM2_OLD" ] && [ -d "$CHECKM2_OLD" ]; then
                echo "        Restaurando o banco anterior / Restoring the previous database: $CHECKM2_OLD -> $CHECKM2_DB_DIR"
                rm -rf "${CHECKM2_DB_DIR:?}"
                mv "$CHECKM2_OLD" "$CHECKM2_DB_DIR"
            fi
            exit 6
        fi
        [ -n "$CHECKM2_OLD" ] && rm -rf "${CHECKM2_OLD:?}"
        CHECKM2_DMND="$(find "$CHECKM2_DB_DIR" -name '*.dmnd' 2>/dev/null | head -1)"
        [ -n "$CHECKM2_DMND" ] && db_version_set checkm2 "$(basename "$CHECKM2_DMND")"
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
    echo "               Tudo pronto para executar o pipeline / Everything is ready to run the pipeline:"
    echo "               $0 --exec --runtime $RUNTIME --container-image $CONTAINER_IMAGE"
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
        if [ -n "$RESUME_SESSION" ]; then
            NF_FLAGS="$NF_FLAGS -resume $RESUME_SESSION"
            echo ">>> [EXEC] Retomando a sessão / Resuming session: $RESUME_SESSION"
        else
            NF_FLAGS="$NF_FLAGS -resume"
            echo ">>> [EXEC] Retomando a ÚLTIMA sessão do histórico / Resuming the LAST session in the history."
            echo "           Se um 'nextflow -preview' rodou depois da execução que você quer retomar,"
            echo "           ele é a última sessão e o cache virá vazio. Confira com:"
            echo "           / If a 'nextflow -preview' ran after the execution you want to resume, it IS"
            echo "           the last session and the cache will come up empty. Check with:"
            echo "               tail -5 .nextflow/history"
        fi
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

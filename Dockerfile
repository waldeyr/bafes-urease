# check=skip=FromPlatformFlagConstDisallowed
# PT-BR: Bioconda distribui a maioria das ferramentas apenas para linux-64;
#        fixar a plataforma evita builds quebrados em hosts ARM (Apple Silicon).
# EN-US: Bioconda ships most tools for linux-64 only; pinning the platform avoids
#        broken builds on ARM hosts (Apple Silicon).
FROM --platform=linux/amd64 mambaorg/micromamba:1.5.8-jammy

LABEL maintainer="CBAFES UnB Urease Mining Team"
LABEL description="Container para pipeline Nextflow de mineracao de urease em genomas bacterianos"
LABEL version="0.3"

USER root

# Instalar dependencias basicas do sistema
RUN apt-get update && apt-get install -y --no-install-recommends \
    procps \
    curl \
    wget \
    git \
    ca-certificates \
    gzip \
    tar \
    unzip \
    && rm -rf /var/lib/apt/lists/*

# PT-BR: Ambiente principal / EN-US: Main environment
COPY environment.yml /tmp/environment.yml
RUN micromamba install -y -n base -f /tmp/environment.yml && \
    micromamba clean --all --yes

# PT-BR: QUAST e CheckM2 em ambientes isolados — suas restrições de versão de
#        BLAST e de Python são incompatíveis com o ambiente principal.
# EN-US: QUAST and CheckM2 in isolated environments — their BLAST and Python
#        version constraints are incompatible with the main environment.
COPY env-quast.yml /tmp/env-quast.yml
COPY env-checkm2.yml /tmp/env-checkm2.yml
RUN micromamba create -y -n quast -f /tmp/env-quast.yml && \
    micromamba create -y -n checkm2 -f /tmp/env-checkm2.yml && \
    micromamba clean --all --yes

# PT-BR: O env principal vem primeiro no PATH; os isolados são anexados ao final
#        para que seus binários auxiliares (ex.: blastn do QUAST) não sombreiem
#        as versões do ambiente principal.
# EN-US: The main env comes first on PATH; isolated envs are appended last so
#        their auxiliary binaries (e.g. QUAST's own blastn) do not shadow the
#        main environment's versions.
ENV PATH="/opt/conda/bin:${PATH}:/opt/conda/envs/quast/bin:/opt/conda/envs/checkm2/bin"
ENV PYTHONUNBUFFERED=1

# PT-BR: Sem diretório de trabalho embutido. O run.sh monta o projeto do host no MESMO
#        caminho absoluto dentro do container e passa '-w' explicitamente. Um /workspace
#        vazio aqui só serviria para mascarar um '-w' errado com um diretório mudo.
# EN-US: No baked-in working directory. run.sh mounts the host project at the SAME
#        absolute path inside the container and passes '-w' explicitly. An empty
#        /workspace here would only mask a wrong '-w' behind a silent empty directory.
WORKDIR /

# Ponto de entrada padrão
CMD ["/bin/bash"]

FROM mambaorg/micromamba:1.5.8-jammy

LABEL maintainer="CBafes UnB Urease Mining Team"
LABEL description="Container para pipeline Nextflow de mineracao de urease em genomas bacterianos"
LABEL version="0.2"

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

# Copiar environment.yml e instalar ambiente no micromamba
COPY environment.yml /tmp/environment.yml
RUN micromamba install -y -n base -f /tmp/environment.yml && \
    micromamba clean --all --yes

# Configurar ambiente
ENV PATH="/opt/conda/bin:${PATH}"
ENV PYTHONUNBUFFERED=1

WORKDIR /workspace

# Ponto de entrada padrão
CMD ["/bin/bash"]

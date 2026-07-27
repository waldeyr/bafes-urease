# Baafes-urease-mining

[![Language: Português](https://img.shields.io/badge/Idioma-Portugu%C3%AAs--BR-blue)](README.md)
[![Language: English](https://img.shields.io/badge/Language-English--US-green)](README.en.md)

> 🌐 **Alternar Idioma / Language Switch**: Leia esta documentação em [English-US](README.en.md).

---

## Resumo do projeto
O **Baafes-urease-mining** é um pipeline automatizado de bioinformática desenvolvido em Nextflow DSL2 e Python 3.11, containerizado via Docker, projetado para realizar a prospecção genômica, identificação funcional e análise evolutiva de genes envolvidos no catabolismo de ureia em genomas bacterianos. 

O foco de aplicação abrange 10 genomas de bactérias aeróbias formadoras de endósporos (AEFB) da coleção CBafes/UnB isoladas de solos do Distrito Federal (acessos GenBank `VKHW00000000.1` a `VKIC00000000.1`). O sistema avalia simultaneamente três frentes metabólicas:
1. **Via Canônica da Urease**: Subunidades catalíticas ($\gamma$ - UreA, $\beta$ - UreB, $\alpha$ - UreC) e proteínas acessórias de maturação (UreD, UreE, UreF, UreG).
2. **Sistema de Transporte Ativo de Ureia**: Complexo ABC $urtABCDE$ (UrtA, UrtB, UrtC, UrtD, UrtE).
3. **Rota Alternativa Independente de Níquel**: Urea carboxilase (`uc`) acoplada à alofanato hidrolase (`ah`).

## Lista de funcionalidades
- **Pre-flight Checks & Validação Prévias de Acesso**: Validação automatizada de acessibilidade e conectividade HTTP/API aos 10 genomas no NCBI, banco Pfam (EMBL-EBI FTP) e banco Bakta (Zenodo) antes do disparo de execuções de alto custo computacional.
- **Aquisição Automática do NCBI**: Download direto dos assemblies FASTA via NCBI Datasets CLI com fallback automático para NCBI Entrez E-utilities.
- **Controle de Qualidade Genômico (QC)**: Avaliação estatística de contiguidade com QUAST (N50, L50, GC%) e completude/contaminação com CheckM2.
- **Anotação Genômica Estrutural e Funcional**: Anotação com Bakta v1.9+ gerando saídas estruturadas JSON, GFF3, FAA e GBFF.
- **Triagem Trifásica de Candidatos**: Consolidação por consenso integrando anotação Bakta (Layer A), BLASTp contra referências curadas (Layer B) e busca de domínios Pfam com HMMER3 (Layer C).
- **Matriz de Presença/Ausência**: Construção automatizada da tabela `urease_presence_absence.tsv` e sumário unificado entre as 10 estirpes.
- **Análise de Vizinhança Gênica (Sintenia)**: Extração da janela genômica de $\pm 10$ kb em torno de $ureC$ e geração de mapa interativo HTML com `clinker`.
- **Reconstrução Filogenética de UreC**: Alinhamento múltiplo com MAFFT, filtragem de posições conservadas com trimAl e inferência de árvore por Máxima Verossimilhança com IQ-TREE 2.
- **Acompanhamento de Progresso por tqdm**: Exibição de barras de progresso interativas gerais e por fase (0 a 100%) através do wrapper `run.sh`.

---

# Arquitetura

## Banco de dados (ou outra forma de armazenamento que estiver sendo usada)

### Diagrama de banco de dados
O pipeline utiliza arquivos tabulares delimitados por tabulação (`.tsv`), documentos estruturados JSON, alinhamentos FASTA e visualizadores HTML persistidos no sistema de arquivos estruturado em diretórios:

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

    ACCESSIONS_TSV ||--|| GENOME_FASTA : "especifica download"
    GENOME_FASTA ||--|| BAKTA_JSON : "anotado em"
    BAKTA_JSON ||--|{ CANDIDATES_TSV : "triado em"
    CANDIDATES_TSV }|--|| PRESENCE_ABSENCE_MATRIX : "consolidado em"
```

### Dicionário de dados

#### Tabela: `data/accessions.tsv`
| Campo | Tipo | Descrição |
| :--- | :--- | :--- |
| `strain` | String (PK) | Identificador único curto da estirpe bacteriana (ex.: S1, S2) |
| `species` | String | Binômio da espécie bacteriana (ex.: Lysinibacillus fusiformis) |
| `genbank_accession` | String | Código de acesso do assembly no GenBank (ex.: VKHW00000000.1) |

#### Tabela: `results/03_candidates/{strain}/{strain}_candidates.tsv`
| Campo | Tipo | Descrição |
| :--- | :--- | :--- |
| `strain` | String (FK) | Identificador da estirpe bacteriana |
| `feature_id` | String (PK) | Código identificador do gene/proteína anotado no Bakta |
| `gene` | String | Símbolo atribuído ao gene (ex.: ureC, urtA, uc) |
| `product` | String | Descrição funcional da proteína |
| `layer_bakta` | Boolean | Sinalização de detecção na Camada A (Bakta) |
| `layer_blast` | Boolean | Sinalização de detecção na Camada B (BLASTp) |
| `layer_hmmer` | Boolean | Sinalização de detecção na Camada C (Pfam HMMER) |
| `confidence` | String | Grau de confiança derivado do consenso (Alta, Media, Baixa) |

#### Tabela: `results/04_summary/urease_presence_absence.tsv`
| Campo | Tipo | Descrição |
| :--- | :--- | :--- |
| `strain` | String (PK) | Identificador único da estirpe |
| `ureC` | Integer (0/1) | Presença (1) ou ausência (0) da subunidade $\alpha$ da urease |
| `ureA` | Integer (0/1) | Presença (1) ou ausência (0) da subunidade $\gamma$ da urease |
| `ureB` | Integer (0/1) | Presença (1) ou ausência (0) da subunidade $\beta$ da urease |
| `ureD` | Integer (0/1) | Presença (1) ou ausência (0) da proteína acessória UreD |
| `ureE` | Integer (0/1) | Presença (1) ou ausência (0) da proteína acessória UreE |
| `ureF` | Integer (0/1) | Presença (1) ou ausência (0) da proteína acessória UreF |
| `ureG` | Integer (0/1) | Presença (1) ou ausência (0) da GTPase acessória UreG |
| `urtA` | Integer (0/1) | Presença (1) ou ausência (0) do transportador de ureia UrtA |
| `uc` | Integer (0/1) | Presença (1) ou ausência (0) da urea carboxilase |
| `ah` | Integer (0/1) | Presença (1) ou ausência (0) da alofanato hidrolase |

---

## Componentes

### Diagrama de componentes

```mermaid
graph TD
    User["Usuário / Operador"] -->|1. Executa| RunWrapper["run.sh Wrapper CLI"]
    RunWrapper -->|2. Valida Recursos| ResourceChecker["resource_checker.py"]
    ResourceChecker -->|3. Checa APIs| NCBI_API["NCBI Datasets / Entrez API"]
    ResourceChecker -->|3. Checa DBs| ExternalDBs["Pfam FTP & Bakta Zenodo"]
    
    RunWrapper -->|4. Dispara com tqdm| ProgressTracker["progress_tracker.py"]
    RunWrapper -->|5. Executa Pipeline| NextflowEngine["Nextflow DSL2 Engine (main.nf)"]
    
    subgraph Container_bafes_urease ["Container Docker: bafes_urease"]
        NextflowEngine --> P0["PREFLIGHT_CHECK"]
        P0 --> P1["DOWNLOAD_GENOME"]
        P1 --> P2["QC (QUAST & CheckM2)"]
        P2 --> P3["BAKTA_ANNOTATE"]
        P3 --> P4["EXTRACT_CANDIDATES (BLASTp + HMMER)"]
        P4 --> P5["MERGE_ALL_STRAINS"]
        P5 --> P6["ANALYZE_SYNTENY (clinker)"]
        P6 --> P7["BUILD_PHYLOGENY (MAFFT + IQ-TREE 2)"]
    end

    NextflowEngine -->|6. Publica Artefatos| ResultsFolder["results/ (Saídas e Relatórios)"]
```

### Tecnologias e suas versões

| Tecnologia / Ferramenta | Versão | Função no Sistema |
| :--- | :--- | :--- |
| **Linux Ubuntu (Jammy)** | 22.04 LTS | Sistema operacional base do container Docker |
| **Python** | 3.11.x | Linguagem para scripts auxiliares, consolidação e tqdm |
| **Java Runtime Environment**| 17.0+ | Pré-requisito para execução do Nextflow |
| **Nextflow** | >= 24.04.0 | Motor de orquestração de workflows bioinformáticos |
| **Docker** | >= 24.0+ | Engine de containerização e isolamento do ambiente |
| **NCBI Datasets CLI / Entrez** | v2 / E-utils | Download e busca de assemblies no GenBank |
| **Bakta** | >= 1.9.0 | Anotação genômica funcional estruturada |
| **QUAST** | >= 5.2.0 | Avaliação de métricas de contiguidade genômica |
| **CheckM2** | >= 1.0.2 | Avaliação de completude e contaminação bacteriana |
| **BLAST+ (blastp)** | >= 2.14.0 | Alinhamento local contra referências proteicas |
| **HMMER3 (hmmscan)** | >= 3.3.2 | Busca de domínios conservados contra banco Pfam-A |
| **clinker** | >= 0.0.28 | Visualização interativa de sintenia e vizinhança gênica |
| **MAFFT** | >= 7.520 | Alinhamento múltiplo de sequências de aminoácidos |
| **trimAl** | >= 1.4.1 | Remoção de posições mal alinhadas ou gappy |
| **IQ-TREE 2** | >= 2.2.0 | Inferência de árvores filogenéticas por Máxima Verossimilhança |
| **pandas / tqdm / Biopython** | Últimas estáveis | Manipulação de matrizes, barras de progresso e FASTA |

---

## Funcionalidades

### Requisitos

| Funcionalidade | Campo de Formulário / Entrada | Campo de Banco / Arquivo de Saída | Regras Aplicadas |
| :--- | :--- | :--- | :--- |
| **Pre-flight Resource Check** | `--accessions`, `--references` | Console log / `preflight_status.txt` | Testa HTTP Status 200/JSON no NCBI, EBI Pfam e Bakta Zenodo. Aborta com Exit Code 4 se insumo crítico falhar. |
| **Download de Genomas** | `genbank_accession` em `data/accessions.tsv` | `results/00_genomes/{strain}.fna` | Resolve o accession WGS master para o assembly (GCA_/GCF_) e baixa via NCBI Datasets CLI, com fallback para o FTP do assembly. Valida contigs, número de bases e alfabeto. **Falha se o genoma real não for obtido.** |
| **QC de Genomas** | FASTA do genoma baixado | `results/01_qc/quast/` e `checkm2/` | QUAST calcula N50, GC% e comprimento total; CheckM2 estima completude/contaminação a partir do banco baixado no bootstrap. |
| **Anotação com Bakta** | FASTA do genoma + `--bakta_db` | `results/02_bakta/{strain}/` (.json, .gff3, .faa, .gbff) | Atribui identificadores RefSeq/UniRef, termos EC (3.5.1.5, 6.3.4.6, 3.5.1.54) e símbolos gênicos. |
| **Triagem de Candidatos** | Bakta JSON, `urease_references.fasta`, `Pfam-A.hmm` | `results/03_candidates/{strain}_candidates.tsv` | Classifica a evidência em 3 camadas: Bakta (A), BLASTp (B, E-value <= 1e-5 + identidade >= 30% + cobertura >= 50%) e HMMER (C, score de domínio >= 20). Confiança Alta (>=2 camadas), Média (1 camada). |
| **Consolidação de Matriz** | Candidatos TSV de todas as estirpes | `results/04_summary/urease_presence_absence.tsv` | Cria matriz binária (0/1) para os 10 genes alvos ($ureA-G$, $urtA$, $uc$, $ah$) por estirpe. |
| **Análise de Sintenia** | Arquivos GBFF das estirpes | `results/05_synteny/synteny.html` | Executa clinker sobre os GenBank anotados (exige >= 2 genomas). Falha explicitamente se não houver GenBank ou se o clinker falhar — nenhum HTML de placeholder é escrito. |
| **Filogenia de UreC** | UreC selecionado a partir dos TSVs de candidatos | `results/06_phylogeny/ureC.treefile`, `ureC_labels.tsv`, `phylogeny_status.txt` | Usa **apenas** as sequências de UreC identificadas na triagem, não o proteoma inteiro. Exige >= 3 UreC; abaixo disso registra o motivo em `phylogeny_status.txt` e não gera árvore. MAFFT -> trimAl `automated1` -> IQ-TREE 2 (`-m MFP`, 1000 ultrafast bootstraps). |
| **Acompanhamento tqdm** | Execução via `run.sh --exec` | Terminal UI com barras interativas | Renderiza o percentual de conclusão por fase e progresso global. |

---

---

## Integridade dos dados

O pipeline **não fabrica dados em nenhuma hipótese**. Toda etapa que não consegue produzir
um resultado real a partir dos genomas do NCBI falha explicitamente (exit != 0) em vez de
substituir a saída por um valor plausível. Em particular:

- **Genomas.** Os accessions do estudo são WGS master (ex.: `VKHW00000000.1`) e esses registros
  não contêm sequência — um `efetch` neles devolve zero bytes. O
  [download_genome.py](bin/download_genome.py) resolve o accession para o assembly
  (`GCA_`/`GCF_`) antes de baixar, e valida contigs, número de bases (>= 500 kb) e alfabeto
  nucleotídico. Sem genoma real, a estirpe não segue no pipeline.
- **QC.** CheckM2 é executado de verdade contra seu banco. Não há valores de
  completude/contaminação escritos no código.
- **Anotação.** Sem o banco do Bakta, [BAKTA_ANNOTATE](main.nf) falha. Não há anotação
  substituta com genes inventados.
- **Referências.** As 22 sequências de [urease_references.fasta](data/urease_references.fasta)
  são baixadas do UniProt a partir de accessions fixadas em
  [urease_reference_accessions.tsv](data/urease_reference_accessions.tsv), com o comprimento
  de cada sequência conferido contra o valor esperado. Regenere com:
  ```bash
  python3 bin/fetch_references.py
  ```
- **Sintenia e filogenia.** Nenhum HTML de placeholder e nenhuma árvore parcial. Quando há
  menos de 3 sequências de UreC, isso é registrado em `phylogeny_status.txt` como resultado
  legítimo — não há árvore possível — e nenhum arquivo de árvore é escrito.

### Conjunto de referências curadas

| Origem | Genes | Observação |
| :--- | :--- | :--- |
| *Sporosarcina pasteurii* (ex-*Bacillus pasteurii*) | ureA, ureB, ureC, ureD, ureE, ureF, ureG | Operon completo revisado (SwissProt); modelo ureolítico de referência |
| *Bacillus subtilis* 168 | ureA, ureB, ureC | Subunidades estruturais revisadas |
| *Bacillus* sp. TB-90 | ureD, ureE, ureF, ureG, ureH | Acessórias revisadas |
| *Prochlorococcus marinus* MED4 | urtA–urtE | Operon de transporte de ureia (PMM0970–0974) |
| *Oleomonas sagaranensis* | uc | Urease carboxilase bacteriana caracterizada (EC 6.3.4.6) |
| *Pseudomonas* sp. ADP | ah | Alofanato hidrolase AtzF revisada (EC 3.5.1.54) |

---

## Instalação e uso

### Pre-requisitos

Todo o pipeline roda **dentro do container** `bafes_urease`. O host precisa apenas de:

- Docker Engine >= 24.0 (com o usuário no grupo `docker`)
- Git, Bash e `curl`
- Espaço em disco: **~160 GB livres** para o banco completo do Bakta (use `--bakta-db-type light`
  para ~10 GB, ou `--skip-bakta-db` se já tiver o banco)

Nextflow, Bakta, BLAST+, HMMER, MAFFT, trimAl, IQ-TREE, QUAST e clinker vêm todos na imagem —
não é preciso instalar nada disso no host. Para rodar sem Docker (com Conda/Mamba já configurado),
use `--runtime local`.

### Passo a Passo de Execução Recomendada

#### 1. Clonar o repositório
```bash
git clone https://github.com/waldeyr/bafes-urease.git
cd bafes-urease
```

#### 2. Tornar o script wrapper executável
```bash
chmod +x run.sh
```

#### 3. Sequência Quick Start
Execute a sequência de 3 passos recomendada:

```bash
# Passo A: Bootstrap — constrói a imagem, baixa Pfam e o banco completo do Bakta
#          Leva horas na primeira vez (imagem ~20 min + banco Bakta full ~75 GB).
./run.sh --bootstrap --runtime docker --container-image bafes_urease

# Passo B: Build e Pre-flight Check (validação de recursos e acessibilidade)
./run.sh --build --runtime docker --container-image bafes_urease

# Passo C: Execução do pipeline com log detalhado e acompanhamento por tqdm
./run.sh --exec --runtime docker --container-image bafes_urease --verbose
```

O `--bootstrap` é idempotente: cada etapa concluída é pulada nas execuções seguintes.

#### 4. Opções úteis do bootstrap

| Opção | Efeito |
| :--- | :--- |
| `--force-rebuild` | Reconstrói a imagem Docker mesmo que ela já exista |
| `--bakta-db-type light` | Baixa o banco reduzido do Bakta (~10 GB em vez de ~75 GB) |
| `--skip-bakta-db` | Não baixa o banco do Bakta (use se já tiver um em `--bakta-db`) |
| `--runtime local` | Executa no host, sem container (exige Conda/Mamba já configurado) |

#### 5. Configuração via `.env`

Copie o modelo e preencha. O `run.sh` carrega o `.env` automaticamente e repassa as
variáveis para dentro do container. O `.env` está no `.gitignore` — **nunca versione sua chave**.

```bash
cp .env.example .env
```

| Variável | Efeito |
| :--- | :--- |
| `NCBI_API_KEY` | Eleva o limite do E-utilities de 3 para 10 req/s. **Recomendado**: sem ela, o pre-check dos 10 accessions e o download dos genomas ficam limitados a 3 req/s e mais sujeitos a HTTP 429. Gere em [account.ncbi.nlm.nih.gov](https://account.ncbi.nlm.nih.gov/) → Account Settings → API Key Management |
| `NCBI_EMAIL` | E-mail de contato exigido pela política do NCBI E-utilities |
| `BAFES_INSECURE_SSL=1` | Desabilita a verificação TLS (só para redes com proxy TLS interceptador) |

#### 6. Execução Direta via Nextflow (Opcional)

Requer Nextflow e as ferramentas no host, ou o perfil `docker` (um container por processo):

```bash
nextflow run main.nf \
   --accessions data/accessions.tsv \
   --bakta_db ./db/db \
   --pfam_hmm ./Pfam-A.hmm \
   --references data/urease_references.fasta \
   --outdir results \
   -profile docker
```

Os relatórios (`nf_report.html`, `nf_trace.txt`, `nf_timeline.html`, `nf_dag.html`) já vêm
habilitados no `nextflow.config` — não é preciso passar as flags `-with-*`.

#### 7. Estrutura de Saídas Geradas
```
results/
├── 00_preflight/    # Status dos testes de validação prévia
├── 00_genomes/      # Assemblies FASTA baixados do NCBI
├── 01_qc/           # Relatórios de QC (QUAST + CheckM2)
├── 02_bakta/        # Anotações em GFF3, FAA, JSON, GBFF por estirpe
├── 03_candidates/   # Candidatos consolidados por estirpe (Bakta + BLASTp + HMMER)
├── 04_summary/      # Matriz urease_presence_absence.tsv e resumo unificado
├── 05_synteny/      # Visualizador interativo de sintenia (synteny.html)
└── 06_phylogeny/    # Árvore filogenética de UreC (ureC.treefile e alinhamentos)
```

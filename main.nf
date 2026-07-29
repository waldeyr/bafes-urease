#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

/*
========================================================================================
   BAFES UREASE MINING PIPELINE (DSL2)
   PT-BR: Pipeline de mineração genômica de urease, transporte de ureia e via alternativa.
   EN-US: Genomic mining pipeline for urease, urea transport, and alternative pathways.

   PT-BR: PRINCÍPIO — nenhum processo fabrica dados. Toda etapa que não puder produzir
          resultado real a partir dos genomas do NCBI falha explicitamente (exit != 0).
   EN-US: PRINCIPLE — no process fabricates data. Any step that cannot produce a real
          result from the NCBI genomes fails explicitly (exit != 0).
========================================================================================
*/

params.accessions = "data/accessions.tsv"
params.references = "data/urease_references.fasta"
params.pfam_hmm   = "Pfam-A.hmm"
params.bakta_db   = "db/db"
params.checkm2_db = "db/checkm2/CheckM2_database/uniref100.KO.1.dmnd"
params.outdir     = "results"

// PT-BR: As sequências de ureC curadas de params.references entram na árvore como contexto
//        filogenético. Sem elas a árvore fica restrita aos UreC encontrados nas estirpes —
//        que podem ser poucos demais para bootstrap, ou até para haver mais de uma
//        topologia possível.
// EN-US: The curated ureC sequences from params.references enter the tree as phylogenetic
//        context. Without them the tree is limited to the UreC found in the strains — which
//        may be too few for bootstrap, or even for more than one topology to exist.
params.phylo_include_refs = true

// PT-BR: Genes flanqueadores de cada lado do locus de urease no recorte para a sintenia.
//        5 dá um contexto de ~11 genes por região — suficiente para ver rearranjo e perda
//        de gene sem arrastar o cromossomo inteiro para dentro do clinker.
// EN-US: Flanking genes on each side of the urease locus in the synteny slice. 5 gives a
//        ~11-gene context per region — enough to see rearrangement and gene loss without
//        dragging the whole chromosome into clinker.
params.synteny_flank_genes = 5

// PT-BR: Por padrão o recorte ancora só em candidatos de alta confiança ou com perfil Pfam
//        diagnóstico da via. Ligar isto ancora em TODO candidato, inclusive os domínios
//        promíscuos (cobW, CT_C_D), que aparecem espalhados longe de qualquer operon de
//        urease e geram regiões espúrias.
// EN-US: By default the slice anchors only on high-confidence candidates or on ones with a
//        Pfam profile diagnostic of the pathway. Turning this on anchors on EVERY candidate,
//        promiscuous domains included (cobW, CT_C_D), which sit scattered far from any
//        urease operon and produce spurious regions.
params.synteny_all_candidates = false

// PT-BR: Caminhos de banco entram no script de cada processo como texto bruto, e cada
//        task roda a partir do seu próprio work dir — um valor relativo como "db/db"
//        precisa ser resolvido contra o diretório de lançamento antes de ser escrito
//        no script, senão nunca é encontrado.
// EN-US: Database paths are interpolated into each process script as raw text, and every
//        task runs from its own work dir — a relative value like "db/db" must be resolved
//        against the launch directory before it lands in the script, or it is never found.
def absPath(p) {
    file(p).toAbsolutePath().toString()
}

// PT-BR: Processo 0 - Validação Prévia de Recursos / EN-US: Process 0 - Pre-flight Resource Check
process PREFLIGHT_CHECK {
    tag "resource_check"
    publishDir "${params.outdir}/00_preflight", mode: 'copy'

    input:
    path accessions_file
    path references_file

    output:
    path "preflight_status.txt", emit: status

    script:
    """
    python3 ${projectDir}/bin/resource_checker.py \\
        --accessions ${accessions_file} \\
        --references ${references_file} \\
        --bakta-db ${absPath(params.bakta_db)} \\
        --pfam ${absPath(params.pfam_hmm)} | tee preflight_status.txt
    """
}

// PT-BR: Processo 1 - Download Genômico NCBI / EN-US: Process 1 - NCBI Genome Download
process DOWNLOAD_GENOME {
    tag "${strain}"
    publishDir "${params.outdir}/00_genomes", mode: 'copy'

    input:
    tuple val(strain), val(species), val(accession)

    output:
    tuple val(strain), path("${strain}.fna"), emit: genome_fasta

    script:
    // PT-BR: O script resolve o WGS master para o assembly (GCA_/GCF_) e valida o
    //        FASTA obtido. Se o genoma real não vier, ele sai com erro — o processo
    //        falha e a estirpe NÃO segue no pipeline com dados inventados.
    // EN-US: The script resolves the WGS master to its assembly (GCA_/GCF_) and
    //        validates the FASTA. If the real genome does not arrive, it exits with an
    //        error — the process fails and the strain does NOT continue with made-up data.
    """
    python3 ${projectDir}/bin/download_genome.py \\
        --accession ${accession} \\
        --strain ${strain} \\
        --out ${strain}.fna
    """
}

// PT-BR: Processo 2 - Controle de Qualidade QUAST / EN-US: Process 2 - Quality Control QUAST
process QC_QUAST {
    tag "${strain}"
    publishDir { "${params.outdir}/01_qc/quast/${strain}" }, mode: 'copy'

    input:
    tuple val(strain), path(genome_fasta)

    output:
    path "*"

    script:
    """
    quast.py ${genome_fasta} -o . --threads ${task.cpus}
    """
}

// PT-BR: Processo 3 - Controle de Qualidade CheckM2 / EN-US: Process 3 - Quality Control CheckM2
process QC_CHECKM2 {
    tag "${strain}"
    publishDir { "${params.outdir}/01_qc/checkm2/${strain}" }, mode: 'copy'

    input:
    tuple val(strain), path(genome_fasta)

    output:
    path "${strain}_checkm2.tsv"

    script:
    // PT-BR: CheckM2 é executado de verdade. Sem banco, o processo falha — não há
    //        mais valores de completude/contaminação fixados no código.
    // EN-US: CheckM2 actually runs. Without its database the process fails — there are
    //        no more hardcoded completeness/contamination values.
    """
    if [ ! -f "${absPath(params.checkm2_db)}" ]; then
        echo "[ERRO/ERROR] Banco do CheckM2 não encontrado / CheckM2 database not found: ${absPath(params.checkm2_db)}" >&2
        echo "             Rode / Run: ./run.sh --bootstrap" >&2
        exit 1
    fi

    checkm2 predict \\
        --input ${genome_fasta} \\
        --output-directory checkm2_out \\
        --database_path ${absPath(params.checkm2_db)} \\
        --threads ${task.cpus} \\
        --extension .fna \\
        --force

    cp checkm2_out/quality_report.tsv ${strain}_checkm2.tsv
    """
}

// PT-BR: Processo 4 - Anotação Bakta / EN-US: Process 4 - Bakta Functional Annotation
process BAKTA_ANNOTATE {
    tag "${strain}"
    publishDir { "${params.outdir}/02_bakta/${strain}" }, mode: 'copy'

    input:
    tuple val(strain), path(genome_fasta)

    output:
    tuple val(strain), path("${strain}.gff3"), path("${strain}.faa"), path("${strain}.json"), path("${strain}.gbff"), emit: bakta_results

    script:
    // PT-BR: Sem banco do Bakta o processo falha. A anotação fictícia com um gene
    //        ureC inventado foi removida — ela produzia candidatos falsos.
    // EN-US: Without the Bakta database the process fails. The fake annotation with a
    //        made-up ureC gene was removed — it produced false candidates.
    """
    if [ ! -d "${absPath(params.bakta_db)}" ]; then
        echo "[ERRO/ERROR] Banco do Bakta não encontrado / Bakta database not found: ${absPath(params.bakta_db)}" >&2
        echo "             Rode / Run: ./run.sh --bootstrap" >&2
        exit 1
    fi

    bakta --db ${absPath(params.bakta_db)} \\
        --prefix ${strain} \\
        --output . \\
        --threads ${task.cpus} \\
        --force \\
        ${genome_fasta}
    """
}

// PT-BR: Processo 5 - Triagem Trifásica de Candidatos / EN-US: Process 5 - 3-Layer Candidate Screening
process EXTRACT_CANDIDATES {
    tag "${strain}"
    publishDir { "${params.outdir}/03_candidates/${strain}" }, mode: 'copy'

    input:
    tuple val(strain), path(gff3), path(faa), path(json_file), path(gbff)
    path ref_fasta
    path pfam_files
    // PT-BR: O script entra como arquivo de entrada, não por caminho absoluto: só assim o
    //        Nextflow inclui o conteúdo dele no hash da tarefa. Chamado por ${projectDir}
    //        ele fica fora do hash, e um `-resume` depois de corrigir a triagem devolveria
    //        silenciosamente os TSVs antigos como se fossem novos.
    // EN-US: The script comes in as an input file, not via an absolute path: only then does
    //        Nextflow include its content in the task hash. Called through ${projectDir} it
    //        stays out of the hash, and a `-resume` after fixing the screening would
    //        silently hand back the old TSVs as if they were fresh.
    path script_py

    output:
    tuple val(strain), path("${strain}_candidates.tsv"), emit: candidates_tsv
    tuple val(strain), path("${strain}.gbff"), emit: gbff_out
    path "${strain}_blast.tsv"
    path "${strain}_hmmer.tbl"

    script:
    def pfam_base = file(params.pfam_hmm).name
    // PT-BR: BLASTp já filtra por E-value; identidade e cobertura são aplicadas no
    //        script Python. Erros de ferramenta não são mais mascarados com `|| true`.
    // EN-US: BLASTp already filters by E-value; identity and coverage are applied in the
    //        Python script. Tool errors are no longer masked with `|| true`.
    """
    makeblastdb -in ${ref_fasta} -dbtype prot -out ref_db

    blastp -query ${faa} \\
        -db ref_db \\
        -out ${strain}_blast.tsv \\
        -outfmt '6 qseqid sseqid pident length qstart qend qlen slen evalue bitscore' \\
        -evalue 1e-5 \\
        -num_threads ${task.cpus}

    if [ ! -f "${pfam_base}.h3i" ]; then
        echo "[ERRO/ERROR] Pfam não indexado / Pfam not pressed: ${pfam_base}.h3i ausente." >&2
        echo "             Rode / Run: ./run.sh --bootstrap" >&2
        exit 1
    fi

    hmmscan --domtblout ${strain}_hmmer.tbl \\
        --cpu ${task.cpus} \\
        ${pfam_base} ${faa} > /dev/null

    python3 ${script_py} \\
        --strain ${strain} \\
        --bakta-json ${json_file} \\
        --blast-tsv ${strain}_blast.tsv \\
        --hmmer-tbl ${strain}_hmmer.tbl \\
        --out-tsv ${strain}_candidates.tsv
    """
}

// PT-BR: Processo 6 - Matriz de Presença/Ausência / EN-US: Process 6 - Presence/Absence Matrix
process MERGE_ALL_STRAINS {
    publishDir "${params.outdir}/04_summary", mode: 'copy'

    input:
    path candidate_tsvs
    // PT-BR: Entra como arquivo para que o conteúdo do script conte no hash da tarefa.
    // EN-US: Comes in as a file so the script's content counts towards the task hash.
    path script_py

    output:
    path "urease_presence_absence.tsv", emit: matrix
    path "summary_unified.tsv"

    script:
    """
    mkdir -p input_tsvs
    cp ${candidate_tsvs} input_tsvs/
    python3 ${script_py} \\
        --candidates-dir input_tsvs \\
        --out-matrix urease_presence_absence.tsv \\
        --out-summary summary_unified.tsv
    """
}

// PT-BR: Processo 7a - Recorte do Locus de Urease / EN-US: Process 7a - Urease Locus Extraction
//
// PT-BR: O clinker alinha todas as proteínas contra todas, em todos os pares de entradas.
//        Passar os 10 genomas inteiros (~46.500 CDS) é da ordem de 10^9 alinhamentos: a
//        etapa estourou as 8 h de `time` e foi morta com SIGTERM (exit 143). O recorte
//        prévio do locus reduz isso a algumas dezenas de genes por estirpe — que é o caso
//        de uso para o qual o clinker foi escrito — e produz um HTML que abre no navegador.
// EN-US: clinker aligns every protein against every other, across every pair of inputs.
//        Feeding it the 10 whole genomes (~46,500 CDS) is on the order of 10^9 alignments:
//        the step blew past its 8 h `time` limit and was killed with SIGTERM (exit 143).
//        Slicing the locus first cuts that to a few dozen genes per strain — the use case
//        clinker was written for — and yields an HTML a browser can actually open.
process EXTRACT_LOCUS {
    tag "${strain}"
    publishDir "${params.outdir}/05_synteny/loci", mode: 'copy'

    input:
    tuple val(strain), path(gbff), path(candidates)
    // PT-BR: Entra como arquivo para que o conteúdo do script conte no hash da tarefa.
    // EN-US: Comes in as a file so the script's content counts towards the task hash.
    path script_py

    output:
    // PT-BR: O GenBank é opcional porque "esta estirpe não tem urease" é um resultado
    //        biológico legítimo — e nesse caso nenhuma região é inventada só para a
    //        estirpe aparecer no mapa. O status é sempre emitido e diz o porquê.
    // EN-US: The GenBank is optional because "this strain has no urease" is a legitimate
    //        biological outcome — and in that case no region is invented just to get the
    //        strain onto the map. The status is always emitted and says why.
    path "${strain}.gbk", emit: locus_gbk, optional: true
    path "${strain}_locus_status.txt", emit: status

    script:
    def all_flag = params.synteny_all_candidates ? "--all-candidates \\\n        " : ""
    """
    python3 ${script_py} \\
        --strain ${strain} \\
        --gbff ${gbff} \\
        --candidates ${candidates} \\
        --flank-genes ${params.synteny_flank_genes} \\
        ${all_flag}--out ${strain}.gbk \\
        --status ${strain}_locus_status.txt
    """
}

// PT-BR: Processo 7b - Análise de Sintenia (clinker) / EN-US: Process 7b - Synteny Analysis (clinker)
process ANALYZE_SYNTENY {
    publishDir "${params.outdir}/05_synteny", mode: 'copy'

    input:
    path locus_gbks

    output:
    path "synteny.html"

    script:
    // PT-BR: Só os .gbk de locus são encenados aqui, então o diretório da tarefa já é o
    //        diretório de entrada — não há mais cópia para gbk_inputs/, que quebrava
    //        quando nenhuma estirpe tinha locus (`cp` sem argumento de origem).
    // EN-US: Only the locus .gbk files are staged here, so the task directory already is
    //        the input directory — no more copy into gbk_inputs/, which broke when no
    //        strain had a locus (`cp` with no source argument).
    """
    python3 ${projectDir}/bin/synteny_analysis.py \\
        --gbk-dir . \\
        --out-dir .
    """
}

// PT-BR: Processo 8 - Filogenia UreC / EN-US: Process 8 - UreC Phylogeny Reconstruction
process BUILD_PHYLOGENY {
    publishDir "${params.outdir}/06_phylogeny", mode: 'copy'

    input:
    path candidate_tsvs
    path faa_files
    path ref_fasta
    // PT-BR: Entra como arquivo para que o conteúdo do script conte no hash da tarefa.
    // EN-US: Comes in as a file so the script's content counts towards the task hash.
    path script_py

    output:
    // PT-BR: O status é sempre emitido; a árvore é opcional porque "menos de 3 UreC"
    //        é um resultado biológico legítimo, não uma falha — e nesse caso nenhuma
    //        árvore é inventada.
    // EN-US: The status is always emitted; the tree is optional because "fewer than 3
    //        UreC" is a legitimate biological outcome, not a failure — and in that case
    //        no tree is invented.
    path "phylogeny_status.txt"
    path "ureC*", optional: true

    script:
    def refs_arg = params.phylo_include_refs ? "--references ${ref_fasta} \\\n        " : ""
    """
    cat ${faa_files} > all_proteins.faa
    python3 ${script_py} \\
        --fasta-in all_proteins.faa \\
        --candidates-dir . \\
        ${refs_arg}--threads ${task.cpus} \\
        --out-dir .
    """
}

workflow {
    // PT-BR: log.info precisa ficar dentro do workflow — o parser estrito do
    //        Nextflow rejeita instruções no nível superior do script.
    // EN-US: log.info must live inside the workflow — Nextflow's strict parser
    //        rejects statements at the script's top level.
    log.info """
========================================================================================
 BAFES-UREASE-MINING PIPELINE
========================================================================================
 Accessions File  : ${params.accessions}
 References FASTA : ${params.references}
 Pfam HMM File    : ${absPath(params.pfam_hmm)}
 Bakta DB Path    : ${absPath(params.bakta_db)}
 CheckM2 DB Path  : ${absPath(params.checkm2_db)}
 Output Directory : ${params.outdir}
========================================================================================
"""

    accessions_ch = Channel.fromPath(params.accessions)
        .splitCsv(header: true, sep: '\t')
        .map { row -> tuple(row.strain, row.species, row.genbank_accession) }

    // PT-BR: O Pfam e seus índices do hmmpress precisam ser encenados juntos.
    // EN-US: Pfam and its hmmpress indices must be staged together.
    pfam_ch = Channel.fromPath("${params.pfam_hmm}*").collect()

    PREFLIGHT_CHECK(file(params.accessions), file(params.references))

    DOWNLOAD_GENOME(accessions_ch)
    QC_QUAST(DOWNLOAD_GENOME.out.genome_fasta)
    QC_CHECKM2(DOWNLOAD_GENOME.out.genome_fasta)

    BAKTA_ANNOTATE(DOWNLOAD_GENOME.out.genome_fasta)

    EXTRACT_CANDIDATES(
        BAKTA_ANNOTATE.out.bakta_results,
        file(params.references),
        pfam_ch,
        file("${projectDir}/bin/extract_urease_candidates.py")
    )

    all_candidates = EXTRACT_CANDIDATES.out.candidates_tsv.map { _strain, tsv -> tsv }.collect()
    all_faas = BAKTA_ANNOTATE.out.bakta_results.map { _strain, _gff, faa, _json, _gbff -> faa }.collect()

    MERGE_ALL_STRAINS(
        all_candidates,
        file("${projectDir}/bin/merge_strain_results.py")
    )

    // PT-BR: join casa o GenBank e o TSV da MESMA estirpe pela chave; sem isso, dois
    //        `collect()` independentes poderiam parear S1.gbff com o TSV de outra estirpe.
    // EN-US: join matches the GenBank and the TSV of the SAME strain by key; without it, two
    //        independent `collect()` calls could pair S1.gbff with another strain's TSV.
    locus_input = EXTRACT_CANDIDATES.out.gbff_out
        .join(EXTRACT_CANDIDATES.out.candidates_tsv)

    EXTRACT_LOCUS(
        locus_input,
        file("${projectDir}/bin/extract_urease_locus.py")
    )

    // PT-BR: ifEmpty([]) garante que o processo rode mesmo sem nenhum locus: aí o
    //        synteny_analysis.py falha explicitamente ("nenhum GenBank"), em vez de o
    //        Nextflow simplesmente pular a etapa e a sintenia sumir do relatório em silêncio.
    // EN-US: ifEmpty([]) makes the process run even with no locus at all: then
    //        synteny_analysis.py fails explicitly ("no GenBank"), rather than Nextflow just
    //        skipping the step and synteny vanishing from the report in silence.
    ANALYZE_SYNTENY(EXTRACT_LOCUS.out.locus_gbk.collect().ifEmpty([]))
    BUILD_PHYLOGENY(
        all_candidates,
        all_faas,
        file(params.references),
        file("${projectDir}/bin/build_phylogeny.py")
    )
}

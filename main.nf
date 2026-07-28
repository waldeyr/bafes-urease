#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

/*
========================================================================================
   BAAFES UREASE MINING PIPELINE (DSL2)
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
        --bakta-db ${params.bakta_db} \\
        --pfam ${params.pfam_hmm} | tee preflight_status.txt
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
    if [ ! -f "${params.checkm2_db}" ]; then
        echo "[ERRO/ERROR] Banco do CheckM2 não encontrado / CheckM2 database not found: ${params.checkm2_db}" >&2
        echo "             Rode / Run: ./run.sh --bootstrap" >&2
        exit 1
    fi

    checkm2 predict \\
        --input ${genome_fasta} \\
        --output-directory checkm2_out \\
        --database_path ${params.checkm2_db} \\
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
    if [ ! -d "${params.bakta_db}" ]; then
        echo "[ERRO/ERROR] Banco do Bakta não encontrado / Bakta database not found: ${params.bakta_db}" >&2
        echo "             Rode / Run: ./run.sh --bootstrap" >&2
        exit 1
    fi

    bakta --db ${params.bakta_db} \\
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

    output:
    tuple val(strain), path("${strain}_candidates.tsv"), emit: candidates_tsv
    path "${strain}.gbff", emit: gbff_out
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

    python3 ${projectDir}/bin/extract_urease_candidates.py \\
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

    output:
    path "urease_presence_absence.tsv", emit: matrix
    path "summary_unified.tsv"

    script:
    """
    mkdir -p input_tsvs
    cp ${candidate_tsvs} input_tsvs/
    python3 ${projectDir}/bin/merge_strain_results.py \\
        --candidates-dir input_tsvs \\
        --out-matrix urease_presence_absence.tsv \\
        --out-summary summary_unified.tsv
    """
}

// PT-BR: Processo 7 - Análise de Sintenia (clinker) / EN-US: Process 7 - Synteny Analysis (clinker)
process ANALYZE_SYNTENY {
    publishDir "${params.outdir}/05_synteny", mode: 'copy'

    input:
    path gbff_files

    output:
    path "synteny.html"

    script:
    """
    mkdir -p gbk_inputs
    cp ${gbff_files} gbk_inputs/
    python3 ${projectDir}/bin/synteny_analysis.py \\
        --gbk-dir gbk_inputs \\
        --out-dir .
    """
}

// PT-BR: Processo 8 - Filogenia UreC / EN-US: Process 8 - UreC Phylogeny Reconstruction
process BUILD_PHYLOGENY {
    publishDir "${params.outdir}/06_phylogeny", mode: 'copy'

    input:
    path candidate_tsvs
    path faa_files

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
    """
    cat ${faa_files} > all_proteins.faa
    python3 ${projectDir}/bin/build_phylogeny.py \\
        --fasta-in all_proteins.faa \\
        --candidates-dir . \\
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
 BAAFES-UREASE-MINING PIPELINE
========================================================================================
 Accessions File  : ${params.accessions}
 References FASTA : ${params.references}
 Pfam HMM File    : ${params.pfam_hmm}
 Bakta DB Path    : ${params.bakta_db}
 CheckM2 DB Path  : ${params.checkm2_db}
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
        pfam_ch
    )

    all_candidates = EXTRACT_CANDIDATES.out.candidates_tsv.map { _strain, tsv -> tsv }.collect()
    all_gbffs = EXTRACT_CANDIDATES.out.gbff_out.collect()
    all_faas = BAKTA_ANNOTATE.out.bakta_results.map { _strain, _gff, faa, _json, _gbff -> faa }.collect()

    MERGE_ALL_STRAINS(all_candidates)
    ANALYZE_SYNTENY(all_gbffs)
    BUILD_PHYLOGENY(all_candidates, all_faas)
}

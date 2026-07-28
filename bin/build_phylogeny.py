#!/usr/bin/env python3
"""
build_phylogeny.py
PT-BR: Seleciona as sequencias de UreC, alinha (MAFFT), apara (trimAl) e infere a arvore (IQ-TREE 2).
EN-US: Selects UreC sequences, aligns (MAFFT), trims (trimAl), and infers the tree (IQ-TREE 2).

PT-BR: Tres garantias — (1) as sequencias das estirpes vem apenas da triagem de UreC, nunca
       do proteoma inteiro; (2) falha de ferramenta aborta com erro, em vez de deixar para
       tras um arquivo parcial que passaria por resultado; (3) o status declara sempre
       quantos taxons sao das estirpes e quantos sao de referencia, para que a arvore nunca
       seja lida como se todas as folhas fossem genomas do estudo.
EN-US: Three guarantees — (1) the strain sequences come only from the UreC screening, never
       from the whole proteome; (2) a tool failure aborts with an error instead of leaving
       behind a partial file that would pass for a result; (3) the status always states how
       many taxa are from strains and how many are references, so the tree is never read as
       if every leaf were a genome from the study.

PT-BR: Sobre a contagem de taxons: o IQ-TREE recusa bootstrap com menos de 4 sequencias, e
       uma arvore nao-enraizada de 3 folhas tem uma unica topologia possivel. Por isso as
       sequencias de UreC de referencia (curadas, com proveniencia registrada) entram como
       contexto filogenetico, e ha uma guarda explicita para os casos degenerados.
EN-US: On taxon counts: IQ-TREE refuses bootstrap with fewer than 4 sequences, and an
       unrooted 3-leaf tree has a single possible topology. Hence the reference UreC
       sequences (curated, with recorded provenance) enter as phylogenetic context, and
       there is an explicit guard for the degenerate cases.
"""

import sys
import os
import glob
import argparse
import subprocess

# PT-BR: Termos que identificam a subunidade alfa da urease / EN-US: Terms identifying the urease alpha subunit
UREC_GENE_TERMS = ("urec",)
UREC_PRODUCT_TERMS = ("urease subunit alpha", "urease alpha", "urea amidohydrolase subunit alpha")

# PT-BR: Minimo de taxons para inferir uma arvore e para que o bootstrap faca sentido.
# EN-US: Minimum taxa to infer a tree at all and for bootstrap to make sense.
MIN_TAXA_FOR_TREE = 3
MIN_TAXA_FOR_BOOTSTRAP = 4
BOOTSTRAP_REPLICATES = 1000

# PT-BR: Tamanho maximo do nome do organismo no rotulo da folha, para a arvore continuar
#        legivel — o mapeamento completo fica em ureC_labels.tsv.
# EN-US: Maximum organism name length in the leaf label, to keep the tree readable — the
#        full mapping lives in ureC_labels.tsv.
MAX_ORGANISM_CHARS = 40


def read_fasta(path):
    """
    PT-BR: Le um FASTA e devolve {id: (cabecalho, sequencia)} usando o 1o token como id.
    EN-US: Reads a FASTA and returns {id: (header, sequence)} using the 1st token as id.
    """
    records = {}
    header = None
    chunks = []
    with open(path, 'r', encoding='utf-8', errors='replace') as handle:
        for line in handle:
            line = line.rstrip('\n')
            if line.startswith('>'):
                if header is not None:
                    records[header.split()[0]] = (header, ''.join(chunks))
                header = line[1:]
                chunks = []
            elif header is not None:
                chunks.append(line.strip())
    if header is not None:
        records[header.split()[0]] = (header, ''.join(chunks))
    return records


def collect_urec_ids(candidates_dir):
    """
    PT-BR: Varre os TSVs de candidatos e devolve {feature_id: estirpe} para os UreC.
    EN-US: Scans the candidate TSVs and returns {feature_id: strain} for the UreC entries.
    """
    wanted = {}
    for tsv in sorted(glob.glob(os.path.join(candidates_dir, "*_candidates.tsv"))):
        with open(tsv, 'r', encoding='utf-8') as handle:
            header = handle.readline().rstrip('\n').split('\t')
            try:
                i_id = header.index('feature_id')
                i_gene = header.index('gene')
                i_product = header.index('product')
            except ValueError:
                continue
            i_strain = header.index('strain') if 'strain' in header else -1
            for line in handle:
                parts = line.rstrip('\n').split('\t')
                if len(parts) <= max(i_id, i_gene, i_product, i_strain):
                    continue
                gene = parts[i_gene].strip().lower()
                product = parts[i_product].strip().lower()
                if any(t in gene for t in UREC_GENE_TERMS) or any(t in product for t in UREC_PRODUCT_TERMS):
                    strain = parts[i_strain].strip() if i_strain >= 0 else ''
                    wanted[parts[i_id].strip()] = strain
    return wanted


def collect_reference_urec(references_fasta):
    """
    PT-BR: Extrai as sequencias de UreC do FASTA de referencia. O cabecalho tem o formato
           ref|<gene>|UniProt:<acesso>|<organismo>|reviewed|<n>aa, entao o gene e comparado
           por igualdade no 2o campo — nao por substring, que casaria com ureC de outros
           genes so por conterem as letras.
    EN-US: Extracts the UreC sequences from the reference FASTA. The header is formatted
           ref|<gene>|UniProt:<accession>|<organism>|reviewed|<n>aa, so the gene is compared
           by equality on the 2nd field — not by substring, which would match other genes
           merely for containing the letters.
    """
    selected = []
    if not (os.path.exists(references_fasta) and os.path.getsize(references_fasta) > 0):
        return selected

    for rid, (header, sequence) in sorted(read_fasta(references_fasta).items()):
        fields = header.split('|')
        if len(fields) < 4 or fields[1].strip().lower() != 'urec':
            continue
        accession = fields[2].split(':')[-1].strip()
        organism = fields[3].strip()[:MAX_ORGANISM_CHARS]
        selected.append({
            'label_seed': f"REF_{accession}_{organism}",
            'source_id': rid,
            'header': header,
            'sequence': sequence,
            'origin': 'referencia/reference',
        })
    return selected


def sanitise(text, used):
    """
    PT-BR: O IQ-TREE trunca nomes em '|' e recusa nomes repetidos, entao os rotulos sao
           saneados e desambiguados aqui.
    EN-US: IQ-TREE truncates names at '|' and rejects duplicates, so labels are sanitised
           and disambiguated here.
    """
    safe = ''.join(ch if (ch.isalnum() or ch in '_-.') else '_' for ch in text)
    candidate = safe
    suffix = 1
    while candidate in used:
        suffix += 1
        candidate = f"{safe}_{suffix}"
    used.add(candidate)
    return candidate


def run(cmd, label):
    """
    PT-BR: Executa um comando externo e aborta se ele falhar.
    EN-US: Runs an external command and aborts if it fails.
    """
    print(f"  $ {' '.join(cmd)}")
    try:
        proc = subprocess.run(cmd, capture_output=True, text=True)
    except FileNotFoundError:
        print(f"[ERRO/ERROR] '{cmd[0]}' nao encontrado no PATH / not found on PATH.")
        return False
    if proc.returncode != 0:
        print(f"[ERRO/ERROR] {label} falhou / failed (exit {proc.returncode}).")
        if proc.stderr:
            print(proc.stderr.strip()[:2000])
        return False
    return True


def main():
    parser = argparse.ArgumentParser(description="UreC Phylogeny Reconstruction (PT-BR / EN-US)")
    parser.add_argument("--fasta-in", required=True, help="FASTA com todas as proteinas / FASTA with all proteins")
    parser.add_argument("--candidates-dir", required=True, help="Diretorio dos TSVs de candidatos / Candidate TSV directory")
    parser.add_argument("--out-dir", required=True, help="Diretorio de saida / Output directory")
    parser.add_argument("--references", default=None,
                        help="FASTA de referencia; suas sequencias de ureC entram como contexto "
                             "/ Reference FASTA; its ureC sequences enter as context")
    parser.add_argument("--threads", type=int, default=1,
                        help="Threads para MAFFT e IQ-TREE / Threads for MAFFT and IQ-TREE")
    args = parser.parse_args()

    threads = max(1, args.threads)

    os.makedirs(args.out_dir, exist_ok=True)
    status_path = os.path.join(args.out_dir, "phylogeny_status.txt")
    urec_fasta = os.path.join(args.out_dir, "ureC.faa")
    aln_out = os.path.join(args.out_dir, "ureC.aln")
    trim_out = os.path.join(args.out_dir, "ureC_trimmed.aln")
    labels_path = os.path.join(args.out_dir, "ureC_labels.tsv")
    tree_prefix = os.path.join(args.out_dir, "ureC")

    def write_status(text):
        with open(status_path, 'w', encoding='utf-8') as handle:
            handle.write(text + '\n')
        print(text)

    if not os.path.exists(args.fasta_in) or os.path.getsize(args.fasta_in) == 0:
        print(f"[ERRO/ERROR] FASTA de proteinas ausente ou vazio / protein FASTA missing or empty: {args.fasta_in}")
        return 1

    proteins = read_fasta(args.fasta_in)
    urec_ids = collect_urec_ids(args.candidates_dir)
    print(f"Candidatos UreC na triagem / UreC candidates in screening: {len(urec_ids)}")

    entries = []
    missing = []
    for pid in sorted(urec_ids):
        if pid not in proteins:
            missing.append(pid)
            continue
        header, sequence = proteins[pid]
        strain = urec_ids[pid]
        entries.append({
            'label_seed': f"{strain}_{pid}" if strain else pid,
            'source_id': pid,
            'header': header,
            'sequence': sequence,
            'origin': 'estirpe/strain',
        })

    strain_count = len(entries)
    if missing:
        print(f"[AVISO/WARNING] {len(missing)} id(s) de UreC sem sequencia correspondente no FASTA "
              f"/ UreC id(s) with no matching sequence in the FASTA: {', '.join(missing)}")

    reference_count = 0
    if args.references:
        references = collect_reference_urec(args.references)
        reference_count = len(references)
        if reference_count:
            print(f"Sequencias UreC de referencia adicionadas como contexto "
                  f"/ Reference UreC sequences added as context: {reference_count}")
            entries.extend(references)
        else:
            print(f"[AVISO/WARNING] Nenhuma sequencia de ureC encontrada em / No ureC sequence "
                  f"found in: {args.references}")

    total = len(entries)
    composition = (f"{strain_count} de estirpe + {reference_count} de referencia "
                   f"/ {strain_count} from strains + {reference_count} references")

    if total < MIN_TAXA_FOR_TREE:
        write_status(
            f"SEM ARVORE / NO TREE: apenas {total} sequencia(s) de UreC ({composition}); "
            f"IQ-TREE exige no minimo {MIN_TAXA_FOR_TREE}. / only {total} UreC sequence(s); "
            f"IQ-TREE requires at least {MIN_TAXA_FOR_TREE}. Nenhuma arvore foi gerada e "
            f"nenhum dado foi fabricado. / No tree was produced and no data was fabricated."
        )
        return 0

    # PT-BR: O mapeamento e gravado para que cada folha continue rastreavel ate o
    #        feature_id (ou a referencia) de origem.
    # EN-US: The mapping is written so every leaf stays traceable back to its source
    #        feature_id (or reference).
    used = set()
    mapping = []
    with open(urec_fasta, 'w', encoding='utf-8') as handle:
        for entry in entries:
            label = sanitise(entry['label_seed'], used)
            mapping.append((label, entry['origin'], entry['source_id'], entry['header']))
            handle.write(f">{label}\n")
            sequence = entry['sequence']
            for i in range(0, len(sequence), 60):
                handle.write(sequence[i:i + 60] + '\n')

    with open(labels_path, 'w', encoding='utf-8') as handle:
        handle.write("tree_label\torigin\tsource_id\toriginal_header\n")
        for label, origin, source_id, header in mapping:
            handle.write(f"{label}\t{origin}\t{source_id}\t{header}\n")

    print(f"Alinhando {total} sequencias UreC com MAFFT / Aligning {total} UreC sequences with MAFFT...")
    try:
        with open(aln_out, 'w', encoding='utf-8') as out:
            proc = subprocess.run(["mafft", "--auto", "--thread", str(threads), urec_fasta],
                                  stdout=out, stderr=subprocess.PIPE, text=True)
    except FileNotFoundError:
        print("[ERRO/ERROR] 'mafft' nao encontrado no PATH / not found on PATH.")
        return 1
    if proc.returncode != 0 or os.path.getsize(aln_out) == 0:
        print(f"[ERRO/ERROR] MAFFT falhou / failed (exit {proc.returncode}).")
        if proc.stderr:
            print(proc.stderr.strip()[:2000])
        return 1

    if not run(["trimal", "-in", aln_out, "-out", trim_out, "-automated1"], "trimAl"):
        return 1

    # PT-BR: Abaixo de 4 taxons o IQ-TREE recusa o bootstrap ("It makes no sense to perform
    #        bootstrap with less than 4 sequences"), entao a arvore e inferida sem suporte
    #        e o status diz isso explicitamente, em vez de a execucao inteira abortar.
    # EN-US: Below 4 taxa IQ-TREE refuses bootstrap ("It makes no sense to perform bootstrap
    #        with less than 4 sequences"), so the tree is inferred without support and the
    #        status says so explicitly, instead of aborting the whole run.
    bootstrap = total >= MIN_TAXA_FOR_BOOTSTRAP
    iqtree_cmd = ["iqtree", "-s", trim_out, "--prefix", tree_prefix, "-m", "MFP",
                  "-T", str(threads)]
    if bootstrap:
        iqtree_cmd += ["-bb", str(BOOTSTRAP_REPLICATES)]
    else:
        print(f"[AVISO/WARNING] {total} taxons (<{MIN_TAXA_FOR_BOOTSTRAP}): bootstrap omitido "
              f"/ taxa (<{MIN_TAXA_FOR_BOOTSTRAP}): bootstrap omitted.")
    iqtree_cmd.append("-redo")

    if not run(iqtree_cmd, "IQ-TREE 2"):
        return 1

    treefile = f"{tree_prefix}.treefile"
    if not os.path.exists(treefile) or os.path.getsize(treefile) == 0:
        print(f"[ERRO/ERROR] IQ-TREE terminou sem gerar / finished without producing: {treefile}")
        return 1

    if bootstrap:
        write_status(
            f"ARVORE GERADA / TREE BUILT: {total} taxons UreC ({composition}), "
            f"{BOOTSTRAP_REPLICATES} ultrafast bootstraps. / {total} UreC taxa, "
            f"{BOOTSTRAP_REPLICATES} ultrafast bootstraps. Arquivo / File: {treefile}"
        )
    else:
        write_status(
            f"ARVORE SEM SUPORTE / TREE WITHOUT SUPPORT: {total} taxons UreC ({composition}); "
            f"abaixo de {MIN_TAXA_FOR_BOOTSTRAP} taxons o bootstrap nao se aplica e com "
            f"{MIN_TAXA_FOR_TREE} existe uma unica topologia nao-enraizada possivel. "
            f"/ below {MIN_TAXA_FOR_BOOTSTRAP} taxa bootstrap does not apply, and with "
            f"{MIN_TAXA_FOR_TREE} there is a single possible unrooted topology. "
            f"Arquivo / File: {treefile}"
        )
    return 0


if __name__ == "__main__":
    sys.exit(main())

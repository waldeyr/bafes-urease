#!/usr/bin/env python3
"""
extract_urease_locus.py
PT-BR: Recorta do GenBank anotado apenas a(s) regiao(oes) do locus de urease, com um
       numero configuravel de genes flanqueadores de cada lado, para alimentar o clinker.
EN-US: Slices from the annotated GenBank only the urease locus region(s), with a
       configurable number of flanking genes on each side, to feed clinker.

PT-BR: MOTIVO — o clinker faz alinhamento all-vs-all entre TODAS as proteinas de TODOS os
       pares de entradas. Com 10 genomas completos (~46.500 CDS) isso e da ordem de 10^9
       alinhamentos par-a-par: a etapa estourou o limite de 8 h e foi morta (SIGTERM/143).
       O clinker foi feito para clusters genicos, nao para cromossomos inteiros; e um HTML
       com 46.500 genes nao e legivel em navegador nenhum.
EN-US: WHY — clinker runs an all-vs-all alignment between EVERY protein of EVERY pair of
       inputs. With 10 complete genomes (~46,500 CDS) that is on the order of 10^9 pairwise
       alignments: the step blew through its 8 h limit and was killed (SIGTERM/143).
       clinker is built for gene clusters, not whole chromosomes; and an HTML holding
       46,500 genes is not readable in any browser.

PT-BR: Este script NUNCA inventa uma regiao. Se a estirpe nao tem nenhum gene ancora da via
       da urease, ele nao escreve GenBank algum e registra isso no arquivo de status — uma
       estirpe sem urease e um resultado biologico legitimo, nao uma falha.
EN-US: This script NEVER invents a region. If the strain has no urease-pathway anchor gene,
       it writes no GenBank at all and records that in the status file — a strain without
       urease is a legitimate biological outcome, not a failure.
"""

import sys
import os
import csv
import argparse

from Bio import SeqIO

# PT-BR: Perfis Pfam especificos da via da urease (subunidades estruturais, chaperonas de
#        maturacao e transporte de niquel). Um acerto em qualquer um deles ancora uma regiao
#        mesmo com confianca "Media", porque o perfil ja e diagnostico por si so.
# EN-US: Pfam profiles specific to the urease pathway (structural subunits, maturation
#        chaperones and nickel transport). A hit on any of them anchors a region even at
#        "Medium" confidence, because the profile is diagnostic on its own.
#
# PT-BR: Repare que 'cobW' e 'CT_C_D' NAO estao aqui de proposito: sao dominios promiscuos
#        (COG0523 e um dominio C-terminal generico) que aparecem espalhados pelo genoma,
#        longe de qualquer operon de urease. Ancorar neles criaria regioes espurias.
# EN-US: Note that 'cobW' and 'CT_C_D' are deliberately absent: they are promiscuous domains
#        (COG0523 and a generic C-terminal domain) scattered across the genome, far from any
#        urease operon. Anchoring on them would create spurious regions.
UREASE_ANCHOR_PROFILES = {
    "Urease_alpha",
    "Urease_beta",
    "Urease_gamma",
    "UreD",
    "UreE_N",
    "UreE_C",
    "UreF",
    "UreG",
    "NicO",
}

# PT-BR: Nomes de gene que ancoram uma regiao independentemente do perfil Pfam.
# EN-US: Gene names that anchor a region regardless of the Pfam profile.
UREASE_ANCHOR_GENES = {
    "urea", "ureb", "urec", "ured", "uree", "uref", "ureg", "ureh", "urei", "urej",
    "utp", "urta", "urtb", "urtc", "urtd", "urte", "yut", "dur3",
}


def is_high_confidence(value):
    """PT-BR: A coluna de confianca e bilingue ("Alta/High"). / EN-US: The confidence column is bilingual ("Alta/High")."""
    return value.strip().lower().startswith("alta")


def read_anchor_ids(candidates_tsv, include_all):
    """
    PT-BR: Devolve (ancoras, todos) — ids de feature que ancoram regioes, e todos os ids
           candidatos (usados so para relatar o que caiu dentro da janela).
    EN-US: Returns (anchors, all) — feature ids that anchor regions, and every candidate id
           (used only to report what landed inside the window).
    """
    anchors = {}
    every = {}
    with open(candidates_tsv, newline="") as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            feature_id = (row.get("feature_id") or "").strip()
            if not feature_id:
                continue
            gene = (row.get("gene") or "").strip()
            profile = (row.get("hmmer_profile") or "").strip()
            label = gene or profile or feature_id
            every[feature_id] = label

            if include_all:
                anchors[feature_id] = label
                continue
            if is_high_confidence(row.get("confidence") or ""):
                anchors[feature_id] = label
            elif profile in UREASE_ANCHOR_PROFILES:
                anchors[feature_id] = label
            elif gene.lower() in UREASE_ANCHOR_GENES:
                anchors[feature_id] = label
    return anchors, every


def feature_locus_tag(feature):
    tags = feature.qualifiers.get("locus_tag")
    return tags[0] if tags else None


def group_indices(indices, max_gap):
    """
    PT-BR: Agrupa indices de CDS ancora: dois grupos se separam quando ha mais de `max_gap`
           genes entre eles. Ancoras vizinhas (um operon) caem num grupo so.
    EN-US: Groups anchor CDS indices: two groups split when more than `max_gap` genes sit
           between them. Neighbouring anchors (an operon) fall into a single group.
    """
    groups = []
    current = [indices[0]]
    for idx in indices[1:]:
        if idx - current[-1] <= max_gap:
            current.append(idx)
        else:
            groups.append(current)
            current = [idx]
    groups.append(current)
    return groups


def main():
    parser = argparse.ArgumentParser(
        description="Urease locus extraction for synteny (PT-BR / EN-US)"
    )
    parser.add_argument("--strain", required=True, help="Identificador da estirpe / Strain identifier")
    parser.add_argument("--gbff", required=True, help="GenBank anotado pelo Bakta / Bakta-annotated GenBank")
    parser.add_argument("--candidates", required=True, help="TSV de candidatos / Candidates TSV")
    parser.add_argument("--out", required=True, help="GenBank de saida / Output GenBank")
    parser.add_argument("--status", required=True, help="Arquivo de status / Status file")
    parser.add_argument("--flank-genes", type=int, default=5,
                        help="Genes flanqueadores de cada lado / Flanking genes on each side")
    parser.add_argument("--all-candidates", action="store_true",
                        help="Ancora em todo candidato, inclusive dominios promiscuos / "
                             "Anchor on every candidate, promiscuous domains included")
    args = parser.parse_args()

    if args.flank_genes < 0:
        print("[ERRO/ERROR] --flank-genes nao pode ser negativo / cannot be negative.")
        return 1

    anchors, every_candidate = read_anchor_ids(args.candidates, args.all_candidates)

    status_lines = [
        f"strain\t{args.strain}",
        f"candidates_total\t{len(every_candidate)}",
        f"anchors_total\t{len(anchors)}",
        f"flank_genes\t{args.flank_genes}",
    ]

    if not anchors:
        # PT-BR: Sem ancora nao ha locus. Nenhum GenBank e escrito — a estirpe simplesmente
        #        nao entra no mapa de sintenia, e o status diz exatamente por que.
        # EN-US: No anchor means no locus. No GenBank is written — the strain simply does not
        #        enter the synteny map, and the status says exactly why.
        status_lines.append("regions_total\t0")
        status_lines.append("result\tNO_UREASE_ANCHOR")
        status_lines.append(
            "note\tNenhum gene ancora da via da urease / No urease-pathway anchor gene"
        )
        with open(args.status, "w") as fh:
            fh.write("\n".join(status_lines) + "\n")
        print(f"[INFO] {args.strain}: nenhuma ancora de urease / no urease anchor. "
              f"Estirpe fora da sintenia / strain excluded from synteny.")
        return 0

    regions = []
    max_gap = 2 * args.flank_genes

    for record in SeqIO.parse(args.gbff, "genbank"):
        cds = [f for f in record.features if f.type == "CDS"]
        cds.sort(key=lambda f: int(f.location.start))

        hit_idx = [i for i, f in enumerate(cds) if feature_locus_tag(f) in anchors]
        if not hit_idx:
            continue

        for group in group_indices(hit_idx, max_gap):
            lo = max(0, group[0] - args.flank_genes)
            hi = min(len(cds) - 1, group[-1] + args.flank_genes)

            start = int(cds[lo].location.start)
            end = int(cds[hi].location.end)

            sub = record[start:end]
            # PT-BR: O fatiamento do Biopython nao carrega as anotacoes do registro; sem
            #        molecule_type o escritor GenBank recusa a gravacao.
            # EN-US: Biopython slicing does not carry the record annotations over; without
            #        molecule_type the GenBank writer refuses to write.
            sub.annotations["molecule_type"] = record.annotations.get("molecule_type", "DNA")
            sub.id = f"{record.id}_{start + 1}_{end}"
            # PT-BR: O nome vira a linha LOCUS e precisa ser curto. / EN-US: The name becomes
            #        the LOCUS line and has to stay short.
            sub.name = f"{args.strain}_{len(regions) + 1}"
            labels = sorted({anchors[feature_locus_tag(cds[i])] for i in group})
            sub.description = (
                f"{args.strain} urease locus {len(regions) + 1} "
                f"[{record.id}:{start + 1}-{end}] anchors={','.join(labels)}"
            )
            regions.append(sub)

            status_lines.append(
                f"region\t{sub.name}\t{record.id}\t{start + 1}\t{end}\t"
                f"{hi - lo + 1}\t{','.join(labels)}"
            )

    if not regions:
        # PT-BR: Havia ancoras no TSV mas nenhuma casou com um locus_tag do GenBank — isso e
        #        inconsistencia entre as duas saidas do Bakta, nao um resultado biologico.
        # EN-US: The TSV had anchors but none matched a GenBank locus_tag — that is an
        #        inconsistency between Bakta's two outputs, not a biological result.
        print(f"[ERRO/ERROR] {args.strain}: {len(anchors)} ancora(s) no TSV, nenhuma encontrada "
              f"no GenBank / anchors in the TSV, none found in the GenBank.")
        print(f"             Os locus_tag de {args.candidates} nao batem com {args.gbff}.")
        print(f"             The locus_tags in {args.candidates} do not match {args.gbff}.")
        return 1

    with open(args.out, "w") as fh:
        SeqIO.write(regions, fh, "genbank")

    total_genes = sum(len([f for f in r.features if f.type == "CDS"]) for r in regions)
    status_lines.append(f"regions_total\t{len(regions)}")
    status_lines.append(f"genes_total\t{total_genes}")
    status_lines.append("result\tOK")
    with open(args.status, "w") as fh:
        fh.write("\n".join(status_lines) + "\n")

    print(f"[OK] {args.strain}: {len(regions)} regiao(oes) / region(s), "
          f"{total_genes} genes -> {args.out}")
    return 0


if __name__ == "__main__":
    sys.exit(main())

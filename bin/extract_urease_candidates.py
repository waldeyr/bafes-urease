#!/usr/bin/env python3
"""
extract_urease_candidates.py
PT-BR: Consolidacao de evidencias em 3 camadas (Bakta, BLASTp, Pfam HMMER) por estirpe.
EN-US: 3-Layer evidence consolidation (Bakta, BLASTp, Pfam HMMER) per strain.

PT-BR: Os limiares abaixo sao os documentados no README. Antes eles nao eram aplicados:
       todo acerto de BLAST contava como evidencia, o que inflava artificialmente a
       confianca dos candidatos.
EN-US: The thresholds below are the ones documented in the README. They used not to be
       applied: every BLAST hit counted as evidence, which artificially inflated
       candidate confidence.

PT-BR: A Camada C tinha o mesmo problema, em escala maior: o perfil Pfam era lido mas
       nunca filtrado, entao qualquer proteina com qualquer dominio conhecido do Pfam-A
       entrava na triagem (3.836 das 3.854 linhas de S1 eram acertos sem relacao alguma
       com urease). Agora so contam os perfis de UREASE_PFAM_ACCESSIONS.
EN-US: Layer C had the same problem on a larger scale: the Pfam profile was read but never
       filtered, so every protein carrying any known Pfam-A domain entered the screening
       (3,836 of S1's 3,854 rows were hits with no relation to urease). Only the profiles
       in UREASE_PFAM_ACCESSIONS count now.

PT-BR: Pontuacao — cada camada vale 1 ponto, exceto o BLASTp cujo melhor acerto e uma
       referencia NAO revisada (TrEMBL): ela e registrada (blast_ref_status = unreviewed) mas
       vale 0, porque a anotacao do TrEMBL e automatica como a do Bakta e nao e evidencia
       independente. >= 2 pontos = Alta; 1 = Media; 0 com acerto TrEMBL = Baixa.
EN-US: Scoring — each layer is worth 1 point, except a BLASTp whose best hit is an UNREVIEWED
       reference (TrEMBL): it is recorded (blast_ref_status = unreviewed) but is worth 0,
       because TrEMBL annotation is automatic like Bakta's and is not independent evidence.
       >= 2 points = High; 1 = Medium; 0 with a TrEMBL hit = Low.

PT-BR: Filtros antes da pontuacao, contra o ruido das referencias TrEMBL — UrtD/UrtE sao a
       ATPase de um transportador ABC, e a 30 % de identidade casavam com toda ATPase ABC do
       genoma: (1) acerto contra TrEMBL so vale com identidade >= --trembl-min-identity;
       (2) UrtB-E so valem a ate --urt-operon-window CDS de uma UrtA no mesmo contig. O que
       cai fica registrado em --filtered-tsv, com o motivo.
EN-US: Filters before scoring, against the TrEMBL reference noise — UrtD/UrtE are the ATPase
       of an ABC transporter, and at 30 % identity they matched every ABC ATPase in the
       genome: (1) a hit on TrEMBL only counts at identity >= --trembl-min-identity;
       (2) UrtB-E only count within --urt-operon-window CDS of a UrtA on the same contig.
       Whatever is dropped is recorded in --filtered-tsv, with the reason.

PT-BR: Familia KipI/PxpABC (5-oxoprolinase) — o perfil CT_C_D (PF02682, antigo AHS1,
       "allophanate hydrolase subunit 1") NAO marca alofanato hidrolase: a AtzF verdadeira e
       da familia amidase (PF01425). CT_C_D aparece sozinho na PxpB/KipI (~250 aa), presente
       em toda estirpe, e junto da biotina carboxilase na urea carboxilase (~1.200 aa). O
       Bakta ainda da as PxpB/PxpC o produto "Allophanate hydrolase" (nome do cluster
       UniRef), o que marcou ah=1 falso em S5. Por isso: (3) genes pxpA/B/C e kipI/A nao
       entram na Camada A (motivo kipI_family); (4) CT_C_D so pontua na Camada C se a mesma
       proteina tiver dominio de biotina carboxilase (motivo ct_c_d_without_bc).
EN-US: KipI/PxpABC family (5-oxoprolinase) — the CT_C_D profile (PF02682, formerly AHS1,
       "allophanate hydrolase subunit 1") does NOT mark allophanate hydrolase: the true AtzF
       belongs to the amidase family (PF01425). CT_C_D sits alone in PxpB/KipI (~250 aa),
       present in every strain, and next to biotin carboxylase in urea carboxylase
       (~1,200 aa). Bakta also gives PxpB/PxpC the product "Allophanate hydrolase" (the
       UniRef cluster name), which scored a false ah=1 in S5. Hence: (3) pxpA/B/C and kipI/A
       genes do not enter Layer A (reason kipI_family); (4) CT_C_D only scores in Layer C if
       the same protein carries a biotin carboxylase domain (reason ct_c_d_without_bc).
"""

import sys
import os
import argparse
import json
import pandas as pd

# PT-BR: Termos EC e genes alvo / EN-US: Target EC numbers and gene symbols
TARGET_EC = ["3.5.1.5", "6.3.4.6", "3.5.1.54"]
TARGET_GENES = ["ureC", "ureA", "ureB", "ureD", "ureE", "ureF", "ureG", "ureH",
                "urtA", "urtB", "urtC", "urtD", "urtE", "uc", "ah"]

# PT-BR: Limiares de evidencia / EN-US: Evidence thresholds
BLAST_MAX_EVALUE = 1e-5
BLAST_MIN_IDENTITY = 30.0
BLAST_MIN_COVERAGE = 50.0
HMMER_MIN_DOMAIN_SCORE = 20.0

# PT-BR: Padroes dos filtros contra o ruido do TrEMBL (ver docstring).
# EN-US: Defaults of the filters against TrEMBL noise (see docstring).
TREMBL_MIN_IDENTITY = 80.0
URT_OPERON_WINDOW = 5
URT_ANCHOR = "urta"
URT_PARTNERS = {"urtb", "urtc", "urtd", "urte"}

# PT-BR: Familia KipI/PxpABC (5-oxoprolinase), fora da via alternativa (ver docstring).
# EN-US: KipI/PxpABC family (5-oxoprolinase), outside the alternative route (see docstring).
KIPI_FAMILY_GENES = {"pxpa", "pxpb", "pxpc", "kipi", "kipa"}

# PT-BR: Perfis Pfam que constituem evidencia de urease (Camada C). Acessoes sem versao.
#        PF01979 (Amidohydro_1) foi deixado de fora de proposito: e o acerto de maior score
#        do ureC (232.2 em S1), mas cobre a superfamilia inteira de amido-hidrolases; o
#        PF00449 marca a mesma proteina com 179.0, entao nada se perde e o ruido cai.
#        Os transportadores urt* nao entram aqui — PF00005 (ABC_tran) e largo demais —
#        e seguem cobertos pelas Camadas A (Bakta) e B (BLASTp, contra as UrtA-E de
#        Firmicutes nao revisadas que o fetch_references.py baixa).
# EN-US: Pfam profiles that constitute urease evidence (Layer C). Accessions without version.
#        PF01979 (Amidohydro_1) is deliberately excluded: it is ureC's top-scoring hit
#        (232.2 in S1) but spans the whole amidohydrolase superfamily; PF00449 marks the
#        same protein at 179.0, so nothing is lost and the noise drops. The urt*
#        transporters are not listed here — PF00005 (ABC_tran) is far too broad — and stay
#        covered by Layers A (Bakta) and B (BLASTp, against the unreviewed Firmicutes UrtA-E
#        that fetch_references.py downloads).
UREASE_PFAM_ACCESSIONS = {
    "PF00449": "Urease_alpha (ureC)",
    "PF00699": "Urease_beta (ureB)",
    "PF00547": "Urease_gamma (ureA)",
    "PF01774": "UreD",
    "PF02814": "UreE_N",
    "PF05194": "UreE_C",
    "PF01730": "UreF",
    "PF02492": "cobW (UreG/COG0523)",
    "PF03824": "NicO (UreH, transporte de niquel / nickel transport)",
    "PF02682": "CT_C_D (carboxiltransferase; so com biotina carboxilase = urea carboxilase "
               "/ carboxyltransferase; only with biotin carboxylase = urea carboxylase)",
}

# PT-BR: CT_C_D so vale com um destes dominios de biotina carboxilase na mesma proteina —
#        e a arquitetura da urea carboxilase. Sozinho, CT_C_D e a PxpB/KipI.
# EN-US: CT_C_D only counts with one of these biotin carboxylase domains on the same
#        protein — the urea carboxylase architecture. Alone, CT_C_D is PxpB/KipI.
CT_C_D_ACCESSION = "PF02682"
BIOTIN_CARBOXYLASE_PFAM = {"PF00289", "PF02786", "PF02785"}
BIOTIN_CARBOXYLASE_NAMES = {"biotin_carb_n", "cpsase_l_d2", "biotin_carb_c"}

# PT-BR: Reserva para o caso de a coluna de acessao vir como '-' no domtblout.
# EN-US: Fallback for when the accession column comes as '-' in the domtblout.
UREASE_PFAM_NAMES = {"urease_alpha", "urease_beta", "urease_gamma", "ured", "uree_n",
                     "uree_c", "uref", "cobw", "nico", "ct_c_d"}

# PT-BR: Status do 5o campo do cabecalho escrito pelo fetch_references.py
#        (ref|gene|UniProt:acesso|organismo|status|Naa). Um FASTA proprio fora desse formato
#        fica "unknown" e pontua como revisado — o comportamento de antes desta regra.
# EN-US: Status from the 5th field of the header written by fetch_references.py
#        (ref|gene|UniProt:accession|organism|status|Naa). A custom FASTA outside that format
#        is "unknown" and scores as reviewed — the behaviour from before this rule.
REF_STATUS_UNREVIEWED = "unreviewed"
REF_STATUS_UNKNOWN = "unknown"

# PT-BR: Colunas do -outfmt 6 usado no main.nf / EN-US: Columns of the -outfmt 6 used in main.nf
BLAST_COLUMNS = ["qseqid", "sseqid", "pident", "length", "qstart", "qend",
                 "qlen", "slen", "evalue", "bitscore"]


def load_bakta_features(json_file):
    """PT-BR: Features do JSON do Bakta. / EN-US: Features from the Bakta JSON."""
    if not os.path.exists(json_file):
        return []
    with open(json_file, 'r', encoding='utf-8') as f:
        return json.load(f).get('features', [])


def cds_positions(features):
    """
    PT-BR: Ordem genica: {locus ou id: (contig, indice da CDS no contig)}, com as CDS de
           cada contig ordenadas pela posicao inicial.
    EN-US: Gene order: {locus or id: (contig, index of the CDS on the contig)}, with each
           contig's CDS sorted by start position.
    """
    by_contig = {}
    for feat in features:
        if feat.get('type') == 'cds':
            by_contig.setdefault(feat.get('sequence') or '', []).append(feat)
    positions = {}
    for contig, feats in by_contig.items():
        for index, feat in enumerate(sorted(feats, key=lambda f: f.get('start') or 0)):
            for key in (feat.get('locus'), feat.get('id')):
                if key:
                    positions[key] = (contig, index)
    return positions


def parse_bakta_json(features, rejected=None):
    """
    PT-BR: Filtra das features do Bakta os genes e termos EC alvos (Camada A). Genes da
           familia KipI/PxpABC saem mesmo com "allophanate" no produto e vao para `rejected`.
    EN-US: Filters target genes and EC numbers out of the Bakta features (Layer A). KipI/PxpABC
           family genes are dropped even with "allophanate" in the product, into `rejected`.
    """
    candidates = []
    for feat in features:
        gene = feat.get('gene') or ''
        product = feat.get('product') or ''
        dbxrefs = feat.get('db_xrefs') or feat.get('dbxrefs') or []
        # PT-BR: O Bakta usa 'locus' no cabecalho do .faa; guardamos os dois para casar
        #        com os IDs vindos do BLASTp e do HMMER.
        # EN-US: Bakta uses 'locus' in the .faa header; we keep both so they match the
        #        IDs coming from BLASTp and HMMER.
        feat_id = feat.get('locus') or feat.get('id') or ''
        alt_id = feat.get('id') or ''

        is_match = False
        if any(g.lower() == str(gene).lower() for g in TARGET_GENES):
            is_match = True
        elif any(ec in str(dbxrefs) for ec in TARGET_EC):
            is_match = True
        elif any(term in product.lower() for term in ("urease", "urea", "allophanate")):
            is_match = True

        if is_match and str(gene).lower() in KIPI_FAMILY_GENES:
            if rejected is not None:
                rejected.append({'feature_id': feat_id, 'layer': 'bakta',
                                 'reference': f"{gene}: {product}", 'identity': '',
                                 'reason': 'kipI_family'})
            continue

        if is_match:
            candidates.append({
                'feature_id': feat_id,
                'alt_id': alt_id,
                'gene': gene or 'unnamed',
                'product': product,
                'bakta_hit': True
            })
    return candidates


def reference_gene(sseqid):
    """PT-BR: Gene da referencia (2o campo). / EN-US: Reference gene (2nd field)."""
    fields = str(sseqid).split('|')
    return fields[1] if len(fields) >= 6 and fields[0] == 'ref' else ''


def reference_status(sseqid):
    """PT-BR: Status de revisao da referencia. / EN-US: Review status of the reference."""
    fields = str(sseqid).split('|')
    if len(fields) >= 6 and fields[0] == 'ref' and fields[4] in ('reviewed', REF_STATUS_UNREVIEWED):
        return fields[4]
    return REF_STATUS_UNKNOWN


def parse_blast_tsv(blast_file, trembl_min_identity=TREMBL_MIN_IDENTITY, rejected=None):
    """
    PT-BR: Extrai IDs de consulta do BLASTp que passam nos limiares (Camada B). Acertos
           contra TrEMBL abaixo de trembl_min_identity saem; quando a proteina fica sem
           acerto algum, o melhor deles vai para `rejected`.
    EN-US: Extracts BLASTp query IDs that pass the thresholds (Layer B). Hits on TrEMBL
           below trembl_min_identity are dropped; when the protein is left with no hit at
           all, the best of them goes into `rejected`.
    """
    hits = {}
    if not (os.path.exists(blast_file) and os.path.getsize(blast_file) > 0):
        return hits

    df = pd.read_csv(blast_file, sep='\t', header=None, names=BLAST_COLUMNS)
    if df.empty:
        return hits

    # PT-BR: Cobertura da consulta coberta pelo alinhamento / EN-US: Query coverage by the alignment
    df['coverage'] = (df['qend'] - df['qstart'] + 1) / df['qlen'] * 100.0

    kept = df[(df['evalue'] <= BLAST_MAX_EVALUE) &
              (df['pident'] >= BLAST_MIN_IDENTITY) &
              (df['coverage'] >= BLAST_MIN_COVERAGE)]

    dropped = len(df) - len(kept)
    if dropped:
        print(f"  [FILTRO/FILTER] BLASTp: {dropped} acerto(s) descartado(s) pelos limiares "
              f"(E<={BLAST_MAX_EVALUE}, id>={BLAST_MIN_IDENTITY}%, cov>={BLAST_MIN_COVERAGE}%) "
              f"/ hit(s) dropped by thresholds.")

    unreviewed = kept['sseqid'].map(reference_status) == REF_STATUS_UNREVIEWED
    weak_trembl = kept[unreviewed & (kept['pident'] < trembl_min_identity)]
    kept = kept[~(unreviewed & (kept['pident'] < trembl_min_identity))]
    if len(weak_trembl):
        print(f"  [FILTRO/FILTER] BLASTp: {len(weak_trembl)} acerto(s) contra TrEMBL abaixo de "
              f"{trembl_min_identity:g}% de identidade descartado(s) / TrEMBL hit(s) below "
              f"{trembl_min_identity:g}% identity dropped.")
        if rejected is not None:
            surviving = set(kept['qseqid'].astype(str))
            seen = set()
            for _, row in weak_trembl.sort_values('bitscore', ascending=False).iterrows():
                qid = str(row['qseqid'])
                if qid in surviving or qid in seen:
                    continue
                seen.add(qid)
                rejected.append({'feature_id': qid, 'layer': 'blast',
                                 'reference': str(row['sseqid']),
                                 'identity': round(float(row['pident']), 1),
                                 'reason': 'trembl_identity'})

    # PT-BR: Guarda o melhor acerto por consulta / EN-US: Keep the best hit per query
    for _, row in kept.sort_values('bitscore', ascending=False).iterrows():
        qid = str(row['qseqid'])
        if qid not in hits:
            hits[qid] = {
                'subject': str(row['sseqid']),
                'ref_gene': reference_gene(row['sseqid']),
                'status': reference_status(row['sseqid']),
                'identity': float(row['pident']),
                'coverage': float(row['coverage']),
                'evalue': float(row['evalue']),
            }
    return hits


def parse_hmmer_tbl(hmmer_file, rejected=None):
    """
    PT-BR: Extrai IDs de consulta do HMMER3 domtblout acima do score de dominio (Camada C).
    EN-US: Extracts HMMER3 domtblout query IDs above the domain score (Layer C).

    PT-BR: No domtblout do hmmscan a consulta (proteina) e a coluna 4, a acessao Pfam do
           perfil e a coluna 2 e o score do dominio e a coluna 14. O codigo anterior lia a
           coluna 1, que e o nome do perfil Pfam — por isso a Camada C nunca casava com nada.
    EN-US: In hmmscan's domtblout the query (protein) is column 4, the profile's Pfam
           accession is column 2 and the domain score is column 14. The previous code read
           column 1, the Pfam profile name — which is why Layer C never matched anything.

    PT-BR: CT_C_D so entra se a proteina tambem tiver biotina carboxilase (ver docstring do
           modulo); senao vai para `rejected`.
    EN-US: CT_C_D only enters if the protein also carries biotin carboxylase (see the module
           docstring); otherwise it goes into `rejected`.
    """
    hits = {}
    ct_c_d_hits = {}
    with_biotin_carboxylase = set()
    if not (os.path.exists(hmmer_file) and os.path.getsize(hmmer_file) > 0):
        return hits

    below = 0
    off_target = 0
    with open(hmmer_file, 'r', encoding='utf-8') as f:
        for line in f:
            if line.startswith('#'):
                continue
            parts = line.split()
            if len(parts) < 14:
                continue
            profile = parts[0]
            accession = parts[1].split('.')[0].upper()
            query = parts[3]
            try:
                domain_score = float(parts[13])
            except ValueError:
                continue

            # PT-BR: A biotina carboxilase esta fora da lista, entao e lida antes do filtro.
            # EN-US: Biotin carboxylase is off the list, so it is read before the filter.
            if domain_score >= HMMER_MIN_DOMAIN_SCORE and (
                    accession in BIOTIN_CARBOXYLASE_PFAM
                    or profile.lower() in BIOTIN_CARBOXYLASE_NAMES):
                with_biotin_carboxylase.add(query)

            # PT-BR: So dominios de urease contam como evidencia — o hmmscan roda contra o
            #        Pfam-A inteiro, entao sem este filtro a triagem admite o proteoma todo.
            # EN-US: Only urease domains count as evidence — hmmscan runs against the whole
            #        Pfam-A, so without this filter the screening admits the entire proteome.
            if accession not in UREASE_PFAM_ACCESSIONS and profile.lower() not in UREASE_PFAM_NAMES:
                off_target += 1
                continue

            if domain_score < HMMER_MIN_DOMAIN_SCORE:
                below += 1
                continue

            target = ct_c_d_hits if (accession == CT_C_D_ACCESSION
                                     or profile.lower() == "ct_c_d") else hits
            if query not in target or domain_score > target[query]['score']:
                target[query] = {'profile': profile, 'score': domain_score}

    lone_ct_c_d = 0
    for query, hit in sorted(ct_c_d_hits.items()):
        if query not in with_biotin_carboxylase:
            lone_ct_c_d += 1
            if rejected is not None:
                rejected.append({'feature_id': query, 'layer': 'hmmer',
                                 'reference': f"{hit['profile']} ({hit['score']:g})",
                                 'identity': '', 'reason': 'ct_c_d_without_bc'})
            continue
        if query not in hits or hit['score'] > hits[query]['score']:
            hits[query] = hit

    if off_target:
        print(f"  [FILTRO/FILTER] HMMER: {off_target} dominio(s) de perfis fora da lista de "
              f"urease descartado(s) / domain(s) from non-urease profiles dropped.")
    if below:
        print(f"  [FILTRO/FILTER] HMMER: {below} dominio(s) de urease abaixo do score "
              f"{HMMER_MIN_DOMAIN_SCORE} descartado(s) / urease domain(s) below score dropped.")
    if lone_ct_c_d:
        print(f"  [FILTRO/FILTER] HMMER: {lone_ct_c_d} proteina(s) com CT_C_D sem biotina "
              f"carboxilase (PxpB/KipI) descartada(s) / protein(s) with CT_C_D but no biotin "
              f"carboxylase (PxpB/KipI) dropped.")
    return hits


def main():
    parser = argparse.ArgumentParser(description="Extract and Consolidate Urease Candidates (PT-BR / EN-US)")
    parser.add_argument("--strain", required=True, help="Identificador da estirpe / Strain ID")
    parser.add_argument("--bakta-json", required=True, help="Bakta JSON output")
    parser.add_argument("--blast-tsv", required=True, help="BLASTp TSV output")
    parser.add_argument("--hmmer-tbl", required=True, help="HMMER domtblout")
    parser.add_argument("--out-tsv", required=True, help="Caminho do TSV de saida / Output TSV path")
    parser.add_argument("--filtered-tsv", default=None,
                        help="TSV com a evidencia descartada pelos filtros / TSV of the evidence "
                             "dropped by the filters (default: <out-tsv sem / without "
                             "_candidates.tsv>_filtered.tsv)")
    parser.add_argument("--trembl-min-identity", type=float, default=TREMBL_MIN_IDENTITY,
                        help="Identidade minima (%%) para acertos contra TrEMBL / Minimum identity "
                             "(%%) for hits on TrEMBL")
    parser.add_argument("--urt-operon-window", type=int, default=URT_OPERON_WINDOW,
                        help="Distancia maxima, em CDS, de UrtB-E ate uma UrtA / Maximum distance, "
                             "in CDS, from UrtB-E to a UrtA")
    args = parser.parse_args()

    features = load_bakta_features(args.bakta_json)
    rejected = []
    bakta_cands = parse_bakta_json(features, rejected)
    blast_hits = parse_blast_tsv(args.blast_tsv, args.trembl_min_identity, rejected)
    hmmer_hits = parse_hmmer_tbl(args.hmmer_tbl, rejected)

    # PT-BR: Indexa os candidatos do Bakta pelos dois identificadores possiveis — o lookup
    #        por 'alt_id' e defensivo, caso BLASTp ou HMMER tragam o id em vez do locus.
    # EN-US: Index Bakta candidates by both possible identifiers — the 'alt_id' lookup is
    #        defensive, in case BLASTp or HMMER bring the id instead of the locus.
    by_id = {}
    for cand in bakta_cands:
        for key in (cand['feature_id'], cand['alt_id']):
            if key:
                by_id[key] = cand

    # PT-BR: So os ids canonicos geram linha. Antes o conjunto era a uniao das chaves de
    #        'by_id', o que emitia DUAS linhas por feature do Bakta sempre que 'locus' e
    #        'id' diferiam (ex.: MHBDHA_02422 e MHBDHAAFLD_2601 sao o mesmo ureC de S1).
    #        Isso duplicava a contagem de candidatos e fazia a filogenia procurar no .faa
    #        um id que o Bakta nunca escreve la.
    # EN-US: Only canonical ids produce a row. The set used to be the union of 'by_id' keys,
    #        which emitted TWO rows per Bakta feature whenever 'locus' and 'id' differed
    #        (e.g. MHBDHA_02422 and MHBDHAAFLD_2601 are the same ureC of S1). That doubled
    #        the candidate count and made the phylogeny look in the .faa for an id Bakta
    #        never writes there.
    canonical_ids = {c['feature_id'] for c in bakta_cands if c['feature_id']}
    all_feat_ids = canonical_ids | set(blast_hits.keys()) | set(hmmer_hits.keys())

    # PT-BR: Contexto de operon: UrtB-E so valem perto de uma UrtA (Bakta ou BLASTp).
    # EN-US: Operon context: UrtB-E only count near a UrtA (Bakta or BLASTp).
    positions = cds_positions(features)

    def urt_roles(fid):
        b_info = by_id.get(fid)
        blast = blast_hits.get(fid)
        return ((b_info['gene'].lower() if b_info else ''),
                (blast['ref_gene'].lower() if blast else ''))

    anchors = {}
    for fid in all_feat_ids:
        if URT_ANCHOR in urt_roles(fid) and fid in positions:
            contig, index = positions[fid]
            anchors.setdefault(contig, []).append(index)

    def near_urta(fid):
        if fid not in positions:
            return False
        contig, index = positions[fid]
        return any(abs(index - a) <= args.urt_operon_window for a in anchors.get(contig, []))

    dropped_bakta = set()
    for fid in sorted(all_feat_ids):
        bakta_role, blast_role = urt_roles(fid)
        if not (bakta_role in URT_PARTNERS or blast_role in URT_PARTNERS) or near_urta(fid):
            continue
        if bakta_role in URT_PARTNERS:
            dropped_bakta.add(fid)
            rejected.append({'feature_id': fid, 'layer': 'bakta',
                             'reference': f"{by_id[fid]['gene']}: {by_id[fid]['product']}",
                             'identity': '', 'reason': 'urt_without_urtA'})
        if blast_role in URT_PARTNERS:
            blast = blast_hits.pop(fid)
            rejected.append({'feature_id': fid, 'layer': 'blast', 'reference': blast['subject'],
                             'identity': round(blast['identity'], 1),
                             'reason': 'urt_without_urtA'})

    results = []
    for fid in sorted(all_feat_ids):
        b_info = by_id.get(fid)
        in_bakta = b_info is not None and fid not in dropped_bakta
        blast = blast_hits.get(fid)
        hmmer = hmmer_hits.get(fid)

        # PT-BR: Pontuacao (ver docstring): o BLASTp contra TrEMBL nao pontua.
        # EN-US: Scoring (see docstring): BLASTp against TrEMBL scores nothing.
        blast_scores = blast is not None and blast['status'] != REF_STATUS_UNREVIEWED
        points = sum([in_bakta, blast_scores, hmmer is not None])
        if points >= 2:
            confidence = "Alta/High"
        elif points == 1:
            confidence = "Media/Medium"
        elif blast is not None:
            confidence = "Baixa/Low"
        else:
            continue  # PT-BR: sem evidencia nenhuma / EN-US: no evidence at all

        # PT-BR: O identificador alternativo do Bakta vira coluna, nao linha propria.
        # EN-US: Bakta's alternate identifier becomes a column, not a row of its own.
        alt_id = b_info['alt_id'] if b_info else ''
        if alt_id == fid:
            alt_id = ''

        results.append({
            'strain': args.strain,
            'feature_id': fid,
            'bakta_alt_id': alt_id,
            'gene': b_info['gene'] if b_info else '',
            'product': b_info['product'] if b_info else '',
            'layer_bakta': in_bakta,
            'layer_blast': blast is not None,
            'layer_hmmer': hmmer is not None,
            'blast_subject': blast['subject'] if blast else '',
            'blast_identity': round(blast['identity'], 1) if blast else '',
            'blast_coverage': round(blast['coverage'], 1) if blast else '',
            'blast_evalue': blast['evalue'] if blast else '',
            'blast_ref_status': blast['status'] if blast else '',
            'hmmer_profile': hmmer['profile'] if hmmer else '',
            'hmmer_score': hmmer['score'] if hmmer else '',
            'evidence_points': points,
            'confidence': confidence,
        })

    columns = ['strain', 'feature_id', 'bakta_alt_id', 'gene', 'product', 'layer_bakta',
               'layer_blast', 'layer_hmmer', 'blast_subject', 'blast_identity',
               'blast_coverage', 'blast_evalue', 'blast_ref_status', 'hmmer_profile',
               'hmmer_score', 'evidence_points', 'confidence']

    df_out = pd.DataFrame(results, columns=columns)
    df_out.to_csv(args.out_tsv, sep='\t', index=False)

    filtered_tsv = args.filtered_tsv or (
        args.out_tsv[:-len('_candidates.tsv')] + '_filtered.tsv'
        if args.out_tsv.endswith('_candidates.tsv') else args.out_tsv + '.filtered.tsv')
    pd.DataFrame([dict(strain=args.strain, **r) for r in rejected],
                 columns=['strain', 'feature_id', 'layer', 'reference', 'identity', 'reason']
                 ).to_csv(filtered_tsv, sep='\t', index=False)

    high = sum(1 for r in results if r['confidence'] == "Alta/High")
    low = sum(1 for r in results if r['confidence'] == "Baixa/Low")
    print(f"[OK] Candidatos consolidados salvos em / Consolidated candidates saved to: {args.out_tsv} "
          f"({len(df_out)} candidatos, {high} de alta confianca, {low} de baixa (so TrEMBL) "
          f"/ candidates, {high} high confidence, {low} low (TrEMBL only))")
    print(f"[OK] Evidencia descartada pelos filtros / Evidence dropped by the filters: "
          f"{len(rejected)} -> {filtered_tsv}")
    return 0


if __name__ == "__main__":
    sys.exit(main())

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
generate_report.py
PT-BR: Gera um relatório Web único e autocontido a partir da árvore de resultados do
       pipeline bafes-urease, mostrando o caminho percorrido pelos dados e os
       principais resultados em gráficos SVG desenhados sem dependências externas.
EN-US: Builds a single self-contained web report from the bafes-urease results tree,
       showing the data journey and the main results as SVG charts drawn without any
       external dependency.

PT-BR: Somente a biblioteca padrao do Python 3.9+ e usada. O script roda no host, fora
       do container, sobre resultados já copiados, e nunca carrega arquivos grandes na
       memória (as tabelas HMMER, os JSON/GBFF do Bakta e os genomas são lidos por
       streaming apenas para contagem, ou não são lidos).
EN-US: Python 3.9+ standard library only. The script runs on the host, outside the
       container, over results that were already copied, and never loads large files
       into memory (HMMER tables, Bakta JSON/GBFF and genomes are streamed for counting
       only, or not read at all).

PT-BR: A taxonomia das estirpes NAO vêm da coluna 'species' de data/accessions.tsv, que
       e um placeholder ('Bacillus sp. Sn' em 9 das 10 linhas). Ela e derivada da linha
       de definição do RefSeq no próprio genoma baixado, que traz gênero, espécie e a
       designação de estirpe do depositante.
EN-US: Strain taxonomy does NOT come from the 'species' column of data/accessions.tsv,
       which is a placeholder ('Bacillus sp. Sn' on 9 of the 10 rows). It is derived from
       the RefSeq definition line of the downloaded genome itself, which carries genus,
       species and the depositor's strain designation.
"""

import argparse
import csv
import glob
import html as html_mod
import math
import os
import re
import sys
from collections import Counter, OrderedDict

# ---------------------------------------------------------------------------
# PT-BR: Constantes de domínio / EN-US: Domain constants
# ---------------------------------------------------------------------------

# PT-BR: Genes da matriz, na ordem das colunas de merge_strain_results.py
# EN-US: Matrix genes, in the column order of merge_strain_results.py
UREASE_GENES = ['ureA', 'ureB', 'ureC', 'ureD', 'ureE', 'ureF', 'ureG', 'ureH',
                'urtA', 'uc', 'ah']

GENE_TOOLTIP = {
    'ureA': ('subunidade gama da urease', 'urease subunit gamma'),
    'ureB': ('subunidade beta da urease', 'urease subunit beta'),
    'ureC': ('subunidade alfa da urease (sítio catalítico)', 'urease subunit alpha (catalytic site)'),
    'ureD': ('proteína acessória UreD', 'accessory protein UreD'),
    'ureE': ('chaperona de níquel UreE', 'nickel chaperone UreE'),
    'ureF': ('proteína acessória UreF', 'accessory protein UreF'),
    'ureG': ('GTPase acessória UreG', 'accessory GTPase UreG'),
    'ureH': ('proteína acessória UreH', 'accessory protein UreH'),
    'urtA': ('transportador de ureia', 'urea transporter'),
    'uc': ('urea carboxilase (via alternativa)', 'urea carboxylase (alternative route)'),
    'ah': ('alofanato hidrolase (via alternativa)', 'allophanate hydrolase (alternative route)'),
}

# PT-BR: Padroes de produto por gene, identicos aos de merge_strain_results.py, para que
#        o destaque nos mapas de sintenia use exatamente o mesmo critério da matriz.
# EN-US: Per-gene product patterns, identical to merge_strain_results.py, so the synteny
#        map highlighting uses exactly the same criterion as the matrix.
PRODUCT_PATTERNS = {
    'ureA': [r'urease\s+subunit\s+gamma', r'urease\s+gamma\s+subunit'],
    'ureB': [r'urease\s+subunit\s+beta', r'urease\s+beta\s+subunit'],
    'ureC': [r'urease\s+subunit\s+alpha', r'urease\s+alpha\s+subunit',
             r'urea\s+amidohydrolase\s+subunit\s+alpha'],
    'ureD': [r'urease\s+accessory\s+protein\s+ured'],
    'ureE': [r'urease\s+accessory\s+protein\s+uree'],
    'ureF': [r'urease\s+accessory\s+protein\s+uref'],
    'ureG': [r'urease\s+accessory\s+protein\s+ureg'],
    'ureH': [r'urease\s+accessory\s+protein\s+ureh'],
    'urtA': [r'urea\s+(abc\s+)?transport', r'urea\s+transporter'],
    'uc': [r'urea\s+carboxylase'],
    'ah': [r'allophanate\s+hydrolase'],
}
COMPILED_PATTERNS = {g: [re.compile(p, re.IGNORECASE) for p in pats]
                     for g, pats in PRODUCT_PATTERNS.items()}

# PT-BR: Perfis Pfam curados da camada HMMER (extract_urease_candidates.py)
# EN-US: Curated Pfam profiles of the HMMER layer (extract_urease_candidates.py)
UREASE_PFAM = {
    'PF00449': 'Urease_alpha', 'PF00699': 'Urease_beta', 'PF00547': 'Urease_gamma',
    'PF01774': 'UreD', 'PF02814': 'UreE_N', 'PF05194': 'UreE_C', 'PF01730': 'UreF',
    'PF02492': 'cobW', 'PF03824': 'NicO', 'PF02682': 'CT_C_D',
}
CORE_UREASE_PROFILES = {'Urease_alpha', 'Urease_beta', 'Urease_gamma',
                        'UreD', 'UreE_N', 'UreE_C', 'UreF'}

PIPELINE_ORDER = ['PREFLIGHT_CHECK', 'DOWNLOAD_GENOME', 'QC_QUAST', 'QC_CHECKM2',
                  'BAKTA_ANNOTATE', 'EXTRACT_CANDIDATES', 'MERGE_ALL_STRAINS',
                  'EXTRACT_LOCUS', 'ANALYZE_SYNTENY', 'BUILD_PHYLOGENY']

NBSP_THIN = ' '  # PT-BR: separador de milhar neutro / EN-US: locale-neutral thousands separator

LANG = 'pt'  # PT-BR: idioma inicial, definido pela CLI / EN-US: initial language, set by the CLI


# ---------------------------------------------------------------------------
# PT-BR: Utilitarios gerais / EN-US: General utilities
# ---------------------------------------------------------------------------

def esc(value):
    """PT-BR: Escapa texto para HTML/SVG. EN-US: Escapes text for HTML/SVG."""
    return html_mod.escape('' if value is None else str(value), quote=True)


def natural_key(name):
    """PT-BR: Ordena S1 < S2 < ... < S10. EN-US: Sorts S1 < S2 < ... < S10."""
    return [int(t) if t.isdigit() else t.lower()
            for t in re.split(r'(\d+)', '' if name is None else str(name))]


def to_float(value, default=None):
    try:
        return float(str(value).strip())
    except (TypeError, ValueError):
        return default


def to_int(value, default=None):
    f = to_float(value, None)
    return default if f is None else int(f)


def num(value, decimals=0):
    """PT-BR: Formata número com separador de milhar neutro. EN-US: Formats a number with a neutral thousands separator."""
    f = to_float(value, None)
    if f is None:
        return '-'
    text = ('%%.%df' % decimals) % f
    if '.' in text:
        head, tail = text.split('.', 1)
        tail = '.' + tail
    else:
        head, tail = text, ''
    sign = ''
    if head.startswith('-'):
        sign, head = '-', head[1:]
    groups = []
    while len(head) > 3:
        groups.insert(0, head[-3:])
        head = head[:-3]
    groups.insert(0, head)
    return sign + NBSP_THIN.join(groups) + tail


def fmt_evalue(value):
    """PT-BR: E-value legível, sem os artefatos de float dos TSV. EN-US: Readable E-value, without the TSV float artifacts."""
    f = to_float(value, None)
    if f is None:
        return '-'
    if f == 0.0:
        return '0'
    if f >= 0.001:
        return '%.3g' % f
    mant, exp = ('%.2e' % f).split('e')
    return '%s×10%s' % (mant.rstrip('0').rstrip('.'), superscript(int(exp)))


def superscript(n):
    table = {'0': '⁰', '1': '¹', '2': '²', '3': '³', '4': '⁴',
             '5': '⁵', '6': '⁶', '7': '⁷', '8': '⁸', '9': '⁹',
             '-': '⁻'}
    return ''.join(table.get(c, c) for c in str(n))


def human_bytes(size):
    step = 1024.0
    for unit in ('B', 'KB', 'MB', 'GB', 'TB'):
        if size < step or unit == 'TB':
            return '%.0f%s%s' % (size, NBSP_THIN, unit) if unit == 'B' else '%.1f%s%s' % (size, NBSP_THIN, unit)
        size /= step
    return '%.1f%sTB' % (size, NBSP_THIN)


def human_bp(value):
    v = to_float(value, None)
    if v is None:
        return '-'
    if v >= 1e6:
        return '%.2f%sMb' % (v / 1e6, NBSP_THIN)
    if v >= 1e3:
        return '%.0f%skb' % (v / 1e3, NBSP_THIN)
    return '%s%sbp' % (num(v), NBSP_THIN)


def human_seconds(seconds):
    if seconds is None:
        return '-'
    if seconds < 60:
        return '%.0f%ss' % (seconds, NBSP_THIN)
    if seconds < 3600:
        return '%d%smin' % (round(seconds / 60.0), NBSP_THIN)
    return '%.1f%sh' % (seconds / 3600.0, NBSP_THIN)


def parse_duration(text):
    """PT-BR: Converte '10.4s', '3m 12s', '1h 4m 3s' em segundos. EN-US: Converts '10.4s', '3m 12s', '1h 4m 3s' into seconds."""
    if not text:
        return None
    total = 0.0
    found = False
    for value, unit in re.findall(r'([\d.]+)\s*(ms|s|m|h|d)', str(text)):
        v = to_float(value, 0.0)
        factor = {'ms': 0.001, 's': 1.0, 'm': 60.0, 'h': 3600.0, 'd': 86400.0}[unit]
        total += v * factor
        found = True
    return total if found else None


# PT-BR: Os TSV do pipeline gravam o português sem acento. A página é texto para leitura
#        humana, então o rótulo é acentuado na exibição — o dado em si nunca é alterado.
# EN-US: The pipeline TSVs write Portuguese without accents. The page is human-facing text,
#        so the label is accented for display — the data itself is never modified.
PT_DISPLAY = {
    'referencia': 'referência', 'Referencia': 'Referência',
    'Media': 'Média', 'media': 'média',
    'ancora': 'âncora', 'Ancora': 'Âncora',
    'ancoras': 'âncoras', 'Ancoras': 'Âncoras',
}
_PT_DISPLAY_RE = re.compile(r'\b(%s)\b' % '|'.join(sorted(PT_DISPLAY, key=len, reverse=True)))


def pt_display(text):
    """PT-BR: Acentua rótulos vindos dos arquivos. EN-US: Accents labels coming from the files."""
    return _PT_DISPLAY_RE.sub(lambda m: PT_DISPLAY[m.group(0)], str(text))


def split_bilingual(value):
    """
    PT-BR: Alguns campos já vêm bilíngues nos TSV ('Alta/High', 'estirpe/strain').
           Devolve tupla (pt, en) quando reconhece o padrão, senão a string original.
    EN-US: Some fields are already bilingual in the TSVs ('Alta/High', 'estirpe/strain').
           Returns a (pt, en) tuple when the pattern is recognised, else the plain string.
    """
    if not value or '/' not in str(value):
        return value
    parts = str(value).split('/')
    if len(parts) == 2 and all(p and ' ' not in p for p in parts):
        return (pt_display(parts[0]), parts[1])
    return value


def nice_num(x, round_it):
    if x <= 0:
        return 0.0
    exp = math.floor(math.log10(x))
    f = x / (10.0 ** exp)
    if round_it:
        nf = 1 if f < 1.5 else (2 if f < 3 else (5 if f < 7 else 10))
    else:
        nf = 1 if f <= 1 else (2 if f <= 2 else (5 if f <= 5 else 10))
    return nf * (10.0 ** exp)


def nice_ticks(vmin, vmax, count=5):
    """PT-BR: Ticks 1-2-5. EN-US: 1-2-5 ticks."""
    if vmax <= vmin:
        vmax = vmin + 1.0
    span = nice_num(vmax - vmin, False) or 1.0
    step = nice_num(span / max(1, count - 1), True) or 1.0
    start = math.floor(vmin / step) * step
    ticks = []
    v = start
    while v <= vmax + step * 0.5:
        ticks.append(round(v, 10))
        v += step
    return ticks, step


# ---------------------------------------------------------------------------
# PT-BR: Camada 1 - leitores / EN-US: Layer 1 - readers
# ---------------------------------------------------------------------------

def read_tsv(path, fieldnames=None):
    """PT-BR: Le um TSV como lista de dicts. EN-US: Reads a TSV as a list of dicts."""
    if not path or not os.path.isfile(path):
        return []
    rows = []
    with open(path, 'r', encoding='utf-8', errors='replace', newline='') as handle:
        reader = csv.DictReader(handle, delimiter='\t', fieldnames=fieldnames)
        for row in reader:
            rows.append({(k or ''): ('' if v is None else v) for k, v in row.items()})
    return rows


def read_text(path):
    if not path or not os.path.isfile(path):
        return ''
    with open(path, 'r', encoding='utf-8', errors='replace') as handle:
        return handle.read()


def strain_dirs(results, stage):
    """PT-BR: Subdiretórios por estirpe de um estágio. EN-US: Per-strain subdirectories of a stage."""
    base = os.path.join(results, stage)
    if not os.path.isdir(base):
        return []
    names = [n for n in os.listdir(base) if os.path.isdir(os.path.join(base, n))]
    return sorted(names, key=natural_key)


def scan_fasta(path):
    """
    PT-BR: Uma única passada por um FASTA: conta registros e guarda o primeiro cabecalho.
           Nunca mantém a sequência em memória.
    EN-US: A single pass over a FASTA: counts records and keeps the first header.
           Never holds the sequence in memory.
    """
    count = 0
    first = ''
    if not os.path.isfile(path):
        return 0, ''
    with open(path, 'r', encoding='utf-8', errors='replace') as handle:
        for line in handle:
            if line.startswith('>'):
                count += 1
                if count == 1:
                    first = line[1:].strip()
    return count, first


ORGANISM_RE = re.compile(
    r'^\S+\s+(?P<genus>[A-Z][a-z]+)\s+(?P<species>[a-z]+\.?)\s+strain\s+(?P<strain>\S+)')


def parse_organism(header):
    """
    PT-BR: Extrai gênero, espécie e designação de estirpe da linha de definição do RefSeq,
           ex.: 'NZ_VKHW01000001.1 Lysinibacillus fusiformis strain SDF0005 scaffold1.1, ...'
    EN-US: Extracts genus, species and strain designation from the RefSeq definition line.
    """
    match = ORGANISM_RE.match(header or '')
    if not match:
        return None
    return {'genus': match.group('genus'),
            'species': match.group('species'),
            'designation': match.group('strain')}


def read_strains(results, accessions_path, warn):
    """
    PT-BR: Monta o registro canonico de cada estirpe: id interno, taxonomia real (do genoma),
           acesso GenBank e número de contigs. A taxonomia do sample sheet e apenas fallback.
    EN-US: Builds the canonical record for each strain: internal id, real taxonomy (from the
           genome), GenBank accession and contig count. The sample sheet taxonomy is only a fallback.
    """
    accessions = {}
    for row in read_tsv(accessions_path):
        accessions[row.get('strain', '').strip()] = row

    genome_dir = os.path.join(results, '00_genomes')
    files = sorted(glob.glob(os.path.join(genome_dir, '*.fna')),
                   key=lambda p: natural_key(os.path.basename(p)))
    strains = OrderedDict()
    for path in files:
        code = os.path.basename(path)[:-4]
        contigs, header = scan_fasta(path)
        organism = parse_organism(header)
        if organism is None:
            fallback = accessions.get(code, {}).get('species', '').strip()
            warn('Taxonomia não extraída do genoma / Taxonomy not parsed from genome: %s' % code)
            parts = fallback.split()
            organism = {'genus': parts[0] if parts else code,
                        'species': parts[1] if len(parts) > 1 else 'sp.',
                        'designation': ''}
        strains[code] = {
            'code': code,
            'genus': organism['genus'],
            'species': organism['species'],
            'designation': organism['designation'],
            'accession': accessions.get(code, {}).get('genbank_accession', '').strip(),
            'sheet_species': accessions.get(code, {}).get('species', '').strip(),
            'contigs': contigs,
            'genome_path': path,
            'genome_size': os.path.getsize(path) if os.path.isfile(path) else 0,
        }
    return strains


def label_full(strain):
    """PT-BR: 'Lysinibacillus fusiformis SDF0005'. EN-US: same."""
    binomial = '%s %s' % (strain['genus'], strain['species'])
    return (binomial + ' ' + strain['designation']).strip()


def label_short(strain):
    """PT-BR: 'L. fusiformis SDF0005' para eixos. EN-US: 'L. fusiformis SDF0005' for axes."""
    binomial = '%s. %s' % (strain['genus'][0], strain['species'])
    return (binomial + ' ' + strain['designation']).strip()


def label_html(strain, short=False):
    """PT-BR: Binomial em itálico + designação em roman. EN-US: Italic binomial + roman designation."""
    if short:
        binomial = '%s. %s' % (strain['genus'][0], strain['species'])
    else:
        binomial = '%s %s' % (strain['genus'], strain['species'])
    out = '<i>%s</i>' % esc(binomial)
    if strain['designation']:
        out += ' ' + esc(strain['designation'])
    return out


def read_checkm2(results, strains):
    for code in strains:
        path = os.path.join(results, '01_qc', 'checkm2', code, '%s_checkm2.tsv' % code)
        rows = read_tsv(path)
        strains[code]['checkm2'] = rows[0] if rows else {}


def read_quast(results, strains):
    for code in strains:
        path = os.path.join(results, '01_qc', 'quast', code, 'transposed_report.tsv')
        rows = read_tsv(path)
        strains[code]['quast'] = rows[0] if rows else {}
        strains[code]['quast_dir'] = os.path.join(results, '01_qc', 'quast', code)


def read_bakta(results, strains):
    for code in strains:
        faa = os.path.join(results, '02_bakta', code, '%s.faa' % code)
        proteins, _ = scan_fasta(faa)
        strains[code]['proteins'] = proteins
        strains[code]['bakta_dir'] = os.path.join(results, '02_bakta', code)


BLAST_COLUMNS = ['qseqid', 'sseqid', 'pident', 'length', 'qstart', 'qend',
                 'qlen', 'slen', 'evalue', 'bitscore']


def parse_reference_id(sseqid):
    """
    PT-BR: 'ref|ureC|UniProt:A7Z9N4|Bacillus_velezensis_(strain_...)|reviewed|569aa'
    EN-US: same, split into gene / accession / organism / length.
    """
    parts = str(sseqid).split('|')
    info = {'gene': '', 'uniprot': '', 'organism': str(sseqid), 'length': ''}
    if len(parts) >= 6 and parts[0] == 'ref':
        info['gene'] = parts[1]
        info['uniprot'] = parts[2].replace('UniProt:', '')
        info['organism'] = parts[3].replace('_', ' ')
        info['length'] = parts[5]
    return info


def read_blast(results, strains):
    """PT-BR: HSPs de BLASTp agrupados por proteína consulta. EN-US: BLASTp HSPs grouped by query protein."""
    by_query = OrderedDict()
    total = 0
    for code in strains:
        path = os.path.join(results, '03_candidates', code, '%s_blast.tsv' % code)
        for row in read_tsv(path, fieldnames=BLAST_COLUMNS):
            total += 1
            by_query.setdefault(row['qseqid'], []).append(row)
    for hits in by_query.values():
        hits.sort(key=lambda r: to_float(r['bitscore'], 0.0), reverse=True)
    return by_query, total


def read_locus_status(results, strains):
    """PT-BR: Le os *_locus_status.txt. EN-US: Reads the *_locus_status.txt files."""
    status = OrderedDict()
    for code in strains:
        path = os.path.join(results, '05_synteny', 'loci', '%s_locus_status.txt' % code)
        entry = {'strain': code, 'regions': [], 'result': '', 'note': '',
                 'candidates_total': None, 'anchors_total': None, 'flank_genes': None,
                 'regions_total': 0, 'genes_total': 0}
        for line in read_text(path).splitlines():
            fields = line.rstrip('\n').split('\t')
            if not fields or not fields[0]:
                continue
            key = fields[0]
            if key == 'region' and len(fields) >= 7:
                entry['regions'].append({
                    'id': fields[1], 'contig': fields[2],
                    'start': to_int(fields[3], 0), 'end': to_int(fields[4], 0),
                    'genes': to_int(fields[5], 0),
                    'anchors': [a for a in fields[6].split(',') if a],
                })
            elif len(fields) >= 2:
                if key in ('candidates_total', 'anchors_total', 'flank_genes',
                           'regions_total', 'genes_total'):
                    entry[key] = to_int(fields[1], 0)
                elif key in ('result', 'note'):
                    entry[key] = fields[1]
        status[code] = entry
    return status


# --- PT-BR: parser GenBank mínimo / EN-US: minimal GenBank parser ----------

GBK_LOCUS_RE = re.compile(r'^LOCUS\s+(\S+)\s+(\d+)\s+bp')
GBK_FEATURE_RE = re.compile(r'^ {5}(\S+)\s+(.+)$')
GBK_QUAL_RE = re.compile(r'^ {21}/(\w+)=?(.*)$')
GBK_CONT_RE = re.compile(r'^ {21}(.+)$')


def parse_gbk(path):
    """
    PT-BR: Extrai regiões e CDS de um GenBank sem Biopython. Ignora /translation e /db_xref,
           que são longos e irrelevantes para o mapa.
    EN-US: Extracts regions and CDS from a GenBank file without Biopython. Skips /translation
           and /db_xref, which are long and irrelevant to the map.
    """
    records = []
    record = None
    feature = None
    qual_key = None
    section = 'header'
    last_header = ''
    if not os.path.isfile(path):
        return records
    with open(path, 'r', encoding='utf-8', errors='replace') as handle:
        for raw in handle:
            line = raw.rstrip('\n')
            match = GBK_LOCUS_RE.match(line)
            if match:
                record = {'id': match.group(1), 'length': to_int(match.group(2), 0),
                          'definition': '', 'contig': '', 'start': 0, 'end': 0,
                          'anchors': [], 'features': []}
                records.append(record)
                feature, qual_key, section, last_header = None, None, 'header', 'LOCUS'
                continue
            if record is None:
                continue
            if line.startswith('FEATURES'):
                section = 'features'
                continue
            if line.startswith('ORIGIN') or line.startswith('//'):
                section = 'origin'
                feature, qual_key = None, None
                continue
            if section == 'header':
                if line.startswith('DEFINITION'):
                    record['definition'] = line[12:].strip()
                    last_header = 'DEFINITION'
                elif line.startswith(' ' * 12) and last_header == 'DEFINITION':
                    record['definition'] += ' ' + line.strip()
                elif line[:1] not in (' ', ''):
                    last_header = line.split()[0] if line.split() else ''
                continue
            if section != 'features':
                continue
            match = GBK_FEATURE_RE.match(line)
            if match:
                key, location = match.group(1), match.group(2)
                if key == 'CDS':
                    coords = [int(c) for c in re.findall(r'\d+', location)]
                    feature = {'start': min(coords) if coords else 0,
                               'end': max(coords) if coords else 0,
                               'strand': -1 if 'complement' in location else 1,
                               'locus_tag': '', 'gene': '', 'product': ''}
                    record['features'].append(feature)
                else:
                    feature = None
                qual_key = None
                continue
            match = GBK_QUAL_RE.match(line)
            if match and feature is not None:
                qual_key = match.group(1)
                value = match.group(2).strip().strip('"')
                if qual_key in ('locus_tag', 'gene', 'product'):
                    feature[qual_key] = value
                continue
            match = GBK_CONT_RE.match(line)
            if match and feature is not None and qual_key in ('product', 'gene', 'locus_tag'):
                feature[qual_key] = (feature[qual_key] + ' ' + match.group(1).strip()).strip().strip('"')

    for record in records:
        loc = re.search(r'\[([^\]:]+):(\d+)-(\d+)\]', record['definition'])
        if loc:
            record['contig'] = loc.group(1)
            record['start'] = to_int(loc.group(2), 0)
            record['end'] = to_int(loc.group(3), 0)
        anchors = re.search(r'anchors=([^.\s]+)', record['definition'])
        if anchors:
            record['anchors'] = [a for a in anchors.group(1).split(',') if a]
    return records


def read_loci(results, strains):
    loci = OrderedDict()
    for code in strains:
        path = os.path.join(results, '05_synteny', 'loci', '%s.gbk' % code)
        records = parse_gbk(path)
        if records:
            loci[code] = records
    return loci


# --- PT-BR: parser Newick / EN-US: Newick parser ---------------------------

def parse_newick(text):
    """PT-BR: Descida recursiva sobre o Newick do IQ-TREE. EN-US: Recursive descent over the IQ-TREE Newick."""
    text = (text or '').strip()
    if not text:
        return None
    pos = [0]
    size = len(text)

    def read_label():
        start = pos[0]
        if pos[0] < size and text[pos[0]] == "'":
            pos[0] += 1
            begin = pos[0]
            while pos[0] < size and text[pos[0]] != "'":
                pos[0] += 1
            label = text[begin:pos[0]]
            pos[0] += 1
            return label
        while pos[0] < size and text[pos[0]] not in '(),:;':
            pos[0] += 1
        return text[start:pos[0]].strip()

    def read_node():
        node = {'name': '', 'length': 0.0, 'support': None, 'children': []}
        if pos[0] < size and text[pos[0]] == '(':
            pos[0] += 1
            while True:
                node['children'].append(read_node())
                if pos[0] < size and text[pos[0]] == ',':
                    pos[0] += 1
                    continue
                if pos[0] < size and text[pos[0]] == ')':
                    pos[0] += 1
                break
        label = read_label()
        if pos[0] < size and text[pos[0]] == ':':
            pos[0] += 1
            start = pos[0]
            while pos[0] < size and text[pos[0]] not in '(),;':
                pos[0] += 1
            node['length'] = to_float(text[start:pos[0]], 0.0) or 0.0
        if node['children']:
            node['support'] = to_float(label, None)
        else:
            node['name'] = label
        return node

    return read_node()


def read_iqtree_stats(path):
    """PT-BR: Extrai as estatísticas do relatório do IQ-TREE. EN-US: Extracts stats from the IQ-TREE report."""
    stats = {}
    patterns = {
        'taxa': r'Input data:\s+(\d+)\s+sequences with\s+(\d+)\s+amino-acid sites',
        'constant': r'Number of constant sites:\s+(\d+)\s+\(=\s*([\d.]+)%',
        'informative': r'Number of parsimony informative sites:\s+(\d+)',
        'patterns': r'Number of distinct site patterns:\s+(\d+)',
        'model': r'Best-fit model according to BIC:\s+(\S+)',
        'loglik': r'Log-likelihood of the tree:\s+(-?[\d.]+)',
    }
    text = read_text(path)
    for key, pattern in patterns.items():
        match = re.search(pattern, text)
        if not match:
            continue
        if key == 'taxa':
            stats['taxa'] = to_int(match.group(1))
            stats['sites'] = to_int(match.group(2))
        elif key == 'constant':
            stats['constant'] = to_int(match.group(1))
            stats['constant_pct'] = to_float(match.group(2))
        elif key in ('informative', 'patterns'):
            stats[key] = to_int(match.group(1))
        else:
            stats[key] = match.group(1)
    return stats


def read_trace(path):
    """PT-BR: Le nf_trace.txt e agrega por processo. EN-US: Reads nf_trace.txt and aggregates per process."""
    rows = read_tsv(path)
    tasks = []
    for row in rows:
        name = row.get('name', '')
        process = name.split('(')[0].strip() or name
        tasks.append({
            'process': process,
            'name': name,
            'status': row.get('status', ''),
            'exit': row.get('exit', ''),
            'submit': row.get('submit', ''),
            'duration_s': parse_duration(row.get('duration', '')),
            'realtime_s': parse_duration(row.get('realtime', '')),
            'peak_rss': row.get('peak_rss', ''),
        })
    return tasks


# ---------------------------------------------------------------------------
# PT-BR: Camada 2 - i18n / EN-US: Layer 2 - i18n
# ---------------------------------------------------------------------------

def bi_attrs(value):
    """
    PT-BR: Recebe str ou (pt, en) e devolve (atributos, texto_padrao). Os elementos
           marcados com .i18n são trocados em tempo real pelo setLang() da página,
           tanto em HTML quanto dentro dos SVG.
    EN-US: Takes a str or a (pt, en) tuple and returns (attributes, default_text).
           Elements carrying .i18n are swapped live by the page's setLang(), both in
           HTML and inside the SVGs.
    """
    if isinstance(value, (tuple, list)) and len(value) == 2:
        pt, en = value
        attrs = ' class="i18n" data-pt="%s" data-en="%s"' % (esc(pt), esc(en))
        return attrs, esc(pt if LANG == 'pt' else en)
    return '', esc(value)


def T(value, tag='span', cls=''):
    """PT-BR: Texto bilíngue em HTML. EN-US: Bilingual text in HTML."""
    attrs, text = bi_attrs(value)
    if cls:
        if ' class="i18n"' in attrs:
            attrs = attrs.replace(' class="i18n"', ' class="i18n %s"' % esc(cls))
        else:
            attrs = ' class="%s"' % esc(cls)
    return '<%s%s>%s</%s>' % (tag, attrs, text, tag)


def plain(value):
    """PT-BR: Versão textual no idioma padrao. EN-US: Plain text in the default language."""
    if isinstance(value, (tuple, list)) and len(value) == 2:
        return value[0] if LANG == 'pt' else value[1]
    return '' if value is None else str(value)


# ---------------------------------------------------------------------------
# PT-BR: Camada 3 - primitivas SVG / EN-US: Layer 3 - SVG primitives
# ---------------------------------------------------------------------------

def svg_open(width, height, cls='chart'):
    return ('<svg class="%s" viewBox="0 0 %g %g" width="100%%" '
            'preserveAspectRatio="xMidYMid meet" role="img">' % (esc(cls), width, height))


SVG_CLOSE = '</svg>'


def svg_text(x, y, value, size=11, fill='var(--text-secondary)', anchor='start',
             weight='normal', style='', extra=''):
    attrs, text = bi_attrs(value)
    style_attr = ' font-style="italic"' if style == 'italic' else ''
    return ('<text x="%g" y="%g" font-size="%g" fill="%s" text-anchor="%s" '
            'font-weight="%s"%s%s%s>%s</text>'
            % (x, y, size, fill, anchor, weight, style_attr, attrs, extra, text))


_NARROW_CHARS = set("iljtfIr.,;:'|!()[] ")
_WIDE_CHARS = set("MWmw@%")


def text_width(text, size):
    """
    PT-BR: Estimativa de largura de texto, não medicao. Serve para decidir truncamento e
           colisão de rótulos em tempo de geração, onde não ha motor de fontes disponível.
    EN-US: Text width estimate, not a measurement. Used to decide truncation and label
           collision at generation time, where no font engine is available.
    """
    total = 0.0
    for char in str(text):
        if char in _NARROW_CHARS:
            total += 0.32
        elif char in _WIDE_CHARS:
            total += 0.90
        elif char.isupper() or char.isdigit():
            total += 0.60
        else:
            total += 0.52
    return total * size


def svg_binomial(x, y, binomial, suffix='', size=11, fill='var(--text-primary)',
                 anchor='start', weight='normal'):
    """PT-BR: Nome científico em itálico + sufixo em roman. EN-US: Italic scientific name + roman suffix."""
    tail = ('<tspan font-style="normal"> %s</tspan>' % esc(suffix)) if suffix else ''
    return ('<text x="%g" y="%g" font-size="%g" fill="%s" text-anchor="%s" font-weight="%s">'
            '<tspan font-style="italic">%s</tspan>%s</text>'
            % (x, y, size, fill, anchor, weight, esc(binomial), tail))


def svg_rect(x, y, width, height, fill, rx=0, extra='', title=None):
    body = ('<rect x="%g" y="%g" width="%g" height="%g" rx="%g" fill="%s"%s>'
            % (x, y, max(0.0, width), max(0.0, height), rx, fill, extra))
    body += ('<title>%s</title>' % esc(title)) if title else ''
    return body + '</rect>'


def svg_line(x1, y1, x2, y2, stroke='var(--grid)', width=1, dash=''):
    dash_attr = ' stroke-dasharray="%s"' % dash if dash else ''
    return ('<line x1="%g" y1="%g" x2="%g" y2="%g" stroke="%s" stroke-width="%g"%s/>'
            % (x1, y1, x2, y2, stroke, width, dash_attr))


def svg_path(d, fill='none', stroke='none', width=1, extra='', title=None):
    body = '<path d="%s" fill="%s" stroke="%s" stroke-width="%g"%s>' % (d, fill, stroke, width, extra)
    body += ('<title>%s</title>' % esc(title)) if title else ''
    return body + '</path>'


def svg_circle(cx, cy, r, fill, extra='', title=None):
    body = '<circle cx="%g" cy="%g" r="%g" fill="%s"%s>' % (cx, cy, r, fill, extra)
    body += ('<title>%s</title>' % esc(title)) if title else ''
    return body + '</circle>'


def bar_end_path(x, y, width, height, radius=4):
    """
    PT-BR: Barra horizontal com a extremidade de dado arredondada (4px) e a base reta,
           ancorada na linha de base.
    EN-US: Horizontal bar with a 4px rounded data end and a flat, baseline-anchored base.
    """
    r = min(radius, max(0.0, width), height / 2.0)
    if width <= 0.5:
        return 'M%g,%g L%g,%g' % (x, y, x, y + height)
    return ('M%g,%g L%g,%g Q%g,%g %g,%g L%g,%g Q%g,%g %g,%g L%g,%g Z'
            % (x, y, x + width - r, y, x + width, y, x + width, y + r,
               x + width, y + height - r, x + width, y + height, x + width - r, y + height,
               x, y + height))


def gene_arrow_path(x, y, width, height, strand):
    """PT-BR: Seta de gene. EN-US: Gene arrow."""
    tip = min(9.0, max(2.0, width * 0.35))
    if strand >= 0:
        return ('M%g,%g L%g,%g L%g,%g L%g,%g L%g,%g Z'
                % (x, y, x + width - tip, y, x + width, y + height / 2.0,
                   x + width - tip, y + height, x, y + height))
    return ('M%g,%g L%g,%g L%g,%g L%g,%g L%g,%g Z'
            % (x + width, y, x + tip, y, x, y + height / 2.0,
               x + tip, y + height, x + width, y + height))


# ---------------------------------------------------------------------------
# PT-BR: Camada 4 - gráficos / EN-US: Layer 4 - charts
# ---------------------------------------------------------------------------

def hbar_chart(items, unit=None, colour='var(--series-1)', value_fmt=None, height_per=26,
               label_width=196, width=920, tick_format=None):
    """
    PT-BR: Barras horizontais, uma serie. Rótulos diretos no fim de cada barra; sem legenda
           (uma serie so e nomeada pelo título do bloco). Grade recessiva.
    EN-US: Horizontal bars, single series. Direct labels at each bar end; no legend
           (a single series is named by the block title). Recessive grid.
    items: list of dicts with 'label' (str or html-safe text), 'value' (float),
           optional 'title' (tooltip), 'colour', 'label_html'.
    """
    if not items:
        return empty_note()
    value_fmt = value_fmt or (lambda v: num(v))
    top, bottom, right = 26, 30, 108
    plot_w = width - label_width - right
    height = top + bottom + height_per * len(items)
    vmax = max([to_float(i['value'], 0.0) or 0.0 for i in items] + [0.0])
    ticks, step = nice_ticks(0, vmax or 1.0, 5)
    scale = (plot_w / ticks[-1]) if ticks[-1] else 0.0
    # PT-BR: Casas decimais pelo passo, senao 0/0,5/1 vira "0, 0, 1" no eixo.
    # EN-US: Decimals follow the step, otherwise 0/0.5/1 renders as "0, 0, 1" on the axis.
    decimals = 0 if step >= 1 else min(3, int(math.ceil(-math.log10(step))))
    tick_fmt = tick_format or (lambda v: num(v, decimals))

    out = [svg_open(width, height)]
    for tick in ticks:
        x = label_width + tick * scale
        out.append(svg_line(x, top - 6, x, height - bottom + 2, 'var(--grid)', 1))
        out.append(svg_text(x, height - bottom + 16, tick_fmt(tick), 10, 'var(--muted)', 'middle'))
    if unit:
        out.append(svg_text(width - 2, 12, unit, 10, 'var(--muted)', 'end'))

    for index, item in enumerate(items):
        y = top + index * height_per
        bar_h = height_per - 8
        value = to_float(item['value'], 0.0) or 0.0
        bar_w = value * scale
        fill = item.get('colour', colour)
        if item.get('label_html'):
            out.append(item['label_html'](label_width - 8, y + bar_h * 0.72))
        else:
            out.append(svg_text(label_width - 8, y + bar_h * 0.72, item['label'], 11,
                                'var(--text-primary)', 'end'))
        out.append(svg_path(bar_end_path(label_width, y, bar_w, bar_h), fill,
                            title=item.get('title')))
        out.append(svg_text(label_width + bar_w + 6, y + bar_h * 0.72,
                            value_fmt(value), 10, 'var(--text-secondary)', 'start'))
    out.append(svg_line(label_width, top - 6, label_width, height - bottom + 2,
                        'var(--baseline)', 1))
    out.append(SVG_CLOSE)
    return ''.join(out)


def empty_note():
    return '<p class="empty">%s</p>' % T(('Sem dados para este bloco.', 'No data for this block.'))


def scatter_chart(points, x_label, y_label, width=920, height=420):
    """PT-BR: Dispersão com rótulos diretos. EN-US: Scatter with direct labels."""
    if not points:
        return empty_note()
    left, right, top, bottom = 62, 24, 24, 46
    plot_w, plot_h = width - left - right, height - top - bottom
    xs = [p['x'] for p in points]
    ys = [p['y'] for p in points]
    xticks, _ = nice_ticks(min(xs) * 0.98, max(xs) * 1.02, 5)
    yticks, _ = nice_ticks(min(ys) * 0.98, max(ys) * 1.02, 5)
    x0, x1 = xticks[0], xticks[-1]
    y0, y1 = yticks[0], yticks[-1]

    def sx(v):
        return left + (v - x0) / (x1 - x0 or 1.0) * plot_w

    def sy(v):
        return top + plot_h - (v - y0) / (y1 - y0 or 1.0) * plot_h

    out = [svg_open(width, height)]
    for tick in yticks:
        out.append(svg_line(left, sy(tick), left + plot_w, sy(tick), 'var(--grid)', 1))
        out.append(svg_text(left - 8, sy(tick) + 3, num(tick, 1), 10, 'var(--muted)', 'end'))
    for tick in xticks:
        out.append(svg_text(sx(tick), top + plot_h + 18, num(tick, 2), 10, 'var(--muted)', 'middle'))
    out.append(svg_line(left, top + plot_h, left + plot_w, top + plot_h, 'var(--baseline)', 1))
    out.append(svg_text(left + plot_w / 2.0, height - 8, x_label, 11, 'var(--text-secondary)', 'middle'))
    attrs, text = bi_attrs(y_label)
    out.append('<text x="14" y="%g" font-size="11" fill="var(--text-secondary)" '
               'text-anchor="middle" transform="rotate(-90 14 %g)"%s>%s</text>'
               % (top + plot_h / 2.0, top + plot_h / 2.0, attrs, text))

    # PT-BR: Rótulos diretos com desvio de colisão. Varios genomas caem quase no mesmo ponto
    #        (as duas Bacillus velezensis diferem em 0,02 Mb), e rótulos empilhados são ilegíveis.
    # EN-US: Direct labels with collision avoidance. Several genomes land on nearly the same
    #        point (the two Bacillus velezensis differ by 0.02 Mb), and stacked labels are unreadable.
    placed = []
    marks = []
    labels = []
    for point in sorted(points, key=lambda p: (p['y'], p['x'])):
        cx, cy = sx(point['x']), sy(point['y'])
        text = ('%s %s' % (point['binomial'], point.get('suffix', ''))).strip()
        label_w = text_width(text, 10) + 12
        for dx, dy, anchor in ((10, -8, 'start'), (10, 14, 'start'), (-10, -8, 'end'),
                               (-10, 14, 'end'), (10, -22, 'start'), (10, 28, 'start')):
            lx = cx + dx
            box = (lx if anchor == 'start' else lx - label_w, cy + dy - 9,
                   (lx + label_w) if anchor == 'start' else lx, cy + dy + 3)
            if all(box[2] < o[0] or box[0] > o[2] or box[3] < o[1] or box[1] > o[3]
                   for o in placed) and box[0] > 4 and box[2] < width - 4:
                placed.append(box)
                labels.append(svg_binomial(lx, cy + dy, point['binomial'],
                                           point.get('suffix', ''), 10,
                                           'var(--text-secondary)', anchor))
                break
        marks.append(svg_circle(cx, cy, 6, point['colour'],
                                extra=' stroke="var(--surface-1)" stroke-width="2"',
                                title=point.get('title')))
    out.extend(marks)
    out.extend(labels)
    out.append(SVG_CLOSE)
    return ''.join(out)


def funnel_chart(stages, width=920):
    """
    PT-BR: Funil da triagem. So entram aqui etapas genuinamente sequenciais, isto e, cada uma
           subconjunto da anterior -- as três camadas de evidência NAO são um funil, são filtros
           paralelos, e aparecem em gráfico próprio. A largura e proporcional a log10(valor)
           porque a queda vai de 46 mil para 16; o valor exato e a fração retida ficam impressos
           ao lado de cada faixa, então a leitura nunca depende da largura.
    EN-US: Screening funnel. Only genuinely sequential stages belong here, i.e. each one a subset
           of the previous -- the three evidence layers are NOT a funnel, they are parallel
           filters, and get their own chart. Width is proportional to log10(value) because the
           drop goes from 46k down to 16; the exact value and the retained fraction are printed
           beside each band, so reading it never depends on width.
    """
    if not stages:
        return empty_note()
    band_h, gap = 52, 12
    height = 20 + len(stages) * (band_h + gap)
    # PT-BR: A faixa nunca invade o espaco do texto a direita, então rótulo e nota ficam fora
    #        dela e nenhum texto some sobre o preenchimento.
    # EN-US: The band never invades the text column on the right, so label and note live outside
    #        it and no text is ever lost against the fill.
    max_w = width * 0.42
    values = [max(1.0, to_float(s['value'], 1.0) or 1.0) for s in stages]
    logs = [math.log10(v + 1.0) for v in values]
    top_log = max(logs) or 1.0
    slots = [1, 3, 5] if len(stages) == 3 else list(range(1, len(stages) + 1))
    out = [svg_open(width, height)]
    for index, stage in enumerate(stages):
        slot = slots[index] if index < len(slots) else len(slots)
        y = 10 + index * (band_h + gap)
        w = max(40.0, max_w * (logs[index] / top_log))
        x = 8
        out.append(svg_path(bar_end_path(x, y, w, band_h, 4),
                            'var(--funnel-%d)' % slot, title=plain(stage['label'])))
        out.append(svg_text(x + w + 18, y + 27, num(stage['value']), 22,
                            'var(--text-primary)', 'start', 'bold'))
        text_x = x + w + 18 + text_width(num(stage['value']), 22) + 16
        out.append(svg_text(text_x, y + 22, stage['label'], 12.5,
                            'var(--text-primary)', 'start', 'bold'))
        note = stage['note']
        if index > 0 and values[index - 1] > 0:
            pct = 100.0 * values[index] / values[index - 1]
            fmt = '%.2f%%' if pct < 1 else '%.1f%%'
            note_pt, note_en = (note if isinstance(note, tuple) else (note, note))
            out.append(svg_text(text_x, y + 39,
                                ('%s da etapa anterior · %s' % (fmt % pct, note_pt),
                                 '%s of the previous stage · %s' % (fmt % pct, note_en)),
                                10, 'var(--text-secondary)', 'start'))
        else:
            out.append(svg_text(text_x, y + 39, note, 10, 'var(--text-secondary)', 'start'))
    out.append(SVG_CLOSE)
    return ''.join(out)


def presence_heatmap(strains, matrix, width=920):
    """PT-BR: Matriz binaria presença/ausência. EN-US: Binary presence/absence matrix."""
    if not matrix:
        return empty_note()
    label_w, cell_w, cell_h, gap = 210, 52, 30, 3
    top = 44
    height = top + len(matrix) * (cell_h + gap) + 26
    out = [svg_open(width, height)]
    for col, gene in enumerate(UREASE_GENES):
        x = label_w + col * (cell_w + gap) + cell_w / 2.0
        out.append(svg_text(x, top - 12, gene, 11, 'var(--text-primary)', 'middle', 'bold'))
    for row, entry in enumerate(matrix):
        strain = entry['strain']
        y = top + row * (cell_h + gap)
        out.append(svg_binomial(label_w - 12, y + 20, entry['binomial'], entry['suffix'],
                                11, 'var(--text-primary)', 'end'))
        for col, gene in enumerate(UREASE_GENES):
            x = label_w + col * (cell_w + gap)
            present = entry['values'].get(gene) == '1'
            tip_pt, tip_en = GENE_TOOLTIP.get(gene, (gene, gene))
            title = '%s — %s: %s' % (entry['full'], gene,
                                          plain(('presente', 'present')) if present
                                          else plain(('ausente', 'absent')))
            title += ' (%s)' % (tip_pt if LANG == 'pt' else tip_en)
            if present:
                out.append(svg_rect(x, y, cell_w, cell_h, 'var(--series-1)', 4, title=title))
                out.append(svg_path('M%g,%g l%g,%g l%g,%g' % (x + cell_w / 2 - 7, y + cell_h / 2,
                                                              5, 5, 9, -11),
                                    'none', 'var(--on-series-1)', 2,
                                    extra=' stroke-linecap="round" stroke-linejoin="round"'))
            else:
                out.append(svg_rect(x, y, cell_w, cell_h, 'var(--cell-empty)', 4,
                                    extra=' stroke="var(--grid)" stroke-width="1"', title=title))
    out.append(SVG_CLOSE)
    return ''.join(out)


def phylogram(root, leaf_meta, width=980, leaf_step=18):
    """
    PT-BR: Filograma retangular. x = distância acumulada até a raiz; y das folhas em ordem de
           travessia (sem cruzamento de ramos); y interno = média dos filhos. Suporte impresso
           apenas quando >= 70 -- com 37 nós internos o resto vira ruído.
    EN-US: Rectangular phylogram. x = accumulated distance to the root; leaf y in traversal
           order (no crossing branches); internal y = mean of children. Support printed only
           when >= 70 -- with 37 internal nodes the rest is noise.
    """
    if root is None:
        return empty_note()
    leaves = []

    def assign_x(node, base):
        node['x'] = base + (node['length'] or 0.0)
        if node['children']:
            for child in node['children']:
                assign_x(child, node['x'])
        else:
            leaves.append(node)

    assign_x(root, 0.0)
    if not leaves:
        return empty_note()
    for index, leaf in enumerate(leaves):
        leaf['y'] = float(index)

    def assign_y(node):
        if node['children']:
            ys = [assign_y(c) for c in node['children']]
            node['y'] = (min(ys) + max(ys)) / 2.0
        return node['y']

    assign_y(root)

    label_w = 330
    left, top = 12, 26
    plot_w = width - label_w - left - 16
    max_x = max(leaf['x'] for leaf in leaves) or 1.0
    height = top + len(leaves) * leaf_step + 54
    scale = plot_w / max_x

    def px(value):
        return left + value * scale

    def py(value):
        return top + value * leaf_step + leaf_step / 2.0

    out = [svg_open(width, height, 'chart tree')]

    def draw(node, parent_x):
        colour = 'var(--branch)'
        out.append(svg_line(px(parent_x), py(node['y']), px(node['x']), py(node['y']),
                            colour, 1.4))
        if node['children']:
            ys = [child['y'] for child in node['children']]
            out.append(svg_line(px(node['x']), py(min(ys)), px(node['x']), py(max(ys)),
                                colour, 1.4))
            if node['support'] is not None and node['support'] >= 70:
                out.append(svg_text(px(node['x']) - 3, py(node['y']) - 3,
                                    '%g' % node['support'], 8, 'var(--muted)', 'end'))
            for child in node['children']:
                draw(child, node['x'])
        else:
            meta = leaf_meta.get(node['name'], {})
            is_strain = meta.get('origin') == 'strain'
            fill = 'var(--series-1)' if is_strain else 'var(--muted)'
            out.append(svg_line(px(node['x']), py(node['y']), px(max_x) + 6, py(node['y']),
                                'var(--grid)', 1, dash='1 3'))
            out.append(svg_circle(px(node['x']), py(node['y']), 4 if is_strain else 2.5, fill,
                                  extra=' stroke="var(--surface-1)" stroke-width="1.5"',
                                  title=meta.get('title', node['name'])))
            out.append(svg_binomial(px(max_x) + 14, py(node['y']) + 3.5,
                                    meta.get('binomial', node['name']), meta.get('suffix', ''),
                                    10.5,
                                    'var(--text-primary)' if is_strain else 'var(--text-secondary)',
                                    'start', 'bold' if is_strain else 'normal'))

    draw(root, 0.0)

    bar_value = nice_num(max_x / 5.0, True) or (max_x / 5.0)
    bar_px = bar_value * scale
    bar_y = height - 22
    out.append(svg_line(left, bar_y, left + bar_px, bar_y, 'var(--text-secondary)', 2))
    out.append(svg_line(left, bar_y - 4, left, bar_y + 4, 'var(--text-secondary)', 2))
    out.append(svg_line(left + bar_px, bar_y - 4, left + bar_px, bar_y + 4, 'var(--text-secondary)', 2))
    out.append(svg_text(left + bar_px + 8, bar_y + 4,
                        ('%g substituições/sítio' % bar_value,
                         '%g substitutions/site' % bar_value),
                        10, 'var(--muted)'))
    out.append(SVG_CLOSE)
    return ''.join(out)


def locus_map(record, strain, candidate_tags, width=920):
    """PT-BR: Mapa de setas de genes de uma região. EN-US: Gene-arrow map of one region."""
    left, right = 14, 14
    track_h, top = 30, 46
    height = top + track_h + 46
    span = max(1, record['length'])
    plot_w = width - left - right
    out = [svg_open(width, height, 'chart locus')]
    out.append(svg_text(left, 16, '%s — %s:%s–%s' % (record['id'], record['contig'],
                                                               num(record['start']), num(record['end'])),
                        11, 'var(--text-primary)', 'start', 'bold'))
    out.append(svg_text(width - right, 16,
                        '%s %s · %s %s' % (num(record['length']), 'bp',
                                                num(len(record['features'])),
                                                plain(('genes', 'genes'))),
                        10, 'var(--muted)', 'end'))
    out.append(svg_line(left, top + track_h / 2.0, left + plot_w, top + track_h / 2.0,
                        'var(--grid)', 2))

    for feature in record['features']:
        start, end = feature['start'], feature['end']
        x = left + (start - 1) / float(span) * plot_w
        w = max(3.0, (end - start + 1) / float(span) * plot_w)
        tier = classify_feature(feature, candidate_tags)
        fill = {'urease': 'var(--series-1)', 'candidate': 'var(--series-2)'}.get(tier, 'var(--gene-other)')
        title = '%s%s\n%s\n%s–%s (%s)' % (
            feature['locus_tag'],
            ' · %s' % feature['gene'] if feature['gene'] else '',
            feature['product'] or '-', num(start), num(end),
            '+' if feature['strand'] >= 0 else '-')
        out.append(svg_path(gene_arrow_path(x, top, w, track_h, feature['strand']), fill,
                            stroke='var(--surface-1)', width=2, title=title))
        if w >= 30 and feature['gene']:
            out.append(svg_text(x + w / 2.0, top - 8, feature['gene'], 9.5,
                                'var(--text-primary)' if tier != 'other' else 'var(--muted)',
                                'middle', 'bold' if tier == 'urease' else 'normal'))
    ticks, _ = nice_ticks(0, span, 6)
    for tick in ticks:
        if tick > span:
            continue
        x = left + tick / float(span) * plot_w
        out.append(svg_line(x, top + track_h + 4, x, top + track_h + 9, 'var(--grid)', 1))
        out.append(svg_text(x, top + track_h + 22, human_bp(record['start'] + tick), 9,
                            'var(--muted)', 'middle'))
    out.append(SVG_CLOSE)
    return ''.join(out)


def classify_feature(feature, candidate_tags):
    gene = (feature.get('gene') or '').strip()
    product = feature.get('product') or ''
    for target, patterns in COMPILED_PATTERNS.items():
        if gene.lower() == target.lower():
            return 'urease'
        for pattern in patterns:
            if pattern.search(product):
                return 'urease'
    if feature.get('locus_tag') in candidate_tags:
        return 'candidate'
    return 'other'


# ---------------------------------------------------------------------------
# PT-BR: Camada 5 - montagem do HTML / EN-US: Layer 5 - HTML assembly
# ---------------------------------------------------------------------------

SECTIONS = [
    ('fluxo', ('Fluxo do pipeline', 'Pipeline flow')),
    ('entrada', ('Entrada e proveniência', 'Input and provenance')),
    ('qualidade', ('Qualidade dos genomas', 'Genome quality')),
    ('anotação', ('Anotação estrutural', 'Structural annotation')),
    ('triagem', ('Triagem em três camadas', 'Three-layer screening')),
    ('matriz', ('Presença e ausência', 'Presence and absence')),
    ('candidatos', ('Candidatos e evidências', 'Candidates and evidence')),
    ('pfam', ('Perfis Pfam', 'Pfam profiles')),
    ('sintenia', ('Sintenia do locus', 'Locus synteny')),
    ('filogenia', ('Filogenia de UreC', 'UreC phylogeny')),
    ('execução', ('Execução', 'Execution')),
    ('arquivos', ('Índice de arquivos', 'File index')),
]


def card(title, body, sub=None, extra_class=''):
    out = ['<section class="card %s">' % esc(extra_class)]
    out.append('<h3>%s</h3>' % T(title))
    if sub:
        out.append('<p class="sub">%s</p>' % T(sub))
    out.append(body)
    out.append('</section>')
    return ''.join(out)


def legend(entries):
    out = ['<ul class="legend">']
    for colour, text in entries:
        out.append('<li><span class="swatch" style="background:%s"></span>%s</li>' % (colour, T(text)))
    out.append('</ul>')
    return ''.join(out)


def table(headers, rows, cls='data', sortable=True, caption=None):
    out = ['<div class="table-wrap"><table class="%s"%s>' % (esc(cls),
                                                             ' data-sortable="1"' if sortable else '')]
    if caption:
        out.append('<caption>%s</caption>' % T(caption))
    out.append('<thead><tr>')
    for header in headers:
        out.append('<th scope="col">%s</th>' % T(header))
    out.append('</tr></thead><tbody>')
    for row in rows:
        out.append('<tr>' + ''.join(row) + '</tr>')
    out.append('</tbody></table></div>')
    return ''.join(out)


def opt_cell(source, key, formatter, cls='right'):
    """
    PT-BR: Celula de metrica opcional. Ausencia vira travessao, nunca zero -- um QUAST que
           não rodou não pode ser lido como um genoma de tamanho zero.
    EN-US: Optional metric cell. A missing value renders as a dash, never as zero -- a QUAST
           that never ran must not read as a zero-length genome.
    """
    raw = (source or {}).get(key, '')
    if raw is None or str(raw).strip() == '':
        return cell('<span class="muted">—</span>', -1, cls)
    return cell(formatter(raw), to_float(raw, -1), cls)


def cell(content, sort_value=None, cls=''):
    attrs = ''
    if sort_value is not None:
        attrs += ' data-sort="%s"' % esc(sort_value)
    if cls:
        attrs += ' class="%s"' % esc(cls)
    return '<td%s>%s</td>' % (attrs, content)


# ---------------------------------------------------------------------------

def build_context(results, accessions_path, out_dir, warn):
    ctx = {'results': results, 'out_dir': out_dir}
    strains = read_strains(results, accessions_path, warn)
    if not strains:
        warn('Nenhum genoma encontrado em 00_genomes / No genome found in 00_genomes')
    read_checkm2(results, strains)
    read_quast(results, strains)
    read_bakta(results, strains)
    ctx['strains'] = strains

    ctx['preflight'] = read_text(os.path.join(results, '00_preflight', 'preflight_status.txt'))
    ctx['candidates'] = read_tsv(os.path.join(results, '04_summary', 'summary_unified.tsv'))
    matrix_rows = read_tsv(os.path.join(results, '04_summary', 'urease_presence_absence.tsv'))
    ctx['matrix'] = sorted(matrix_rows, key=lambda r: natural_key(r.get('strain', '')))
    ctx['blast'], ctx['blast_total'] = read_blast(results, strains)
    ctx['locus_status'] = read_locus_status(results, strains)
    ctx['loci'] = read_loci(results, strains)
    ctx['phylo_status'] = read_text(os.path.join(results, '06_phylogeny', 'phylogeny_status.txt')).strip()
    ctx['iqtree'] = read_iqtree_stats(os.path.join(results, '06_phylogeny', 'ureC.iqtree'))
    ctx['tree'] = parse_newick(read_text(os.path.join(results, '06_phylogeny', 'ureC.treefile')))
    ctx['tree_labels'] = read_tsv(os.path.join(results, '06_phylogeny', 'ureC_labels.tsv'))
    ctx['tasks'] = read_trace(os.path.join(results, 'nf_trace.txt'))

    refs_count = 0
    match = re.search(r'urease_references\.fasta\s*\((\d+)\s*seqs\)', ctx['preflight'])
    if match:
        refs_count = to_int(match.group(1), 0)
    ctx['references_count'] = refs_count

    ctx['candidate_tags'] = {}
    for row in ctx['candidates']:
        ctx['candidate_tags'].setdefault(row.get('strain', ''), set()).add(row.get('feature_id', ''))
    return ctx


def rel(ctx, *parts):
    """PT-BR: Caminho relativo do relatório até um artefato. EN-US: Relative path from the report to an artifact."""
    target = os.path.join(ctx['results'], *parts)
    try:
        return esc(os.path.relpath(target, ctx['out_dir']))
    except ValueError:
        return esc(target)


def link(ctx, label, *parts):
    target = os.path.join(ctx['results'], *parts)
    if not os.path.exists(target):
        return '<span class="missing">%s</span>' % T(label)
    return '<a class="artifact" href="%s">%s</a>' % (rel(ctx, *parts), T(label))


# --- PT-BR: seções / EN-US: sections ---------------------------------------

def section_hero(ctx):
    strains = ctx['strains']
    proteins = sum(s.get('proteins', 0) for s in strains.values())
    contigs = sum(s.get('contigs', 0) for s in strains.values())
    high = sum(1 for r in ctx['candidates'] if r.get('confidence', '').startswith('Alta'))
    positive = sum(1 for row in ctx['matrix']
                   if any(row.get(g) == '1' for g in UREASE_GENES))
    dates = sorted(t['submit'] for t in ctx['tasks'] if t['submit'])
    window = ''
    if dates:
        window = '%s → %s' % (dates[0][:16], dates[-1][:16])

    tiles = [
        (num(len(strains)), ('genomas analisados', 'genomes analysed')),
        (num(contigs), ('contigs', 'contigs')),
        (num(proteins), ('proteínas preditas (Bakta)', 'predicted proteins (Bakta)')),
        (num(len(ctx['candidates'])), ('candidatos a urease', 'urease candidates')),
        (num(high), ('candidatos de alta confiança', 'high-confidence candidates')),
        (num(positive), ('estirpes urease-positivas', 'urease-positive strains')),
    ]
    out = ['<header class="hero" id="topo">']
    out.append('<p class="eyebrow">%s</p>' % T(('Relatório do pipeline BAFES-Urease',
                                                'BAFES-Urease pipeline report')))
    out.append('<h1>%s</h1>' % T(ctx['title']))
    out.append('<p class="lead">%s</p>' % T((
        'Mineração genômica da via da urease em %s genomas, do acesso GenBank até o locus '
        'sintênico e a filogenia de UreC.' % num(len(strains)),
        'Genomic mining of the urease pathway across %s genomes, from the GenBank accession '
        'to the syntenic locus and the UreC phylogeny.' % num(len(strains)))))
    if window:
        out.append('<p class="meta">%s <span class="mono">%s</span></p>'
                   % (T(('Janela de execução:', 'Execution window:')), esc(window)))
    out.append('<div class="tiles">')
    for value, label in tiles:
        out.append('<div class="tile"><span class="tile-value">%s</span>%s</div>'
                   % (esc(value), T(label, tag='span', cls='tile-label')))
    out.append('</div></header>')
    return ''.join(out)


PIPELINE_NODES = [
    # id, column, row, is_input
    ('accessions', 1, 0, True),
    ('PREFLIGHT_CHECK', 0, 1, False),
    ('DOWNLOAD_GENOME', 1, 1, False),
    ('QC_QUAST', 0, 2, False),
    ('QC_CHECKM2', 1, 2, False),
    ('BAKTA_ANNOTATE', 2, 2, False),
    ('EXTRACT_CANDIDATES', 2, 3, False),
    ('MERGE_ALL_STRAINS', 1, 4, False),
    ('EXTRACT_LOCUS', 2, 4, False),
    ('BUILD_PHYLOGENY', 3, 4, False),
    ('ANALYZE_SYNTENY', 2, 5, False),
]

PIPELINE_EDGES = [
    ('accessions', 'PREFLIGHT_CHECK'), ('accessions', 'DOWNLOAD_GENOME'),
    ('DOWNLOAD_GENOME', 'QC_QUAST'), ('DOWNLOAD_GENOME', 'QC_CHECKM2'),
    ('DOWNLOAD_GENOME', 'BAKTA_ANNOTATE'), ('BAKTA_ANNOTATE', 'EXTRACT_CANDIDATES'),
    ('EXTRACT_CANDIDATES', 'MERGE_ALL_STRAINS'), ('EXTRACT_CANDIDATES', 'EXTRACT_LOCUS'),
    ('EXTRACT_CANDIDATES', 'BUILD_PHYLOGENY'), ('EXTRACT_LOCUS', 'ANALYZE_SYNTENY'),
]


def section_flow(ctx):
    strains = ctx['strains']
    by_process = {}
    for task in ctx['tasks']:
        by_process.setdefault(task['process'], []).append(task)

    proteins = sum(s.get('proteins', 0) for s in strains.values())
    loci_regions = sum(len(v) for v in ctx['loci'].values())
    details = {
        'accessions': (('data/accessions.tsv', 'data/accessions.tsv'),
                       ('%s acessos GenBank' % num(len(strains)),
                        '%s GenBank accessions' % num(len(strains))), 'entrada'),
        'PREFLIGHT_CHECK': (None, ('5 pre-checks', '5 pre-checks'), 'entrada'),
        'DOWNLOAD_GENOME': (None, ('%s genomas · %s contigs' % (num(len(strains)),
                                                                     num(sum(s['contigs'] for s in strains.values()))),
                                   '%s genomes · %s contigs' % (num(len(strains)),
                                                                     num(sum(s['contigs'] for s in strains.values())))),
                            'entrada'),
        'QC_QUAST': (None, ('métricas de montagem', 'assembly metrics'), 'qualidade'),
        'QC_CHECKM2': (None, ('completude e contaminação', 'completeness and contamination'), 'qualidade'),
        'BAKTA_ANNOTATE': (None, ('%s proteínas' % num(proteins), '%s proteins' % num(proteins)), 'anotação'),
        'EXTRACT_CANDIDATES': (None, ('3 camadas → %s candidatos' % num(len(ctx['candidates'])),
                                      '3 layers → %s candidates' % num(len(ctx['candidates']))),
                               'triagem'),
        'MERGE_ALL_STRAINS': (None, ('matriz 10×11', '10×11 matrix'), 'matriz'),
        'EXTRACT_LOCUS': (None, ('%s regiões' % num(loci_regions), '%s regions' % num(loci_regions)), 'sintenia'),
        'BUILD_PHYLOGENY': (None, ('%s táxons UreC' % num(ctx['iqtree'].get('taxa', 0)),
                                   '%s UreC taxa' % num(ctx['iqtree'].get('taxa', 0))), 'filogenia'),
        'ANALYZE_SYNTENY': (None, ('clinker', 'clinker'), 'sintenia'),
    }

    col_w, col_gap, box_h, row_gap = 214, 22, 62, 34
    left, top = 12, 14
    width = left * 2 + 4 * col_w + 3 * col_gap
    rows = max(node[2] for node in PIPELINE_NODES) + 1
    height = top * 2 + rows * box_h + (rows - 1) * row_gap

    def box_rect(node_id):
        for nid, col, row, is_input in PIPELINE_NODES:
            if nid == node_id:
                x = left + col * (col_w + col_gap)
                y = top + row * (box_h + row_gap)
                return x, y, col_w, box_h, is_input
        return None

    out = [svg_open(width, height, 'chart flow')]
    out.append('<defs><marker id="arrow" viewBox="0 0 10 10" refX="9" refY="5" '
               'markerWidth="6" markerHeight="6" orient="auto-start-reverse">'
               '<path d="M0,0 L10,5 L0,10 z" fill="var(--baseline)"/></marker></defs>')

    for src, dst in PIPELINE_EDGES:
        sx, sy, sw, sh, _ = box_rect(src)
        dx, dy, dw, dh, _ = box_rect(dst)
        x1, y1 = sx + sw / 2.0, sy + sh
        x2, y2 = dx + dw / 2.0, dy
        my = (y1 + y2) / 2.0
        path = 'M%g,%g L%g,%g L%g,%g L%g,%g' % (x1, y1, x1, my, x2, my, x2, y2 - 2)
        out.append(svg_path(path, 'none', 'var(--baseline)', 1.4,
                            extra=' marker-end="url(#arrow)"'))

    for node_id, col, row, is_input in PIPELINE_NODES:
        x, y, w, h, _ = box_rect(node_id)
        title, subtitle, anchor = details[node_id]
        tasks = by_process.get(node_id, [])
        ok = tasks and all(t['status'] in ('COMPLETED', 'CACHED') for t in tasks)
        fill = 'var(--flow-input)' if is_input else 'var(--flow-box)'
        stroke = 'var(--good)' if ok else 'var(--grid)'
        out.append('<a href="#%s">' % esc(anchor))
        out.append(svg_rect(x, y, w, h, fill, 8,
                            extra=' stroke="%s" stroke-width="1.5"%s'
                                  % (stroke, ' stroke-dasharray="4 3"' if is_input else '')))
        label = title if title else (node_id, node_id)
        attrs, text = bi_attrs(label)
        out.append('<text x="%g" y="%g" font-size="11.5" font-family="ui-monospace,monospace" '
                   'font-weight="700" fill="var(--text-primary)"%s>%s</text>'
                   % (x + 12, y + 24, attrs, text))
        out.append(svg_text(x + 12, y + 43, subtitle, 10, 'var(--text-secondary)'))
        if tasks:
            out.append(svg_text(x + w - 12, y + 24,
                                '%s×' % num(len(tasks)), 10, 'var(--muted)', 'end'))
        out.append('</a>')
    out.append(SVG_CLOSE)

    body = ''.join(out)
    body += legend([('var(--good)', ('processo concluído sem falha', 'process completed with no failure')),
                    ('var(--flow-input)', ('entrada do usuário', 'user input'))])
    body += ('<p class="note">%s</p>'
             % T(('Cada caixa é um processo do main.nf; clique para ir a seção correspondente. '
                  'As bases de dados (Pfam-A, banco do Bakta v6.0 e as referências curadas do UniProt) '
                  'entram nos processos que as consomem e estão detalhadas em Entrada e proveniência.',
                  'Each box is a main.nf process; click to jump to the matching section. '
                  'The databases (Pfam-A, the Bakta v6.0 database and the curated UniProt references) '
                  'feed the processes that consume them and are detailed under Input and provenance.')))
    return card(('Do acesso GenBank a árvore de UreC', 'From GenBank accession to the UreC tree'),
                body,
                ('%s processos do Nextflow, do download do genoma aos resultados agregados.'
                 % num(len(PIPELINE_ORDER)),
                 '%s Nextflow processes, from genome download to the aggregated results.'
                 % num(len(PIPELINE_ORDER))))


def section_input(ctx):
    strains = ctx['strains']
    rows = []
    for code, strain in strains.items():
        rows.append([
            cell(label_html(strain), plain(label_full(strain))),
            cell('<span class="mono muted">%s</span>' % esc(code), code),
            cell('<span class="mono">%s</span>' % esc(strain['accession'] or '-'), strain['accession']),
            cell(num(strain['contigs']), strain['contigs'], 'right'),
            opt_cell(strain.get('quast'), 'Total length', human_bp),
            cell(link(ctx, ('genoma', 'genome'), '00_genomes', '%s.fna' % code)),
        ])
    body = table([('Organismo', 'Organism'), ('ID interno', 'Internal ID'),
                  ('Acesso GenBank', 'GenBank accession'), ('Contigs', 'Contigs'),
                  ('Tamanho', 'Size'), ('Arquivo', 'File')], rows)
    # PT-BR: O texto i18n e trocado via textContent, então marcacao inline não sobreviveria a
    #        troca de idioma: os trechos em código ficam fora dos nós traduzidos.
    # EN-US: i18n text is swapped via textContent, so inline markup would not survive a language
    #        switch: the code fragments stay outside the translated nodes.
    body += ('<p class="note warn">%s <code>species</code> %s <code>data/accessions.tsv</code>%s '
             '<code>Bacillus sp. Sn</code>%s</p>' % (
                 T(('A taxonomia acima foi lida da linha de definição do RefSeq de cada genoma '
                    'baixado, e não da coluna',
                    'The taxonomy above was read from the RefSeq definition line of each '
                    'downloaded genome, not from the column')),
                 T(('de', 'of')),
                 T((': naquele arquivo, nove das dez linhas trazem o rótulo genérico',
                    ': in that file, nine of the ten rows carry the generic label')),
                 T(('. Quatro das estirpes nem sequer são Bacillus: Heyndrickxia, Peribacillus, '
                    'Paenibacillus e Brevibacillus.',
                    '. Four of the strains are not even Bacillus: Heyndrickxia, Peribacillus, '
                    'Paenibacillus and Brevibacillus.'))))

    checks = []
    for line in ctx['preflight'].splitlines():
        stripped = line.strip()
        if stripped.startswith('>>>'):
            checks.append({'title': stripped.lstrip('> ').strip(), 'ok': 0, 'retry': 0, 'items': []})
        elif checks and stripped.startswith('[OK]'):
            checks[-1]['ok'] += 1
            checks[-1]['items'].append(stripped[4:].strip())
        elif checks and stripped.startswith('[RETRY'):
            checks[-1]['retry'] += 1
    pre_body = ['<ul class="checks">']
    for check in checks:
        # PT-BR: O título já vem bilíngue no arquivo, separado por ' / ' depois do prefixo.
        # EN-US: The title is already bilingual in the file, separated by ' / ' after the prefix.
        prefix, _, rest = check['title'].partition(']')
        halves = rest.strip().rstrip('.').split(' / ')
        title = (('%s] %s' % (prefix, pt_display(halves[0])), '%s] %s' % (prefix, halves[1]))
                 if len(halves) == 2 else check['title'])
        pre_body.append('<li><span class="dot ok"></span>%s'
                        '<span class="muted"> · %s [OK]%s</span></li>'
                        % (T(title, cls='mono'), num(check['ok']),
                           ('· %s [RETRY]' % num(check['retry'])) if check['retry'] else ''))
    pre_body.append('</ul>')
    success = '[SUCESSO/SUCCESS]' in ctx['preflight']
    pre_body.append('<p class="note">%s %s</p>'
                    % (T(('Referências curadas do UniProt (Firmicutes, revisadas):',
                          'Curated UniProt references (Firmicutes, reviewed):')),
                       num(ctx['references_count'])))
    if success:
        pre_body.append('<p class="ok-line"><span class="dot ok"></span>%s</p>'
                        % T(('Todos os pre-checks passaram.', 'All pre-checks passed.')))
    pre_body.append('<p class="links">%s</p>'
                    % link(ctx, ('preflight_status.txt', 'preflight_status.txt'),
                           '00_preflight', 'preflight_status.txt'))

    return (card(('Os %s organismos' % num(len(strains)),
                  'The %s organisms' % num(len(strains))), body,
                 ('Estirpes depositadas pelo grupo, recuperadas do NCBI pelo acesso WGS master.',
                  'Strains deposited by the group, retrieved from NCBI through the WGS master accession.'))
            + card(('Verificação previa de recursos', 'Resource pre-flight check'), ''.join(pre_body),
                   ('Executado antes do pipeline: acessos no NCBI, Pfam-A, banco do Bakta, '
                    'referências e ferramentas.',
                    'Run before the pipeline: NCBI accessions, Pfam-A, Bakta database, references and tools.')))


def section_quality(ctx):
    strains = ctx['strains']
    rows = []
    for code, strain in strains.items():
        quast = strain.get('quast', {})
        checkm2 = strain.get('checkm2', {})
        rows.append([
            cell(label_html(strain, short=True), plain(label_full(strain)), 'nowrap'),
            opt_cell(quast, '# contigs', lambda v: num(v)),
            opt_cell(quast, 'Total length', human_bp),
            opt_cell(quast, 'N50', human_bp),
            opt_cell(quast, 'Largest contig', human_bp),
            opt_cell(quast, 'GC (%)', lambda v: '%s%%' % num(v, 2)),
            opt_cell(checkm2, 'Completeness', lambda v: '%s%%' % num(v, 1)),
            opt_cell(checkm2, 'Contamination', lambda v: '%s%%' % num(v, 2)),
            cell(link(ctx, ('QUAST', 'QUAST'), '01_qc', 'quast', code, 'report.html')),
        ])
    tbl = table([('Organismo', 'Organism'), ('Contigs', 'Contigs'), ('Tamanho', 'Size'),
                 ('N50', 'N50'), ('Maior contig', 'Largest contig'), ('GC', 'GC'),
                 ('Completude', 'Completeness'), ('Contaminação', 'Contamination'),
                 ('Relatório', 'Report')], rows)

    positive = urease_positive_codes(ctx)
    points = []
    for code, strain in strains.items():
        quast = strain.get('quast', {})
        size = to_float(quast.get('Total length'), 0.0) or 0.0
        gc = to_float(quast.get('GC (%)'), 0.0) or 0.0
        is_pos = code in positive
        points.append({
            'x': size / 1e6, 'y': gc,
            'colour': 'var(--series-1)' if is_pos else 'var(--gene-other)',
            'binomial': '%s. %s' % (strain['genus'][0], strain['species']),
            'suffix': strain['designation'],
            'title': '%s · %s · GC %.2f%%' % (label_full(strain), human_bp(size), gc),
        })
    scatter = scatter_chart(points, ('Tamanho do genoma (Mb)', 'Genome size (Mb)'),
                            ('Conteúdo GC (%)', 'GC content (%)'))
    scatter += legend([('var(--series-1)', ('estirpe urease-positiva', 'urease-positive strain')),
                       ('var(--gene-other)', ('estirpe urease-negativa', 'urease-negative strain'))])

    n50_items = []
    for code, strain in strains.items():
        n50 = to_float(strain.get('quast', {}).get('N50'), 0.0) or 0.0
        n50_items.append({
            'value': n50,
            'label': label_short(strain),
            'label_html': make_binomial_label(strain),
            'title': '%s · N50 %s' % (label_full(strain), human_bp(n50)),
        })
    n50_items.sort(key=lambda i: i['value'], reverse=True)
    n50 = hbar_chart(n50_items, ('pares de base', 'base pairs'), value_fmt=human_bp,
                     tick_format=human_bp)

    contam_items = []
    for code, strain in strains.items():
        value = to_float(strain.get('checkm2', {}).get('Contamination'), 0.0) or 0.0
        contam_items.append({
            'value': value,
            'label': label_short(strain),
            'label_html': make_binomial_label(strain),
            'title': '%s · %.2f%%' % (label_full(strain), value),
        })
    contam_items.sort(key=lambda i: i['value'], reverse=True)
    contam = hbar_chart(contam_items, ('% de contaminação', '% contamination'),
                        value_fmt=lambda v: '%.2f%%' % v)

    # PT-BR: A frase de resumo e derivada dos valores lidos, nunca afirmada de antemao.
    # EN-US: The summary sentence is derived from the values read, never asserted up front.
    completeness = [to_float(s.get('checkm2', {}).get('Completeness'), None) for s in strains.values()]
    completeness = [c for c in completeness if c is not None]
    contamination = [to_float(s.get('checkm2', {}).get('Contamination'), None) for s in strains.values()]
    contamination = [c for c in contamination if c is not None]
    if completeness:
        lo, hi, n = num(min(completeness), 1), num(max(completeness), 1), num(len(completeness))
        span_pt = 'de %s%%' % lo if lo == hi else 'entre %s%% e %s%%' % (lo, hi)
        span_en = 'of %s%%' % lo if lo == hi else 'between %s%% and %s%%' % (lo, hi)
        qc_sub = ('QUAST e CheckM2 por estirpe. Completude %s nas %s montagens avaliadas.'
                  % (span_pt, n),
                  'QUAST and CheckM2 per strain. Completeness %s across the %s assemblies evaluated.'
                  % (span_en, n))
    else:
        qc_sub = ('Sem resultados de CheckM2 nesta árvore de resultados.',
                  'No CheckM2 results in this results tree.')
    contam_sub = None
    if contamination:
        contam_sub = ('Nenhuma montagem passa de %s%% de contaminação.' % num(max(contamination), 2),
                      'No assembly exceeds %s%% contamination.' % num(max(contamination), 2))

    return (card(('Métricas de montagem e completude', 'Assembly and completeness metrics'), tbl,
                 qc_sub)
            + card(('Tamanho do genoma e conteúdo GC', 'Genome size and GC content'), scatter,
                   ('Cada ponto é um genoma; a cor marca o resultado final da triagem de urease.',
                    'Each dot is a genome; colour marks the final urease screening outcome.'))
            + card(('Contiguidade (N50)', 'Contiguity (N50)'), n50)
            + card(('Contaminação estimada (CheckM2)', 'Estimated contamination (CheckM2)'),
                   contam, contam_sub))


def make_binomial_label(strain):
    """PT-BR: Rótulo de eixo com binomial em itálico. EN-US: Axis label with italic binomial."""
    def render(x, y):
        return svg_binomial(x, y, '%s. %s' % (strain['genus'][0], strain['species']),
                            strain['designation'], 11, 'var(--text-primary)', 'end')
    return render


def urease_positive_codes(ctx):
    codes = set()
    for row in ctx['matrix']:
        if any(row.get(g) == '1' for g in UREASE_GENES):
            codes.add(row.get('strain', ''))
    return codes


def section_annotation(ctx):
    strains = ctx['strains']
    items = []
    for code, strain in strains.items():
        items.append({
            'value': strain.get('proteins', 0),
            'label': label_short(strain),
            'label_html': make_binomial_label(strain),
            'title': '%s · %s %s' % (label_full(strain), num(strain.get('proteins', 0)),
                                          plain(('proteínas', 'proteins'))),
        })
    items.sort(key=lambda i: i['value'], reverse=True)
    chart = hbar_chart(items, ('proteínas preditas', 'predicted proteins'))
    total = sum(s.get('proteins', 0) for s in strains.values())
    body = chart + ('<p class="note">%s <strong>%s</strong>. %s</p>'
                    % (T(('Total de proteínas anotadas pelo Bakta:',
                          'Total proteins annotated by Bakta:')), esc(num(total)),
                       T(('Este é o universo sobre o qual a triagem de urease foi aplicada.',
                          'This is the universe the urease screening was applied to.'))))
    files = []
    for code in strains:
        files.append('<li>%s %s %s %s</li>' % (
            link(ctx, ('%s.faa' % code, '%s.faa' % code), '02_bakta', code, '%s.faa' % code),
            link(ctx, ('%s.gff3' % code, '%s.gff3' % code), '02_bakta', code, '%s.gff3' % code),
            link(ctx, ('%s.gbff' % code, '%s.gbff' % code), '02_bakta', code, '%s.gbff' % code),
            link(ctx, ('%s.json' % code, '%s.json' % code), '02_bakta', code, '%s.json' % code)))
    body += '<details class="files"><summary>%s</summary><ul>%s</ul></details>' % (
        T(('Arquivos de anotação por estirpe', 'Per-strain annotation files')), ''.join(files))
    return card(('Proteínas preditas por genoma (Bakta v1.12.0, DB v6.0)',
                 'Predicted proteins per genome (Bakta v1.12.0, DB v6.0)'), body)


def section_screening(ctx):
    proteins = sum(s.get('proteins', 0) for s in ctx['strains'].values())
    blast_queries = len(ctx['blast'])
    candidates = ctx['candidates']
    bakta_layer = sum(1 for r in candidates if r.get('layer_bakta') == 'True')
    blast_layer = sum(1 for r in candidates if r.get('layer_blast') == 'True')
    hmm_layer = sum(1 for r in candidates if r.get('layer_hmmer') == 'True')
    high = sum(1 for r in candidates if r.get('confidence', '').startswith('Alta'))

    # PT-BR: So etapas realmente aninhadas entram no funil. As três camadas não são aninhadas
    #        (38 dos 58 candidatos vieram apenas do HMMER), então viram um gráfico próprio.
    # EN-US: Only genuinely nested stages go into the funnel. The three layers are not nested
    #        (38 of the 58 candidates came from HMMER alone), so they get their own chart.
    stages = [
        {'value': proteins, 'label': ('Proteínas preditas', 'Predicted proteins'),
         'note': ('todas as CDS anotadas pelo Bakta nos dez genomas',
                  'every CDS annotated by Bakta across the ten genomes')},
        {'value': len(candidates), 'label': ('Candidatos retidos', 'Retained candidates'),
         'note': ('união das três camadas de evidência', 'union of the three evidence layers')},
        {'value': high, 'label': ('Alta confiança', 'High confidence'),
         'note': ('duas ou mais camadas concordantes', 'two or more agreeing layers')},
    ]
    funnel = funnel_chart(stages)
    funnel += ('<p class="note">%s</p>' % T((
        'A largura das faixas usa escala logarítmica: a queda vai de dezenas de milhares para '
        'dezesseis, e em escala linear as duas últimas faixas seriam invisíveis. O valor exato e '
        'a fração retida estão impressos ao lado de cada faixa, então a leitura não depende da '
        'largura. Cada etapa é subconjunto estrito da anterior.',
        'Band width uses a logarithmic scale: the drop goes from tens of thousands down to sixteen, '
        'and on a linear scale the last two bands would be invisible. The exact value and the '
        'retained fraction are printed beside each band, so reading it does not depend on width. '
        'Each stage is a strict subset of the previous one.')))

    lane_items = [
        {'value': hmm_layer, 'label': 'HMMER',
         'title': '%s: %s %s' % ('HMMER', num(hmm_layer), plain(('candidatos', 'candidates')))},
        {'value': bakta_layer, 'label': 'Bakta',
         'title': '%s: %s %s' % ('Bakta', num(bakta_layer), plain(('candidatos', 'candidates')))},
        {'value': blast_layer, 'label': 'BLASTp',
         'title': '%s: %s %s' % ('BLASTp', num(blast_layer), plain(('candidatos', 'candidates')))},
    ]
    lanes = hbar_chart(lane_items, ('candidatos apoiados', 'supported candidates'), label_width=120)
    lanes += ('<p class="note">%s</p>' % T((
        'As três camadas são filtros independentes aplicados ao mesmo conjunto de proteínas, '
        'não etapas encadeadas: as barras se sobrepõem e somam mais que os %s candidatos, porque '
        'um mesmo candidato pode ser apoiado por duas ou três camadas. O BLASTp produziu %s '
        'alinhamentos sobre %s proteínas distintas antes dos limiares de identidade e cobertura.'
        % (num(len(candidates)), num(ctx['blast_total']), num(blast_queries)),
        'The three layers are independent filters applied to the same protein set, not chained '
        'stages: the bars overlap and add up to more than the %s candidates, because one candidate '
        'may be supported by two or three layers. BLASTp produced %s alignments over %s distinct '
        'proteins before the identity and coverage thresholds.'
        % (num(len(candidates)), num(ctx['blast_total']), num(blast_queries)))))

    combos = Counter()
    for row in candidates:
        key = (row.get('layer_bakta') == 'True', row.get('layer_blast') == 'True',
               row.get('layer_hmmer') == 'True')
        combos[key] += 1
    combo_items = []
    for key, count in combos.most_common():
        names = []
        for flag, name in zip(key, ('Bakta', 'BLASTp', 'HMMER')):
            if flag:
                names.append(name)
        joined = ' + '.join(names) if names else '-'
        combo_items.append({
            'value': count,
            'label': joined,
            'colour': 'var(--series-1)' if sum(key) >= 2 else 'var(--gene-other)',
            'title': '%s: %s %s' % (joined, count, plain(('candidatos', 'candidates'))),
        })
    combo = hbar_chart(combo_items, ('candidatos', 'candidates'), label_width=170)
    combo += legend([('var(--series-1)', ('alta confiança (≥ 2 camadas)', 'high confidence (≥ 2 layers)')),
                     ('var(--gene-other)', ('confiança média (1 camada)', 'medium confidence (1 layer)'))])

    criteria = ('<div class="criteria">'
                '<div><h4>%s</h4><p>%s</p></div>'
                '<div><h4>%s</h4><p>%s</p></div>'
                '<div><h4>%s</h4><p>%s</p></div></div>') % (
        T(('Camada 1 · Bakta', 'Layer 1 · Bakta')),
        T(('Símbolo do gene na lista alvo (ureA–ureH, urtA, uc, ah), número EC 3.5.1.5 / '
           '6.3.4.6 / 3.5.1.54, ou produto contendo urease, urea ou allophanate.',
           'Gene symbol in the target list (ureA–ureH, urtA, uc, ah), EC number 3.5.1.5 / '
           '6.3.4.6 / 3.5.1.54, or product containing urease, urea or allophanate.')),
        T(('Camada 2 · BLASTp', 'Layer 2 · BLASTp')),
        T(('Alinhamento contra as referências curadas do UniProt com E ≤ 1e-5, identidade '
           '≥ 30% e cobertura da consulta ≥ 50%.',
           'Alignment against the curated UniProt references with E ≤ 1e-5, identity '
           '≥ 30% and query coverage ≥ 50%.')),
        T(('Camada 3 · HMMER', 'Layer 3 · HMMER')),
        T(('hmmscan contra o Pfam-A completo, retendo apenas domínios com score ≥ 20 cujo '
           'perfil esta na lista curada de dez acessos da via da urease.',
           'hmmscan against the full Pfam-A, keeping only domains with score ≥ 20 whose '
           'profile is in the curated list of ten urease-pathway accessions.')))

    return (card(('Como um genoma inteiro vira dezesseis proteínas',
                  'How a whole genome becomes sixteen proteins'), criteria,
                 ('Três critérios independentes; a confiança é alta quando pelo menos dois concordam.',
                  'Three independent criteria; confidence is high when at least two agree.'))
            + card(('Funil da triagem', 'Screening funnel'), funnel)
            + card(('Rendimento de cada camada, isoladamente',
                    'Yield of each layer, on its own'), lanes)
            + card(('Quais camadas sustentam cada candidato', 'Which layers support each candidate'),
                   combo))


def section_matrix(ctx):
    strains = ctx['strains']
    entries = []
    for row in ctx['matrix']:
        code = row.get('strain', '')
        strain = strains.get(code)
        if strain is None:
            continue
        entries.append({
            'strain': code,
            'binomial': '%s. %s' % (strain['genus'][0], strain['species']),
            'suffix': strain['designation'],
            'full': label_full(strain),
            'values': {g: row.get(g, '0') for g in UREASE_GENES},
        })
    heat = presence_heatmap(strains, entries)
    heat += legend([('var(--series-1)', ('gene presente', 'gene present')),
                    ('var(--cell-empty)', ('gene ausente', 'gene absent'))])

    rows = []
    for entry in entries:
        cells = [cell('<i>%s</i> %s' % (esc(entry['binomial']), esc(entry['suffix'])), entry['full']),
                 cell('<span class="mono muted">%s</span>' % esc(entry['strain']), entry['strain'])]
        total = 0
        for gene in UREASE_GENES:
            present = entry['values'].get(gene) == '1'
            total += 1 if present else 0
            cells.append(cell('<span class="bit %s">%s</span>' % ('on' if present else 'off',
                                                                  '1' if present else '0'),
                              '1' if present else '0', 'center'))
        cells.append(cell('<strong>%s</strong>' % num(total), total, 'right'))
        rows.append(cells)
    headers = [('Organismo', 'Organism'), ('ID', 'ID')] + [(g, g) for g in UREASE_GENES] + \
              [('Total', 'Total')]
    tbl = table(headers, rows)

    readings = []
    for entry in entries:
        present = [g for g in UREASE_GENES if entry['values'].get(g) == '1']
        if not present:
            continue
        readings.append('<li><i>%s</i> %s — <span class="mono">%s</span></li>'
                        % (esc(entry['binomial']), esc(entry['suffix']), esc(', '.join(present))))
    negatives = [e for e in entries if not any(e['values'].get(g) == '1' for g in UREASE_GENES)]
    reading = '<div class="reading"><h4>%s</h4><ul>%s</ul><p>%s <strong>%s</strong> %s</p></div>' % (
        T(('Leitura biológica', 'Biological reading')),
        ''.join(readings),
        T(('Urease-negativas:', 'Urease-negative:')),
        esc(num(len(negatives))),
        '· ' + ', '.join('<i>%s</i> %s' % (esc(e['binomial']), esc(e['suffix']))
                              for e in negatives))

    return (card(('Matriz de presença e ausência dos genes da via da urease',
                  'Presence/absence matrix of the urease pathway genes'),
                 heat + reading,
                 ('Presença exige igualdade do símbolo do gene ou um produto casando os padrões '
                  'curados — nunca substring, que marcaria "hydroxyurea kinase" como ureA.',
                  'Presence requires gene-symbol equality or a product matching the curated patterns '
                  '— never substring matching, which would score "hydroxyurea kinase" as ureA.'))
            + card(('Tabela da matriz', 'Matrix table'), tbl))


def layer_badges(row):
    out = []
    for field, name in (('layer_bakta', 'B'), ('layer_blast', 'L'), ('layer_hmmer', 'H')):
        on = row.get(field) == 'True'
        out.append('<span class="badge %s" title="%s">%s</span>'
                   % ('on' if on else 'off',
                      esc({'B': 'Bakta', 'L': 'BLASTp', 'H': 'HMMER'}[name]), name))
    return ''.join(out)


def section_candidates(ctx):
    strains = ctx['strains']
    blast = ctx['blast']
    rows_html = []
    for index, row in enumerate(ctx['candidates']):
        code = row.get('strain', '')
        strain = strains.get(code)
        organism = label_html(strain, short=True) if strain else esc(code)
        sort_org = label_full(strain) if strain else code
        conf = split_bilingual(row.get('confidence', ''))
        conf_cls = 'high' if str(row.get('confidence', '')).startswith('Alta') else 'medium'
        hits = blast.get(row.get('feature_id', ''), [])
        expandable = bool(hits)
        rows_html.append(
            '<tr class="cand %s"%s>' % ('expandable' if expandable else '',
                                        ' data-target="evi-%d"' % index if expandable else '')
            + cell(organism, sort_org, 'nowrap')
            + cell('<span class="mono">%s</span>' % esc(row.get('feature_id', '')),
                   row.get('feature_id', ''))
            + cell('<span class="mono gene">%s</span>' % esc(row.get('gene', '') or '—'),
                   row.get('gene', ''))
            + cell(esc(row.get('product', '') or '—'), row.get('product', ''))
            + cell(layer_badges(row), row.get('layer_bakta', '') + row.get('layer_blast', '')
                   + row.get('layer_hmmer', ''), 'center')
            + cell(('%s%%' % num(to_float(row.get('blast_identity'), 0.0), 1))
                   if row.get('blast_identity') else '—',
                   to_float(row.get('blast_identity'), -1), 'right')
            + cell(('%s%%' % num(to_float(row.get('blast_coverage'), 0.0), 0))
                   if row.get('blast_coverage') else '—',
                   to_float(row.get('blast_coverage'), -1), 'right')
            + cell('<span class="mono">%s</span>' % esc(fmt_evalue(row.get('blast_evalue')))
                   if row.get('blast_evalue') else '—',
                   to_float(row.get('blast_evalue'), 1e9), 'right')
            + cell('<span class="mono">%s</span>' % esc(row.get('hmmer_profile', '') or '—'),
                   row.get('hmmer_profile', ''))
            + cell(num(to_float(row.get('hmmer_score'), 0.0), 1) if row.get('hmmer_score') else '—',
                   to_float(row.get('hmmer_score'), -1), 'right')
            + cell('<span class="conf %s">%s</span>' % (conf_cls, T(conf)),
                   row.get('confidence', ''))
            + '</tr>')
        if expandable:
            shown = hits[:10]
            hit_rows = []
            for hit in shown:
                info = parse_reference_id(hit['sseqid'])
                hit_rows.append(
                    '<tr><td><span class="mono gene">%s</span></td>'
                    '<td><i>%s</i></td><td class="mono">%s</td>'
                    '<td class="right">%s%%</td><td class="right">%s</td>'
                    '<td class="right mono">%s</td><td class="right">%s</td></tr>'
                    % (esc(info['gene'] or '-'), esc(info['organism']), esc(info['uniprot']),
                       num(to_float(hit['pident'], 0.0), 1), num(to_int(hit['length'], 0)),
                       esc(fmt_evalue(hit['evalue'])), num(to_float(hit['bitscore'], 0.0), 0)))
            more = ''
            if len(hits) > len(shown):
                more = '<p class="note">%s</p>' % T(
                    ('Mostrando os %d melhores de %d alinhamentos.' % (len(shown), len(hits)),
                     'Showing the top %d of %d alignments.' % (len(shown), len(hits))))
            rows_html.append(
                '<tr class="evidence" id="evi-%d" hidden><td colspan="11">'
                '<h5>%s</h5><table class="inner"><thead><tr>'
                '<th>%s</th><th>%s</th><th>UniProt</th><th>%s</th><th>%s</th>'
                '<th>E-value</th><th>%s</th></tr></thead><tbody>%s</tbody></table>%s</td></tr>'
                % (index, T(('Alinhamentos BLASTp contra as referências curadas',
                             'BLASTp alignments against the curated references')),
                   T(('Gene', 'Gene')), T(('Organismo de referência', 'Reference organism')),
                   T(('Identidade', 'Identity')), T(('Extensao', 'Length')), T(('Bitscore', 'Bitscore')),
                   ''.join(hit_rows), more))

    ph_pt, ph_en = ('Filtrar por organismo, gene, produto, perfil...',
                    'Filter by organism, gene, product, profile...')
    lb_pt, lb_en = ('Filtrar candidatos', 'Filter candidates')
    controls = ('<div class="controls">'
                '<input type="search" class="filter" data-table="candidates" '
                'data-ph-pt="%s" data-ph-en="%s" data-al-pt="%s" data-al-en="%s" '
                'placeholder="%s" aria-label="%s">'
                '<label class="chk"><input type="checkbox" class="only-high" data-table="candidates"> %s</label>'
                '</div>') % (esc(ph_pt), esc(ph_en), esc(lb_pt), esc(lb_en),
                             esc(ph_pt if LANG == 'pt' else ph_en),
                             esc(lb_pt if LANG == 'pt' else lb_en),
                             T(('somente alta confiança', 'high confidence only')))

    headers = [('Organismo', 'Organism'), ('Locus tag', 'Locus tag'), ('Gene', 'Gene'),
               ('Produto', 'Product'), ('Camadas', 'Layers'), ('Identidade', 'Identity'),
               ('Cobertura', 'Coverage'), ('E-value', 'E-value'), ('Perfil Pfam', 'Pfam profile'),
               ('Score', 'Score'), ('Confiança', 'Confidence')]
    tbl = ['<div class="table-wrap"><table class="data" id="candidates" data-sortable="1"><thead><tr>']
    for header in headers:
        tbl.append('<th scope="col">%s</th>' % T(header))
    tbl.append('</tr></thead><tbody>')
    tbl.extend(rows_html)
    tbl.append('</tbody></table></div>')

    body = controls + ''.join(tbl)
    body += ('<p class="note">%s %s · %s</p>'
             % (T(('Clique em uma linha para ver as evidências de BLASTp.',
                   'Click a row to see the BLASTp evidence.')),
                T(('Total de alinhamentos carregados:', 'Total alignments loaded:')),
                esc(num(ctx['blast_total']))))
    body += '<p class="links">%s</p>' % link(ctx, ('summary_unified.tsv', 'summary_unified.tsv'),
                                             '04_summary', 'summary_unified.tsv')
    return card(('Os %s candidatos retidos' % num(len(ctx['candidates'])),
                 'The %s retained candidates' % num(len(ctx['candidates']))), body,
                ('Cada linha é uma proteína predita que passou por pelo menos uma das três camadas.',
                 'Each row is a predicted protein that cleared at least one of the three layers.'))


def section_pfam(ctx):
    counts = Counter()
    for row in ctx['candidates']:
        profile = (row.get('hmmer_profile') or '').strip()
        if profile:
            counts[profile] += 1
    items = []
    for profile, count in counts.most_common():
        core = profile in CORE_UREASE_PROFILES
        accession = next((acc for acc, name in UREASE_PFAM.items() if name == profile), '')
        items.append({
            'value': count, 'label': profile,
            'colour': 'var(--series-1)' if core else 'var(--series-2)',
            'title': '%s%s: %s' % (profile, ' (%s)' % accession if accession else '',
                                   num(count)),
        })
    chart = hbar_chart(items, ('candidatos', 'candidates'), label_width=150)
    chart += legend([('var(--series-1)', ('domínio estrutural da urease', 'urease structural domain')),
                     ('var(--series-2)', ('domínio acessório ou de metalochaperona',
                                          'accessory or metallochaperone domain'))])
    note = ('<p class="note">%s</p>' % T((
        'cobW e CT_C_D são perfis de metalochaperonas ligadas ao cobalto e ao níquel. Eles estão '
        'na lista curada porque a maturação da urease depende da inserção de níquel, mas sozinhos '
        'não provam função de urease — por isso quase todo candidato apoiado apenas por eles '
        'fica com confiança média. Nenhum deles é âncora de locus em extract_urease_locus.py.',
        'cobW and CT_C_D are cobalt- and nickel-linked metallochaperone profiles. They are in the '
        'curated list because urease maturation depends on nickel insertion, but on their own they '
        'do not prove urease function — which is why nearly every candidate supported only by '
        'them ends up at medium confidence. Neither is a locus anchor in extract_urease_locus.py.')))
    return card(('Distribuição dos perfis Pfam entre os candidatos',
                 'Distribution of Pfam profiles across the candidates'), chart + note)


def section_synteny(ctx):
    strains = ctx['strains']
    blocks = []
    for code, records in ctx['loci'].items():
        strain = strains.get(code)
        tags = ctx['candidate_tags'].get(code, set())
        header = ('<h4 class="locus-head">%s <span class="mono muted">%s</span></h4>'
                  % (label_html(strain) if strain else esc(code), esc(code)))
        maps = []
        for record in records:
            anchors = ', '.join(record['anchors'])
            maps.append('<div class="locus-block">%s<p class="note">%s <span class="mono">%s</span></p></div>'
                        % (locus_map(record, strain, tags),
                           T(('Âncoras:', 'Anchors:')), esc(anchors or '-')))
        blocks.append('<div class="locus">%s%s<p class="links">%s</p></div>'
                      % (header, ''.join(maps),
                         link(ctx, ('%s.gbk' % code, '%s.gbk' % code), '05_synteny', 'loci', '%s.gbk' % code)))

    missing = []
    for code, status in ctx['locus_status'].items():
        if status['result'] == 'OK':
            continue
        strain = strains.get(code)
        # PT-BR: O 'note' já vem bilíngue no arquivo, separado por ' / '.
        # EN-US: The 'note' is already bilingual in the file, separated by ' / '.
        note_text = status['note'] or plain(('Nenhum gene âncora da via da urease.',
                                             'No urease-pathway anchor gene.'))
        parts = note_text.split('/')
        note = ((pt_display(parts[0].strip()), parts[1].strip())
                if len(parts) == 2 else note_text)
        missing.append('<li><span class="dot warn"></span>%s <span class="mono muted">%s</span>'
                       '<span class="muted"> · %s %s · %s</span></li>'
                       % (label_html(strain, short=True) if strain else esc(code), esc(code),
                          T(('candidatos:', 'candidates:')), esc(num(status['candidates_total'] or 0)),
                          T(note)))

    body = legend([('var(--series-1)', ('gene da via da urease', 'urease pathway gene')),
                   ('var(--series-2)', ('outro candidato da triagem', 'other screening candidate')),
                   ('var(--gene-other)', ('gene vizinho', 'neighbouring gene'))])
    body += ''.join(blocks) if blocks else empty_note()
    body += ('<p class="links">%s</p>' % link(ctx, ('Comparação interativa (clinker)',
                                                    'Interactive comparison (clinker)'),
                                              '05_synteny', 'synteny.html'))
    if missing:
        body += ('<div class="missing-box"><h4>%s</h4><ul class="checks">%s</ul></div>'
                 % (T(('Estirpes sem locus de urease', 'Strains with no urease locus')),
                    ''.join(missing)))

    regions = sum(len(v) for v in ctx['loci'].values())
    return card(('Contexto genômico dos genes de urease', 'Genomic context of the urease genes'),
                body,
                ('%s regiões extraídas com 5 genes de flanco em cada lado do bloco de âncoras.'
                 % num(regions),
                 '%s regions extracted with 5 flanking genes on each side of the anchor block.'
                 % num(regions)))


def section_phylogeny(ctx):
    labels = {}
    strains = ctx['strains']
    for row in ctx['tree_labels']:
        tree_label = row.get('tree_label', '')
        origin = 'strain' if row.get('origin', '').startswith('estirpe') else 'reference'
        if origin == 'strain':
            code = tree_label.split('_')[0]
            strain = strains.get(code)
            if strain:
                binomial = '%s %s' % (strain['genus'], strain['species'])
                suffix = '%s · %s' % (strain['designation'], row.get('source_id', ''))
                title = '%s — %s' % (label_full(strain), row.get('original_header', ''))
            else:
                binomial, suffix, title = tree_label, '', tree_label
        else:
            info = parse_reference_id(row.get('source_id', ''))
            organism = info['organism']
            binomial = re.split(r'\s+\(strain', organism)[0].strip() or organism
            words = binomial.split()
            if len(words) > 2:
                binomial = ' '.join(words[:2])
            suffix = '(%s)' % info['uniprot'] if info['uniprot'] else ''
            title = organism
        labels[tree_label] = {'origin': origin, 'binomial': binomial,
                              'suffix': suffix, 'title': title}

    tree_svg = phylogram(ctx['tree'], labels)
    tree_svg += legend([('var(--series-1)', ('UreC das estirpes deste estudo',
                                             'UreC from this study’s strains')),
                        ('var(--muted)', ('UreC de referência (UniProt, revisada)',
                                          'reference UreC (UniProt, reviewed)'))])

    stats = ctx['iqtree']
    facts = [
        (('Taxons', 'Taxa'), num(stats.get('taxa', 0))),
        (('Sítios alinhados (após trimAl)', 'Aligned sites (after trimAl)'), num(stats.get('sites', 0))),
        (('Sítios informativos por parcimônia', 'Parsimony informative sites'),
         num(stats.get('informative', 0))),
        (('Sítios constantes', 'Constant sites'),
         '%s (%s%%)' % (num(stats.get('constant', 0)), num(stats.get('constant_pct', 0.0), 1))),
        (('Modelo (BIC)', 'Model (BIC)'), stats.get('model', '-')),
        (('Log-verossimilhança', 'Log-likelihood'), num(to_float(stats.get('loglik'), 0.0), 2)),
    ]
    facts_html = ['<dl class="facts">']
    for label, value in facts:
        facts_html.append('<div><dt>%s</dt><dd>%s</dd></div>' % (T(label), esc(value)))
    facts_html.append('</dl>')

    status = ctx['phylo_status']
    method = ('<p class="note">%s</p>' % T((
        'Sequências de UreC das estirpes e das referências curadas alinhadas com MAFFT --auto, '
        'aparadas com trimAl -automated1 e reconstruídas com IQ-TREE 2 (-m MFP) e 1000 ultrafast '
        'bootstraps. Os valores de suporte aparecem nos nós internos apenas quando ≥ 70.',
        'UreC sequences from the strains and the curated references aligned with MAFFT --auto, '
        'trimmed with trimAl -automated1 and reconstructed with IQ-TREE 2 (-m MFP) and 1000 ultrafast '
        'bootstraps. Support values are shown on internal nodes only when ≥ 70.')))

    files = '<p class="links">%s %s %s %s</p>' % (
        link(ctx, ('ureC.treefile', 'ureC.treefile'), '06_phylogeny', 'ureC.treefile'),
        link(ctx, ('ureC_trimmed.aln', 'ureC_trimmed.aln'), '06_phylogeny', 'ureC_trimmed.aln'),
        link(ctx, ('ureC.iqtree', 'ureC.iqtree'), '06_phylogeny', 'ureC.iqtree'),
        link(ctx, ('ureC_labels.tsv', 'ureC_labels.tsv'), '06_phylogeny', 'ureC_labels.tsv'))

    rows = []
    for row in ctx['tree_labels']:
        origin = split_bilingual(row.get('origin', ''))
        rows.append([
            cell('<span class="mono">%s</span>' % esc(row.get('tree_label', '')),
                 row.get('tree_label', '')),
            cell(T(origin), plain(origin)),
            cell('<span class="mono">%s</span>' % esc(row.get('source_id', '')),
                 row.get('source_id', '')),
        ])
    decoder = table([('Rótulo na árvore', 'Tree label'), ('Origem', 'Origin'),
                     ('Identificador de origem', 'Source identifier')], rows)

    body = ''.join(facts_html) + method + tree_svg + files
    if status:
        body += '<p class="note">%s</p>' % esc(status)
    return (card(('Árvore de máxima verossimilhança da subunidade alfa (UreC)',
                  'Maximum-likelihood tree of the alpha subunit (UreC)'), body,
                 ('UreC é a subunidade catalítica; é o marcador usual para posicionar a urease '
                  'de um genoma frente a ureases revisadas.',
                  'UreC is the catalytic subunit; it is the usual marker for placing a genome’s '
                  'urease against reviewed ureases.'))
            + card(('Decodificador de rótulos', 'Label decoder'), decoder))


def section_execution(ctx):
    tasks = ctx['tasks']
    if not tasks:
        return card(('Execução do Nextflow', 'Nextflow execution'), empty_note())
    statuses = Counter(t['status'] for t in tasks)
    failed = sum(count for status, count in statuses.items()
                 if status not in ('COMPLETED', 'CACHED'))

    by_process = {}
    for task in tasks:
        by_process.setdefault(task['process'], []).append(task)
    items = []
    for process in PIPELINE_ORDER:
        entries = by_process.get(process)
        if not entries:
            continue
        durations = [t['realtime_s'] or 0.0 for t in entries]
        total = sum(durations)
        items.append({
            'value': total, 'label': process,
            'title': '%s: %s %s · %s %s' % (process, num(len(entries)),
                                                 plain(('tarefas', 'tasks')),
                                                 human_seconds(total),
                                                 plain(('de CPU acumulada', 'accumulated')))})
    items.sort(key=lambda i: i['value'], reverse=True)
    chart = hbar_chart(items, ('tempo acumulado', 'accumulated time'), label_width=200,
                       value_fmt=human_seconds, tick_format=human_seconds)

    tiles = ('<div class="tiles small">'
             '<div class="tile"><span class="tile-value">%s</span>%s</div>'
             '<div class="tile"><span class="tile-value">%s</span>%s</div>'
             '<div class="tile"><span class="tile-value">%s</span>%s</div>'
             '<div class="tile"><span class="tile-value %s">%s</span>%s</div></div>') % (
        esc(num(len(tasks))), T(('tarefas', 'tasks'), tag='span', cls='tile-label'),
        esc(num(statuses.get('COMPLETED', 0))), T(('concluídas', 'completed'), tag='span', cls='tile-label'),
        esc(num(statuses.get('CACHED', 0))), T(('reaproveitadas do cache', 'reused from cache'),
                                               tag='span', cls='tile-label'),
        'good' if failed == 0 else 'bad', esc(num(failed)),
        T(('falhas', 'failures'), tag='span', cls='tile-label'))

    links = '<p class="links">%s %s %s %s</p>' % (
        link(ctx, ('Relatório de execução', 'Execution report'), 'nf_report.html'),
        link(ctx, ('Linha do tempo', 'Timeline'), 'nf_timeline.html'),
        link(ctx, ('Grafo de dependências', 'Dependency graph'), 'nf_dag.html'),
        link(ctx, ('nf_trace.txt', 'nf_trace.txt'), 'nf_trace.txt'))

    return card(('Execução do Nextflow', 'Nextflow execution'), tiles + chart + links,
                ('Tempo acumulado por processo, somando todas as tarefas do processo.',
                 'Accumulated time per process, summing every task of that process.'))


def section_files(ctx):
    results = ctx['results']
    groups = OrderedDict()
    for root, dirs, files in os.walk(results):
        dirs.sort(key=natural_key)
        for name in sorted(files, key=natural_key):
            if name.startswith('.'):
                continue
            full = os.path.join(root, name)
            relative = os.path.relpath(full, results)
            top = relative.split(os.sep)[0] if os.sep in relative else '.'
            groups.setdefault(top, []).append((relative, os.path.getsize(full)))

    out = []
    for top in sorted(groups, key=natural_key):
        entries = groups[top]
        total = sum(size for _, size in entries)
        rows = []
        for relative, size in entries:
            rows.append('<li><a href="%s"><span class="mono">%s</span></a>'
                        '<span class="muted"> %s</span></li>'
                        % (esc(os.path.relpath(os.path.join(results, relative), ctx['out_dir'])),
                           esc(relative), esc(human_bytes(size))))
        name = ('<span class="mono">%s</span>' % esc(top)) if top != '.' else \
            T(('raiz dos resultados', 'results root'))
        out.append('<details class="files"><summary>%s'
                   '<span class="muted"> · %s %s · %s</span></summary><ul>%s</ul></details>'
                   % (name, esc(num(len(entries))), T(('arquivos', 'files')),
                      esc(human_bytes(total)), ''.join(rows)))
    return card(('Todos os artefatos do run', 'Every artifact of the run'), ''.join(out),
                ('Caminhos relativos ao próprio relatório; os links abrem os arquivos originais.',
                 'Paths relative to the report itself; the links open the original files.'))


# ---------------------------------------------------------------------------
# PT-BR: CSS e JS embutidos / EN-US: Inline CSS and JS
# ---------------------------------------------------------------------------

CSS = """
:root{
  color-scheme:light;
  --surface-1:#fcfcfb; --plane:#f9f9f7;
  --text-primary:#0b0b0b; --text-secondary:#52514e; --muted:#898781;
  --grid:#e1e0d9; --baseline:#c3c2b7; --border:rgba(11,11,11,.10);
  --series-1:#2a78d6; --series-2:#eb6834; --series-3:#1baf7a;
  --on-series-1:#ffffff;
  --gene-other:#c3c2b7; --cell-empty:#f2f1ed; --branch:#a9a7a0;
  --flow-box:#ffffff; --flow-input:#f0efec;
  --good:#0ca30c; --warning:#fab219; --critical:#d03b3b;
  --funnel-1:#86b6ef; --funnel-2:#5598e7; --funnel-3:#2a78d6;
  --funnel-4:#1c5cab; --funnel-5:#104281;
  --funnel-ink-1:#0b0b0b; --funnel-ink-2:#0b0b0b; --funnel-ink-3:#ffffff;
  --funnel-ink-4:#ffffff; --funnel-ink-5:#ffffff;
  --shadow:0 1px 2px rgba(11,11,11,.05);
}
@media (prefers-color-scheme:dark){
  :root:where(:not([data-theme="light"])){
    color-scheme:dark;
    --surface-1:#1a1a19; --plane:#0d0d0d;
    --text-primary:#ffffff; --text-secondary:#c3c2b7; --muted:#898781;
    --grid:#2c2c2a; --baseline:#383835; --border:rgba(255,255,255,.10);
    --series-1:#3987e5; --series-2:#d95926; --series-3:#199e70;
    --on-series-1:#0b0b0b;
    --gene-other:#4a4a46; --cell-empty:#232322; --branch:#5c5b56;
    --flow-box:#232322; --flow-input:#2c2c2a;
    --funnel-1:#b7d3f6; --funnel-2:#86b6ef; --funnel-3:#5598e7;
    --funnel-4:#2a78d6; --funnel-5:#184f95;
    --funnel-ink-1:#0b0b0b; --funnel-ink-2:#0b0b0b; --funnel-ink-3:#0b0b0b;
    --funnel-ink-4:#ffffff; --funnel-ink-5:#ffffff;
    --shadow:none;
  }
}
:root[data-theme="dark"]{
  color-scheme:dark;
  --surface-1:#1a1a19; --plane:#0d0d0d;
  --text-primary:#ffffff; --text-secondary:#c3c2b7; --muted:#898781;
  --grid:#2c2c2a; --baseline:#383835; --border:rgba(255,255,255,.10);
  --series-1:#3987e5; --series-2:#d95926; --series-3:#199e70;
  --on-series-1:#0b0b0b;
  --gene-other:#4a4a46; --cell-empty:#232322; --branch:#5c5b56;
  --flow-box:#232322; --flow-input:#2c2c2a;
  --funnel-1:#b7d3f6; --funnel-2:#86b6ef; --funnel-3:#5598e7;
  --funnel-4:#2a78d6; --funnel-5:#184f95;
  --funnel-ink-1:#0b0b0b; --funnel-ink-2:#0b0b0b; --funnel-ink-3:#0b0b0b;
  --funnel-ink-4:#ffffff; --funnel-ink-5:#ffffff;
  --shadow:none;
}
*{box-sizing:border-box}
body{margin:0;background:var(--plane);color:var(--text-primary);
  font-family:system-ui,-apple-system,"Segoe UI",sans-serif;font-size:15px;line-height:1.55}
.mono{font-family:ui-monospace,SFMono-Regular,Menlo,Consolas,monospace;font-size:.92em}
.muted{color:var(--muted)}
.layout{display:grid;grid-template-columns:246px minmax(0,1fr);gap:28px;
  max-width:1420px;margin:0 auto;padding:0 24px 80px}
nav.side{position:sticky;top:0;align-self:start;height:100vh;overflow:auto;padding:22px 0;
  border-right:1px solid var(--border)}
nav.side .brand{font-weight:700;font-size:13px;letter-spacing:.06em;text-transform:uppercase;
  color:var(--muted);margin:0 0 14px}
nav.side ol{list-style:none;margin:0;padding:0;counter-reset:s}
nav.side li{counter-increment:s}
nav.side a{display:block;padding:6px 12px 6px 10px;border-left:2px solid transparent;
  color:var(--text-secondary);text-decoration:none;font-size:13.5px;border-radius:0 6px 6px 0}
nav.side a::before{content:counter(s) ".";color:var(--muted);margin-right:7px;
  font-variant-numeric:tabular-nums}
nav.side a:hover{background:var(--surface-1);color:var(--text-primary)}
nav.side a.active{border-left-color:var(--series-1);color:var(--text-primary);font-weight:600}
.toolbar{display:flex;gap:8px;margin:18px 0 0;padding:0 10px}
.toolbar button{flex:1;padding:6px 8px;font:inherit;font-size:12.5px;cursor:pointer;
  background:var(--surface-1);color:var(--text-secondary);
  border:1px solid var(--border);border-radius:7px}
.toolbar button:hover{color:var(--text-primary)}
.toolbar button[aria-pressed="true"]{background:var(--series-1);border-color:var(--series-1);
  color:var(--on-series-1);font-weight:600}
main{min-width:0;padding-top:26px}
.hero{padding:26px 0 8px}
.eyebrow{margin:0;font-size:12.5px;letter-spacing:.09em;text-transform:uppercase;color:var(--muted)}
h1{margin:6px 0 10px;font-size:31px;line-height:1.2;letter-spacing:-.01em}
.lead{margin:0 0 6px;max-width:70ch;color:var(--text-secondary);font-size:16.5px}
.meta{margin:2px 0 0;font-size:13px;color:var(--muted)}
.tiles{display:grid;grid-template-columns:repeat(auto-fit,minmax(158px,1fr));gap:12px;margin:22px 0 6px}
.tiles.small{grid-template-columns:repeat(auto-fit,minmax(132px,1fr));margin:0 0 18px}
.tile{background:var(--surface-1);border:1px solid var(--border);border-radius:11px;
  padding:13px 15px;box-shadow:var(--shadow)}
.tile-value{display:block;font-size:26px;font-weight:700;line-height:1.15;letter-spacing:-.01em}
.tile-value.good{color:var(--good)} .tile-value.bad{color:var(--critical)}
.tile-label{display:block;font-size:12.5px;color:var(--text-secondary);margin-top:2px}
h2.sec{margin:44px 0 4px;font-size:13px;letter-spacing:.09em;text-transform:uppercase;
  color:var(--muted);padding-top:14px;border-top:1px solid var(--border)}
.card{background:var(--surface-1);border:1px solid var(--border);border-radius:13px;
  padding:20px 22px;margin:16px 0;box-shadow:var(--shadow)}
.card h3{margin:0 0 4px;font-size:18px;letter-spacing:-.005em}
.card h4{margin:16px 0 6px;font-size:14.5px}
.card h5{margin:0 0 8px;font-size:13px;color:var(--text-secondary)}
.sub{margin:0 0 14px;color:var(--text-secondary);font-size:14px;max-width:82ch}
.note{margin:12px 0 0;font-size:13px;color:var(--text-secondary);max-width:88ch}
.note code{background:var(--plane);padding:1px 5px;border-radius:4px;font-size:.92em}
.note.warn{border-left:3px solid var(--warning);padding-left:12px}
.empty{color:var(--muted);font-style:italic;margin:8px 0}
svg.chart{display:block;max-width:100%;height:auto;margin:6px 0 2px;overflow:visible}
svg.chart a{cursor:pointer}
svg.chart a:hover rect{filter:brightness(.97)}
.legend{display:flex;flex-wrap:wrap;gap:6px 20px;list-style:none;margin:10px 0 0;padding:0;
  font-size:12.5px;color:var(--text-secondary)}
.legend li{display:flex;align-items:center;gap:7px}
.swatch{width:12px;height:12px;border-radius:3px;border:1px solid var(--border);flex:none}
.table-wrap{overflow-x:auto;margin:12px 0 0;border:1px solid var(--border);border-radius:10px}
table{border-collapse:collapse;width:100%;font-size:13.5px}
table.inner{border:0;margin:0}
caption{text-align:left;padding:8px 12px;color:var(--muted);font-size:12.5px}
th,td{padding:7px 12px;text-align:left;border-bottom:1px solid var(--border);vertical-align:top}
thead th{position:sticky;top:0;background:var(--surface-1);color:var(--text-secondary);
  font-weight:600;font-size:12.5px;white-space:nowrap;z-index:1}
table[data-sortable] thead th{cursor:pointer;user-select:none}
table[data-sortable] thead th::after{content:"";opacity:.35;margin-left:5px}
table[data-sortable] thead th[aria-sort="ascending"]::after{content:"\\2191";opacity:1}
table[data-sortable] thead th[aria-sort="descending"]::after{content:"\\2193";opacity:1}
tbody tr:last-child td{border-bottom:0}
td.right,th.right{text-align:right;font-variant-numeric:tabular-nums}
td.center{text-align:center}
td.nowrap{white-space:nowrap}
tr.expandable{cursor:pointer}
tr.expandable:hover td{background:var(--plane)}
tr.evidence>td{background:var(--plane);padding:14px 16px}
.gene{font-weight:600;color:var(--text-primary)}
.badge{display:inline-block;min-width:19px;padding:1px 0;margin-right:3px;border-radius:4px;
  font-size:11px;font-weight:700;text-align:center;font-family:ui-monospace,monospace}
.badge.on{background:var(--series-1);color:var(--on-series-1)}
.badge.off{background:var(--cell-empty);color:var(--muted)}
.bit{display:inline-block;min-width:20px;font-family:ui-monospace,monospace;font-weight:600}
.bit.on{color:var(--series-1)} .bit.off{color:var(--muted)}
.conf{font-size:12.5px;font-weight:600;white-space:nowrap}
.conf.high{color:var(--series-1)} .conf.medium{color:var(--text-secondary);font-weight:400}
.controls{display:flex;flex-wrap:wrap;gap:12px;align-items:center;margin:12px 0 0}
.controls input[type=search]{flex:1;min-width:220px;padding:7px 11px;font:inherit;font-size:13.5px;
  background:var(--plane);color:var(--text-primary);
  border:1px solid var(--border);border-radius:8px}
.chk{display:flex;align-items:center;gap:6px;font-size:13.5px;color:var(--text-secondary);cursor:pointer}
.checks{list-style:none;margin:6px 0 0;padding:0;font-size:13.5px}
.checks li{padding:4px 0;display:flex;align-items:center;gap:8px;flex-wrap:wrap}
.dot{width:8px;height:8px;border-radius:50%;flex:none}
.dot.ok{background:var(--good)} .dot.warn{background:var(--warning)}
.ok-line{display:flex;align-items:center;gap:8px;margin:10px 0 0;font-size:13.5px;color:var(--good)}
.links{margin:12px 0 0;font-size:13px;display:flex;flex-wrap:wrap;gap:8px 16px}
a.artifact{color:var(--series-1);text-decoration:none;border-bottom:1px solid var(--border)}
a.artifact:hover{border-bottom-color:var(--series-1)}
.missing{color:var(--muted);text-decoration:line-through}
.criteria{display:grid;grid-template-columns:repeat(auto-fit,minmax(238px,1fr));gap:16px;margin-top:8px}
.criteria h4{margin:0 0 4px;font-size:13.5px;color:var(--series-1)}
.criteria p{margin:0;font-size:13px;color:var(--text-secondary)}
.reading{margin-top:18px;padding:14px 16px;background:var(--plane);border-radius:10px}
.reading h4{margin:0 0 8px;font-size:14px}
.reading ul{margin:0;padding-left:18px;font-size:13.5px}
.reading p{margin:10px 0 0;font-size:13.5px;color:var(--text-secondary)}
.facts{display:grid;grid-template-columns:repeat(auto-fit,minmax(178px,1fr));gap:12px;margin:6px 0 4px}
.facts div{background:var(--plane);border-radius:9px;padding:10px 12px}
.facts dt{font-size:12px;color:var(--muted)}
.facts dd{margin:2px 0 0;font-size:16px;font-weight:600;font-variant-numeric:tabular-nums}
.locus{margin:22px 0 0;padding-top:14px;border-top:1px solid var(--border)}
.locus-head{margin:0 0 2px;font-size:15px}
.locus-block{margin:10px 0 0}
.missing-box{margin-top:22px;padding:14px 16px;background:var(--plane);border-radius:10px}
.missing-box h4{margin:0 0 4px;font-size:14px}
details.files{margin:10px 0 0;font-size:13px}
details.files summary{cursor:pointer;color:var(--text-secondary);padding:5px 0}
details.files ul{list-style:none;margin:6px 0 0;padding:0 0 0 14px;
  columns:2;column-gap:26px}
details.files li{padding:2px 0;break-inside:avoid}
details.files a{color:var(--series-1);text-decoration:none}
footer{max-width:1420px;margin:0 auto;padding:26px 24px 48px;color:var(--muted);font-size:12.5px;
  border-top:1px solid var(--border)}
@media (max-width:960px){
  .layout{grid-template-columns:1fr;gap:0}
  nav.side{position:static;height:auto;border-right:0;border-bottom:1px solid var(--border)}
  nav.side ol{display:flex;flex-wrap:wrap;gap:2px}
  nav.side a{border-left:0;padding:5px 9px;border-radius:6px}
  nav.side a.active{background:var(--surface-1)}
  details.files ul{columns:1}
}
@media print{
  nav.side,.toolbar{display:none}
  .layout{grid-template-columns:1fr;max-width:none}
  .card{break-inside:avoid;box-shadow:none}
}
"""

JS = """
(function(){
  var root=document.documentElement;

  function setLang(lang){
    var nodes=document.querySelectorAll('.i18n');
    for(var i=0;i<nodes.length;i++){
      var v=nodes[i].getAttribute('data-'+lang);
      if(v!==null) nodes[i].textContent=v;
    }
    var ph=document.querySelectorAll('[data-ph-'+lang+']');
    for(var p=0;p<ph.length;p++){
      ph[p].setAttribute('placeholder',ph[p].getAttribute('data-ph-'+lang));
      var al=ph[p].getAttribute('data-al-'+lang);
      if(al) ph[p].setAttribute('aria-label',al);
    }
    root.setAttribute('lang',lang==='pt'?'pt-BR':'en');
    var btns=document.querySelectorAll('[data-lang]');
    for(var j=0;j<btns.length;j++){
      btns[j].setAttribute('aria-pressed',btns[j].getAttribute('data-lang')===lang?'true':'false');
    }
    try{localStorage.setItem('bafes-lang',lang);}catch(e){}
  }

  function setTheme(theme){
    if(theme) root.setAttribute('data-theme',theme); else root.removeAttribute('data-theme');
    var btns=document.querySelectorAll('[data-theme-set]');
    for(var i=0;i<btns.length;i++){
      btns[i].setAttribute('aria-pressed',btns[i].getAttribute('data-theme-set')===(theme||'auto')?'true':'false');
    }
    try{localStorage.setItem('bafes-theme',theme||'auto');}catch(e){}
  }

  document.addEventListener('click',function(ev){
    var el=ev.target.closest ? ev.target.closest('[data-lang],[data-theme-set]') : null;
    if(!el) return;
    if(el.hasAttribute('data-lang')) setLang(el.getAttribute('data-lang'));
    else { var t=el.getAttribute('data-theme-set'); setTheme(t==='auto'?null:t); }
  });

  // PT-BR: expansão das evidências / EN-US: evidence drill-down
  document.addEventListener('click',function(ev){
    var row=ev.target.closest ? ev.target.closest('tr.expandable') : null;
    if(!row) return;
    var target=document.getElementById(row.getAttribute('data-target'));
    if(target) target.hidden=!target.hidden;
  });

  // PT-BR: ordenação das tabelas / EN-US: table sorting
  document.addEventListener('click',function(ev){
    var th=ev.target.closest ? ev.target.closest('table[data-sortable] thead th') : null;
    if(!th) return;
    var table=th.closest('table'), tbody=table.tBodies[0];
    var index=Array.prototype.indexOf.call(th.parentNode.children,th);
    var asc=th.getAttribute('aria-sort')!=='ascending';
    var heads=th.parentNode.children;
    for(var h=0;h<heads.length;h++) heads[h].removeAttribute('aria-sort');
    th.setAttribute('aria-sort',asc?'ascending':'descending');
    var groups=[],rows=tbody.rows,i;
    for(i=0;i<rows.length;i++){
      if(rows[i].classList.contains('evidence')){
        if(groups.length) groups[groups.length-1].extra.push(rows[i]);
      } else groups.push({main:rows[i],extra:[]});
    }
    function key(tr){
      var cell=tr.children[index];
      if(!cell) return '';
      var raw=cell.hasAttribute('data-sort')?cell.getAttribute('data-sort'):cell.textContent;
      var n=parseFloat(String(raw).replace(/[\\u202f\\s,%]/g,''));
      return isNaN(n)?String(raw).toLowerCase():n;
    }
    groups.sort(function(a,b){
      var ka=key(a.main),kb=key(b.main);
      if(typeof ka==='number'&&typeof kb==='number') return asc?ka-kb:kb-ka;
      ka=String(ka);kb=String(kb);
      return asc?ka.localeCompare(kb):kb.localeCompare(ka);
    });
    var frag=document.createDocumentFragment();
    for(i=0;i<groups.length;i++){
      frag.appendChild(groups[i].main);
      for(var j=0;j<groups[i].extra.length;j++) frag.appendChild(groups[i].extra[j]);
    }
    tbody.appendChild(frag);
  });

  // PT-BR: filtros da tabela de candidatos / EN-US: candidate table filters
  // PT-BR: A pagina esta acentuada, entao a busca ignora acentos: quem digita
  //        "media" ou "regiao" encontra "Media" e "regiao".
  // EN-US: The page is accented, so the search is accent-insensitive: typing
  //        "media" or "regiao" still finds "Media" and "regiao".
  function fold(s){
    s=String(s).toLowerCase();
    return s.normalize ? s.normalize('NFD').replace(/[\u0300-\u036f]/g,'') : s;
  }
  function applyFilter(tableId){
    var table=document.getElementById(tableId); if(!table) return;
    var input=document.querySelector('.filter[data-table="'+tableId+'"]');
    var only=document.querySelector('.only-high[data-table="'+tableId+'"]');
    var term=fold((input&&input.value||'').trim());
    var high=only&&only.checked;
    var rows=table.tBodies[0].rows,shown=null;
    for(var i=0;i<rows.length;i++){
      var row=rows[i];
      if(row.classList.contains('evidence')){ if(shown===false) row.hidden=true; continue; }
      var conf=row.querySelector('.conf'); var isHigh=conf&&conf.classList.contains('high');
      var text=fold(row.textContent);
      var ok=(!term||text.indexOf(term)>-1)&&(!high||isHigh);
      row.hidden=!ok; shown=ok;
    }
  }
  document.addEventListener('input',function(ev){
    if(ev.target.classList.contains('filter')) applyFilter(ev.target.getAttribute('data-table'));
  });
  document.addEventListener('change',function(ev){
    if(ev.target.classList.contains('only-high')) applyFilter(ev.target.getAttribute('data-table'));
  });

  // PT-BR: navegação ativa / EN-US: scrollspy
  var links={},sections=[];
  var anchors=document.querySelectorAll('nav.side a[href^="#"]');
  for(var a=0;a<anchors.length;a++){
    var id=anchors[a].getAttribute('href').slice(1);
    links[id]=anchors[a];
    var sec=document.getElementById(id); if(sec) sections.push(sec);
  }
  if('IntersectionObserver' in window && sections.length){
    var obs=new IntersectionObserver(function(entries){
      entries.forEach(function(entry){
        var link=links[entry.target.id]; if(!link) return;
        if(entry.isIntersecting){
          for(var k in links) links[k].classList.remove('active');
          link.classList.add('active');
        }
      });
    },{rootMargin:'-10% 0px -75% 0px'});
    for(var s=0;s<sections.length;s++) obs.observe(sections[s]);
  }

  try{
    var savedLang=localStorage.getItem('bafes-lang');
    if(savedLang) setLang(savedLang); else setLang(root.getAttribute('data-default-lang')||'pt');
    var savedTheme=localStorage.getItem('bafes-theme');
    setTheme(savedTheme&&savedTheme!=='auto'?savedTheme:null);
  }catch(e){}
})();
"""


def build_html(ctx):
    sections = OrderedDict([
        ('fluxo', section_flow),
        ('entrada', section_input),
        ('qualidade', section_quality),
        ('anotação', section_annotation),
        ('triagem', section_screening),
        ('matriz', section_matrix),
        ('candidatos', section_candidates),
        ('pfam', section_pfam),
        ('sintenia', section_synteny),
        ('filogenia', section_phylogeny),
        ('execução', section_execution),
        ('arquivos', section_files),
    ])

    nav = ['<nav class="side"><p class="brand">BAFES‑Urease</p><ol>']
    for anchor, title in SECTIONS:
        nav.append('<li><a href="#%s">%s</a></li>' % (esc(anchor), T(title)))
    nav.append('</ol>')
    nav.append('<div class="toolbar" role="group">'
               '<button type="button" data-lang="pt" aria-pressed="false">PT</button>'
               '<button type="button" data-lang="en" aria-pressed="false">EN</button>'
               '</div>')
    nav.append('<div class="toolbar" role="group">'
               '<button type="button" data-theme-set="auto" aria-pressed="true">%s</button>'
               '<button type="button" data-theme-set="light" aria-pressed="false">%s</button>'
               '<button type="button" data-theme-set="dark" aria-pressed="false">%s</button>'
               '</div>' % (T(('Auto', 'Auto')), T(('Claro', 'Light')), T(('Escuro', 'Dark'))))
    nav.append('</nav>')

    main = [section_hero(ctx)]
    for anchor, title in SECTIONS:
        main.append('<h2 class="sec" id="%s">%s</h2>' % (esc(anchor), T(title)))
        try:
            main.append(sections[anchor](ctx))
        except Exception as error:  # PT-BR: uma seção quebrada não derruba o relatório
            main.append(card(title, '<p class="empty">%s %s</p>'
                             % (T(('Falha ao montar esta seção:', 'Failed to build this section:')),
                                esc(error))))

    doc = ['<!doctype html><html lang="%s" data-default-lang="%s"><head>' % (
        'pt-BR' if LANG == 'pt' else 'en', LANG)]
    doc.append('<meta charset="utf-8">')
    doc.append('<meta name="viewport" content="width=device-width,initial-scale=1">')
    doc.append('<title>%s</title>' % esc(plain(ctx['title'])))
    doc.append('<style>%s</style>' % CSS)
    doc.append('</head><body>')
    doc.append('<div class="layout">%s<main>%s</main></div>' % (''.join(nav), ''.join(main)))
    doc.append('<footer>%s</footer>' % T((
        'Gerado por bin/generate_report.py a partir da árvore de resultados do pipeline '
        'bafes-urease. Todos os gráficos foram desenhados diretamente em SVG; nenhum dado foi '
        'reprocessado ou estimado.',
        'Generated by bin/generate_report.py from the bafes-urease results tree. Every chart is '
        'drawn directly in SVG; no data was reprocessed or estimated.')))
    doc.append('<script>%s</script>' % JS)
    doc.append('</body></html>')
    return ''.join(doc)


# ---------------------------------------------------------------------------

def main():
    global LANG
    parser = argparse.ArgumentParser(
        description='Relatório Web dos resultados do bafes-urease / '
                    'Web report of the bafes-urease results (PT-BR / EN-US)')
    parser.add_argument('--results', required=True,
                        help='Diretório de resultados do pipeline / Pipeline results directory')
    parser.add_argument('--out', default=None,
                        help='Arquivo HTML de saída / Output HTML file '
                             '(padrao / default: <results>/../report.html)')
    parser.add_argument('--accessions', default=None,
                        help='data/accessions.tsv (opcional, só para o acesso GenBank / '
                             'optional, only for the GenBank accession)')
    parser.add_argument('--lang-default', default='pt', choices=['pt', 'en'],
                        help='Idioma inicial da página / Initial page language')
    parser.add_argument('--title', default=None, help='Título do relatório / Report title')
    args = parser.parse_args()

    results = os.path.abspath(args.results)
    if not os.path.isdir(results):
        print('[ERRO/ERROR] Diretório de resultados inexistente / Results directory not found: %s'
              % results)
        return 1

    out_path = os.path.abspath(args.out) if args.out else \
        os.path.join(os.path.dirname(results), 'report.html')
    out_dir = os.path.dirname(out_path) or '.'
    if not os.path.isdir(out_dir):
        print('[ERRO/ERROR] Diretório de saída inexistente / Output directory not found: %s' % out_dir)
        return 1

    LANG = args.lang_default

    accessions = args.accessions
    if accessions is None:
        guess = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
                             'data', 'accessions.tsv')
        accessions = guess if os.path.isfile(guess) else None

    warnings = []

    def warn(message):
        warnings.append(message)
        print('[AVISO/WARNING] %s' % message)

    ctx = build_context(results, accessions, out_dir, warn)
    ctx['title'] = (args.title if args.title else
                    ('Urease em Bacillales: mineração genômica de dez estirpes',
                     'Urease in Bacillales: genomic mining of ten strains'))

    document = build_html(ctx)
    with open(out_path, 'w', encoding='utf-8') as handle:
        handle.write(document)

    size = os.path.getsize(out_path)
    print('[OK] Relatório gerado / Report generated: %s (%s, %d %s)'
          % (out_path, human_bytes(size).replace(NBSP_THIN, ' '), len(SECTIONS),
             'seções / sections'))
    print('[OK] %s estirpes, %s candidatos, %s táxons na árvore / strains, candidates, tree taxa'
          % (len(ctx['strains']), len(ctx['candidates']), ctx['iqtree'].get('taxa', 0)))
    if warnings:
        print('[AVISO/WARNING] %d aviso(s) durante a geração / warning(s) during generation'
              % len(warnings))
    return 0


if __name__ == '__main__':
    sys.exit(main())

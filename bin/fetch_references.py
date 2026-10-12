#!/usr/bin/env python3
"""
fetch_references.py
PT-BR: Baixa do UniProt as sequencias de referencia definidas por uma consulta REST, sempre
       restrita a Firmicutes (taxonomy_id:1239), em dois grupos:
         1. urease (protein_name:urease) — so revisadas (Swiss-Prot);
         2. transporte de ureia e via alternativa — genes urtA-E, uca (urea carboxilase) e
            atzF (alofanato hidrolase), por nome exato de gene, revisadas ou nao.
EN-US: Downloads the reference sequences from UniProt as defined by a REST query, always
       restricted to Firmicutes (taxonomy_id:1239), in two groups:
         1. urease (protein_name:urease) — reviewed only (Swiss-Prot);
         2. urea transport and the alternative pathway — genes urtA-E, uca (urea carboxylase)
            and atzF (allophanate hydrolase), by exact gene name, reviewed or not.

PT-BR: Por que o grupo 2 aceita TrEMBL: o Swiss-Prot nao tem nenhuma UrtA-E em organismo
       algum, e as unicas urea carboxilases / alofanato hidrolases revisadas sao de fungos
       (DUR1,2, lamA) e de Pseudomonas (AtzF). Sem TrEMBL a Camada B nunca cobriria esses
       genes. O custo e que a anotacao dessas entradas e automatica, como a do Bakta — o
       acerto de BLASTp para elas confirma similaridade, nao funcao verificada.
EN-US: Why group 2 accepts TrEMBL: Swiss-Prot holds no UrtA-E in any organism, and the only
       reviewed urea carboxylases / allophanate hydrolases are fungal (DUR1,2, lamA) or from
       Pseudomonas (AtzF). Without TrEMBL, Layer B would never cover these genes. The cost is
       that these entries are annotated automatically, like Bakta's — a BLASTp hit on them
       confirms similarity, not verified function.

PT-BR: Por que nome de gene e nao numero EC: ec:3.5.1.54 em Firmicutes tambem traz as
       "alofanato hidrolases" de duas subunidades (KipI/KipA), que em Bacillus sao a
       5-oxoprolinase PxpB/PxpC; e ec:6.3.4.6 traz fragmentos "urea amidolyase related
       protein". gene_exact:atzF e gene_exact:uca ficam so com as enzimas canonicas.
EN-US: Why gene name and not EC number: ec:3.5.1.54 in Firmicutes also brings the two-subunit
       "allophanate hydrolases" (KipI/KipA), which in Bacillus are the 5-oxoprolinase
       PxpB/PxpC; and ec:6.3.4.6 brings "urea amidolyase related protein" fragments.
       gene_exact:atzF and gene_exact:uca keep only the canonical enzymes.

PT-BR: Toda sequencia vem do endpoint /stream do UniProt. Nada e sintetizado localmente:
       falha de rede, resposta truncada ou conjunto vazio abortam sem escrever o FASTA.
EN-US: Every sequence comes from UniProt's /stream endpoint. Nothing is synthesised locally:
       a network failure, a truncated response or an empty result set aborts without
       writing the FASTA.
"""

import sys
import os
import gzip
import time
import hashlib
import argparse
import urllib.parse
import urllib.request
import urllib.error

# PT-BR: Genes do grupo 2 (ver docstring), buscados por nome exato de gene.
# EN-US: Group 2 genes (see docstring), searched by exact gene name.
TRANSPORT_AND_ALT_GENES = ("urtA", "urtB", "urtC", "urtD", "urtE", "uca", "atzF")

# PT-BR: Consulta oficial das referencias deste estudo.
# EN-US: Official reference query for this study.
UNIPROT_QUERY = (
    "(taxonomy_id:1239) AND ("
    "((protein_name:urease) AND (reviewed:true)) OR ("
    + " OR ".join(f"gene_exact:{g}" for g in TRANSPORT_AND_ALT_GENES)
    + "))"
)
USER_AGENT = "bafes-urease-mining/0.4 (https://github.com/cbafes-unb/bafes-urease)"

# PT-BR: Piso de sanidade — a consulta devolve da ordem de mil registros; um numero muito
#        baixo indica resposta truncada ou consulta alterada, e nao um resultado real.
# EN-US: Sanity floor — the query returns on the order of a thousand records; a very low
#        number means a truncated response or a changed query, not a real result.
MIN_SEQUENCES = 100

# PT-BR: Grupo externo da arvore de UreC: as UreC revisadas de Klebsiella aerogenes e de
#        Proteus mirabilis — ureases-modelo com estrutura resolvida, de Gammaproteobacteria,
#        claramente fora de Firmicutes. Vao para um FASTA proprio e NAO entram nas
#        referencias do BLASTp, que continuam restritas a Firmicutes.
# EN-US: Outgroup for the UreC tree: the reviewed UreC of Klebsiella aerogenes and Proteus
#        mirabilis — model ureases with solved structures, from Gammaproteobacteria, clearly
#        outside Firmicutes. They go into a FASTA of their own and do NOT enter the BLASTp
#        references, which stay restricted to Firmicutes.
OUTGROUP_ACCESSIONS = ("P18314", "P17086")
OUTGROUP_QUERY = (
    "(" + " OR ".join(f"accession:{a}" for a in OUTGROUP_ACCESSIONS) + ") AND (reviewed:true)"
)


def build_url(query):
    """
    PT-BR: Monta a URL do endpoint /stream para uma consulta arbitraria.
    EN-US: Builds the /stream endpoint URL for an arbitrary query.
    """
    params = urllib.parse.urlencode({
        'compressed': 'true',
        'format': 'fasta',
        'query': query,
    })
    return f"https://rest.uniprot.org/uniprotkb/stream?{params}"


UNIPROT_STREAM_URL = build_url(UNIPROT_QUERY)


def fetch_stream(url, attempts=4, timeout=300):
    """
    PT-BR: Baixa o FASTA completo da consulta, com retry, backoff e gunzip.
    EN-US: Downloads the full FASTA for the query, with retry, backoff and gunzip.

    PT-BR: Devolve (texto, cabecalhos_http). O gzip so descomprime por inteiro se o corpo
           chegou completo — resposta cortada no meio levanta erro em vez de virar FASTA parcial.
    EN-US: Returns (text, http_headers). gzip only inflates if the body arrived complete —
           a cut-off response raises instead of becoming a partial FASTA.
    """
    backoff = 2.0
    last_error = None

    for attempt in range(1, attempts + 1):
        try:
            req = urllib.request.Request(url, headers={'User-Agent': USER_AGENT})
            with urllib.request.urlopen(req, timeout=timeout) as resp:
                payload = resp.read()
                headers = dict(resp.headers.items())
            if payload[:2] == b'\x1f\x8b':
                payload = gzip.decompress(payload)
            return payload.decode('utf-8'), headers
        except Exception as exc:  # HTTPError, URLError, timeout, BadGzipFile, EOFError
            last_error = exc
            if attempt == attempts:
                break
            print(f"    [RETRY {attempt}/{attempts - 1}] {exc} — aguardando / waiting {backoff:.0f}s")
            time.sleep(backoff)
            backoff *= 2

    raise RuntimeError(f"Falha ao baixar do UniProt / Failed to download from UniProt: {last_error}")


def parse_uniprot_header(header):
    """
    PT-BR: Extrai (accession, gene, organismo, status) de um cabecalho FASTA do UniProt,
           no formato db|accession|entry Nome OS=... OX=... GN=... PE=... SV=...
    EN-US: Extracts (accession, gene, organism, status) from a UniProt FASTA header in the
           db|accession|entry Name OS=... OX=... GN=... PE=... SV=... format.
    """
    ident, _, description = header.partition(' ')
    parts = ident.split('|')
    if len(parts) != 3:
        raise ValueError(f"Cabecalho UniProt inesperado / Unexpected UniProt header: {header!r}")

    db, accession, entry_name = parts
    status = 'reviewed' if db == 'sp' else 'unreviewed'

    fields = {}
    key = None
    for token in description.split():
        if len(token) > 3 and token[2] == '=' and token[:2].isalpha() and token[:2].isupper():
            key = token[:2]
            fields[key] = token[3:]
        elif key is not None:
            fields[key] += '_' + token

    gene = fields.get('GN') or entry_name
    organism = fields.get('OS') or 'unknown_organism'
    return accession, gene, organism, status


def wrap(sequence, width=60):
    return '\n'.join(sequence[i:i + width] for i in range(0, len(sequence), width))


def parse_fasta(text):
    """
    PT-BR: Percorre o FASTA multi-registro devolvendo pares (cabecalho, sequencia).
    EN-US: Walks the multi-record FASTA yielding (header, sequence) pairs.
    """
    header = None
    chunks = []
    for raw in text.splitlines():
        line = raw.strip()
        if not line:
            continue
        if line.startswith('>'):
            if header is not None:
                yield header, ''.join(chunks)
            header = line[1:]
            chunks = []
        else:
            if header is None:
                raise ValueError("FASTA invalido: sequencia antes do cabecalho / "
                                 "Invalid FASTA: sequence before header")
            chunks.append(line)
    if header is not None:
        yield header, ''.join(chunks)


def write_provenance(path, url, headers, n_seqs, fasta_path, digest, extra=()):
    """
    PT-BR: Registra a procedencia do download. A consulta e viva: o conjunto muda a cada
           release do UniProt, entao o que garante reprodutibilidade e este registro
           (release, data, contagem e hash) junto do FASTA versionado.
    EN-US: Records the download provenance. The query is live: the set changes with every
           UniProt release, so reproducibility rests on this record (release, date, count
           and hash) alongside the versioned FASTA.
    """
    lines = [
        "# PT-BR: Procedencia das referencias de urease / EN-US: Urease reference provenance",
        f"query_url\t{url}",
        f"uniprot_release\t{headers.get('X-UniProt-Release', 'unknown')}",
        f"uniprot_release_date\t{headers.get('X-UniProt-Release-Date', 'unknown')}",
        f"downloaded_at\t{headers.get('Date', 'unknown')}",
        f"fasta\t{os.path.basename(fasta_path)}",
        f"sequences\t{n_seqs}",
        *(f"{key}\t{value}" for key, value in extra),
        f"sha256\t{digest}",
    ]
    with open(path, 'w', encoding='utf-8') as handle:
        handle.write('\n'.join(lines) + '\n')


def download_records(url):
    """
    PT-BR: Baixa e analisa a consulta. Devolve (registros, cabecalhos_http), sendo cada
           registro um dict com accession, gene, organism, status e sequence. Levanta
           RuntimeError (rede) ou ValueError (resposta invalida) — nunca devolve parcial.
    EN-US: Downloads and parses the query. Returns (records, http_headers), each record a
           dict with accession, gene, organism, status and sequence. Raises RuntimeError
           (network) or ValueError (invalid response) — never returns a partial set.
    """
    text, headers = fetch_stream(url)
    records = []
    seen = set()
    for header, sequence in parse_fasta(text):
        accession, gene, organism, status = parse_uniprot_header(header)
        if not sequence:
            raise ValueError(f"Sequencia vazia / Empty sequence: {accession}")
        if accession in seen:
            continue
        seen.add(accession)
        records.append({
            'accession': accession,
            'gene': gene,
            'organism': '_'.join(organism.split()),
            'status': status,
            'sequence': sequence,
        })
    return records, headers


def species_key(organism):
    """
    PT-BR: Reduz o organismo ao binomio, descartando cepa e subespecie:
           'Staphylococcus_aureus_(strain_JH9)' -> 'Staphylococcus_aureus'. Organismos sem
           especie definida ('Bacillus_sp._(strain_TB-90)') ficam com o nome inteiro, para
           que isolados distintos nao colapsem num so.
    EN-US: Reduces the organism to its binomial, dropping strain and subspecies:
           'Staphylococcus_aureus_(strain_JH9)' -> 'Staphylococcus_aureus'. Organisms with
           no species assigned ('Bacillus_sp._(strain_TB-90)') keep the full name, so that
           distinct isolates do not collapse into one.
    """
    words = organism.split('_(')[0].split('_')
    if len(words) < 2 or words[1] == 'sp.':
        return organism
    return '_'.join(words[:2])


def deduplicate(records):
    """
    PT-BR: Mantem uma sequencia por (gene, especie). Sem isso, as dezenas de cepas de
           S. aureus — quase identicas entre si — dominavam o conjunto (94 de 218 na
           release 2026_03) e inflavam esse ramo da arvore sem acrescentar diversidade.
           Escolha deterministica: revisada primeiro, depois a menor acessao.
    EN-US: Keeps one sequence per (gene, species). Without it, the dozens of S. aureus
           strains — near-identical to one another — dominated the set (94 of 218 on
           release 2026_03) and inflated that branch of the tree without adding diversity.
           Deterministic choice: reviewed first, then the lowest accession.
    """
    kept = {}
    for rec in sorted(records, key=lambda r: (r['status'] != 'reviewed', r['accession'])):
        key = (rec['gene'].lower(), species_key(rec['organism']))
        kept.setdefault(key, rec)
    return sorted(kept.values(), key=lambda r: (r['gene'].lower(), r['organism'], r['accession']))


def to_fasta(records):
    return ''.join(
        f">ref|{r['gene']}|UniProt:{r['accession']}|{r['organism']}|{r['status']}|"
        f"{len(r['sequence'])}aa\n{wrap(r['sequence'])}\n"
        for r in records
    )


def write_set(records, out, provenance, url, headers, extra=()):
    """
    PT-BR: Grava o FASTA e a procedencia correspondente.
    EN-US: Writes the FASTA and its matching provenance.
    """
    body = to_fasta(records)
    digest = hashlib.sha256(body.encode('utf-8')).hexdigest()
    os.makedirs(os.path.dirname(os.path.abspath(out)), exist_ok=True)
    with open(out, 'w', encoding='utf-8') as handle:
        handle.write(body)
    write_provenance(provenance, url, headers, len(records), out, digest, extra)
    return digest


def main():
    parser = argparse.ArgumentParser(
        description="Fetch urease, urea transport and alternative-pathway reference sequences "
                    "from UniProt (Firmicutes), plus the UreC outgroup")
    parser.add_argument("--url", default=None,
                        help="URL completa do endpoint /stream / Full /stream endpoint URL")
    parser.add_argument("--query", default=None,
                        help=f"Consulta UniProt / UniProt query (default: {UNIPROT_QUERY!r})")
    parser.add_argument("--out", default="data/urease_references.fasta",
                        help="FASTA de saida / Output FASTA")
    parser.add_argument("--provenance", default=None,
                        help="Arquivo de procedencia / Provenance file "
                             "(default: <out>.provenance.txt)")
    parser.add_argument("--min-sequences", type=int, default=MIN_SEQUENCES,
                        help="Minimo aceitavel de sequencias / Minimum acceptable sequence count")
    parser.add_argument("--keep-redundant", action="store_true",
                        help="Nao deduplicar por (gene, especie) / Do not deduplicate by (gene, species)")
    parser.add_argument("--outgroup-out", default="data/urec_outgroup.fasta",
                        help="FASTA do grupo externo da arvore de UreC / UreC tree outgroup FASTA")
    parser.add_argument("--skip-outgroup", action="store_true",
                        help="Nao baixar o grupo externo / Do not download the outgroup")
    args = parser.parse_args()

    if args.url and args.query:
        print("[ERRO/ERROR] Use --url ou --query, nao os dois / Use --url or --query, not both.")
        return 2

    if args.url:
        url = args.url
    elif args.query:
        url = build_url(args.query)
    else:
        url = UNIPROT_STREAM_URL

    print(">>> Baixando referencias do UniProt / Downloading references from UniProt...")
    print(f"    {url}")

    # PT-BR: Tudo e baixado e validado antes de qualquer escrita: uma falha no grupo externo
    #        nao pode deixar referencias novas ao lado de um grupo externo de outra release.
    # EN-US: Everything is downloaded and validated before anything is written: an outgroup
    #        failure must not leave new references next to an outgroup from another release.
    try:
        records, headers = download_records(url)
        if not args.skip_outgroup:
            outgroup_url = build_url(OUTGROUP_QUERY)
            print(f"    {outgroup_url}")
            outgroup, outgroup_headers = download_records(outgroup_url)
    except RuntimeError as exc:
        print(f"[ERRO/ERROR] {exc}")
        print("             Nenhuma sequencia sera inventada / No sequence will be invented.")
        return 1
    except ValueError as exc:
        print(f"[ERRO/ERROR] Resposta do UniProt invalida / Invalid UniProt response: {exc}")
        print("             Nenhum arquivo foi escrito / No file was written.")
        return 1

    # PT-BR: O piso vale para a resposta bruta, antes da deduplicacao: e ela que denuncia
    #        uma resposta truncada.
    # EN-US: The floor applies to the raw response, before deduplication: that is what
    #        gives away a truncated response.
    if len(records) < args.min_sequences:
        print(f"[ERRO/ERROR] Apenas {len(records)} sequencia(s) obtida(s), minimo {args.min_sequences} "
              f"/ only {len(records)} sequence(s) retrieved, minimum {args.min_sequences}.")
        print("             Resposta provavelmente truncada / Response likely truncated.")
        print("             Nenhum arquivo foi escrito / No file was written.")
        return 1

    if not args.skip_outgroup:
        got = {r['accession'] for r in outgroup}
        missing = [a for a in OUTGROUP_ACCESSIONS if a not in got]
        if missing:
            print(f"[ERRO/ERROR] Grupo externo incompleto / Incomplete outgroup: faltam / missing "
                  f"{', '.join(missing)}")
            print("             Nenhum arquivo foi escrito / No file was written.")
            return 1

    raw_count = len(records)
    if not args.keep_redundant:
        records = deduplicate(records)
    reviewed = sum(1 for r in records if r['status'] == 'reviewed')

    provenance = args.provenance or (os.path.splitext(args.out)[0] + '.provenance.txt')
    digest = write_set(records, args.out, provenance, url, headers, extra=[
        ("sequences_raw", raw_count),
        ("deduplication", "none" if args.keep_redundant
         else "one per (gene, species); reviewed first, then lowest accession"),
        ("reviewed", reviewed),
        ("unreviewed", len(records) - reviewed),
    ])

    print(f"\n[OK] {len(records)} referencias reais escritas em / real references written to: {args.out}")
    print(f"     {raw_count} na resposta / in the response; {reviewed} revisadas / reviewed, "
          f"{len(records) - reviewed} nao revisadas / unreviewed")
    print(f"[OK] Procedencia registrada em / Provenance recorded in: {provenance}")
    print(f"     UniProt release {headers.get('X-UniProt-Release', 'unknown')} | sha256 {digest[:16]}...")

    if not args.skip_outgroup:
        outgroup_prov = os.path.splitext(args.outgroup_out)[0] + '.provenance.txt'
        write_set(outgroup, args.outgroup_out, outgroup_prov, outgroup_url, outgroup_headers)
        print(f"[OK] Grupo externo ({len(outgroup)} UreC) escrito em / Outgroup written to: "
              f"{args.outgroup_out}")
    return 0


if __name__ == "__main__":
    sys.exit(main())

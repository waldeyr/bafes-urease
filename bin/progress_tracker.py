#!/usr/bin/env python3
"""
progress_tracker.py
PT-BR: Monitor de progresso do pipeline, alimentado pelo arquivo de trace REAL do Nextflow.
EN-US: Pipeline progress monitor, driven by Nextflow's REAL trace file.

PT-BR: Este monitor nao simula nada. Ele le `results/nf_trace.txt` (habilitado em
       nextflow.config) e conta tarefas concluidas de verdade. Ele tambem NUNCA anuncia
       sucesso: quem sabe se a execucao deu certo e o codigo de saida do Nextflow, e quem
       reporta isso e o run.sh. A versao anterior era um `time.sleep()` de 8 segundos que
       imprimia "Pipeline concluido com sucesso!" enquanto o Bakta ainda rodava — em
       .logs/exec-20260728-143927.log ela declarou sucesso 8 horas antes de a execucao
       falhar de fato na etapa de sintenia.
EN-US: This monitor simulates nothing. It reads `results/nf_trace.txt` (enabled in
       nextflow.config) and counts tasks that actually completed. It also NEVER announces
       success: the authority on whether a run succeeded is Nextflow's exit code, and run.sh
       is what reports it. The previous version was an 8-second `time.sleep()` loop that
       printed "Pipeline completed successfully!" while Bakta was still running — in
       .logs/exec-20260728-143927.log it declared success 8 hours before the run actually
       failed at the synteny step.
"""

import sys
import os
import csv
import time
import signal
import argparse

# PT-BR: O tracker roda no host, fora do container, onde o tqdm pode nao existir.
#        Ele e puramente cosmetico, entao a ausencia da lib degrada para texto simples
#        em vez de derrubar o monitor com um traceback.
# EN-US: The tracker runs on the host, outside the container, where tqdm may be absent.
#        It is purely cosmetic, so a missing library degrades to plain text instead of
#        killing the monitor with a traceback.
try:
    from tqdm import tqdm
except ImportError:
    tqdm = None

# PT-BR: Fases do pipeline: rotulo, processos do Nextflow que a compoem, e quantas tarefas
#        por estirpe cada uma gera (0 = processo agregado, roda uma vez so).
# EN-US: Pipeline phases: label, the Nextflow processes making it up, and how many tasks per
#        strain each one spawns (0 = aggregate process, runs exactly once).
PHASES = [
    ("Pre-flight Check",                  ["PREFLIGHT_CHECK"],                    0),
    ("Download Genomes (NCBI)",           ["DOWNLOAD_GENOME"],                    1),
    ("Quality Control (QC)",              ["QC_QUAST", "QC_CHECKM2"],             1),
    ("Functional Annotation (Bakta)",     ["BAKTA_ANNOTATE"],                     1),
    ("Candidate Screening (BLAST/HMMER)", ["EXTRACT_CANDIDATES"],                 1),
    ("Presence/Absence Matrix",           ["MERGE_ALL_STRAINS"],                  0),
    ("Urease Locus Extraction",           ["EXTRACT_LOCUS"],                      1),
    ("Synteny Analysis (clinker)",        ["ANALYZE_SYNTENY"],                    0),
    ("UreC Phylogeny (IQ-TREE)",          ["BUILD_PHYLOGENY"],                    0),
]

# PT-BR: Status do trace que contam como tarefa encerrada com exito.
# EN-US: Trace statuses that count as a task finished successfully.
DONE_STATUSES = {"COMPLETED", "CACHED"}
FAILED_STATUSES = {"FAILED", "ABORTED"}
# PT-BR: Tarefas ainda em voo. Contar isto e o que da vida ao monitor durante as horas de
#        Bakta, em que nada conclui — e, no fim, uma tarefa que ficou RUNNING e justamente a
#        assinatura de uma morta por estouro de `time`: o Nextflow aborta na
#        IllegalThreadStateException antes de marcar a linha como FAILED (foi o que houve com
#        ANALYZE_SYNTENY em .logs/exec-20260728-143927.log).
# EN-US: Tasks still in flight. Counting these is what keeps the monitor alive during the
#        hours of Bakta, when nothing completes — and, at the end, a task left RUNNING is
#        precisely the signature of one killed by a blown `time`: Nextflow aborts on the
#        IllegalThreadStateException before marking the row FAILED (which is what happened to
#        ANALYZE_SYNTENY in .logs/exec-20260728-143927.log).
RUNNING_STATUSES = {"RUNNING", "SUBMITTED", "NEW"}

_stop = False


def _handle_stop(signum, frame):
    global _stop
    _stop = True


def count_strains(accessions_path):
    """
    PT-BR: Numero de estirpes = linhas de dados do TSV de accessions. Sem ele, o total
           esperado e desconhecido e o monitor mostra so contagens absolutas.
    EN-US: Number of strains = data rows in the accessions TSV. Without it the expected
           total is unknown and the monitor shows plain counts only.
    """
    if not accessions_path or not os.path.isfile(accessions_path):
        return None
    try:
        with open(accessions_path, newline="") as fh:
            rows = [r for r in csv.DictReader(fh, delimiter="\t") if any(v.strip() for v in r.values() if v)]
        return len(rows) or None
    except OSError:
        return None


def process_name(task_name):
    """PT-BR: "BAKTA_ANNOTATE (S1)" -> "BAKTA_ANNOTATE". / EN-US: same."""
    return task_name.split(" (", 1)[0].strip()


def read_trace(trace_path):
    """
    PT-BR: Le o trace e devolve {processo: {"done": n, "failed": n, "running": n}}. O Nextflow
           reescreve o arquivo durante a execucao, entao uma leitura parcial pode pegar linha
           truncada — linhas malformadas sao ignoradas, nao estimadas.
    EN-US: Reads the trace and returns {process: {"done": n, "failed": n, "running": n}}.
           Nextflow rewrites the file during the run, so a partial read may catch a truncated
           line — malformed lines are skipped, never guessed at.
    """
    counts = {}
    try:
        with open(trace_path, newline="") as fh:
            reader = csv.DictReader(fh, delimiter="\t")
            if not reader.fieldnames or "name" not in reader.fieldnames:
                return counts
            for row in reader:
                name = row.get("name")
                status = (row.get("status") or "").strip().upper()
                if not name:
                    continue
                bucket = counts.setdefault(
                    process_name(name), {"done": 0, "failed": 0, "running": 0}
                )
                if status in DONE_STATUSES:
                    bucket["done"] += 1
                elif status in FAILED_STATUSES:
                    bucket["failed"] += 1
                elif status in RUNNING_STATUSES:
                    bucket["running"] += 1
    except (OSError, csv.Error):
        return counts
    return counts


def trace_is_fresh(trace_path, started):
    """
    PT-BR: `trace.overwrite = true` só trunca o arquivo quando o Nextflow sobe. Entre o
           lançamento do monitor e esse momento, o arquivo em disco ainda é o da execução
           ANTERIOR — lê-lo mostraria "Bakta 10/10" antes de a execução atual anotar nada.
           Só contamos o trace depois que ele for modificado após o início do monitor.
    EN-US: `trace.overwrite = true` only truncates the file once Nextflow starts. Between
           launching the monitor and that moment, the file on disk is still the PREVIOUS
           run's — reading it would show "Bakta 10/10" before the current run recorded
           anything. We only count the trace once it has been modified after the monitor started.
    """
    try:
        return os.path.getmtime(trace_path) >= started
    except OSError:
        return False


def phase_totals(counts, n_strains):
    """
    PT-BR: Agrega por fase. Devolve lista de (rotulo, feitas, falhas, rodando, esperadas|None).
    EN-US: Aggregates per phase. Returns (label, done, failed, running, expected|None) tuples.
    """
    out = []
    for label, procs, per_strain in PHASES:
        done = sum(counts.get(p, {}).get("done", 0) for p in procs)
        failed = sum(counts.get(p, {}).get("failed", 0) for p in procs)
        running = sum(counts.get(p, {}).get("running", 0) for p in procs)
        if per_strain and n_strains:
            expected = len(procs) * n_strains
        elif per_strain:
            expected = None
        else:
            expected = len(procs)
        out.append((label, done, failed, running, expected))
    return out


def render_plain(totals, elapsed, last_render):
    """
    PT-BR: Sem TTY (o caso do `screen`, onde a saida do Nextflow se mistura), imprime uma
           linha so quando algo muda de verdade — nao a cada poll.
    EN-US: With no TTY (the `screen` case, where Nextflow's own output interleaves), prints a
           line only when something actually changes — not on every poll.
    """
    snapshot = tuple((d, f, r) for _, d, f, r, _ in totals)
    if snapshot == last_render:
        return last_render

    parts = []
    for label, done, failed, running, expected in totals:
        if done == 0 and failed == 0 and running == 0:
            continue
        target = f"/{expected}" if expected else ""
        flag = ""
        if running:
            flag += f" +{running} rodando/running"
        if failed:
            flag += f" ({failed} FAILED)"
        parts.append(f"{label}: {done}{target}{flag}")

    mins, secs = divmod(int(elapsed), 60)
    hours, mins = divmod(mins, 60)
    stamp = f"{hours:02d}:{mins:02d}:{secs:02d}"
    if parts:
        print(f"  [{stamp}] " + " | ".join(parts), flush=True)
    else:
        print(f"  [{stamp}] Nenhuma tarefa concluída ainda / No task completed yet.", flush=True)
    return snapshot


def run_tracker(trace_path, accessions_path, poll_seconds):
    n_strains = count_strains(accessions_path)

    print("\n[PT-BR] Monitorando o trace real do Nextflow / [EN-US] Monitoring Nextflow's real trace")
    print(f"        {trace_path}")
    if n_strains:
        print(f"        {n_strains} estirpe(s) / strain(s)\n", flush=True)
    else:
        print("        [AVISO/WARNING] Total de estirpes desconhecido; sem percentual."
              " / Strain total unknown; no percentage.\n", flush=True)

    interactive = sys.stdout.isatty() and tqdm is not None
    started = time.time()
    last_render = None

    bars = {}
    if interactive:
        for position, (label, _, _) in enumerate(PHASES):
            bars[label] = tqdm(total=100, desc=f"{label:<34}", position=position,
                               leave=True, bar_format="{desc} {n:3.0f}% |{bar}| {postfix}")

    while not _stop:
        counts = read_trace(trace_path) if trace_is_fresh(trace_path, started) else {}
        totals = phase_totals(counts, n_strains)

        if interactive:
            for label, done, failed, running, expected in totals:
                bar = bars[label]
                bar.n = min(100, int(100 * done / expected)) if expected else 0
                note = f"{done}/{expected}" if expected else f"{done}"
                if running:
                    note += f"  +{running} running"
                if failed:
                    note += f"  {failed} FAILED"
                bar.set_postfix_str(note, refresh=False)
                bar.refresh()
        else:
            last_render = render_plain(totals, time.time() - started, last_render)

        # PT-BR: Dorme em fatias curtas para responder ao sinal de parada sem atraso.
        # EN-US: Sleeps in short slices so the stop signal is answered without delay.
        waited = 0.0
        while waited < poll_seconds and not _stop:
            time.sleep(0.5)
            waited += 0.5

    if interactive:
        for bar in bars.values():
            bar.close()

    # PT-BR: Resumo final factual. Nenhuma afirmacao sobre sucesso — o veredito e do
    #        codigo de saida do Nextflow, impresso pelo run.sh.
    # EN-US: Factual closing summary. No claim about success — the verdict belongs to
    #        Nextflow's exit code, printed by run.sh.
    counts = read_trace(trace_path) if os.path.isfile(trace_path) else {}
    totals = phase_totals(counts, n_strains)
    total_done = sum(d for _, d, _, _, _ in totals)
    total_failed = sum(f for _, _, f, _, _ in totals)
    total_running = sum(r for _, _, _, r, _ in totals)

    print("\n--- Tarefas registradas no trace / Tasks recorded in the trace ---", flush=True)
    for label, done, failed, running, expected in totals:
        target = f"/{expected}" if expected else ""
        flag = ""
        if running:
            flag += f"   {running} ainda RUNNING/still RUNNING"
        if failed:
            flag += f"   {failed} FALHOU/FAILED"
        print(f"  {label:<34} {done}{target}{flag}", flush=True)
    print(f"  {'TOTAL':<34} {total_done} concluída(s)/done, "
          f"{total_failed} falha(s)/failed, {total_running} em voo/in flight", flush=True)

    if total_running:
        # PT-BR: Uma tarefa que ficou RUNNING depois do fim da execução quase sempre foi
        #        morta por estouro do `time` do processo: o Nextflow aborta com
        #        "process hasn't exited" antes de marcar a linha como FAILED. Confira o
        #        `time` do processo em nextflow.config e o .exitcode no work dir (143 = SIGTERM).
        # EN-US: A task left RUNNING after the run ended was almost always killed by a blown
        #        process `time` limit: Nextflow aborts with "process hasn't exited" before it
        #        can mark the row FAILED. Check that process's `time` in nextflow.config and
        #        the .exitcode in its work dir (143 = SIGTERM).
        print("  [PT-BR] Tarefa(s) presas em RUNNING: quase sempre estouro do `time` do")
        print("          processo — o Nextflow aborta com \"process hasn't exited\" antes de")
        print("          marcar FAILED. Veja o `time` em nextflow.config e o .exitcode no")
        print("          work dir (143 = SIGTERM).")
        print("  [EN-US] Task(s) stuck in RUNNING: almost always a blown process `time` limit —")
        print("          Nextflow aborts with \"process hasn't exited\" before marking FAILED.")
        print("          Check `time` in nextflow.config and the .exitcode in its work dir")
        print("          (143 = SIGTERM).", flush=True)

    print("  [PT-BR] O resultado da execução vem do código de saída do Nextflow, abaixo.")
    print("  [EN-US] The run's outcome comes from Nextflow's exit code, below.\n", flush=True)
    return 0


def main():
    parser = argparse.ArgumentParser(
        description="Nextflow progress monitor driven by the trace file (PT-BR / EN-US)"
    )
    parser.add_argument("--trace", default="results/nf_trace.txt",
                        help="Arquivo de trace do Nextflow / Nextflow trace file")
    parser.add_argument("--accessions", default="data/accessions.tsv",
                        help="TSV de accessions, para saber o total de estirpes / "
                             "Accessions TSV, to learn the strain total")
    parser.add_argument("--poll", type=float, default=15.0,
                        help="Intervalo entre leituras, em segundos / Seconds between reads")
    parser.add_argument("--verbose", action="store_true",
                        help="Mantido por compatibilidade com o run.sh / Kept for run.sh compatibility")
    args = parser.parse_args()

    signal.signal(signal.SIGTERM, _handle_stop)
    signal.signal(signal.SIGINT, _handle_stop)

    return run_tracker(args.trace, args.accessions, max(1.0, args.poll))


if __name__ == "__main__":
    sys.exit(main())

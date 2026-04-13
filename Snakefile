configfile: "config/config.yaml"

from pathlib import Path

shell.executable("/bin/bash")

OUTDIR = config.get("output_dir", "results")
ASSAY = config.get("assay_type", "rt_primer_25")

INPUT_GENOMES = config["input"]["genomes_fasta"]
REF_FASTA = config["input"]["ref_fasta"]
REF_ID = config["input"]["ref_id"]

QC_DIR = f"{OUTDIR}/01_qc"
ALN_DIR = f"{OUTDIR}/02_alignment"
CAND_DIR = f"{OUTDIR}/03_candidates"
CONS_DIR = f"{OUTDIR}/04_conservation"
ACC_DIR = f"{OUTDIR}/05_accessibility"
SPEC_DIR = f"{OUTDIR}/06_specificity"
PANEL_DIR = f"{OUTDIR}/07_panel"
LOG_DIR = f"{OUTDIR}/logs"

CLEAN_FASTA = f"{QC_DIR}/rsvb_clean.fasta"
QC_TSV = f"{QC_DIR}/rsvb_qc.tsv"
PLUS_REF_FASTA = f"{ALN_DIR}/rsvb_plus_ref.fasta"
PLUS_REF_ALN = f"{ALN_DIR}/rsvb_plus_ref.aln.fa"
CANDIDATES_TSV = f"{CAND_DIR}/candidates.tsv"
CONSERVATION_TSV = f"{CONS_DIR}/candidates.conservation.tsv"
LUNP_FILE = f"{ACC_DIR}/{REF_ID}_lunp"
ACCESS_TSV = f"{ACC_DIR}/candidates.access.tsv"
SCORED_TSV = f"{SPEC_DIR}/candidates.scored.tsv"
RANKED_TSV = f"{PANEL_DIR}/ranked.strict.tsv"
PANEL_TSV = f"{PANEL_DIR}/panel.strict.tsv"


def bool_flag(value, flag):
    return flag if bool(value) else ""


def optional_flag(value, flag):
    if value is None:
        return ""
    sval = str(value).strip()
    if sval == "" or sval.lower() == "none":
        return ""
    return f"{flag} {sval}"


def specificity_optional_args():
    sp = config.get("specificity", {})
    pieces = []
    pieces.append(optional_flag(sp.get("human_blast_db_prefix", ""), "--human-blast-db"))
    pieces.append(optional_flag(sp.get("virus_blast_db_prefix", ""), "--virus-blast-db"))
    pieces.append(optional_flag(sp.get("human_blast_tsv", ""), "--human-blast-tsv"))
    pieces.append(optional_flag(sp.get("virus_blast_tsv", ""), "--virus-blast-tsv"))
    pieces.append(optional_flag(sp.get("mfeprimer_tsv", ""), "--mfeprimer-tsv"))
    pieces.append(optional_flag(sp.get("cache_dir", f"{SPEC_DIR}/.blast_cache"), "--cache-dir"))
    return " ".join([p for p in pieces if p])


rule all:
    input:
        PANEL_TSV,
        RANKED_TSV,
        SCORED_TSV,
        ACCESS_TSV,
        CONSERVATION_TSV,
        CANDIDATES_TSV,
        PLUS_REF_ALN,
        LUNP_FILE,
        QC_TSV,


rule qc_rsvb:
    input:
        fasta=INPUT_GENOMES,
    output:
        clean=CLEAN_FASTA,
        qc=QC_TSV,
    log:
        f"{LOG_DIR}/qc_rsvb.log",
    threads: 1
    resources:
        mem_mb=4000,
        time_min=30,
    conda:
        "envs/python_bio.yaml"
    shell:
        r"""
        mkdir -p {QC_DIR} {LOG_DIR}
        python scripts/qc_rsvb.py \
          --input-fasta {input.fasta} \
          --output-fasta {output.clean} \
          --qc-tsv {output.qc} \
          --min-length {config[qc][min_length]} \
          --max-n-fraction {config[qc][max_n_fraction]} \
          > {log} 2>&1
        """


rule combine_ref_and_clean:
    input:
        ref=REF_FASTA,
        clean=CLEAN_FASTA,
    output:
        PLUS_REF_FASTA,
    log:
        f"{LOG_DIR}/combine_ref_and_clean.log",
    run:
        Path(ALN_DIR).mkdir(parents=True, exist_ok=True)
        Path(LOG_DIR).mkdir(parents=True, exist_ok=True)
        with open(output[0], "wb") as out:
            for fp in [input.ref, input.clean]:
                with open(fp, "rb") as fh:
                    out.write(fh.read())
        with open(log[0], "w") as lg:
            lg.write(f"Wrote combined FASTA -> {output[0]}\n")


rule mafft_align:
    input:
        PLUS_REF_FASTA,
    output:
        PLUS_REF_ALN,
    log:
        f"{LOG_DIR}/mafft_align.log",
    threads: lambda wc: int(config["threads"].get("mafft", 4))
    resources:
        mem_mb=16000,
        time_min=240,
    conda:
        "envs/mafft.yaml"
    shell:
        r"""
        mkdir -p {ALN_DIR} {LOG_DIR}
        {config[alignment][mafft_bin]} {config[alignment][mafft_args]} --thread {threads} {input} > {output} 2> {log}
        """


rule generate_candidates:
    input:
        ref=REF_FASTA,
    output:
        CANDIDATES_TSV,
    log:
        f"{LOG_DIR}/generate_candidates.log",
    threads: 1
    resources:
        mem_mb=4000,
        time_min=30,
    conda:
        "envs/python_bio.yaml"
    params:
        require_terminal=lambda wc: bool_flag(config["candidates"].get("require_terminal_matchable_base", True), "--require-terminal-matchable-base"),
    shell:
        r"""
        mkdir -p {CAND_DIR} {LOG_DIR}
        python scripts/generate_candidates.py \
          --ref-fasta {input.ref} \
          --output-tsv {output} \
          --assay-type {ASSAY} \
          --window {config[candidates][window]} \
          --step {config[candidates][step]} \
          --min-gc {config[candidates][min_gc]} \
          --max-gc {config[candidates][max_gc]} \
          --max-homopolymer {config[candidates][max_homopolymer]} \
          --max-low-complexity {config[candidates][max_low_complexity]} \
          --low-complexity-k {config[candidates][low_complexity_k]} \
          --min-tm {config[candidates][min_tm]} \
          --max-tm {config[candidates][max_tm]} \
          --min-gc-3p5 {config[candidates][min_gc_3p5]} \
          --max-gc-3p5 {config[candidates][max_gc_3p5]} \
          --max-3p-homopolymer {config[candidates][max_3p_homopolymer]} \
          {params.require_terminal} \
          > {log} 2>&1
        """


rule score_conservation:
    input:
        aln=PLUS_REF_ALN,
        cand=CANDIDATES_TSV,
    output:
        CONSERVATION_TSV,
    log:
        f"{LOG_DIR}/score_conservation.log",
    threads: 1
    resources:
        mem_mb=8000,
        time_min=60,
    conda:
        "envs/python_bio.yaml"
    shell:
        r"""
        mkdir -p {CONS_DIR} {LOG_DIR}
        python scripts/score_conservation.py \
          --alignment-fasta {input.aln} \
          --ref-id {REF_ID} \
          --candidates-tsv {input.cand} \
          --output-tsv {output} \
          --assay-type {ASSAY} \
          --three-prime-nt {config[conservation][three_prime_nt]} \
          > {log} 2>&1
        """


rule rnaplfold:
    input:
        ref=REF_FASTA,
    output:
        lunp=LUNP_FILE,
    log:
        f"{LOG_DIR}/rnaplfold.log",
    threads: 1
    resources:
        mem_mb=8000,
        time_min=60,
    conda:
        "envs/viennarna.yaml"
    params:
        no_lp=lambda wc: "--noLP" if bool(config["accessibility"].get("plfold_no_lp", False)) else "",
        bin=lambda wc: config["accessibility"].get("rnaplfold_bin", "RNAplfold"),
        ulen=lambda wc: config["accessibility"].get("u_length", 8),
        win=lambda wc: config["accessibility"].get("plfold_window", 150),
        span=lambda wc: config["accessibility"].get("plfold_max_span", 100),
    shell:
        r"""
        set -euo pipefail

        mkdir -p {ACC_DIR} {LOG_DIR}

        out_lunp=$(python -c 'import os,sys; print(os.path.abspath(sys.argv[1]))' "{output.lunp}")
        out_log=$(python -c 'import os,sys; print(os.path.abspath(sys.argv[1]))' "{log}")
        in_ref=$(python -c 'import os,sys; print(os.path.abspath(sys.argv[1]))' "{input.ref}")

        tmpdir=$(mktemp -d)
        trap 'rm -rf "$tmpdir"' EXIT

        cp "$in_ref" "$tmpdir/ref.fa"
        cd "$tmpdir"

        "{params.bin}" -W {params.win} -L {params.span} -u {params.ulen} {params.no_lp} \
          < ref.fa > /dev/null 2> "$out_log"

        lunp=$(find . -maxdepth 1 -name '*_lunp' | head -n 1)
        test -n "$lunp"
        cp "$lunp" "$out_lunp"
        """


rule score_accessibility:
    input:
        cand=CONSERVATION_TSV,
        lunp=LUNP_FILE,
    output:
        ACCESS_TSV,
    log:
        f"{LOG_DIR}/score_accessibility.log",
    threads: 1
    resources:
        mem_mb=8000,
        time_min=60,
    conda:
        "envs/python_bio.yaml"
    shell:
        r"""
        mkdir -p {ACC_DIR} {LOG_DIR}
        python scripts/score_accessibility.py \
          --candidates-tsv {input.cand} \
          --lunp-file {input.lunp} \
          --u-length {config[accessibility][u_length]} \
          --three-prime-nt {config[accessibility][three_prime_nt]} \
          --output-tsv {output} \
          --assay-type {ASSAY} \
          > {log} 2>&1
        """


rule score_specificity:
    input:
        ACCESS_TSV,
    output:
        SCORED_TSV,
    log:
        f"{LOG_DIR}/score_specificity.log",
    threads: lambda wc: int(config["threads"].get("specificity", 8))
    resources:
        mem_mb=32000,
        time_min=360,
    conda:
        "envs/blast.yaml"
    params:
        optional=lambda wc: specificity_optional_args(),
    shell:
        r"""
        set -euo pipefail
        mkdir -p {SPEC_DIR} {LOG_DIR}

        python scripts/score_specificity.py \
          --candidates-tsv {input} \
          --output-tsv {output} \
          --assay-type {ASSAY} \
          --blastn-bin {config[specificity][blastn_bin]} \
          --task {config[specificity][task]} \
          --word-size {config[specificity][word_size]} \
          --evalue {config[specificity][evalue]} \
          --max-target-seqs {config[specificity][max_target_seqs]} \
          --max-hsps {config[specificity][max_hsps]} \
          --num-threads {threads} \
          --min-pident {config[specificity][min_pident]} \
          --min-query-cover {config[specificity][min_query_cover]} \
          --three-prime-anchor-nt {config[specificity][three_prime_anchor_nt]} \
          --ungapped {config[specificity][ungapped]} \
          --blast-perc-identity {config[specificity][blast_perc_identity]} \
          {params.optional} \
          > {log} 2>&1
        """


rule select_panel:
    input:
        SCORED_TSV,
    output:
        ranked=RANKED_TSV,
        panel=PANEL_TSV,
    log:
        f"{LOG_DIR}/select_panel.log",
    threads: 1
    resources:
        mem_mb=8000,
        time_min=60,
    conda:
        "envs/python_bio.yaml"
    shell:
        r"""
        mkdir -p {PANEL_DIR} {LOG_DIR}
        python scripts/select_panel.py \
          --scored-tsv {input} \
          --ranked-output {output.ranked} \
          --panel-output {output.panel} \
          --assay-type {ASSAY} \
          --max-sites {config[panel][max_sites]} \
          --min-site-score {config[panel][min_site_score]} \
          --min-gap {config[panel][min_gap]} \
          --hard-min-gap {config[panel][hard_min_gap]} \
          --max-overlap-fraction {config[panel][max_overlap_fraction]} \
          --min-new-bp {config[panel][min_new_bp]} \
          --cross-penalty-weight {config[panel][cross_penalty_weight]} \
          --cross-penalty-sum-weight {config[panel][cross_penalty_sum_weight]} \
          --safe-run {config[panel][safe_run]} \
          --risky-run {config[panel][risky_run]} \
          --safe-terminal-run {config[panel][safe_terminal_run]} \
          --risky-terminal-run {config[panel][risky_terminal_run]} \
          --safe-match-fraction {config[panel][safe_match_fraction]} \
          --risky-match-fraction {config[panel][risky_match_fraction]} \
          --hard-max-terminal-run {config[panel][hard_max_terminal_run]} \
          --hard-max-run {config[panel][hard_max_run]} \
          --min-accessibility-score {config[panel][min_accessibility_score]} \
          --min-access-3p-terminal {config[panel][min_access_3p_terminal]} \
          --min-robustness-score {config[panel][min_robustness_score]} \
          --min-cov-3p-1mm {config[panel][min_cov_3p_1mm]} \
          --min-cov-terminal-match {config[panel][min_cov_terminal_match]} \
          --max-human-hits {config[panel][max_human_hits]} \
          --max-human-anchored-hits {config[panel][max_human_anchored_hits]} \
          --max-human-nearperfect-hits {config[panel][max_human_nearperfect_hits]} \
          --max-virus-anchored-hits {config[panel][max_virus_anchored_hits]} \
          > {log} 2>&1
        """

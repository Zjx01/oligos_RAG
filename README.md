# Snakemake workflow manual for RT-primer panel design

This document describes how to run and maintain the Snakemake workflow for **25 nt biotinylated RT-primer panel design**.

---

## 1. Purpose of the workflow

This workflow is designed to build a **short RT-primer panel** for viral genomes such as RSV or HAV.

The workflow performs these stages:

1. clean the viral genome set
2. align the genomes together with a reference
3. generate candidate 25 nt primer-binding sites
4. score each candidate for:
   - conservation / robustness
   - RNA accessibility
   - host off-target risk
   - thermo / synthesis quality
5. select a final primer panel, with emphasis on:
   - good biological quality
   - low off-target risk
   - better genome-wide spacing

This workflow is intended for a **biotinylated 25 nt RT-primer strategy**, not a classical long 80–120 nt capture-bait workflow.

---

## 2. Workflow overview

### Stage 1: QC
Input viral genomes are filtered to remove poor-quality or duplicate entries.

Outputs:
- cleaned FASTA
- QC summary TSV

### Stage 2: Alignment
The cleaned sequences are combined with the reference FASTA and aligned with MAFFT.

This is important because candidate coordinates are defined on the reference, while conservation is estimated from the whole alignment.

Outputs:
- combined FASTA
- multiple sequence alignment

### Stage 3: Candidate generation
Sliding windows are generated on the reference genome.

For RT-primer mode, this stage focuses on:
- 25 nt length
- GC range
- Tm range
- 3′-end behavior
- homopolymers
- low complexity

Outputs:
- `candidates.tsv`

### Stage 4: Conservation scoring
Each candidate is checked against the alignment.

Important outputs include:
- whole-primer mismatch tolerance
- 3′-end mismatch tolerance
- terminal match behavior
- robustness score

Outputs:
- `candidates.conservation.tsv`

### Stage 5: Accessibility scoring
RNAplfold is used to estimate whether the **target RNA region** is structurally exposed enough for primer binding.

For RT primers, this step gives more weight to:
- local unpaired probability
- 3′-proximal accessibility

Outputs:
- `_lunp` file from RNAplfold
- `candidates.access.tsv`

### Stage 6: Specificity scoring
Each candidate is checked against a human BLAST database.

For RT-primer mode, this stage penalizes:
- general host hits
- 3′-anchored host hits
- near-perfect host hits

Outputs:
- `candidates.scored.tsv`

### Stage 7: Panel selection
The final selector chooses a primer panel that tries to:
- keep high-quality primers
- avoid dangerous primers
- spread primers more evenly across the genome

Outputs:
- ranked candidate table
- final panel TSV

---

## 3. Expected directory structure

A typical repo layout is:

```text
oligos_RAG/
├── Snakefile
├── config/
│   └── config.yaml
├── envs/
│   ├── python_bio.yaml
│   ├── blast.yaml
│   ├── viennarna.yaml
│   └── mafft.yaml
├── scripts/
│   ├── qc_rsvb.py
│   ├── generate_candidates.py
│   ├── score_conservation.py
│   ├── score_accessibility.py
│   ├── score_specificity.py
│   └── select_panel.py
├── data/
│   ├── genomes.fa
│   └── ref.fa
└── results/
```

---

## 4. Required inputs

You need at least these files.

### `genomes.fa`
A FASTA containing the viral genomes to be used for design.

These should ideally be:
- near full-length
- low N content
- relevant to the strain diversity you want to capture

### `ref.fa`
A FASTA containing **exactly one reference sequence**.

This reference is used for:
- candidate coordinate generation
- alignment coordinate mapping
- RNA accessibility scoring

### Human BLAST DB prefix
This is not a directory and not a FASTA file.

It should be a BLAST database **prefix**, for example:

```text
/lustre/home/zhaoj11/databases/blastdb_human/human_bg
```

not:

```text
/lustre/home/zhaoj11/databases/blastdb_human/
```

---

## 5. Example config structure

Your `config.yaml` should typically contain these sections.

### `input`
Defines core input files.

```yaml
input:
  genomes_fasta: data/rsvb_genomes.fasta
  ref_fasta: data/PP_002W6UA.1.fa
  ref_id: PP_002W6UA.1
```

### `qc`
Controls sequence filtering.

```yaml
qc:
  min_length: 14000
  max_n_fraction: 0.01
```

### `alignment`
Controls MAFFT.

```yaml
alignment:
  mafft_bin: mafft
  mafft_args: --auto
```

### `candidates`
Controls candidate generation.

```yaml
candidates:
  window: 25
  step: 1
  min_gc: 0.32
  max_gc: 0.64
  min_tm: 54
  max_tm: 64
  max_homopolymer: 5
  max_low_complexity: 0.85
  low_complexity_k: 3
  min_gc_3p5: 1
  max_gc_3p5: 4
  max_3p_homopolymer: 3
```

### `conservation`
Controls 3′-end evaluation.

```yaml
conservation:
  three_prime_nt: 8
```

### `accessibility`
Controls RNAplfold.

```yaml
accessibility:
  rnaplfold_bin: RNAplfold
  u_length: 8
  three_prime_nt: 8
  plfold_window: 150
  plfold_max_span: 100
  plfold_no_lp: false
```

### `specificity`
Controls BLAST specificity scoring.

```yaml
specificity:
  blastn_bin: blastn
  task: auto
  word_size: 7
  evalue: 1000.0
  max_target_seqs: 10
  max_hsps: 1
  min_pident: 80
  min_query_cover: 0.80
  three_prime_anchor_nt: 8
  ungapped: auto
  blast_perc_identity: auto
  human_blast_db_prefix: /lustre/home/zhaoj11/databases/blastdb_human/human_bg
  cache_dir: results/06_specificity/.blast_cache
```

### `panel`
Controls final panel stringency.

```yaml
panel:
  max_sites: 24
  min_site_score: 0.56
  min_gap: 250
  hard_min_gap: 40
  max_overlap_fraction: 0.20
  min_new_bp: 12
  cross_penalty_weight: 0.22
  cross_penalty_sum_weight: 0.08
  safe_run: 6
  risky_run: 9
  safe_terminal_run: 3
  risky_terminal_run: 5
  safe_match_fraction: 0.50
  risky_match_fraction: 0.75
  hard_max_terminal_run: 7
  hard_max_run: 10
  min_accessibility_score: 0.10
  min_access_3p_terminal: 0.15
  min_robustness_score: 0.90
  min_cov_3p_1mm: 0.95
  min_cov_terminal_match: 0.95
  max_human_hits: 10
  max_human_anchored_hits: 0
  max_human_nearperfect_hits: 0
  max_virus_anchored_hits: 0
```

---

## 6. Main outputs

### `results/01_qc/`
- cleaned FASTA
- QC TSV

### `results/02_alignment/`
- reference + cleaned FASTA
- MAFFT alignment

### `results/03_candidates/`
- candidate primer table

### `results/04_conservation/`
- conservation-scored candidate table

### `results/05_accessibility/`
- RNAplfold `_lunp`
- accessibility-scored candidate table

### `results/06_specificity/`
- specificity-scored candidate table
- BLAST cache directory

### `results/07_panel/`
- ranked candidate list
- final panel TSV

---

## 7. How to run the workflow

### Dry run

```bash
snakemake -n -p
```

### Local execution

```bash
snakemake --use-conda --cores 8 -p
```

### SLURM execution

```bash
snakemake \
  --use-conda \
  --jobs 20 \
  --latency-wait 60 \
  --rerun-incomplete \
  --default-resources mem_mb=8000 time_min=120 \
  --cluster "sbatch -c {threads} --mem={resources.mem_mb}M -t {resources.time_min}" \
  -p
```

### Force rerun of final steps

```bash
snakemake \
  --use-conda \
  --jobs 20 \
  --latency-wait 60 \
  --rerun-incomplete \
  --default-resources mem_mb=8000 time_min=120 \
  --cluster "sbatch -c {threads} --mem={resources.mem_mb}M -t {resources.time_min}" \
  -R score_specificity select_panel \
  -p
```

---

## 8. How to rerun after editing scripts

Changing a script does not always force Snakemake to rerun the rule automatically.

The safest options are:

### Option A: force rerun selected rules

```bash
snakemake \
  --use-conda \
  --jobs 20 \
  --latency-wait 60 \
  --rerun-incomplete \
  --default-resources mem_mb=8000 time_min=120 \
  --cluster "sbatch -c {threads} --mem={resources.mem_mb}M -t {resources.time_min}" \
  -R score_specificity select_panel \
  -p
```

### Option B: delete outputs and rerun

```bash
rm -f results/06_specificity/candidates.scored.tsv
rm -f results/07_panel/ranked.strict.tsv results/07_panel/panel.strict.tsv
snakemake --use-conda --cores 8 -p
```

---

## 9. Common problems and fixes

### A. `Reference ID ... not found in alignment`
Cause:
- the reference sequence used for candidate generation is not present in the alignment

Fix:
- include the reference FASTA together with cleaned sequences before MAFFT
- make sure `ref_id` exactly matches the FASTA header in the alignment

### B. `RNAplfold binary not found`
Cause:
- ViennaRNA is not available in the conda env

Fix:
- use an env that includes `viennarna`
- verify with:
```bash
which RNAplfold
RNAplfold --version
```

### C. `BLAST Database error: No alias or index file found`
Cause:
- `--human-blast-db` points to a directory, not a database prefix

Fix:
- use a BLAST prefix, for example:
```text
/lustre/home/zhaoj11/databases/blastdb_human/human_bg
```

### D. `score_specificity.py: error: unrecognized arguments: --cache-dir`
Cause:
- malformed optional-argument string in `Snakefile`

Fix:
- ensure optional arguments are joined with plain spaces, not backslash escapes
- `specificity_optional_args()` should return:
```python
" ".join(args)
```

### E. `AttributeError: Namespace object has no attribute 'second_pass'`
Cause:
- old `select_panel.py` still in use

Fix:
- replace it with the corrected selector
- rerun `select_panel`

### F. `rnaplfold` rule fails after `cd "$tmpdir"`
Cause:
- output and log paths become relative to the temp dir

Fix:
- resolve paths to absolute paths before `cd`
- the corrected `Snakefile` already does this

---

## 10. Primer orientation

Your final `panel.tsv` or `panel.strict.tsv` contains **target-site sequences**, not necessarily the oligos you order.

For RT primers, the primer you order is typically the **reverse complement** of `sequence_ref`.

So:

- `sequence_ref` = target sequence on the reference
- `rt_primer_seq` = reverse complement to order

For a final production workflow, it is best to include both columns.

---

## 11. Recommended operating procedure

Use this workflow in three phases.

### Phase 1: technical validation
Goal:
- make sure all rules run
- verify paths, envs, BLAST DB, RNAplfold

### Phase 2: biological filtering
Goal:
- inspect `candidates.scored.tsv`
- confirm accessibility and specificity are reasonable
- tune thresholds if needed

### Phase 3: panel optimization
Goal:
- produce a final panel with:
  - good specificity
  - good accessibility
  - stronger 3′ robustness
  - more even genome spacing

---

## 12. Best-practice notes

For **25 nt RT primers**:
- give more weight to:
  - 3′ exactness
  - 3′ accessibility
  - anchored host off-targets
- do not use long-bait assumptions from 80–120 nt capture workflows
- prefer evenly spaced primers if the goal is whole-genome recovery

For **whole-genome coverage**:
- avoid heavy clustering in one genomic region
- inspect the largest inter-primer gaps
- the selector should be spacing-aware, not only score-greedy

---

## 13. Summary

This Snakemake workflow is a **short RT-primer design pipeline** that:

- starts from a viral genome collection
- generates candidate primers on a fixed reference
- scores them for conservation, structure, and specificity
- selects a final primer panel with strong biological constraints

The most important things to remember are:

- the reference must be included in the alignment
- RNAplfold and BLAST must be correctly configured
- BLAST DB arguments must be **prefixes**
- RT primers should usually be **reverse complements** of `sequence_ref`
- rerun changed rules explicitly when you update scripts

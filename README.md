# RSVB biotinylated 25 nt RT-primer workflow

This package modifies the original long-capture-bait workflow for a new assay concept:

- **25 nt biotinylated oligos**
- oligo first binds RNA
- **reverse transcription elongation** is initiated from the oligo
- the elongated product is then **pulled down** for sequencing

Because of that assay change, the scoring logic is no longer long-bait-like.
The updated code now treats the candidates as **short RT primers** rather than 80–120 nt capture baits.

## Main changes

1. **Candidate generation is primer-aware**
   - default `--window 25`
   - default `--step 1`
   - stronger filters for GC, Tm, homopolymers, and low complexity
   - new 3' features:
     - `gc_3p5`
     - `terminal_base`
     - `terminal_is_gc`
     - `three_prime_max_homopolymer`

2. **Conservation scoring now cares about the 3' end**
   - still reports `cov_0mm`, `cov_1mm`, `cov_2mm`
   - adds:
     - `cov_3p_0mm`
     - `cov_3p_1mm`
     - `cov_terminal_match`
     - `mean_exact_3p_suffix_frac`
   - `robustness_score` is reweighted for short RT primers

3. **Accessibility scoring is 3'-weighted**
   - default `--u-length 8`
   - adds:
     - `access_3p_terminal`
     - `access_3p_near_mean`
   - `accessibility_score` now weights the primer 3' end more heavily

4. **Specificity scoring is short-primer-aware**
   - default BLAST behavior remains `auto`, which resolves to `blastn-short` for 25-mers
   - adds short-primer host-hit summaries:
     - `human_offtarget_anchored_hits`
     - `human_offtarget_nearperfect_hits`
     - same for virus background if used
   - `specificity_score` penalizes 3'-anchored host hits more strongly

5. **Panel selection is primer-aware**
   - denser default spacing than long-bait mode
   - more weight on thermo / specificity / 3'-sensitive robustness
   - stricter cross-dimer / terminal complementarity defaults for short oligos

```bash 
python scripts/qc_rsvb.py \
  --input-fasta data/rsvb_genomes.fasta \
  --output-fasta results/rsvb_clean.fasta \
  --qc-tsv results/rsvb_qc.tsv
```

```bash
mafft --auto results/rsvb_clean.fasta > results/rsvb_clean.aln.fa
```

### 1. Generate 25-mer RT-primer candidates
```bash
python scripts/generate_candidates.py \
  --ref-fasta data/PP_002W6UA.1.fa \
  --output-tsv results/candidates.tsv \
  --assay-type rt_primer_25
```

### 2. Score conservation using an alignment that includes the reference sequence
```bash
python scripts/score_conservation.py \
  --alignment-fasta results/rsvb_plus_ref.aln.fa \
  --ref-id PP_002W6UA.1 \
  --candidates-tsv results/candidates.tsv \
  --output-tsv results/candidates.conservation.tsv \
  --assay-type rt_primer_25 \
  --three-prime-nt 8
```

### 3. Score accessibility
```bash
python scripts/score_accessibility.py \
  --candidates-tsv results/candidates.conservation.tsv \
  --lunp-file results/RSVB_REF_lunp \
  --u-length 8 \
  --output-tsv results/candidates.access.tsv \
  --assay-type rt_primer_25 \
  --three-prime-nt 8
```

if you don't have a lunp file, you can input the ref.fasta and the system will handle it 
```bash
  python scripts/score_accessibility.py \
  --candidates-tsv results/candidates.conservation.tsv \
  --ref-fasta data/PP_002W6UA.1.fa \
  --u-length 8 \
  --output-tsv results/candidates.access.tsv \
  --assay-type rt_primer_25 \
  --three-prime-nt 8
```

### 4. Score specificity
```bash
python scripts/score_specificity.py \
  --candidates-tsv results/candidates.access.tsv \
  --human-blast-db ~/databases/blastdb_human/human_bg \
  --min-pident 80 \
  --min-query-cover 0.80 \
  --num-threads 8 \
  --output-tsv results/candidates.scored.tsv \
  --assay-type rt_primer_25 \
  --three-prime-anchor-nt 8
```

### 5. Select the final RT-primer panel
```bash
python scripts/select_panel.py \
  --scored-tsv results/candidates.scored.tsv \
  --ranked-output results/ranked.tsv \
  --panel-output results/panel.tsv \
  --assay-type rt_primer_25 \
  --max-sites 24
```


## Practical note
This code now matches the **short RT-primer assay concept** much better than the original long-bait workflow, but it is still a computational prioritization framework.
Final success will still depend on:
- RT chemistry
- hybridization conditions
- pull-down efficiency
- RNA fragmentation state
- whether your 25-mer can tolerate real 3' mismatches in the wet lab


# RSVB panel MVP framework

This is a minimal first-pass scaffold for RSVB whole-genome binding-site panel design.

## Included scripts

1. `scripts/qc_rsvb.py`
2. `scripts/generate_candidates.py`
3. `scripts/score_conservation.py`
4. `scripts/score_accessibility.py`
5. `scripts/score_specificity.py`
6. `scripts/select_panel.py`

## Suggested execution order

```bash
python scripts/qc_rsvb.py \
  --input-fasta data/rsvb_genomes.fasta \
  --output-fasta results/rsvb_clean.fasta \
  --qc-tsv results/rsvb_qc.tsv

mafft --auto results/rsvb_clean.fasta > results/rsvb_clean.aln.fa

python scripts/generate_candidates.py \
  --ref-fasta data/rsvb_ref.fa \
  --output-tsv results/candidate_sites.tsv \
  --window 80 \
  --step 10

python scripts/score_conservation.py \
  --alignment-fasta results/rsvb_clean.aln.fa \
  --ref-id RSVB_REF \
  --candidates-tsv results/candidate_sites.tsv \
  --output-tsv results/candidates.conservation.tsv

python scripts/score_accessibility.py \
  --candidates-tsv results/candidates.conservation.tsv \
  --output-tsv results/candidates.access.tsv

python scripts/score_specificity.py \
  --candidates-tsv results/candidates.access.tsv \
  --output-tsv results/candidates.scored.tsv

python scripts/select_panel.py \
  --scored-tsv results/candidates.scored.tsv \
  --ranked-output results/rsvb_panel_candidates_ranked.tsv \
  --panel-output results/rsvb_panel_final.tsv \
  --max-sites 12
```

## Notes

- `score_accessibility.py` currently supports an external per-position unpaired-probability profile, or a placeholder heuristic when none is supplied.
- `score_specificity.py` supports external summarized off-target hits; otherwise it fills placeholder values.
- `select_panel.py` currently uses a greedy set-cover style selection with overlap and spacing penalties.
- For production use, the next upgrades should be:
  - real RNAplfold/NUPACK integration
  - real BLAST/MFEprimer integration
  - explicit cross-dimer scoring between selected oligos
  - subtype/time-split validation

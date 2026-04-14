### 1. Robustness score

This asks:

“Will this primer still work across sequence diversity in the alignment?”

It comes from score_conservation.py and mainly reflects how conserved the target site is across all aligned genomes.

For RT primers, it especially cares about:

overall mismatch tolerance across the whole 25 nt site
conservation near the primer 3′-critical region
exactness of the terminal base match

Typical components include:

cov_0mm
cov_1mm
cov_2mm
cov_3p_0mm
cov_3p_1mm
cov_terminal_match
Interpretation
high robustness = this target region is conserved across the alignment, especially where the primer 3′ end matters
low robustness = this primer may fail on many strains / variants

### 2. Accessibility score

This asks:

“Is the target RNA region physically exposed enough for the primer to bind?”

It comes from score_accessibility.py, usually using RNAplfold.

This is a property of the target binding region, not the primer itself.

For RT primers, it gives extra weight to the part of the target where the primer 3′ end initiates binding.

Typical components include:

access_mean
access_min
access_3p_terminal
access_3p_near_mean
Interpretation
high accessibility = the target region is often unpaired / exposed
low accessibility = the target may be buried in RNA secondary structure, making priming harder


### 3. Thermo score

This asks:

“Does the actual RT primer oligo have reasonable thermodynamic properties?”

This is now primer-centric, not just target-region-centric.

It uses properties of rt_primer_seq, such as:

primer GC fraction
primer Tm
3′ clamp strength
whether the terminal base is G/C
whether the 3′ end has problematic homopolymer behavior

Typical fields include:

rt_primer_gc
rt_primer_tm
rt_primer_terminal_base
rt_primer_terminal_is_gc
rt_primer_gc_3p5
Interpretation
high thermo score = Tm is in the desired range and the 3′ end looks good for priming
low thermo score = primer may be too weak, too strong, or have a poor 3′ end

### 4. Specificity score

This asks:

“Will this RT primer also bind or initiate on the wrong background sequences?”

It comes from score_specificity.py.

For RT-primer mode, it should be based on the actual primer sequence (rt_primer_seq), because that is the oligo that can misprime on background.

It usually considers:

total human off-target hits
3′-anchored off-target hits
near-perfect off-target hits
optionally virus-background hits if you include them later

Typical fields include:

human_offtarget_hits
human_offtarget_anchored_hits
human_offtarget_nearperfect_hits
specificity_score
Interpretation
high specificity = little evidence this primer will bind or initiate on human/background templates
low specificity = more risk of non-specific priming

For RT primers, anchored 3′ off-targets are especially bad.

### 5. Synthesis score

This asks:

“Is this primer sequence easy and safe to synthesize and use as an oligo?”

This is more of a sequence hygiene / manufacturability score.

It penalizes things like:

long homopolymers
low-complexity sequence
repetitive composition
problematic 3′ sequence simplicity

Typical fields include:

max_homopolymer
three_prime_max_homopolymer
low_complexity_fraction
Interpretation
high synthesis score = sequence is cleaner and less troublesome as an oligo
low synthesis score = sequence may be repetitive, unstable, or harder to work with




*** Simple summary of each score ***
robustness = how well the target site survives viral variation
accessibility = how open the target RNA site is
thermo = whether the actual primer has good Tm / GC / 3′ properties
specificity = how unlikely the primer is to misprime on background
synthesis = whether the oligo sequence itself is clean and manufacturable

site_score =
0.28 * robustness
+ 0.24 * accessibility
+ 0.18 * thermo
+ 0.22 * specificity
+ 0.08 * synthesis
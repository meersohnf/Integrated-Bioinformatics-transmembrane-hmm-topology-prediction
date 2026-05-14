# Transmembrane Protein Topology Prediction with a Hidden Markov Model

A Python implementation of a discrete Hidden Markov Model (HMM) for predicting
whether each amino acid in a membrane protein is inside the cell, outside the
cell, or spanning the lipid bilayer.

This project was originally made for my final graduate integrated bioinformatics
assignment. It was later refactored into a cleaner code sample
to show HMM implementation, dynamic programming, log-space probability handling,
and biological interpretation of model output.

Some support files were provided as part of the coursework environment. The main
focus of this repository is the HMM implementation and evaluation workflow.

## Background

Membrane proteins account for a large fraction of proteins in most genomes and
are important in transport, signaling, energy production, and drug targeting.
Knowing the topology of a transmembrane protein helps explain which parts of the
protein cross the membrane and which regions face the inside or outside of the
cell.


## Approach

### States

The model uses four hidden topology states:

| State | Meaning |
|---|---|
| `I` | Inside the cell |
| `O` | Outside the cell |
| `M` | Membrane-spanning, inside to outside |
| `N` | Membrane-spanning, outside to inside |

Separating membrane crossings into `M` and `N` helps enforce a biologically
important constraint: a predicted transmembrane segment should enter and exit
on opposite sides of the bilayer. Without this distinction, the model can
generate topologically impossible paths.

### Model Source and Evaluation

Transition and emission probabilities were derived from annotated non-*E. coli*
membrane proteins. The *E. coli* sequences are held out from model construction
and used only for evaluation.

The learned probabilities were then manually adjusted to zero out forbidden
transitions, such as direct `I -> O` movement or invalid membrane-state switches.
The probability mass from those forbidden transitions was redistributed to the
topologically valid alternative.

### Decoding

Two decoding strategies are applied to each evaluation sequence:

- **Viterbi decoding** finds the single globally most probable state path using
  dynamic programming.
- **Posterior decoding** selects the most probable state at each position by
  combining the forward and backward algorithms.

These two methods can produce different paths because Viterbi optimizes the
whole path, while posterior decoding optimizes each position independently.

### Null Comparison

Each sequence is also scored against an all-inside baseline path. The log-odds
ratio between the Viterbi path and this baseline measures how much more
consistent the sequence is with a transmembrane topology than with a fully
intracellular protein.

### Numerical Stability

Protein sequences can be hundreds of residues long. Multiplying many small
transition and emission probabilities can cause floating-point underflow.

This project uses a `LogFloat` class to keep probability calculations in log
space while still allowing readable arithmetic operations in the HMM methods.

## Files

| File | Purpose |
|---|---|
| `transmembrane_hmm.py` | Main HMM implementation and analysis script |
| `FastA_V2.py` | FASTA parser |
| `log_float.py` | Log-space arithmetic class |
| `160_membrane_prots.txt` | Annotated protein sequences and topology paths |

## Usage

Place all required files in the same folder, then run:

```bash
python transmembrane_hmm.py
```

No external dependencies beyond the Python standard library are required.

## Output

For each held-out *E. coli* sequence, the script prints the known topology path,
Viterbi-decoded path, posterior-decoded path, baseline all-inside path, and the
log-odds ratio against that baseline.

Example output:

```text
>ATPL_ECOLI
Sequence length: 79 amino acids
Actual:    SOOOOOOOOOOMMMMMMMMMMMMMMMMMMMMMIIIIIIIIIIIIIIIIIIIIIMMMMMMMMMMMMMMMMMMMMMOOOOOO
Viterbi:   SOOOOOOONNNNNNNNNNNNNNNNNNNNNNNNNNIIIIIIIIIIIIIIIIIIIMMMMMMMMMMMMMMMMMMMMMMMMMMM
Posterior: SOOOOOOONNNNNNNNNNNNNNNNNNNNNNNNNIIIIIIIIIIIIIIIIIIIIMMMMMMMMMMMMMMMMMMMMMMMMMMM
Baseline:  SIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII

  Viterbi log-probability:           -220.2458
  Baseline log-probability:          -233.4767
  Log-odds ratio (Viterbi/baseline): 13.2308
  Odds ratio (real space):           557278.00
```

This example shows that the Viterbi-decoded membrane topology explains the
sequence much better than the all-inside baseline path.

## Interpretation

The final model improves one of the main issues from the simpler one-membrane
state model: membrane segments are less likely to enter and exit on the same
side of the bilayer.

However, this is still a simplified educational model. It can identify
membrane-like regions, but it may still struggle with:

- choosing the correct inside/outside orientation
- identifying exact membrane boundaries
- producing clean posterior paths
- modeling realistic transmembrane helix lengths

More advanced tools, such as TMHMM-style models, use richer state structures
with helix cores, helix caps, loop states, and more detailed validation.

## Author

Francisco Meersohn

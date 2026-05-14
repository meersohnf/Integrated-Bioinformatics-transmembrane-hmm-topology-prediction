# Transmembrane Protein Topology Prediction with a Hidden Markov Model

This repository contains a Python implementation of a discrete Hidden Markov Model (HMM) for predicting transmembrane protein topology from amino acid sequences.

The model predicts whether each residue is:

- `I`: inside / cytoplasmic
- `O`: outside / non-cytoplasmic
- `M`: membrane-spanning, inside to outside
- `N`: membrane-spanning, outside to inside

The project was originally developed for a graduate bioinformatics assignment and later cleaned into a more readable code sample. The goal is to demonstrate HMM implementation, dynamic programming, log-space probability handling, and biological interpretation of model output.

## Project Overview

The input file contains annotated membrane protein sequences with known topology paths. The model parameters used in the final script were derived from annotated non-*E. coli* membrane proteins and then manually adjusted to enforce more biologically realistic topology.

Held-out *E. coli* sequences are used for evaluation only.

The final model uses two membrane-spanning states, `M` and `N`, instead of a single membrane state. This prevents the model from predicting membrane segments that enter and exit on the same side of the bilayer.

## Main Features

- Parses pseudo-FASTA sequence/path records
- Estimates transition and emission probabilities from labeled data
- Implements Viterbi decoding
- Implements forward and backward algorithms
- Implements posterior decoding
- Scores known or proposed hidden paths
- Uses log-space arithmetic to avoid numerical underflow
- Compares predicted transmembrane paths against an all-inside baseline

## Repository Contents

| File | Purpose |
|---|---|
| `transmembrane_hmm.py` | Main HMM implementation |
| `FastA_V2.py` | FASTA-style parser |
| `log_float.py` | Log-space probability helper |
| `160_membrane_prots.txt` | Annotated protein sequences |

## Methods Implemented

The `HMM` class includes the following methods:

### `learn(seqs, paths)`

Estimates transition and emission probabilities from paired observed sequences and known hidden state paths.

### `evaluate(path)`

Computes the joint probability of an observed sequence and a specified hidden state path.

### `forward()`

Computes the total probability of the observed sequence summed over all possible hidden paths.

### `backward()`

Computes the same total sequence probability using the backward dynamic programming formulation.

### `posterior()`

Uses the forward and backward tables to select the most probable state at each sequence position.

### `viterbi()`

Finds the single most probable hidden state path for the observed sequence.

### `get_sequences_and_paths(filename)`

Reads the input file and separates each record into an annotation, amino acid sequence, and known topology path.

## Why Log-Space Arithmetic Is Used

Protein sequences can be hundreds of residues long. Multiplying many small transition and emission probabilities can quickly underflow in standard floating-point arithmetic.

This project uses the provided `LogFloat` class to keep probability calculations numerically stable. `LogFloat` stores values internally in log space while still supporting readable arithmetic operations such as multiplication, division, addition, and comparison.

## Model Design

The final topology model uses separate membrane-crossing states:

- `M`: inside → membrane → outside
- `N`: outside → membrane → inside

Forbidden transitions are assigned probability zero. For example, the model does not allow a direct `I → O` transition or a membrane segment that exits on the same side it entered.

This adjustment makes the predicted paths more biologically realistic than the earlier one-membrane-state model.

## Usage

Place all required files in the same folder:

```text
transmembrane_hmm.py
FastA_V2.py
log_float.py
160_membrane_prots.txt
```

Then run:

```bash
python transmembrane_hmm.py
```

The script uses relative file paths, so no hard-coded system paths are required.

## Example Output

The script evaluates held-out *E. coli* sequences. For each sequence, it reports the known topology path, Viterbi-decoded path, posterior-decoded path, all-inside baseline path, and log-odds comparison.

```text
>EXBD_ECOLI
Sequence length: 141 amino acids
Actual:    SIIIIIIIIIIIIIMMMMMMMMMMMMMMMMMMMMMOOOOOOOO...
Viterbi:   SOOOOOOOOOOOOONNNNNNNNNNNNNNNNNNNNIIIIIIII...
Posterior: SOOOOOOOOOOOOONNNNNNNNNNNNNIIIIIIIIIIIIIII...
Baseline:  SIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII...

  Viterbi log-probability:           -412.4743
  Baseline log-probability:          -416.5502
  Log-odds ratio (Viterbi/baseline): 4.0759
  Odds ratio (real space):           59.14
```

A positive log-odds ratio means the Viterbi-decoded transmembrane path explains the sequence better than forcing the whole protein to be intracellular.

## Interpretation

The final model improves one major issue from the simpler model: membrane segments are less likely to enter and exit on the same side of the bilayer.

However, this is still a simplified educational model. It can identify membrane-like regions, but it may still struggle with:

- choosing the correct inside/outside orientation
- identifying exact membrane boundaries
- producing clean posterior paths
- modeling realistic transmembrane helix lengths

More advanced tools, such as TMHMM-style models, use richer state structures with helix cores, helix caps, loop states, and more detailed validation.

## Notes

This project was developed as part of a graduate bioinformatics HMM assignment and refactored into a cleaner code sample. Some support files were provided as part of the coursework environment.

## Requirements

No external Python packages are required.

The code uses only the Python standard library plus local support files included in the repository.

## Author

Francisco Meersohn

"""
Hidden Markov Model for Transmembrane Domain Prediction
Author: Francisco Meersohn

Predicts whether each amino acid position in a membrane protein is inside the
cell (I), outside the cell (O), or spanning the bilayer (M or N). Two membrane
states enforce biologically valid topology: M for inside->outside crossings,
N for outside->inside crossings.

Transition and emission probabilities were derived from annotated non-E. coli
membrane proteins, then manually adjusted to forbid topologically impossible
transitions (e.g. direct I->O, same-side membrane exits). E. coli sequences
are held out for evaluation only.

Each evaluation sequence is decoded via Viterbi and posterior algorithms, then
scored against an all-inside null path as a log-odds ratio.

Dependencies: FastA_V2.py, log_float.py, 160_membrane_prots.txt
"""

import math
import re

from FastA_V2 import FastA
from log_float import LogFloat


def normalize_rows(probability_dict):
    """Normalize each row of a nested probability dict so values sum to 1."""
    normalized = {}
    for state, row in probability_dict.items():
        total = sum(row.values())
        if total == 0:
            raise ValueError(f"Probability row for state '{state}' sums to zero.")
        normalized[state] = {symbol: value / total for symbol, value in row.items()}
    return normalized


class HMM:

    def __init__(self, transitions=None, emissions=None, sequence=None):
        self.sequence          = sequence
        self.sequence_len      = len(sequence) if sequence else 0
        self.observed_sequence = sequence.upper() if sequence else ""

        self.transitions = transitions or {
            "S": {"F": 0.5,  "L": 0.5},
            "F": {"F": 0.95, "L": 0.05},
            "L": {"L": 0.90, "F": 0.10},
        }
        self.emissions = emissions or {
            "S": {"_": 1},
            "F": {str(i): 1 / 6 for i in range(1, 7)},
            "L": {"1": 0.1, "2": 0.1, "3": 0.1, "4": 0.1, "5": 0.1, "6": 0.5},
        }

        self.forward_table   = []
        self.backward_table  = []
        self.viterbi_table   = []
        self.posterior_table = []

    def learn(self, seqs, paths):
        """Estimate transition and emission probabilities from labelled sequences."""
        transitions = {}
        emissions   = {}

        for seq, path in zip(seqs, paths):
            for position in range(1, len(seq)):
                old_state = path[position - 1]
                new_state = path[position]
                emission  = seq[position]

                if old_state not in transitions:
                    transitions[old_state] = {}
                transitions[old_state][new_state] = transitions[old_state].get(new_state, 0) + 1

                if new_state not in emissions:
                    emissions[new_state] = {}
                emissions[new_state][emission] = emissions[new_state].get(emission, 0) + 1

        for state, next_states in transitions.items():
            total = sum(next_states.values())
            transitions[state] = {k: v / total for k, v in next_states.items()}

        for state, symbol_counts in emissions.items():
            total = sum(symbol_counts.values())
            emissions[state] = {k: v / total for k, v in symbol_counts.items()}

        emissions["S"] = {"_": 1}

        self.transitions = transitions
        self.emissions   = emissions
        return self.transitions, self.emissions

    def get_sequences_and_paths(self, filename):
        """Parse a pseudo-FASTA file with sequences and topology paths separated by '#'."""
        annotations, sequences, paths = [], [], []

        for annotation, combined in FastA(filename):
            combined = re.sub(r"\s+", "", combined)
            sequence, actual_path = combined.split("#")
            sequence    = "_" + sequence.upper()
            actual_path = "S" + actual_path.upper()

            if len(sequence) != len(actual_path):
                raise ValueError(
                    f"Length mismatch for {annotation}: "
                    f"sequence={len(sequence)}, path={len(actual_path)}"
                )
            annotations.append(annotation.upper())
            sequences.append(sequence)
            paths.append(actual_path)

        return annotations, sequences, paths

    def evaluate(self, path):
        """Compute joint probability P(sequence, path) for a given path."""
        joint_probability = LogFloat(1)

        for position in range(self.sequence_len):
            current_state   = path[position]
            observed_symbol = self.observed_sequence[position]

            if observed_symbol not in self.emissions[current_state]:
                return LogFloat(0)

            joint_probability *= self.emissions[current_state][observed_symbol]

            if position > 0:
                previous_state = path[position - 1]
                if current_state not in self.transitions[previous_state]:
                    return LogFloat(0)
                joint_probability *= self.transitions[previous_state][current_state]

        return joint_probability

    def forward(self):
        """Forward algorithm: compute P(sequence) summed over all paths."""
        self.forward_table = []

        first_column = {state: LogFloat(0) for state in self.transitions}
        first_column["S"] = LogFloat(1)
        self.forward_table.append(first_column)

        for position in range(1, self.sequence_len + 1):
            observed_symbol = self.observed_sequence[position - 1]
            self.forward_table.append({})

            for current_state in self.transitions:
                if observed_symbol not in self.emissions[current_state]:
                    self.forward_table[position][current_state] = LogFloat(0)
                    continue

                if position == 1:
                    self.forward_table[position][current_state] = (
                        self.forward_table[0][current_state]
                        * self.emissions[current_state][observed_symbol]
                    )
                else:
                    total = LogFloat(0)
                    for previous_state in self.transitions:
                        if current_state in self.transitions[previous_state]:
                            total += (
                                self.forward_table[position - 1][previous_state]
                                * self.transitions[previous_state][current_state]
                            )
                    self.forward_table[position][current_state] = (
                        total * self.emissions[current_state][observed_symbol]
                    )

        return sum(self.forward_table[self.sequence_len].values())

    def backward(self):
        """Backward algorithm: compute P(sequence) summed over all paths."""
        self.backward_table = [{} for _ in range(self.sequence_len + 1)]

        for state in self.transitions:
            self.backward_table[self.sequence_len][state] = LogFloat(1)

        for position in range(self.sequence_len - 1, 0, -1):
            next_observed = self.observed_sequence[position]
            for current_state in self.transitions:
                total = LogFloat(0)
                for next_state in self.transitions:
                    if (next_state in self.transitions[current_state]
                            and next_observed in self.emissions[next_state]):
                        total += (
                            self.transitions[current_state][next_state]
                            * self.emissions[next_state][next_observed]
                            * self.backward_table[position + 1][next_state]
                        )
                self.backward_table[position][current_state] = total

        start_symbol = self.observed_sequence[0]
        for state in self.transitions:
            if start_symbol in self.emissions[state]:
                self.backward_table[0][state] = (
                    self.emissions[state][start_symbol]
                    * self.backward_table[1][state]
                )
            else:
                self.backward_table[0][state] = LogFloat(0)

        return sum(self.backward_table[0].values())

    def posterior(self):
        """Posterior decoding: per-position most probable state using forward * backward."""
        probability_of_sequence = self.forward()
        self.backward()

        self.posterior_table = []
        posterior_path       = []

        for position in range(self.sequence_len):
            self.posterior_table.append({})
            table_position   = position + 1
            best_state       = None
            best_probability = LogFloat(0)

            for state in self.transitions:
                if probability_of_sequence.evaluate_log() is None:
                    posterior_prob = LogFloat(0)
                else:
                    posterior_prob = (
                        self.forward_table[table_position][state]
                        * self.backward_table[table_position][state]
                        / probability_of_sequence
                    )
                self.posterior_table[position][state] = posterior_prob

                if posterior_prob > best_probability:
                    best_probability = posterior_prob
                    best_state       = state

            posterior_path.append(best_state)

        return posterior_path

    def viterbi(self):
        """Viterbi decoding: most probable single path via dynamic programming."""
        vtable         = [{state: LogFloat(0) for state in self.transitions}]
        possible_paths = {state: [state] for state in self.transitions}
        vtable[0]["S"] = LogFloat(1)

        for position in range(1, self.sequence_len):
            vtable.append({})
            new_paths       = {}
            observed_symbol = self.observed_sequence[position]

            for current_state in self.transitions:
                candidates = []
                for previous_state in self.transitions:
                    if current_state not in self.transitions[previous_state]:
                        continue
                    if self.transitions[previous_state][current_state] == 0:
                        continue  # skip forbidden transitions
                    prob = (
                        vtable[position - 1][previous_state]
                        * self.transitions[previous_state][current_state]
                        * self.emissions[current_state][observed_symbol]
                    )
                    candidates.append((prob, previous_state))

                if candidates:
                    best_prob, best_prev = max(candidates)
                    vtable[position][current_state] = best_prob
                    new_paths[current_state]        = possible_paths[best_prev] + [current_state]
                else:
                    vtable[position][current_state] = LogFloat(0)
                    new_paths[current_state]        = possible_paths[current_state]

            possible_paths     = new_paths
            self.viterbi_table = vtable

        best_prob, best_path = max(
            (vtable[self.sequence_len - 1][state], possible_paths[state])
            for state in self.transitions
        )
        return best_prob, best_path


# Hand-tuned two-state membrane topology model (2h)
# Derived from non-E. coli sequences in 160_membrane_prots.txt, then edited to enforce:
#   M: I -> membrane -> O
#   N: O -> membrane -> I
# Forbidden transitions set to 0; mass redistributed to the one valid alternative.

learned_transitions_2h = {
    "S": {"O": 0.5714285714285714, "I": 0.42857142857142855, "M": 0,                    "N": 0},
    "O": {"O": 0.9910156722130454, "M": 0,                   "N": 0.008984327786954596, "I": 0},
    "M": {"M": 0.9559774671395785, "I": 0,                   "O": 0.04402253286042145,  "N": 0},
    "I": {"I": 0.9813914344492238, "M": 0.018608565550776137,"N": 0,                    "O": 0},
    "N": {"N": 0.9559774671395785, "I": 0.04402253286042145, "O": 0,                    "M": 0},
}

learned_emissions_2h = {
    "O": {
        "M": 0.0214, "K": 0.0481, "R": 0.0502, "L": 0.0892, "V": 0.0623,
        "C": 0.0282, "D": 0.0560, "W": 0.0199, "A": 0.0622, "S": 0.0722,
        "N": 0.0512, "T": 0.0642, "E": 0.0607, "Y": 0.0382, "G": 0.0673,
        "Q": 0.0421, "F": 0.0413, "I": 0.0457, "P": 0.0556, "H": 0.0240,
    },
    "M": {
        "Y": 0.0442, "F": 0.0922, "V": 0.1189, "I": 0.1193, "Q": 0.0100,
        "T": 0.0512, "L": 0.1663, "P": 0.0275, "C": 0.0211, "M": 0.0379,
        "S": 0.0558, "W": 0.0282, "A": 0.1038, "R": 0.0058, "G": 0.0716,
        "D": 0.0074, "E": 0.0083, "N": 0.0174, "H": 0.0081, "K": 0.0048,
    },
    "N": {  # same amino-acid bias as M; direction enforced by transitions
        "Y": 0.0442, "F": 0.0922, "V": 0.1189, "I": 0.1193, "Q": 0.0100,
        "T": 0.0512, "L": 0.1663, "P": 0.0275, "C": 0.0211, "M": 0.0379,
        "S": 0.0558, "W": 0.0282, "A": 0.1038, "R": 0.0058, "G": 0.0716,
        "D": 0.0074, "E": 0.0083, "N": 0.0174, "H": 0.0081, "K": 0.0048,
    },
    "I": {
        "L": 0.0794, "N": 0.0407, "R": 0.0825, "E": 0.0673, "S": 0.0789,
        "Y": 0.0271, "F": 0.0317, "T": 0.0539, "K": 0.0729, "G": 0.0677,
        "A": 0.0721, "W": 0.0115, "D": 0.0510, "V": 0.0565, "P": 0.0588,
        "I": 0.0397, "Q": 0.0402, "M": 0.0254, "H": 0.0229, "C": 0.0199,
    },
    "S": {"_": 1},
}


if __name__ == "__main__":

    reader = HMM(sequence="_")
    annotations, sequences, paths = reader.get_sequences_and_paths("160_membrane_prots.txt")

    n_eval  = sum(1 for ann in annotations if "ECOLI" in ann)
    n_model = len(annotations) - n_eval
    print(f"Total sequences:                      {len(annotations)}")
    print(f"Model source sequences (non-E. coli): {n_model}")
    print(f"Evaluation sequences (E. coli):       {n_eval}\n")

    learned_transitions = learned_transitions_2h
    learned_emissions   = normalize_rows(learned_emissions_2h)

    # Verify transition rows sum to 1 before running
    for state, row in learned_transitions.items():
        row_sum = sum(row.values())
        if not math.isclose(row_sum, 1.0, rel_tol=1e-6):
            raise ValueError(f"Transition row for '{state}' sums to {row_sum:.6f}")

    for annotation, sequence, actual_path in zip(annotations, sequences, paths):

        if "ECOLI" not in annotation:
            continue

        model = HMM(
            transitions=learned_transitions,
            emissions=learned_emissions,
            sequence=sequence,
        )

        viterbi_probability, viterbi_path = model.viterbi()
        posterior_path = model.posterior()

        baseline_path        = "S" + "I" * (len(sequence) - 1)
        baseline_probability = model.evaluate(baseline_path)

        viterbi_log    = viterbi_probability.evaluate_log()
        baseline_log   = baseline_probability.evaluate_log()
        log_odds_ratio = viterbi_log - baseline_log

        try:
            odds_ratio = math.exp(log_odds_ratio)
        except OverflowError:
            odds_ratio = float("inf")

        print(f">{annotation}")
        print(f"Sequence length: {len(sequence) - 1} amino acids")
        print(f"Actual:    {actual_path}")
        print(f"Viterbi:   {''.join(viterbi_path)}")
        print(f"Posterior: {''.join(posterior_path)}")
        print(f"Baseline:  {baseline_path}")
        print()
        print(f"  Viterbi log-probability:           {viterbi_log:.4f}")
        print(f"  Baseline log-probability:          {baseline_log:.4f}")
        print(f"  Log-odds ratio (Viterbi/baseline): {log_odds_ratio:.4f}")
        print(f"  Odds ratio (real space):           {odds_ratio:.2f}")
        print()

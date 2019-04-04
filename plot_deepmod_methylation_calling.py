#!/usr/bin/env python
"""Plot data for deepmod methylation calling"""
########################################################################
# File: plot_deepmod_methylation_calling.py
#  executable: plot_deepmod_methylation_calling.py
#
# Author: Andrew Bailey
# History: Created 04/02/19
########################################################################


import pandas as pd
import numpy as np
import os
import matplotlib.pyplot as plt


from py3helpers.classification import ClassificationMetrics
from py3helpers.utils import list_dir, time_it
from py3helpers.seq_tools import ReverseComplement


class CustomAmbiguityPositions(object):
    def __init__(self, ambig_filepath):
        """Deal with ambiguous positions from a tsv ambiguity position file with the format of
        contig  position            strand  change_from change_to
        'name'  0 indexed position   +/-    C           E


        :param ambig_filepath: path to ambiguity position file"""

        self.ambig_df = self.parseAmbiguityFile(ambig_filepath)

    @staticmethod
    def parseAmbiguityFile(ambig_filepath):
        """Parses a 'ambiguity position file' that should have the format:
            contig  position    strand  change_from change_to

        :param ambig_filepath: path to ambiguity position file
        """
        return pd.read_csv(ambig_filepath, sep='\t',
                           usecols=(0, 1, 2, 3, 4),
                           names=["contig", "position", "strand", "change_from", "change_to"],
                           dtype={"contig": np.str,
                                  "position": np.int,
                                  "strand": np.str,
                                  "change_from": np.str,
                                  "change_to": np.str})

    def getForwardSequence(self, contig, raw_sequence):
        """Edit 'raw_sequence' given a ambiguity positions file. Assumes raw_sequence is forward direction( 5'-3')
        :param contig: which contig the sequence belongs (aka header)
        :param raw_sequence: raw nucleotide sequence
        :return: edited nucleotide sequence
        """
        return self._get_substituted_sequence(contig, raw_sequence, "+")

    def getBackwardSequence(self, contig, raw_sequence):
        """Edit 'raw_sequence' given a ambiguity positions file, Assumes raw_sequence is forward direction( 5'-3')
        :param contig: which contig the sequence belongs (aka header)
        :param raw_sequence: raw nucleotide sequence
        :return: edited nucleotide sequence
        """
        rc = ReverseComplement()
        raw_sequence = rc.complement(raw_sequence)
        return self._get_substituted_sequence(contig, raw_sequence, "-")

    def _get_substituted_sequence(self, contig, raw_sequence, strand):
        """Change the given raw nucleotide sequence using the edits defined in the positions file

        :param contig: name of contig to find
        :param raw_sequence: nucleotide sequence (note: this is note edited in this function)
        :param strand: '+' or '-' to indicate strand
        """
        contif_df = self._get_contig_positions(contig, strand)
        raw_sequence = list(raw_sequence)
        for _, row in contif_df.iterrows():
            if raw_sequence[row["position"]] != row["change_from"]:
                raise RuntimeError(
                    "[CustomAmbiguityPositions._get_substituted_sequence]Illegal substitution requesting "
                    "change from %s to %s, row: %s" % (raw_sequence[row["position"]], row["change_to"], row))
            raw_sequence[row["position"]] = row["change_to"]
        return "".join(raw_sequence)

    def _get_contig_positions(self, contig, strand):
        """Get all unique locations within the positions file

        :param contig: name of contig to find
        :param strand: '+' or '-' to indicate strand
        """
        df = self.ambig_df.loc[
            (self.ambig_df["contig"] == contig) & (self.ambig_df["strand"] == strand)].drop_duplicates()
        assert len(df['position']) == len(set(df['position'])), "Multiple different changes for a single position. {}" \
            .format(df['position'])
        return df


def parse_deepmod_bed(deepmod_bed_path):
    """Parse the summary chromosome bed file output from deepmod

    eg: Chromosome 2 3 C 8 + 2 3 0,0,0 8 0 0
    """

    return pd.read_csv(deepmod_bed_path, sep=" ", header=None,
                       usecols=(0, 1, 3, 5, 9, 10, 11),
                       names=["contig", "start_position", "base", "strand", "n_reads", "modification_percentage",
                              "n_mod_calls"],
                       dtype={"contig": np.str,
                              "start_position": np.int,
                              "base": np.str,
                              "strand": np.str,
                              "n_reads": np.int,
                              "modification_percentage": np.int,
                              "n_mod_calls": np.int})


def aggregate_deepmod_data(deepmod_output_dir):
    deepmod_beds = list_dir(deepmod_output_dir, ext="bed")
    deepmod_bed_data = []
    for bed in deepmod_beds:
        data = parse_deepmod_bed(bed)
        data["E"] = (data["modification_percentage"] / 100)
        data["C"] = 1 - (data["modification_percentage"] / 100)
        deepmod_bed_data.append(data)

    return pd.concat(deepmod_bed_data)


def print_confusion_matrix(tp, fp, fn, tn):
    precision = tp / (tp + fp)
    false_discovery_rate = fp / (tp + fp)

    false_omission_rate = fn / (tn + fn)
    negative_predictive_value = tn / (tn + fn)

    true_positive_rate_recall = tp / (tp + fn)
    false_negative_rate = fn / (tp + fn)

    false_positive_rate = fp / (fp + tn)
    true_negative_rate_specificity = tn / (tn + fp)

    positive_likelihood_ratio = true_positive_rate_recall / false_positive_rate
    negative_likelihood_ratio = false_negative_rate / true_negative_rate_specificity

    diagnostic_odds_ratio = positive_likelihood_ratio / negative_likelihood_ratio

    f1_score = 2 * ((precision * true_positive_rate_recall) / (precision + true_positive_rate_recall))

    return (np.asarray(
        [[tp, fp, precision, false_discovery_rate],
         [fn, tn, false_omission_rate, negative_predictive_value],
         [true_positive_rate_recall, false_positive_rate, positive_likelihood_ratio, diagnostic_odds_ratio],
         [false_negative_rate, true_negative_rate_specificity, negative_likelihood_ratio, f1_score]]))


def plot_confusion_matrix(tp, fp, fn, tn, classes=("E", "C"), title="Confusion Matrix",
                          output_path=None, normalize=False):
    """Plot the confusion matrix with the information of each box explained
    :param classes: classes to label axis
    :param normalize: option to normalize output of confusion matrix
    :param tp: true positives
    :param fp: false positives
    :param fn: false negatives
    :param tn: true negatives
    :param title: title of the plot
    :param output_path: place to save plot
    """
    data = np.asarray([[tp, fp], [fn, tn]])
    # total = sum([tp, fp, fn, tn])
    # precision = tp / (tp + fp)
    # false_discovery_rate = fp / (tp + fp)
    #
    # false_omission_rate = fn / (tn + fn)
    # negative_predictive_value = tn / (tn + fn)
    #
    # true_positive_rate_recall = tp / (tp + fn)
    # false_negative_rate = fn / (tp + fn)
    #
    # false_positive_rate = fp / (fp + tn)
    # true_negative_rate_specificity = tn / (tn + fp)
    #
    # positive_likelihood_ratio = true_positive_rate_recall / false_positive_rate
    # negative_likelihood_ratio = false_negative_rate / true_negative_rate_specificity
    #
    # diagnostic_odds_ratio = positive_likelihood_ratio / negative_likelihood_ratio
    #
    # f1_score = 2 * ((precision * true_positive_rate_recall) / (precision + true_positive_rate_recall))

    if not title:
        if normalize:
            title = 'Normalized confusion matrix'
        else:
            title = 'Confusion matrix, without normalization'

    # Compute confusion matrix
    n_data = data.astype('float') / data.sum(axis=1)[:, np.newaxis]
    fig, ax = plt.subplots()
    im = ax.imshow(data, interpolation='nearest', cmap=plt.cm.Blues)
    ax.figure.colorbar(im, ax=ax)
    # We want to show all ticks...
    ax.set(xticks=np.arange(data.shape[1]),
           yticks=np.arange(data.shape[0]),
           xticklabels=classes, yticklabels=classes,
           title=title,
           ylabel='True label',
           xlabel='Predicted label')

    # Rotate the tick labels and set their alignment.
    plt.setp(ax.get_xticklabels(), rotation=45, ha="right",
             rotation_mode="anchor")

    # Loop over data dimensions and create text annotations.
    fmt = '.2f' if normalize else 'd'
    thresh = data.max() / 2.
    for i in range(data.shape[0]):
        for j in range(data.shape[1]):
            ax.text(j, i, format(data[i, j], 'd') + format(n_data[i, j], '.2f'),
                    ha="center", va="center",
                    color="white" if data[i, j] > thresh else "black")

    fig.tight_layout()

    if output_path is not None:
        plt.savefig(output_path)
    else:
        plt.show()
    return True


def main():
    cpg_positions_file = "/Users/andrewbailey/data/references/ecoli/CG_ecoli_k12_mg1655_C_E.positions"
    modified_deepmod_output_dir = "/Users/andrewbailey/CLionProjects/DeepMod/ecoli_pcr_MSssI_R9"
    canonical_deepmod_output_dir = "/Users/andrewbailey/CLionProjects/DeepMod/ecoli_pcr_MSssI_R9"

    output_dir = "/Users/andrewbailey/CLionProjects/modification_detection_pipeline/output_dir/plotting_output/"
    log_file_path = os.path.join(output_dir, "confusion_matrices_file.txt")
    cpg_positions = CustomAmbiguityPositions(cpg_positions_file)

    canonical_data = aggregate_deepmod_data(canonical_deepmod_output_dir)
    canonical_data["E_label"] = 0
    canonical_data["C_label"] = 1

    modified_data = aggregate_deepmod_data(modified_deepmod_output_dir)
    modified_data["E_label"] = 1
    modified_data["C_label"] = 0

    tps = 0
    fps = 0
    tns = 0
    fns = 0
    all_data = []
    with open(log_file_path, "w") as log_file:
        chromosomes = set(modified_data["contig"]) | set(canonical_data["contig"])
        strands = set(modified_data["strand"]) | set(canonical_data["strand"])
        for chromosome in chromosomes:
            for strand in strands:
                # get positions for strand and contig
                sc_positions = cpg_positions.ambig_df.loc[(cpg_positions.ambig_df["strand"] == strand) &
                                                          (cpg_positions.ambig_df["contig"] == chromosome)]
                # pare data sets for specific contig and strand to get positions that are cpgs
                mod_sc_data = modified_data.loc[(modified_data["contig"] == chromosome) &
                                                (modified_data["strand"] == strand)]
                mod_methylation_calls = mod_sc_data.loc[mod_sc_data["start_position"].isin(sc_positions["position"])]
                canon_sc_data = canonical_data.loc[(canonical_data["contig"] == chromosome) &
                                                   (canonical_data["strand"] == strand)]
                canon_methylation_calls = canon_sc_data.loc[
                    canon_sc_data["start_position"].isin(sc_positions["position"])]

                # per site
                n_negative_calls = sum(canon_methylation_calls["n_reads"])
                n_false_negatives = sum(canon_methylation_calls["n_mod_calls"])
                n_true_negatives = n_negative_calls - n_false_negatives

                n_positive_calls = sum(mod_methylation_calls["n_reads"])
                n_true_positives = sum(mod_methylation_calls["n_mod_calls"])
                n_false_positives = n_positive_calls - n_true_positives

                tps += n_true_positives
                fps += n_false_positives
                tns += n_true_negatives
                fns += n_false_negatives
                print("Chromosome {} strand {}:".format(chromosome, strand), file=log_file)
                print("Per-call confusion matrix", file=log_file)
                print(print_confusion_matrix(n_true_positives, n_false_positives, n_false_negatives, n_true_negatives),
                      file=log_file)
                plot_confusion_matrix(n_true_positives, n_false_positives, n_false_negatives, n_true_negatives,
                                      normalize=True,
                                      output_path=os.path.join(output_dir, "per_call_{}_{}_confusion_matrix.png".format(strand, chromosome)),
                                      title="Per call CpG Normalized Confusion Matrix {}{}".format(strand, chromosome))

                # per genomic position
                chr_strand_data = pd.concat([canon_methylation_calls, mod_methylation_calls])
                label_data = chr_strand_data.loc[:, ['C_label', "E_label"]]
                prediction_data = chr_strand_data.loc[:, ['C', "E"]]
                label_data.rename(columns={'C_label': 'C', "E_label": "E"}, inplace=True)
                cmh = ClassificationMetrics(label_data, prediction_data)
                cmh.plot_roc("E", os.path.join(output_dir, "per_genomic_site_{}_{}_roc.png".format(chromosome, strand)))
                cmh.plot_precision_recall("E", os.path.join(output_dir,
                                                            "per_genomic_site_{}_{}_"
                                                            "precision_recall.png".format(chromosome, strand)))
                print("Per-genomic-site confusion matrix", file=log_file)
                print(cmh.confusion_matrix(), file=log_file)

                all_data.append(chr_strand_data)

        print("All Chromosomes both strands:", file=log_file)
        print("Per-call confusion matrix", file=log_file)
        print(print_confusion_matrix(tps, fps, fns, tns),
              file=log_file)
        plot_confusion_matrix(tps, fps, fns, tns,
                              normalize=True,
                              output_path=os.path.join(output_dir,
                                                       "all_calls_confusion_matrix.png"),
                              title="All calls CpG "
                                    "Normalized Confusion Matrix")
        all_data = pd.concat(all_data)
        label_data = all_data.loc[:, ['C_label', "E_label"]]
        prediction_data = all_data.loc[:, ['C', "E"]]
        label_data.rename(columns={'C_label': 'C', "E_label": "E"}, inplace=True)
        cmh = ClassificationMetrics(label_data, prediction_data)
        cmh.plot_roc("E", os.path.join(output_dir, "per_genomic_site_all_chromosomes_roc.png"))
        cmh.plot_precision_recall("E", os.path.join(output_dir,
                                                    "per_genomic_site_all_chromosomes"
                                                    "precision_recall.png"))
        print("Per-genomic-site confusion matrix", file=log_file)
        print(cmh.confusion_matrix(), file=log_file)


if __name__ == '__main__':
    print(time_it(main))

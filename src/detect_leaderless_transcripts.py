#!/usr/bin/env python3

import configparser
import csv
import os
from typing import Dict, List


class LeaderlessTranscriptDetector:

    def __init__(self, reads_path: str, threshold: int):
        """
        Init for LeaderlessTranscriptDetector class.
        For each gene identifier count number of corresponding matched overlapping reads.
        A gene is detected as leaderless transcript if the number of matched reads is below threshold.
        :param mapped_reads_path: str contains path to input tsv file
        :param threshold: maximum number of allowed mapped reads for a gene identifier
        """
        self.reads_path = reads_path
        self.threshold = threshold

    def count_mapped_reads_for_gene(self) -> Dict:
        """
        Count number of mapped overlapping reads for each gene identifier.
        :return: Dict contains gene identifier as key and number of mapped overlapping reads as value
        """
        count_reads_dict = dict()
        with open(self.reads_path, 'r') as in_file:
            tsv_reader = csv.reader(in_file, delimiter='\t')
            next(tsv_reader)  # skip header
            for row in tsv_reader:
                gene_id = row[0]
                if gene_id not in count_reads_dict:
                    count_reads_dict[gene_id] = 0
                else:
                    count_reads_dict[gene_id] += 1
        return count_reads_dict

    def get_leaderless_transcripts(self, count_reads_dict: Dict) -> List:
        """
        Find all leaderless transcripts after reads were mapped to a gene identifier.
        Only choose gene identifiers where the number of mapped reads is below threshold.
        :param count_reads_dict: dict contains gene identifier as key and number of mapped overlapping reads as value
        :return: List contains gene identifiers of leaderless transcripts
        """
        leaderless_transcripts_list = []
        for gene_id, number_reads in count_reads_dict.items():
            if number_reads <= self.threshold:
                leaderless_transcripts_list.append(gene_id)
        return leaderless_transcripts_list


def main():
    # adjust data_path for your system in config.ini
    config = configparser.ConfigParser()
    config.read(os.path.join(os.getcwd().split("src")[0], "config/config.ini"))
    threshold = int(config["DETECTION"]["threshold"])
    debug_mode = eval(config["DEBUG_MODE"]["debug"])

    in_path = os.path.join(os.getcwd().split("src")[0], "output_data/reads_before_first_exon.tsv")
    count_reads_out_path = os.path.join(os.getcwd().split("src")[0], "output_data/number_mapped_reads.tsv")
    leaderless_transcripts_out_path = os.path.join(os.getcwd().split("src")[0],
                                                   f"output_data/leaderless_transcripts_threshold_{threshold}.tsv")

    ltd = LeaderlessTranscriptDetector(reads_path=in_path, threshold=threshold)
    count_reads_dict = ltd.count_mapped_reads_for_gene()
    leaderless_transcripts_list = ltd.get_leaderless_transcripts(count_reads_dict=count_reads_dict)

    # write all genes and number of overlapping reads into file
    with open(count_reads_out_path, 'w',
              newline='') as count_out_file:
        # write gene identifier and number of corresponding mapped reads into file
        count_fieldnames = ["GeneIdentifier", "Number_mapped_reads"]
        count_writer = csv.DictWriter(count_out_file, fieldnames=count_fieldnames, delimiter='\t')
        count_writer.writeheader()
        for identifier, number_reads in count_reads_dict.items():
            if debug_mode:
                print(f"GeneIdentifier: {identifier}, number of mapped reads: {number_reads}")
            count_writer.writerow({"GeneIdentifier": identifier, "Number_mapped_reads": number_reads})

    # write leaderless transcripts for given threshold into files
    with open(leaderless_transcripts_out_path, 'w', newline='') as out_file:
        fieldnames = ["GeneIdentifier"]
        writer = csv.DictWriter(out_file, fieldnames=fieldnames, delimiter='\t')
        writer.writeheader()
        for gene_id in leaderless_transcripts_list:
            if debug_mode:
                print(f"GeneIdentifier: {gene_id}")
            writer.writerow({"GeneIdentifier": gene_id})


if __name__ == "__main__":
    main()

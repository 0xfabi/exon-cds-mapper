#!/usr/bin/env python3

import configparser
import csv
import os
from typing import List


class LeaderlessTranscriptDetector:

    def __init__(self, mapped_reads_path: str, threshold: int):
        """
        Init for LeaderlessTranscriptDetector class. Init with path for a given tsv file containing gene identifiers
        and number of corresponding mapped reads and threshold value which defines a leaderless transcript
        :param mapped_reads_path: str contains path to input tsv file
        :param threshold: maximum number of allowed mapped reads for a gene identifier
        """
        self.path = mapped_reads_path
        self.threshold = threshold

    def get_leaderless_transcripts(self) -> List:
        """
        Find all leaderless transcripts after reads were mapped to a gene identifier.
        Only choose gene identifiers where the number of mapped reads is beyond threshold.
        :return: List contains gene identifiers of leaderless transcripts
        """
        leaderless_transcripts_list = []
        with open(self.path, 'r') as in_file:
            tsv_reader = csv.reader(in_file, delimiter='\t')
            next(tsv_reader)  # skip header
            for row in tsv_reader:
                gene_id = row[0]
                read_count = int(row[1])
                # add gene identifier to list if number of mapped reads is beyond threshold
                if read_count <= self.threshold:
                    leaderless_transcripts_list.append((gene_id, read_count))
        return leaderless_transcripts_list


def main():
    # adjust data_path for your system in config.ini
    config = configparser.ConfigParser()
    config.read(os.path.join(os.getcwd().split("src")[0], "config.ini"))
    threshold = int(config["DETECTION"]["threshold"])
    debug_mode = eval(config["DEBUG_MODE"]["debug"])

    in_path = os.path.join(os.getcwd().split("src")[0], "output_data/number_mapped_reads.tsv")
    out_path = os.path.join(os.getcwd().split("src")[0], "output_data/leaderless_transcripts.tsv")

    leaderless_transcript_detection = LeaderlessTranscriptDetector(mapped_reads_path=in_path, threshold=threshold)
    leaderless_transcripts_list = leaderless_transcript_detection.get_leaderless_transcripts()

    with open(out_path, 'w', newline='') as out_file:
        fieldnames = ["GeneIdentifier"]
        writer = csv.DictWriter(out_file, fieldnames=fieldnames, delimiter='\t')
        writer.writeheader()

        for leaderless_transcript in leaderless_transcripts_list:
            gene_id = leaderless_transcript[0]
            if debug_mode:
                print(f"GeneIdentifier: {gene_id}")
            writer.writerow({"GeneIdentifier": gene_id})


if __name__ == "__main__":
    main()

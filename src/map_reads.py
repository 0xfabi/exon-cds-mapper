#!/usr/bin/env python3

import configparser
import csv
import os
from collections import defaultdict
from typing import Dict, List

import pysam


class MapReads:

    def __init__(self):
        """
        Init for MapReads class.
        """
        pass

    def create_gene_lookup(self, path: str) -> Dict:
        """
        Create a look up table for mapping read to a gene identifier of corresponding chromosome.
        :param path: str contains path to matching exon and cds file
        :return: Dict contains chromosome as key and gene specific values as value tuple.
        """
        gene_matching_dict = defaultdict(list)
        with (open(path, mode='r', newline='')) as infile:
            reader = csv.DictReader(infile, delimiter="\t")
            for row in reader:
                gene_matching_dict[row["Chromosome"]].append(
                    (row["GeneIdentifier"], row["Start"], row["Stop"], row["Strand"]))
        # sort each chromosome group by start value
        for k, v in gene_matching_dict.items():
            gene_matching_dict[k] = sorted(v, key=lambda tup: tup[1])
        return gene_matching_dict

    def get_len_from_cigar(self, cigar: List) -> int:
        """
        Calculate sequence length from a given CIGAR string.
        CIGAR string is already preformated by pysam (e.g. 129M1S is represented as [(0, 129), (4, 1)])
        :param cigar: List contains tuples (operation code, number operations)
        :return: length of sequence as int
        """
        seq_len = 0
        operations = {0: "M", 1: "I", 4: "S", 7: "=", 8: "X"}
        for tup in cigar:
            if tup[0] in operations:
                seq_len += tup[1]
        return seq_len


def main():
    # adjust data_path for your system in config.ini
    config = configparser.ConfigParser()
    config.read(os.path.join(os.getcwd().split("src")[0], "config/config.ini"))
    data_path = config["DATA_PATH"]["path"]
    debug_mode = eval(config["DEBUG_MODE"]["debug"])

    bam_path = os.path.join(data_path, "ath.bam")

    mapping_in_path = os.path.join(os.getcwd().split("src")[0], "output_data/matching_exon_cds.tsv")
    read_out_path = os.path.join(os.getcwd().split("src")[0], "output_data/reads_before_first_exon.tsv")

    read_mapper = MapReads()
    gene_matching_dict = read_mapper.create_gene_lookup(path=mapping_in_path)
    sam_file = pysam.AlignmentFile(bam_path, "rb")

    with open(read_out_path, 'w', newline='') as read_out_file:
        read_fieldnames = ["GeneIdentifier", "Start_exon_1", "Read_QNAME", "Read_START"]
        # write all reads into file which can be mapped to a gene identifier
        read_writer = csv.DictWriter(read_out_file, fieldnames=read_fieldnames, delimiter='\t')
        read_writer.writeheader()
        # iterate over each row in bam file
        for row in sam_file:
            if debug_mode:
                print(row)
            # get necessary information for data point from BAM file
            q_name = row.query_name  # contains unique identifier
            chr_name = row.reference_name  # chromosome name
            start_align = row.reference_start  # starting position of read
            seq_len = read_mapper.get_len_from_cigar(row.cigartuples)  # length of read sequence
            if row.is_reverse:  # get the strand orientation
                strand = "-"
                # stopping position of read should be before start because of reverse strand orientation
                stop_align = start_align - seq_len
            else:
                strand = "+"
                stop_align = start_align + seq_len  # stopping position of read

            # for given chromosome try to find corresponding gene identifier for the given read
            matching_gene_identifier = ""
            for genes in gene_matching_dict[chr_name]:
                if genes[3] == strand:
                    # validate the strand orientation and validate if the read overlaps which exon:1 from gene
                    if strand == "+":
                        exon_start = int(genes[1])
                        exon_stop = int(genes[2])
                        if start_align <= exon_start <= stop_align <= exon_stop:
                            matching_gene_identifier = genes[0]
                            if debug_mode:
                                print(f"Read is before gene {matching_gene_identifier}.")
                            break
                    else:
                        if start_align >= exon_start >= stop_align >= exon_stop:
                            matching_gene_identifier = genes[0]
                            if debug_mode:
                                print(f"Read is before gene {matching_gene_identifier}.")
                            break
                else:
                    continue  # try next gene if strand orientation does not match
            # if read could be matched to exon:1 of a gene write read to file
            if matching_gene_identifier:
                if debug_mode:
                    print(
                        f"For GeneIdentifier {genes[0]} read {q_name} ({start_align}) starts before exon:1 ({genes[2]})")
                read_writer.writerow({"GeneIdentifier": genes[0], "Start_exon_1": genes[2], "Read_QNAME": q_name,
                                      "Read_START": start_align})


if __name__ == "__main__":
    main()

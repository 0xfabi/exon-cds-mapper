#!/usr/bin/env python3

import configparser
import csv
import os
from typing import Dict, List

import pandas as pd
import pysam
from map_first_exon_cds import ExonCdsMapper


class MapReads():
    def __init__(self):
        """
        Init for MapReads.
        """
        pass

    def get_gene_lookup(self, path: str) -> Dict:
        """
        Create a dictionary from a given look up file.
        File contains columns for GeneIdentifier, start position of exon-1 and strand orientation.
        :param path: str contains path to gff file
        :return: Dict contains GeneIdentifier as key and tuple of start position and strand orientation as value
        """
        matching_dict = {}
        with (open(path, mode='r', newline='')) as infile:
            reader = csv.DictReader(infile, delimiter="\t")
            for row in reader:
                matching_dict[row["GeneIdentifier"]] = (row["Start"], row["Strand"])
        return matching_dict

    def get_identifier_for_datapoint(self, dataframe: pd.DataFrame, chromosome: str, strand: str, stop: int) -> str:
        """
        Find corresponding GeneIdentifier For a given datapoint from BAM file.
        :param dataframe: table of all entries from given gff file
        :param chromosome: str contains name of chromosome (e.g. Chr1)
        :param strand: str contains strand orientation
        :param stop: int contains stop position for given sequence
        :return: matching GeneIdentifier for datapoint as str
        """
        # select subset for ChrX by grouping dataframe by chromosome name (Chr1, ...) and strand orientation (+, -)
        chr_groups = dataframe.groupby(["Chromosome", "Strand"])
        is_group_found = False
        for group in chr_groups:
            if group[0] == (chromosome, strand):
                is_group_found = True
                # sort subset by position
                if strand == "+":
                    sort_val = "Start"
                else:
                    sort_val = "Stop"
                subset = group[1].sort_values(by=sort_val)
                break
        if not is_group_found:  # raise error if datapoint does not match with given mapping data
            raise ValueError(f"Chromosome {chromosome} and strand {strand} do not exist in given gff mapping file.\n")
        # select possible subgroup identifiers where read is in range of unique sequence
        filtered_subset_group = list(subset.query(f"Start <= '{stop}' & Stop >= '{str(stop)}'").groupby(["Domain"]))
        if len(filtered_subset_group) == 1:
            gene_identifier = filtered_subset_group[0][0]
            return str(gene_identifier)
        else:
            matching_gene_identifiers = ', '.join(g[0] for g in filtered_subset_group)
            raise ValueError(
                f"Unique gene identification is not possible due to multiple matches: {matching_gene_identifiers}\n")

    def get_len_from_cigar(self, cigar: List) -> int:
        """
        Calculate sequence length from a given CIGAR string.
        :param cigar: str contains
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
    config.read(os.path.join(os.getcwd().split("src")[0], "config.ini"))
    data_path = config["DATA_PATH"]["path"]
    debug_mode = config["DEBUG_MODE"]["debug"]

    gff_path = os.path.join(data_path, "Araport11_GFF3_genes_transposons.201606.gff")
    bam_path = os.path.join(data_path, "ath.bam")
    mapping_path = os.path.join(os.getcwd().split("src")[0], "output_data/all_exon.tsv")
    out_path = os.path.join(os.getcwd().split("src")[0], "output_data/read_before_first_exon.tsv")

    read_mapper = MapReads()
    gene_matching_dict = read_mapper.get_gene_lookup(path=mapping_path)
    exon_cds_mapper = ExonCdsMapper(path=gff_path)
    df = pd.DataFrame(data=exon_cds_mapper.extract_input_from_file())
    samfile = pysam.AlignmentFile(bam_path, "rb")

    with open(out_path, 'w', newline='') as out_file:
        fieldnames = ["GeneIdentifier", "Start_exon_1", "Read_QNAME", "Read_START"]
        # write all positions of first exon into file
        writer = csv.DictWriter(out_file, fieldnames=fieldnames, delimiter='\t')
        writer.writeheader()

        for row in samfile:
            if debug_mode: print(row)
            # get necessary information for datapoint from BAM file
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

            # extract identifier from gff dataframe
            try:
                gene_identifier = read_mapper.get_identifier_for_datapoint(dataframe=df, chromosome=chr_name,
                                                                           strand=strand,
                                                                           stop=stop_align)
            except Exception as e:
                if debug_mode: print("Could not find identifier for datapoint: ", e)
                continue  # skip row if error occurs

            # validate if identifier is in mapping file
            try:
                exon_start, exon_strand = gene_matching_dict[gene_identifier]
                exon_start = int(exon_start)
            except Exception as e:
                if debug_mode: print("GeneIdentifier is not part of mapping file: ", e)
                continue  # skip row if error occurs

            # check if starting read position is before exon-1 position of given gene.
            if start_align < exon_start and strand == exon_strand:
                if debug_mode: print(
                    f"For GeneIdentifier {gene_identifier} read {q_name} ({start_align}) starts before exon-1 ({exon_start})")
                writer.writerow({"GeneIdentifier": gene_identifier, "Start_exon_1": exon_start, "Read_QNAME": q_name,
                                 "Read_START": start_align})


if __name__ == "__main__":
    main()

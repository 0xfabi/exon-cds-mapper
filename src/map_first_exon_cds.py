#!/usr/bin/env python3
import configparser
import csv
import os
from collections import OrderedDict
from typing import Dict, List

import pandas as pd


class ExonCdsMapper:

    def __init__(self, path):
        """
        Init for ExonCdsMapper class. Init with path for a given gff file.
        :param path: str contains path to gff file
        """
        self.path = path

    def _extract_id_domain(self, add_info_str: str) -> (str, str):
        """
        Extract specific gene identifier is ID (contains type like gene, exon, cds, ...). and
        global gene identifier (ATXGXXXXX) as domain name.
        :param add_info_str: additional information from input file contains gene identifier
        :return: tuple contains gene identifier ID and domain name
        """
        _id = add_info_str.split(";")[0].replace("ID=", "")
        for split_char in [":", "."]:
            if split_char in _id:
                domain = _id.split(split_char)[0]
                break  # break if domain could be converted
            else:
                domain = _id  # for gene domain is equal to its id
        return _id, domain

    def extract_input_from_file(self) -> List:
        """
        Extract from each line in gff file relevant column values and store them as dict.
        Each dict is appended to list which is returned.
        :return: List contains each row of gff file as dict
        """
        column_names = ["ID", "Domain", "Chromosome", "Type", "Start", "Stop", "Strand"]
        row_list = []
        with open(self.path, "r") as infile:
            csv_reader = csv.reader(infile, delimiter='\t')
            for row in csv_reader:
                if "".join(row).startswith("#"):  # ignore comments.
                    continue
                # extract necessary information from given gff file row
                extracted_values = [*self._extract_id_domain(row[-1]), row[0], row[2], row[3], row[4], row[6]]
                # set e.g. exon:1 as type name instead of exon
                if ":" in extracted_values[0]:
                    extracted_values[3] = extracted_values[0].split(":", 1)[1]
                # create dict where keys are entries from column_names list
                data_dict = OrderedDict((col, extracted_values[idx]) for idx, col in enumerate(column_names))
                row_list.append(data_dict)
        return row_list

    def extract_values_from_df(self, dataframe: pd.DataFrame, _type: str) -> (int, str):
        """
        Extract start, stop and strand position from given sub table containing all rows of a corresponding gene
        identifier for a specific type.
        :param dataframe: table of all row entries for a given gene identifier
        :param _type: type of each row entry (e.g. exon:1, cds:1, gene, ...)
        :return: Dict contains chromosome, start, stop position and strand orientation of a given entry
        """
        if not dataframe.loc[dataframe["Type"] == _type].empty:
            entry_values = dataframe.loc[dataframe["Type"] == _type].values[0]
            chromosome = entry_values[2]
            start = entry_values[4]
            stop = entry_values[5]
            strand = entry_values[6]
            return {"Chromosome": chromosome, "Start": start, "Stop": stop, "Strand": strand}

    def get_first_position(self, gene_part: Dict) -> (str, str, str):
        """
        Depending on the strand orientation return the dict key for start and stop position of a given gene dict.
        :param gene_part: dict of a gene
        :return: three strings containing strand orientation (+, -) and start/stop position of gene
        """
        # check on which strand the gene sequence is
        strand = gene_part["Strand"]
        if strand == "+":  # define lookup direction
            start_pos = "Start"  # lookup from left side
            stop_pos = "Stop"
        else:
            start_pos = "Stop"  # lookup from right side
            stop_pos = "Start"
        return strand, start_pos, stop_pos

    def compare_start_position(self, exon: Dict, cds: Dict) -> (bool, int, str):
        """
        Evaluate if given exon and cds start at the same position.
        :param exon: Dict contains information about strand orientation, start and stop position of a first exon
        :param cds: Dict contains information about strand orientation, start and stop position of a first cds
        :return: quadtruple of
            first value: boolean if start position of first exon/cds are equal
            second value: start position of first exon
            third value: stop position of first exon
            fourth value: char contains strand orientation (forward: +, reverse: -)
        """
        # check on which strand the gene sequence is
        strand, start_pos, stop_pos = self.get_first_position(exon)
        if exon[start_pos] == cds[start_pos]:
            return True, exon[start_pos], exon[stop_pos], strand
        else:
            return False, None, None, None


def print_statistics(statistics: List):
    """
    Print statistics on command line.
    :param statistics: List contains some statistics
    :return:
    """
    for line in statistics:
        print(line)


def main():
    # adjust data_path for your system in config.ini
    config = configparser.ConfigParser()
    config.read(os.path.join(os.getcwd().split("src")[0], "config/config.ini"))
    data_path = config["DATA_PATH"]["path"]
    debug_mode = eval(config["DEBUG_MODE"]["debug"])
    statistics = []  # list contains some statistics which are printed after successful execution

    gff_path = os.path.join(data_path, "Araport11_GFF3_genes_transposons.201606.gff")
    matching_out_path = os.path.join(os.getcwd().split("src")[0], "output_data/matching_exon_cds.tsv")

    exon_cds_mapper = ExonCdsMapper(path=gff_path)
    df = pd.DataFrame(data=exon_cds_mapper.extract_input_from_file())
    # df.set_index("ID", inplace=True)
    statistics.append(f"Total number of entries in gff file: {len(df.index)}")
    number_matches = 0  # information about number of genes where start positions are equal (exon:1 == cds:1)

    with open(matching_out_path, 'w', newline='') as match_file:
        fieldnames = ["GeneIdentifier", "Chromosome", "Start", "Stop", "Strand"]
        # write all positions of first exon into file
        match_writer = csv.DictWriter(match_file, fieldnames=fieldnames, delimiter='\t')
        match_writer.writeheader()

        df_grouped = df.groupby(["Domain"])  # group by domain name to get dataframe for each gene
        statistics.append(f"Total number of genes: {df_grouped.ngroups}")

        for group in df_grouped:
            gene_identifier = group[0]
            gene_df = group[1]  # single gene dataframe
            if debug_mode:
                print(f"gene identifier: {gene_identifier}")
                print(f"table of gene: {gene_df}")
            first_exon = exon_cds_mapper.extract_values_from_df(gene_df, "exon:1")
            first_cds = exon_cds_mapper.extract_values_from_df(gene_df, "CDS:1")

            # append gene identifiers to list, if an exon:1 and cds:1 exist for the specific gene sequence
            if first_exon and first_cds:
                is_equal, start, stop, strand = exon_cds_mapper.compare_start_position(first_exon, first_cds)
                if is_equal:
                    if debug_mode:
                        print(f"GeneIdentifier: {gene_identifier}, Chromosome: {first_exon['Chromosome']}, Start: {start}, Stop: {stop}, Strand: {strand}")
                    number_matches += 1
                    match_writer.writerow({"GeneIdentifier": gene_identifier, "Chromosome": first_exon['Chromosome'], "Start": start, "Stop": stop, "Strand": strand})

        statistics.append(f"Number of matches where start position are equal (exon:1 == cds:1): {number_matches}")
        statistics.append(
            f"Proportion of genes with matching exon:1 and cds:1: {round(((number_matches / df_grouped.ngroups) * 100), 2)}%")
        print_statistics(statistics)


if __name__ == "__main__":
    main()

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
        Init for ExonCdsMapper. Init with path for a given gff file.
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

                data_dict = OrderedDict((col, extracted_values[idx]) for idx, col in enumerate(column_names))
                row_list.append(data_dict)
        return row_list

    def extract_values_from_df(self, dataframe: pd.DataFrame, _type: str) -> (int, str):
        """
        Extract start, stop and strand position from given sub table containing all rows of a corresponding gene
        identifier for a specific type.
        :param dataframe: table of all row entries for a given gene identifier
        :param _type: type of each row entry (e.g. exon:1, cds:1, gene, ...)
        :return: Dict contains start, stop position and strand orientation of a given entry
        """
        if not dataframe.loc[dataframe["Type"] == _type].empty:
            entry_values = dataframe.loc[dataframe["Type"] == _type].values[0]
            start = entry_values[4]
            stop = entry_values[5]
            strand = entry_values[6]
            return {"Start": start, "Stop": stop, "Strand": strand}

    def get_first_position(self, gene_part: Dict) -> (str, str):
        # check on which strand the gene sequence is
        strand = gene_part["Strand"]
        if strand == "+":  # define lookup direction
            pos = "Start"  # lookup from left side
        else:
            pos = "Stop"  # lookup from right side
        return strand, pos

    def compare_start_position(self, exon: Dict, cds: Dict) -> (bool, int, str):
        """
        Evaluate if given exon and cds start at the same position.
        :param exon: Dict contains information about strand orientation, start and stop position of a first exon
        :param cds: Dict contains information about strand orientation, start and stop position of a first cds
        :return: triple of
            first value: boolean if start position of first exon/cds are equal
            second value: start position of first exon
            third value: char contains strand orientation (forward: +, reverse: -)
        """
        # check on which strand the gene sequence is
        strand, pos = self.get_first_position(exon)
        if exon[pos] == cds[pos]:
            return True, exon[pos], strand
        else:
            return False, None, None


def main():
    # adjust data_path for your system in config.ini
    config = configparser.ConfigParser()
    config.read(os.path.join(os.getcwd().split("src")[0], "config.ini"))
    data_path = config["DATA_PATH"]["path"]
    debug_mode = eval(config["DEBUG_MODE"]["debug"])

    gff_path = os.path.join(data_path, "Araport11_GFF3_genes_transposons.201606.gff")
    all_out_path = os.path.join(os.getcwd().split("src")[0], "output_data/all_exon.tsv")
    matching_out_path = os.path.join(os.getcwd().split("src")[0], "output_data/matching_exon_cds.tsv")

    exon_cds_mapper = ExonCdsMapper(path=gff_path)
    df = pd.DataFrame(data=exon_cds_mapper.extract_input_from_file())
    # df.set_index("ID", inplace=True)

    with open(matching_out_path, 'w', newline='') as match_file, open(all_out_path, 'w', newline='') as all_file:
        fieldnames = ["GeneIdentifier", "Start", "Strand"]
        # write all positions of first exon into file
        all_writer = csv.DictWriter(all_file, fieldnames=fieldnames, delimiter='\t')
        all_writer.writeheader()
        match_writer = csv.DictWriter(match_file, fieldnames=fieldnames, delimiter='\t')
        match_writer.writeheader()

        # group by domain name to get dataframe for each gene
        for group in df.groupby(["Domain"]):
            gene_identifier = group[0]
            group_df = group[1]  # single gene dataframe
            first_exon = exon_cds_mapper.extract_values_from_df(group_df, "exon:1")
            first_cds = exon_cds_mapper.extract_values_from_df(group_df, "CDS:1")

            if first_exon:
                strand, pos = exon_cds_mapper.get_first_position(first_exon)
                all_writer.writerow({"GeneIdentifier": gene_identifier, "Start": first_exon[pos], "Strand": strand})

            # append gene identifiers to list, if an exon:1 and cds:1 exist for the specific gene sequence
            if first_exon and first_cds:
                is_equal, pos, strand = exon_cds_mapper.compare_start_position(first_exon, first_cds)
                if is_equal:
                    if debug_mode:
                        print(f"GeneIdentifier: {gene_identifier}, Start: {pos}, Strand: {strand}")
                    match_writer.writerow({"GeneIdentifier": gene_identifier, "Start": pos, "Strand": strand})


if __name__ == "__main__":
    main()

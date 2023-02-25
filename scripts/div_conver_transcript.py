#!/usr/bin/env python
# coding: utf-8

import os
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature, FeatureLocation
from Bio import SeqIO
import argparse as arg
import bz2
import pandas as pd


def cross_info_tus(data, MatrixTus):

    for count,line_predict in enumerate(data):
        
        for id_prot_tu in MatrixTus:
            
            TU_whole_lst = id_prot_tu[1].split(",")
            TU_whole_str = str(id_prot_tu[1])
            
            if len(TU_whole_lst) > 1:
                
                if TU_whole_lst[0].find(line_predict[0]) > -1:
                    data[count][0] = TU_whole_str
                
                if TU_whole_lst[0].find(line_predict[2]) > -1:
                    data[count][2] = TU_whole_str
    
    return data


if __name__ == "__main__":

    parser = arg.ArgumentParser(description='Get convergent/divergent genes. Also generate table of ids with type data [gene-TF]')
    parser.add_argument("--inputPathGBK", dest="inputPathGBK")
    parser.add_argument("--inputPathTUS", dest="inputPathTUS")
    parser.add_argument("--inputPathFasta", dest="inputPathFasta")
    parser.add_argument("--inputPathNetwork", dest="inputPathNetwork")
    parser.add_argument("--outputPath", dest="outputPath")

    args = parser.parse_args()

    print("\n#### PARAMETERS ####\n")
    print("inputPathGBK: ",args.inputPathGBK)
    print("inputPathTUS: ",args.inputPathTUS)
    print("inputPathFasta: ",args.inputPathFasta)
    print("inputPathNetwork: ",args.inputPathNetwork)
    print("outputPath: ",args.outputPath,"\n")
    print("####################\n")

    print("Reading genbank file...")

    # get all sequence records for the specified genbank file
    recs = [rec for rec in SeqIO.parse(args.inputPathGBK, "genbank")]

    # get the sequence record
    for rec in recs:
        feats = [feat for feat in rec.features if feat.type == "CDS"]

    # get protein id, genomic ranges (location) and strand 
    Range = {}
    for count,feat in enumerate(feats):
        try:
        
            key = feat.qualifiers["protein_id"][0]
            value = feat.location
            Range[str(key)] = value

        except KeyError as error:
            print(f"Key not found: {error} in {feat}")

    print(f"Total elements: {count}")

    # Convergent transcription

    print("Getting transcrip type...")

    keys_iterate = list(Range.keys())
    before_id, before_range = keys_iterate[0],Range[keys_iterate[0]]

    conver_ids = []

    for k in keys_iterate[1:]:

        if Range[k].strand != before_range.strand:

            if Range[k].start.real <= before_range.end.real:

                data = [before_id, str(before_range), k, str(Range[k]), str(before_range.end.real - Range[k].start.real + 1)]
                conver_ids.append(data)

        before_id, before_range = k, Range[k]


    # Divergent transcription

    keys_iterate = list(Range.keys())
    before_id, before_range = keys_iterate[0],Range[keys_iterate[0]]

    diver_ids = []

    for k in keys_iterate[1:]:

        if Range[k].strand != before_range.strand:

            if Range[k].strand == 1:
            
                if Range[k].start > before_range.end.real:
                
                    data = [before_id, str(before_range), k, str(Range[k]), str(before_range.end.real - Range[k].start)]
                    diver_ids.append(data)

        before_id, before_range = k, Range[k]


    print("Convergent pairs found: {}".format(len(conver_ids)))
    print("Divergent pairs found: {}".format(len(diver_ids)))

    MatrixTus = []

    print("Reading TUs file...")

    with bz2.open(args.inputPathTUS,"rt") as tFile:

        for line in tFile:
            if line.find("#") == -1:
                line = line.strip("\n").strip().split("\t")
                MatrixTus.append(line)


    print("Cross-information TUs and pairs found...")

    # Divergent transcription
    diver_ids = cross_info_tus(diver_ids, MatrixTus)

    # Convergent transcription
    conver_ids = cross_info_tus(conver_ids, MatrixTus)


    print("Saving files...")

    # Saving divergent
    colmn = ["#element_1","location_firts_gene/firts_element-strand","element_2","location_firts_gene/second_element-strand","distance_bp"]

    diver_ids_df = pd.DataFrame(diver_ids, columns=colmn)

    diver_ids_df.to_csv(args.outputPath + "divergent_genes", index=False, sep="\t")


    # Saving convergent
    colmn = ["#element_1","location_firts_gene/firts_element-strand","element_2","location_firts_gene/second_element-strand","overlapping_bp"]

    diver_ids_df = pd.DataFrame(conver_ids, columns=colmn)

    diver_ids_df.to_csv(args.outputPath + "convergent_genes", index=False, sep="\t")

    # Reading fasta file
    print("Reading fasta file...")

    headers_fasta = [] 
    with open(args.inputPathFasta,"r") as iFilef:
        for line in iFilef:
            if line.find(">") > -1:
                line = line.strip("\n").strip()
                line = line.replace(">","").replace("|","\t")
                line = line.split("\t")[0:3]
                headers_fasta.append(line)

    headers_fasta_df = pd.DataFrame(headers_fasta, columns=["id_1","id_2","common_name"])
    headers_fasta_df["type"] = "gene"

    # Reading network file
    print("Reading network file...")

    transcription_factors = pd.read_csv(args.inputPathNetwork, sep="\t", header=None, usecols=[1])
    transcription_factors = transcription_factors[1].apply("unique")

    # Assigning type [gene-TF]
    print("Creating table of ids...")

    for count,id in enumerate(headers_fasta_df["id_1"]):
    
        for TF in transcription_factors:
    
            if id.find(TF) > -1:
                headers_fasta_df.loc[count]["type"] = "TF"

    print("Saving file...")

    headers_fasta_df.to_csv(args.outputPath + "eq_ids", header=True, index=False, sep="\t")

    print("All done")

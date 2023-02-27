#!/usr/bin/env python
# coding: utf-8


import bz2
import time
import os
import pandas as pd
import argparse as arg

tiempo = time.time()

if __name__ == "__main__":

    # Retrieving bash parameters and printing

    argp = arg.ArgumentParser(description="Cross modified network and transcription unit (TU) files")
    argp.add_argument("--inputNet", dest="inputNet", help="Input file of the extended network")
    argp.add_argument("--inputTus", dest="inputTus", help="Input file of transcriptional units")
    argp.add_argument("--outputPath", dest="outputPath", help="Path to save files")
    args = argp.parse_args()

    print("\n##### PARAMETERS #####")
    print(f"inputNet: {args.inputNet}")
    print(f"inputTus: {args.inputTus}")
    print(f"outputPath: {args.outputPath}")
    print("######################","\n")

    MatrixNet = []
    MatrixTus = []
    FinalOutp = []

    # Opening and parsing each line of both files
    print("Reading files...")
    MatrixNet = pd.read_csv(args.inputNet, header=None,sep="\t", usecols=[1,2]).rename({1:"TF",2:"TG"}, axis=1)

    with bz2.open(args.inputTus,"rt") as tFile:

        for line in tFile:
            # Ignore headers
            if line.find("#") < 0:
                line = line.strip("\n").strip().split("\t")
                MatrixTus.append(line)

    MatrixTus = pd.DataFrame(MatrixTus).rename({0:"strain", 1:"TU", 2:"common_name"}, axis=1)

    print("Cross network-tus info...")
    
    # Assign an id to the different TUS (only ids mapped, see if inside the second for)
    id_tus = 1
    FinalOutp = []
    
    for ind_net,TG in enumerate(MatrixNet["TG"]):

        for ind_tus,TU in enumerate(MatrixTus["TU"]):
            
            # Create a list with the tus
            TU_as_list = TU.replace(" ","").split(",")
            # If the firts element of the TU is equal to the target described in the network, go ahead
            
            if TU_as_list[0].find(TG) > -1:
                # Iterating over the genes inside the TU
                
                for gene in TU_as_list:
                    # Giving a format to the output: STRAIN-TF-TG-TU_id
                    FinalOutp.append([MatrixTus.iloc[ind_tus]["strain"], MatrixNet.iloc[ind_net]["TF"], gene, str(id_tus) + "_tu"])
                
                id_tus+=1

    FinalOutp = pd.DataFrame(FinalOutp)

    # Saving final output into the specified directory
    Prefix = os.path.basename(args.inputNet).replace("_extended_network.txt","")
    FileNameOutput = args.outputPath + Prefix +"_TF_TG_predictionsTU.txt"
    FinalOutp.to_csv(FileNameOutput, header=False, index=False, sep="\t")

    # Message about final file and time
    print("Making {}\t{}".format(FileNameOutput,time.time()-tiempo))
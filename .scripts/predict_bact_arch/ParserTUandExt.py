#!/usr/bin/env python
# coding: utf-8

import sys
import bz2
import time
import os

tiempo = time.time()

# Retrieving bash parameters 

ExtenNet = sys.argv[1]
TusFiles = sys.argv[2]
PathOutpName = sys.argv[3]
TargetCol= sys.argv[4]

MatrixNet = []
MatrixTus = []
FinalOutp = []

# Opening and parsing each line of both files

with open(ExtenNet,"r") as eFile:
    
    for line in eFile:
        
        line = line.strip("\n").strip().split("\t")
        MatrixNet.append(line)


with bz2.open(TusFiles,"rt") as tFile:

	for line in tFile:

		if line.find("#") == -1:
			line = line.strip("\n").strip().split("\t")
			MatrixTus.append(line)

count = 1

# Retrieving number column where are TFs and targets

TargetCol = int(TargetCol.strip()) - 1


for index_Net in range(len(MatrixNet)):
    
    # Retrieving the target of each interaction in the extended network

    TG = MatrixNet[index_Net][TargetCol]

    for index_Tus in range(len(MatrixTus)):

	    # Retrieving trancription units from TU file and making an array based on        
        #Tus = MatrixTus[index_Tus][2].replace(" ","").split(",")
        
        Tus = MatrixTus[index_Tus][1].replace(" ","").split(",")

	    # If TG is found even as a substring in the first transcription unit, then everything else remaining of the TU is considered as regulated by the transcription factor of TG

        if Tus[0].find(TG) > -1:
        
            for genes in Tus:

                # Giving a format to the output: LOCUS-STRAIN-TF-TG-TU_id

                FinalOutp.append("NULL" + "\t" + MatrixTus[index_Tus][0] + "\t" + MatrixNet[index_Net][TargetCol - 1] + "\t" + genes + "\t" + str(count)+"_tu"+"\n")
    
            count+=1

# Saving final output into the specified directory and name

with open(PathOutpName, "w") as oFile:
    
    for line in FinalOutp:
        oFile.write(line)

        
print("(background) Making {}\t{}".format(PathOutpName,time.time()-tiempo))




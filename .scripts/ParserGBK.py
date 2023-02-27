
import os
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature, FeatureLocation
from Bio import SeqIO
import argparse

parser = argparse.ArgumentParser(description="Getting faa file from genbak file")
parser.add_argument("--file", dest="file")
parser.add_argument("--org", dest="org", default="",
                     help="Organism name")
parser.add_argument("--get", dest="get",
                     help="Specify information to get separate them by comma.")
parser.add_argument("--sep", dest="sep", default="\t",
                     help="Character to separate head elements")
args = parser.parse_args()

# get all sequence records for the specified genbank file
recs = [rec for rec in SeqIO.parse(args.file, "genbank")]

# print the CDS sequence feature summary information for each feature in each
# sequence record
for rec in recs:
    feats = [feat for feat in rec.features if feat.type == "CDS"]
    for feat in feats:
        try:
            feat.qualifiers["translation"]
            
            retrieve = str(args.get).split(",")

            print(">", end="")
            for data in retrieve:
                try:
                    
                    if data != "translation": 
                        print("{}{}".format(feat.qualifiers[str(data)][0], args.sep), end="")
                
                except KeyError as error:

                    print("{}{}".format("NODATA", args.sep), end="" )

            print(args.org)

            if "translation" in retrieve:
                print(feat.qualifiers["translation"][0])

        except KeyError as error:
            pass

        
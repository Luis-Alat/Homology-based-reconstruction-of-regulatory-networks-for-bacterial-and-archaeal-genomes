# %%
import pandas as pd
import argparse
import os
import ast
from time import time

# %%
def ParserArguments() -> argparse.ArgumentParser:
    
    '''
    This functions parser the command line arguments supplied by the user
    
    '''

    argp = argparse.ArgumentParser(description="Cross modified network and transcription unit (TU) files")

    argp.add_argument("--input_net", dest="input_net", help="Input file of the extended network", type=str)
    argp.add_argument("--arguments_read_net", dest="arguments_read_net", help="Arguments to be passed to the pandas read csv method while reading --inputNet")
    argp.add_argument("--input_tus", dest="input_tus", help="Input file of transcriptional units", type=str)
    argp.add_argument("--arguments_read_tus", dest="arguments_read_tus", help="Arguments to be passed to the pandas read csv method while reading --inputTus")
    argp.add_argument("--use_first_match", dest="use_first_match", help="True or False", type=bool, default=False)
    argp.add_argument("--output_path_name", dest="output_path_name", help="Path to save files", type=str)

    args = argp.parse_args()

    # Parsing arguments for pandas
    args.arguments_read_net = ast.literal_eval(args.arguments_read_net)
    args.arguments_read_tus = ast.literal_eval(args.arguments_read_tus)

    return args

def ShowArguments(arguments:argparse.ArgumentParser):

    print("\n##### PARAMETERS #####")
    
    for arg in vars(arguments):
        print(f"{arg}: {getattr(arguments, arg)}")
    
    print("######################","\n")

# %%
if __name__ == "__main__":

    start = time()

    args = ParserArguments()
    ShowArguments(args)
    
    network_df = pd.read_csv(args.input_net, **args.arguments_read_net)
    network_df.rename({args.arguments_read_net["usecols"][0]: "TF",
                       args.arguments_read_net["usecols"][1]: "TG"}, axis=1, inplace=True)
    network_df = network_df.astype(str)
    
    tus_df = pd.read_csv(args.input_tus, **args.arguments_read_tus)
    tus_df.rename({0:"id", 1:"tu", 2:"strand"}, axis=1, inplace=True)
    tus_df = tus_df.astype(str)
    
    print("Cross network-tus info...")
    
    final_output = pd.DataFrame(columns=["strand","TF","TG","id"])

    for TF, TG in zip(network_df["TF"], network_df["TG"]):
        for index_tus, current_unit in enumerate(tus_df["tu"]):
        
            if args.use_first_match:
                
                tgs_splitted = current_unit.split(",")

                if TG in tgs_splitted[0]:
                    expected_format = {"strand": tus_df["strand"].iloc[index_tus],
                                        "TF": TF,
                                        "TG": tgs_splitted,
                                        "id": tus_df["id"].iloc[index_tus] + "_tu"}

                    final_output = pd.concat([final_output, pd.DataFrame(expected_format)])

            else:

                if TG in current_unit:
                    expected_format = {"strand": tus_df["strand"].iloc[index_tus],
                                        "TF": TF,
                                        "TG": current_unit.split(","),
                                        "id": tus_df["id"].iloc[index_tus] + "_tu"}

                    final_output = pd.concat([final_output, pd.DataFrame(expected_format)])

    final_output.drop_duplicates(inplace=True)

    print(f"Making {args.output_path_name}")
    final_output.to_csv(args.output_path_name, header=False, index=False, sep="\t")
    print(f"All donde in {time() - start}")



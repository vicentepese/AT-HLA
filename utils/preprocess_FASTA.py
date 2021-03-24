import numpy 
import json
from collections import defaultdict
import pandas as pd

def main():
    
    # Read file and parse codons
    prot_dict = defaultdict(str)
    with open("Resources/hla_prot.fasta") as inFile:
        for line in inFile:
            if ">" in line:
                codon = line.split(" ")[1]
            else:
                prot_dict[codon] += line.split("\n")[0]

    # Convert to pandas dataframe and and modify keys
    df = pd.DataFrame({ "locus": [k.split("*")[0] for k in prot_dict.keys()], \
        "allele": [":".join(k.split(":")[0:2]) for k in prot_dict.keys()],\
        "sequence":[v for v in prot_dict.values()]})

    # Drop duplicated keys but keep first instance 
    df.drop_duplicates(subset = "allele", keep = "first", inplace = True)

    # Write 
    df.to_csv("Resources/AminoAcid_alignment_2021.csv", sep = ",", index=False)
    


if __name__ == "__main__":
    main()
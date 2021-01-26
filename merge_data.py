import pandas as pd 
import os 

def main():

    # Import settings 
    settings = pd.read_json('settings.json')

    # Load each dataset and merge 
    data_merged = []
    for dataset in os.listdir(settings.directory.Indiv_Data):
        data_merged.append(pd.read_csv(settings.directory.Indiv_Data + dataset))
    data_merged = pd.concat(data_merged)

    # Write 
    data_merged.to_csv(settings.file.HLA_Data, index = False)

if __name__ == "__main__":
    main()
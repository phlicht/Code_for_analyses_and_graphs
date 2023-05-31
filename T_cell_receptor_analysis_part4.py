                                            ################################
                                            ### T cell receptor analysis ###
                                            ################################

                                                      ##############
                                                      ### Part 4 ###
                                                      ##############

# Perform Binding Affinity assessment with ERGO-II
# language: python
import ERGOII

directory = "path/to/ERGOII/epitope_pairs/Saureus_IEDB_predicted/"
dataset = "vdjdb"

for file_name in os.listdir(directory):
    if file_name.endswith(".csv"):
        datafile = os.path.join(directory, file_name)
        sample = file_name.split(".")[0]
        print(sample)
        df = pd.read_csv(datafile)

        if __name__ == '__main__':
            df = ERGOII.predict(dataset, datafile)
            filepath = "path/to/ERGOII/results/Saureus_IEDB_predicted/" + sample + "_ERGOII.csv"
            df.to_csv(filepath, index=False)
            pass

#####################################
###       PROCEED WITH PART 5     ###
### to extract meaningful results ###
#####################################
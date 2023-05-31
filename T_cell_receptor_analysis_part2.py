                                            ################################
                                            ### T cell receptor analysis ###
                                            ################################

                                                      ##############
                                                      ### Part 2 ###
                                                      ##############


###################################################
###   Assess TCR Binding Affinity to epitopes   ###
### of S. aureus virulence factors with ERGO-II ###
###################################################

# Clean the data.frame containing TCR information from
# possible asterisks and slash because ERGO-II cannot handle them
import pandas as pd

df = pd.read_csv("/path/to/TCR/Pat1_blood.csv", keep_default_na=False)
df = df.astype(str)
mask = df.apply(lambda x: x.str.contains('\*')).any(axis=1)
clean_df = df.loc[~mask]
clean_df = clean_df.fillna('')
clean_df.to_csv("/path/to/TCR_clean/Pat1_blood.csv", index=False, na_rep='')

###########################
### PROCEED WITH PART 3 ###
###########################
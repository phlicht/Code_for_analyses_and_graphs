                                            ################################
                                            ### T cell receptor analysis ###
                                            ################################

                                                      ##############
                                                      ### Part 1 ###
                                                      ##############

#####################################################
### TCR repertoire analysis with  immunarch 0.7.0 ###
#####################################################
library(immunarch)
library(tidyverse)

# Number of unique clonotypes
repExplore(
  .data = TCR_skin$data,
  .method = "volume",
  .col = "aa",
  .coding = TRUE
)

# Diversity estimation
repDiversity(
  .data = TCR_skin$data,
  .method = "chao1",
  .col = "aa"
)

# Tracking clonotypes across tissues and stages
# Shown is a representitive example for Pat1
Pat1 <- TCR_blood_skin$data[c("Pat1_Blood", "Pat1_Plaque")]
Pat1_TCR_tracking <- trackClonotypes(
    .data = Pat1,
    .which = list("Pat1_Plaque", 10),
    .col = "aa",
    .norm = FALSE
) %>%
    vis(
        .plot = "smooth"
    )
Pat1_TCR_tracking +
    theme(
    axis.text.x = element_text(angle = 0, hjust = .35),
    legend.position = "none"
    )

# TRBV Gene Usage Plot
geneUsage(
  .data = TCR_skin$data,
  .gene = "hs.trbv",
  .quant = "count",
  .ambig = "inc",
  .type = "segment",
  .norm = FALSE
) %>%
    vis(
        .meta = TCR_skin$meta,
        .by = c("Stage", "Staphylococcus_aureus"),
        .plot = "hist"
    )


###################################################
###         PROCEED WITH PART 2 TO              ###
### Assess TCR Binding Affinity to epitopes     ###
### of S. aureus virulence factors with ERGO-II ###
###################################################

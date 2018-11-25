#This is necessary to remove the single exon transcripts with TPM<2

setwd("~/research/polyA_ribozero/equcab3")

require(tidyr)
require(dplyr)
require(stringr)

AH1_Liver_Ribozero_remove <- read.table("683610_Liver_Ribozero_S51_remove.bed", 
                                        header=FALSE, stringsAsFactors = FALSE, sep = "\t")
AH1_Liver_Ribozero_remove <- AH1_Liver_Ribozero_remove[with(AH1_Liver_Ribozero_remove, order(V1,V2)),]

AH1_Liver_Ribozero_all <- read.table("683610_Liver_Ribozero_S51.bed",
                                     header=FALSE, stringsAsFactors = FALSE, sep = "\t")
AH1_Liver_Ribozero_all <- AH1_Liver_Ribozero_all[with(AH1_Liver_Ribozero_all, order(V1,V2)),]

AH1_Liver_Ribozero_keep <- anti_join(AH1_Liver_Ribozero_all, AH1_Liver_Ribozero_remove, by="V4")

AH1_Liver_Ribozero_keep$V7 = AH1_Liver_Ribozero_keep$V3 - AH1_Liver_Ribozero_keep$V2
  
write.table(AH1_Liver_Ribozero_keep, "AH1_Liver_Ribozero_length.bed", 
            row.names = FALSE, col.names = FALSE, quote = FALSE, sep = "\t")


AH1_Liver_remove <- read.table("683610_Liver_S61_remove.bed",
                               header=FALSE, stringsAsFactors = FALSE, sep = "\t")
AH1_Liver_remove <- AH1_Liver_remove[with(AH1_Liver_remove, order(V1,V2)),]

AH1_Liver_all <- read.table("683610_Liver_S61.bed",
                            header=FALSE, stringsAsFactors = FALSE, sep = "\t")
AH1_Liver_all <- AH1_Liver_all[with(AH1_Liver_all, order(V1,V2)),]

AH1_Liver_keep <- anti_join(AH1_Liver_all, AH1_Liver_remove, by="V4")

AH1_Liver_keep$V7 = AH1_Liver_keep$V3 - AH1_Liver_keep$V2

write.table(AH1_Liver_keep, "AH1_Liver_length.bed",
            row.names = FALSE, col.names = FALSE, quote = FALSE, sep = "\t")

AH1_Parietal_Cortex_Ribozero_remove <- read.table("683610_Parietal_Cortex_Ribozero_S53_remove.bed",
                                                  header=FALSE, stringsAsFactors = FALSE, sep = "\t")
AH1_Parietal_Cortex_Ribozero_remove <- AH1_Parietal_Cortex_Ribozero_remove[with(AH1_Parietal_Cortex_Ribozero_remove, order(V1,V2)),]

AH1_Parietal_Cortex_Ribozero_all <- read.table("683610_Parietal_Cortex_Ribozero_S53.bed",
                                               header=FALSE, stringsAsFactors = FALSE, sep = "\t")
AH1_Parietal_Cortex_Ribozero_all <- AH1_Parietal_Cortex_Ribozero_all[with(AH1_Parietal_Cortex_Ribozero_all, order(V1,V2)),]

AH1_Parietal_Cortex_Ribozero_keep <- anti_join(AH1_Parietal_Cortex_Ribozero_all, AH1_Parietal_Cortex_Ribozero_remove, by="V4")

AH1_Parietal_Cortex_Ribozero_keep$V7 = AH1_Parietal_Cortex_Ribozero_keep$V3 - AH1_Parietal_Cortex_Ribozero_keep$V2

write.table(AH1_Parietal_Cortex_Ribozero_keep, "AH1_Parietal_Cortex_Ribozero_length.bed",
            row.names = FALSE, col.names = FALSE, quote = FALSE, sep = "\t")


AH1_Parietal_Cortex_remove <- read.table("683610_Parietal_Cortex_S55_remove.bed",
                                         header=FALSE, stringsAsFactors = FALSE, sep = "\t")
AH1_Parietal_Cortex_remove <- AH1_Parietal_Cortex_remove[with(AH1_Parietal_Cortex_remove, order(V1,V2)),]

AH1_Parietal_Cortex_all <- read.table("683610_Parietal_Cortex_S55.bed",
                                      header=FALSE, stringsAsFactors = FALSE, sep = "\t")
AH1_Parietal_Cortex_all <- AH1_Parietal_Cortex_all[with(AH1_Parietal_Cortex_all, order(V1,V2)),]

AH1_Parietal_Cortex_keep <- anti_join(AH1_Parietal_Cortex_all, AH1_Parietal_Cortex_remove, by="V4")

AH1_Parietal_Cortex_keep$V7 = AH1_Parietal_Cortex_keep$V3 - AH1_Parietal_Cortex_keep$V2

write.table(AH1_Parietal_Cortex_keep, "AH1_Parietal_Cortex_length.bed",
            row.names = FALSE, col.names = FALSE, quote = FALSE, sep = "\t")



AH2_Liver_Ribozero_remove <- read.table("686521_Liver_Ribozero_S50_remove.bed",
                                        header=FALSE, stringsAsFactors = FALSE, sep = "\t")
AH2_Liver_Ribozero_remove <- AH2_Liver_Ribozero_remove[with(AH2_Liver_Ribozero_remove, order(V1,V2)),]

AH2_Liver_Ribozero_all <- read.table("686521_Liver_Ribozero_S50.bed",
                                     header=FALSE, stringsAsFactors = FALSE, sep = "\t")
AH2_Liver_Ribozero_all <- AH2_Liver_Ribozero_all[with(AH2_Liver_Ribozero_all, order(V1,V2)),]

AH2_Liver_Ribozero_keep <- anti_join(AH2_Liver_Ribozero_all, AH2_Liver_Ribozero_remove, by="V4")

AH2_Liver_Ribozero_keep$V7 = AH2_Liver_Ribozero_keep$V3 - AH2_Liver_Ribozero_keep$V2

write.table(AH2_Liver_Ribozero_keep, "AH2_Liver_Ribozero_length.bed",
            row.names = FALSE, col.names = FALSE, quote = FALSE, sep = "\t")


AH2_Liver_remove <- read.table("686521_Liver_S60_remove.bed",
                               header=FALSE, stringsAsFactors = FALSE, sep = "\t")
AH2_Liver_remove <- AH2_Liver_remove[with(AH2_Liver_remove, order(V1,V2)),]

AH2_Liver_all <- read.table("686521_Liver_S60.bed",
                            header=FALSE, stringsAsFactors = FALSE, sep = "\t")
AH2_Liver_all <- AH2_Liver_all[with(AH2_Liver_all, order(V1,V2)),]

AH2_Liver_keep <- anti_join(AH2_Liver_all, AH2_Liver_remove, by="V4")

AH2_Liver_keep$V7 = AH2_Liver_keep$V3 - AH2_Liver_keep$V2

write.table(AH2_Liver_keep, "AH2_Liver_length.bed",
            row.names = FALSE, col.names = FALSE, quote = FALSE, sep = "\t")


AH2_Parietal_Cortex_Ribozero_remove <- read.table("686521_Parietal_Cortex_Ribozero_S52_remove.bed",
                                                  header=FALSE, stringsAsFactors = FALSE, sep = "\t")
AH2_Parietal_Cortex_Ribozero_remove <- AH2_Parietal_Cortex_Ribozero_remove[with(AH2_Parietal_Cortex_Ribozero_remove, order(V1,V2)),]

AH2_Parietal_Cortex_Ribozero_all <- read.table("686521_Parietal_Cortex_Ribozero_S52.bed",
                                               header=FALSE, stringsAsFactors = FALSE, sep = "\t")
AH2_Parietal_Cortex_Ribozero_all <- AH2_Parietal_Cortex_Ribozero_all[with(AH2_Parietal_Cortex_Ribozero_all, order(V1,V2)),]

AH2_Parietal_Cortex_Ribozero_keep <- anti_join(AH2_Parietal_Cortex_Ribozero_all, AH2_Parietal_Cortex_Ribozero_remove, by="V4")

AH2_Parietal_Cortex_Ribozero_keep$V7 = AH2_Parietal_Cortex_Ribozero_keep$V3 - AH2_Parietal_Cortex_Ribozero_keep$V2

write.table(AH2_Parietal_Cortex_Ribozero_keep, "AH2_Parietal_Cortex_Ribozero_length.bed",
            row.names = FALSE, col.names = FALSE, quote = FALSE, sep = "\t")


AH2_Parietal_Cortex_remove <- read.table("686521_Parietal_Cortex_S54_remove.bed",
                                         header=FALSE, stringsAsFactors = FALSE, sep = "\t")
AH2_Parietal_Cortex_remove <- AH2_Parietal_Cortex_remove[with(AH2_Parietal_Cortex_remove, order(V1,V2)),]

AH2_Parietal_Cortex_all <- read.table("686521_Parietal_Cortex_S54.bed",
                                      header=FALSE, stringsAsFactors = FALSE, sep = "\t")
AH2_Parietal_Cortex_all <- AH2_Parietal_Cortex_all[with(AH2_Parietal_Cortex_all, order(V1,V2)),]

AH2_Parietal_Cortex_keep <- anti_join(AH2_Parietal_Cortex_all, AH2_Parietal_Cortex_remove, by="V4")

AH2_Parietal_Cortex_keep$V7 = AH2_Parietal_Cortex_keep$V3 - AH2_Parietal_Cortex_keep$V2

write.table(AH2_Parietal_Cortex_keep, "AH2_Parietal_Cortex_length.bed",
            row.names = FALSE, col.names = FALSE, quote = FALSE, sep = "\t")

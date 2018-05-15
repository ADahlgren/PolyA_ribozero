setwd("~/research/polyA_ribozero/Filtering")

library(dplyr)

#on server, used bedops convert2bed to convert refined gtf from E.Scott and T.Mansour to bed format
#downloaded bed file

refined_bed <- read.table("refined.bed", header = FALSE, stringsAsFactors = FALSE, sep = "\t")
#There was the following step in pipeline, but all files have same (as in no) header
#names(AH1_Liver_Ribozero_f3_bed) <- names(refined_edited_bed)

#Read in result of filter 3
AH1_Liver_Ribozero_f3_bed <- read.table("AH1_Liver_Ribozero_f3.bed", 
                                        header = FALSE, stringsAsFactors = FALSE, sep = "\t")
AH1_Liver_f3_bed <- read.table("AH1_Liver_f3.bed", 
                               header = FALSE, stringsAsFactors = FALSE, sep = "\t")
AH1_PC_Ribozero_f3_bed <- read.table("AH1_PC_Ribozero_f3.bed",
                                     header = FALSE, stringsAsFactors = FALSE, sep = "\t")
AH1_PC_f3_bed <- read.table("AH1_PC_f3.bed",
                            header = FALSE, stringsAsFactors = FALSE, sep = "\t")
AH2_Liver_Ribozero_f3_bed <- read.table("AH2_Liver_Ribozero_f3.bed",
                                        header = FALSE, stringsAsFactors = FALSE, sep = "\t")
AH2_Liver_f3_bed <- read.table("AH2_Liver_f3.bed",
                               header = FALSE, stringsAsFactors = FALSE, sep = "\t")
AH2_PC_Ribozero_f3_bed <- read.table("AH2_PC_Ribozero_f3.bed",
                                     header = FALSE, stringsAsFactors = FALSE, sep = "\t")
AH2_PC_f3_bed <- read.table("AH2_PC_f3.bed",
                            header = FALSE, stringsAsFactors = FALSE, sep = "\t")

#Removing any overlap from custom transcriptome, which should only have protein coding, and these files
AH1_Liver_Ribozero_refined_nolncRNA_bed <- anti_join(refined_bed, AH1_Liver_Ribozero_f3_bed, by="V4")
AH1_Liver_refined_nolncRNA_bed <- anti_join(refined_bed, AH1_Liver_f3_bed, by="V4")
AH1_PC_Ribozero_refined_nolncRNA_bed <- anti_join(refined_bed, AH1_PC_Ribozero_f3_bed, by="V4")
AH1_PC_refined_nolncRNA_bed <- anti_join(refined_bed, AH1_PC_f3_bed, by="V4")
AH2_Liver_Ribozero_refined_nolncRNA_bed <- anti_join(refined_bed, AH2_Liver_Ribozero_f3_bed, by="V4")
AH2_Liver_refined_nolncRNA_bed <- anti_join(refined_bed, AH2_Liver_f3_bed, by="V4")
AH2_PC_Ribozero_refined_nolncRNA_bed <- anti_join(refined_bed, AH2_PC_Ribozero_f3_bed, by="V4")
AH2_PC_refined_nolncRNA_bed <- anti_join(refined_bed, AH2_PC_f3_bed, by="V4")

#order chr
AH1_Liver_Ribozero_refined_nolncRNA_bed <- AH1_Liver_Ribozero_refined_nolncRNA_bed[with
                                           (AH1_Liver_Ribozero_refined_nolncRNA_bed, order(V1,V2)), ]
AH1_Liver_refined_nolncRNA_bed <- AH1_Liver_refined_nolncRNA_bed[with
                                  (AH1_Liver_refined_nolncRNA_bed, order(V1,V2)), ]
AH1_PC_Ribozero_refined_nolncRNA_bed <- AH1_PC_Ribozero_refined_nolncRNA_bed[with
                                        (AH1_PC_Ribozero_refined_nolncRNA_bed, order(V1,V2)), ]
AH1_PC_refined_nolncRNA_bed <- AH1_PC_refined_nolncRNA_bed[with
                               (AH1_PC_refined_nolncRNA_bed, order(V1,V2)), ]
AH2_Liver_Ribozero_refined_nolncRNA_bed <- AH2_Liver_Ribozero_refined_nolncRNA_bed[with
                                           (AH2_Liver_Ribozero_refined_nolncRNA_bed, order(V1,V2)), ]
AH2_Liver_refined_nolncRNA_bed <- AH2_Liver_refined_nolncRNA_bed[with
                                  (AH2_Liver_refined_nolncRNA_bed, order(V1,V2)), ]
AH2_PC_Ribozero_refined_nolncRNA_bed <- AH2_PC_Ribozero_refined_nolncRNA_bed[with
                                        (AH2_PC_Ribozero_refined_nolncRNA_bed, order(V1,V2)), ]
AH2_PC_refined_nolncRNA_bed <- AH2_PC_refined_nolncRNA_bed[with
                               (AH2_PC_refined_nolncRNA_bed, order(V1,V2)), ]


write.table(AH1_Liver_Ribozero_refined_nolncRNA_bed, "AH1_Liver_Ribozero_refined_nolncRNA.bed",
            row.names=F, col.names=F, quote=F, sep = "\t")
write.table(AH1_Liver_refined_nolncRNA_bed, "AH1_Liver_refined_nolncRNA.bed",
            row.names=F, col.names=F, quote=F, sep = "\t")
write.table(AH1_PC_Ribozero_refined_nolncRNA_bed, "AH1_PC_Ribozero_refined_nolncRNA.bed",
            row.names=F, col.names=F, quote=F, sep = "\t")
write.table(AH1_PC_refined_nolncRNA_bed, "AH1_PC_refined_nolncRNA.bed",
            row.names=F, col.names=F, quote=F, sep = "\t")
write.table(AH2_Liver_Ribozero_refined_nolncRNA_bed, "AH2_Liver_Ribozero_refined_nolncRNA.bed",
            row.names=F, col.names=F, quote=F, sep = "\t")
write.table(AH2_Liver_refined_nolncRNA_bed, "AH2_Liver_refined_nolncRNA.bed",
            row.names=F, col.names=F, quote=F, sep = "\t")
write.table(AH2_PC_Ribozero_refined_nolncRNA_bed, "AH2_PC_Ribozero_refined_nolncRNA.bed",
            row.names=F, col.names=F, quote=F, sep = "\t")
write.table(AH2_PC_refined_nolncRNA_bed, "AH2_PC_refined_nolncRNA.bed",
            row.names=F, col.names=F, quote=F, sep = "\t")

##observe overlap of lncRNA with transcriptome
#add 1 kb to end of each gene for merged transcriptome for UTRs
AH1_Liver_Ribozero_refined_nolncRNA_bed["minus"] <- (AH1_Liver_Ribozero_refined_nolncRNA_bed[ ,2]-1000)
AH1_Liver_refined_nolncRNA_bed["minus"] <- (AH1_Liver_refined_nolncRNA_bed[ ,2]-1000)
AH1_PC_Ribozero_refined_nolncRNA_bed["minus"] <- (AH1_PC_Ribozero_refined_nolncRNA_bed[ ,2]-1000)
AH1_PC_refined_nolncRNA_bed["minus"] <- (AH1_PC_refined_nolncRNA_bed[ ,2]-1000)
AH2_Liver_Ribozero_refined_nolncRNA_bed["minus"] <- (AH2_Liver_Ribozero_refined_nolncRNA_bed[ ,2]-1000)
AH2_Liver_refined_nolncRNA_bed["minus"] <- (AH2_Liver_refined_nolncRNA_bed[ ,2]-1000)
AH2_PC_Ribozero_refined_nolncRNA_bed["minus"] <- (AH2_PC_Ribozero_refined_nolncRNA_bed[ ,2]-1000)
AH2_PC_refined_nolncRNA_bed["minus"] <- (AH2_PC_refined_nolncRNA_bed[ ,2]-1000)

AH1_Liver_Ribozero_refined_nolncRNA_bed["add"] <- (AH1_Liver_Ribozero_refined_nolncRNA_bed[ ,3]+1000)
AH1_Liver_refined_nolncRNA_bed["add"] <- (AH1_Liver_refined_nolncRNA_bed[ ,3]+1000)
AH1_PC_Ribozero_refined_nolncRNA_bed["add"] <- (AH1_PC_Ribozero_refined_nolncRNA_bed[ ,3]+1000)
AH1_PC_refined_nolncRNA_bed["add"] <- (AH1_PC_refined_nolncRNA_bed[ ,3]+1000)
AH2_Liver_Ribozero_refined_nolncRNA_bed["add"] <- (AH2_Liver_Ribozero_refined_nolncRNA_bed[ ,3]+1000)
AH2_Liver_refined_nolncRNA_bed["add"] <- (AH2_Liver_refined_nolncRNA_bed[ ,3]+1000)
AH2_PC_Ribozero_refined_nolncRNA_bed["add"] <- (AH2_PC_Ribozero_refined_nolncRNA_bed[ ,3]+1000)
AH2_PC_refined_nolncRNA_bed["add"] <- (AH2_PC_refined_nolncRNA_bed[ ,3]+1000)

#re-order
AH1_Liver_Ribozero_refined_nolncRNA_bed_extended <- AH1_Liver_Ribozero_refined_nolncRNA_bed[ ,
                                                    c("V1", "minus", "add", "V4", "V5", "V6", "V7", "V8", "V9", "V10")]
AH1_Liver_refined_nolncRNA_bed_extended <- AH1_Liver_refined_nolncRNA_bed[ ,
                                            c("V1", "minus", "add", "V4", "V5", "V6", "V7", "V8", "V9", "V10")]
AH1_PC_Ribozero_refined_nolncRNA_bed_extended <- AH1_PC_Ribozero_refined_nolncRNA_bed[ ,
                                                  c("V1", "minus", "add", "V4", "V5", "V6", "V7", "V8", "V9", "V10")]
AH1_PC_refined_nolncRNA_bed_extended <- AH1_PC_refined_nolncRNA_bed[ ,
                                         c("V1", "minus", "add", "V4", "V5", "V6", "V7", "V8", "V9", "V10")]
AH2_Liver_Ribozero_refined_nolncRNA_bed_extended <- AH2_Liver_Ribozero_refined_nolncRNA_bed[ ,
                                                     c("V1", "minus", "add", "V4", "V5", "V6", "V7", "V8", "V9", "V10")]
AH2_Liver_refined_nolncRNA_bed_extended <- AH2_Liver_refined_nolncRNA_bed[ ,
                                            c("V1", "minus", "add", "V4", "V5", "V6", "V7", "V8", "V9", "V10")]
AH2_PC_Ribozero_refined_nolncRNA_bed_extended <- AH2_PC_Ribozero_refined_nolncRNA_bed[ ,
                                                  c("V1", "minus", "add", "V4", "V5", "V6", "V7", "V8", "V9", "V10")]
AH2_PC_refined_nolncRNA_bed_extended <- AH2_PC_refined_nolncRNA_bed[ ,
                                         c("V1", "minus", "add", "V4", "V5", "V6", "V7", "V8", "V9", "V10")]

#re-naming "minus" to V2 and "add" to V3
names(AH1_Liver_Ribozero_refined_nolncRNA_bed_extended)[2] <- paste("V2")
names(AH1_Liver_refined_nolncRNA_bed_extended)[2] <- paste("V2")
names(AH1_PC_Ribozero_refined_nolncRNA_bed_extended)[2] <- paste("V2")
names(AH1_PC_refined_nolncRNA_bed_extended)[2] <- paste("V2")
names(AH2_Liver_Ribozero_refined_nolncRNA_bed_extended)[2] <- paste("V2")
names(AH2_Liver_refined_nolncRNA_bed_extended)[2] <- paste("V2")
names(AH2_PC_Ribozero_refined_nolncRNA_bed_extended)[2] <- paste("V2")
names(AH2_PC_refined_nolncRNA_bed_extended)[2] <- paste("V2")

names(AH1_Liver_Ribozero_refined_nolncRNA_bed_extended)[3] <- paste("V3")
names(AH1_Liver_refined_nolncRNA_bed_extended)[3] <- paste("V3")
names(AH1_PC_Ribozero_refined_nolncRNA_bed_extended)[3] <- paste("V3")
names(AH1_PC_refined_nolncRNA_bed_extended)[3] <- paste("V3")
names(AH2_Liver_Ribozero_refined_nolncRNA_bed_extended)[3] <- paste("V3")
names(AH2_Liver_refined_nolncRNA_bed_extended)[3] <- paste("V3")
names(AH2_PC_Ribozero_refined_nolncRNA_bed_extended)[3] <- paste("V3")
names(AH2_PC_refined_nolncRNA_bed_extended)[3] <- paste("V3")

#convert all negative numbers to 0
AH1_Liver_Ribozero_refined_nolncRNA_bed_extended$V2[AH1_Liver_Ribozero_refined_nolncRNA_bed_extended$V2 < 0] <- 0
AH1_Liver_refined_nolncRNA_bed_extended$V2[AH1_Liver_refined_nolncRNA_bed_extended$V2 < 0] <- 0
AH1_PC_Ribozero_refined_nolncRNA_bed_extended$V2[AH1_PC_Ribozero_refined_nolncRNA_bed_extended$V2 < 0] <- 0
AH1_PC_refined_nolncRNA_bed_extended$V2[AH1_PC_refined_nolncRNA_bed_extended$V2 < 0] <- 0
AH2_Liver_Ribozero_refined_nolncRNA_bed_extended$V2[AH2_Liver_Ribozero_refined_nolncRNA_bed_extended$V2 < 0] <- 0
AH2_Liver_refined_nolncRNA_bed_extended$V2[AH2_Liver_refined_nolncRNA_bed_extended$V2 < 0] <- 0
AH2_PC_Ribozero_refined_nolncRNA_bed_extended$V2[AH2_PC_Ribozero_refined_nolncRNA_bed_extended$V2 < 0] <- 0
AH2_PC_refined_nolncRNA_bed_extended$V2[AH2_PC_refined_nolncRNA_bed_extended$V2 < 0] <- 0

#order chr again
AH1_Liver_Ribozero_refined_nolncRNA_bed_extended <- AH1_Liver_Ribozero_refined_nolncRNA_bed_extended[with(AH1_Liver_Ribozero_refined_nolncRNA_bed_extended,
                                                                                                          order(V1, V2)), ]
AH1_Liver_refined_nolncRNA_bed_extended <- AH1_Liver_refined_nolncRNA_bed_extended[with(AH1_Liver_refined_nolncRNA_bed_extended, 
                                                                                        order(V1, V2)), ]
AH1_PC_Ribozero_refined_nolncRNA_bed_extended <- AH1_PC_Ribozero_refined_nolncRNA_bed_extended[with(AH1_PC_Ribozero_refined_nolncRNA_bed_extended,
                                                                                                    order(V1, V2)), ]
AH1_PC_refined_nolncRNA_bed_extended <- AH1_PC_refined_nolncRNA_bed_extended[with(AH1_PC_refined_nolncRNA_bed_extended,
                                                                                  order(V1, V2)), ]
AH2_Liver_Ribozero_refined_nolncRNA_bed_extended <- AH2_Liver_Ribozero_refined_nolncRNA_bed_extended[with(AH2_Liver_Ribozero_refined_nolncRNA_bed_extended,
                                                                                                          order(V1, V2)), ]
AH2_Liver_refined_nolncRNA_bed_extended <- AH2_Liver_refined_nolncRNA_bed_extended[with(AH2_Liver_refined_nolncRNA_bed_extended,
                                                                                        order(V1, V2)), ]
AH2_PC_Ribozero_refined_nolncRNA_bed_extended <- AH2_PC_Ribozero_refined_nolncRNA_bed_extended[with(AH2_PC_Ribozero_refined_nolncRNA_bed_extended,
                                                                                                    order(V1, V2)), ]
AH2_PC_refined_nolncRNA_bed_extended <- AH2_PC_refined_nolncRNA_bed_extended[with(AH2_PC_refined_nolncRNA_bed_extended,
                                                                                  order(V1, V2)), ]

write.table(AH1_Liver_Ribozero_refined_nolncRNA_bed_extended, "AH1_Liver_Ribozero_refined_nolncRNA_extended.bed", 
            row.names=FALSE, col.names=FALSE, quote=FALSE, sep="\t")
write.table(AH1_Liver_refined_nolncRNA_bed_extended, "AH1_Liver_refined_nolncRNA_extended.bed", 
            row.names=FALSE, col.names=FALSE, quote=FALSE, sep="\t")
write.table(AH1_PC_Ribozero_refined_nolncRNA_bed_extended, "AH1_PC_Ribozero_refined_nolncRNA_extended.bed",
            row.names=FALSE, col.names=FALSE, quote=FALSE, sep="\t")
write.table(AH1_PC_refined_nolncRNA_bed_extended, "AH1_PC_refined_nolncRNA_extended.bed",
            row.names=FALSE, col.names=FALSE, quote=FALSE, sep="\t")
write.table(AH2_Liver_Ribozero_refined_nolncRNA_bed_extended, "AH2_Liver_Ribozero_refined_nolncRNA_extended.bed",
            row.names=FALSE, col.names=FALSE, quote=FALSE, sep="\t")
write.table(AH2_Liver_refined_nolncRNA_bed_extended, "AH2_Liver_refined_nolncRNA_extended.bed",
            row.names=FALSE, col.names=FALSE, quote=FALSE, sep="\t")
write.table(AH2_PC_Ribozero_refined_nolncRNA_bed_extended, "AH2_PC_Ribozero_refined_nolncRNA_extended.bed",
            row.names=FALSE, col.names=FALSE, quote=FALSE, sep="\t")
write.table(AH2_PC_refined_nolncRNA_bed_extended, "AH2_PC_refined_nolncRNA_extended.bed",
            row.names=FALSE, col.names=FALSE, quote=FALSE, sep="\t")

#Isolate -1kb as 5' UTR
AH1_Liver_Ribozero_refined_nolncRNA_bed_5 <- AH1_Liver_Ribozero_refined_nolncRNA_bed[ ,
                                             c("V1","minus","V2","V4","V5","V6","V7","V8","V9","V10")]
AH1_Liver_refined_nolncRNA_bed_5 <- AH1_Liver_refined_nolncRNA_bed[ ,
                                             c("V1","minus","V2","V4","V5","V6","V7","V8","V9","V10")]
AH1_PC_Ribozero_refined_nolncRNA_bed_5 <- AH1_PC_Ribozero_refined_nolncRNA_bed[ ,
                                           c("V1","minus","V2","V4","V5","V6","V7","V8","V9","V10")]
AH1_PC_refined_nolncRNA_bed_5 <- AH1_PC_refined_nolncRNA_bed[ ,
                                 c("V1","minus","V2","V4","V5","V6","V7","V8","V9","V10")]
AH2_Liver_Ribozero_refined_nolncRNA_bed_5 <- AH2_Liver_Ribozero_refined_nolncRNA_bed[ ,
                                             c("V1","minus","V2","V4","V5","V6","V7","V8","V9","V10")]
AH2_Liver_refined_nolncRNA_bed_5 <- AH2_Liver_refined_nolncRNA_bed[ ,
                                    c("V1","minus","V2","V4","V5","V6","V7","V8","V9","V10")]
AH2_PC_Ribozero_refined_nolncRNA_bed_5 <- AH2_PC_Ribozero_refined_nolncRNA_bed[ ,
                                           c("V1","minus","V2","V4","V5","V6","V7","V8","V9","V10")]
AH2_PC_refined_nolncRNA_bed_5 <- AH2_PC_refined_nolncRNA_bed[ ,
                                 c("V1","minus","V2","V4","V5","V6","V7","V8","V9","V10")]
#Make all negatives, 0
AH1_Liver_Ribozero_refined_nolncRNA_bed_5$minus[AH1_Liver_Ribozero_refined_nolncRNA_bed_5$minus < 0] <- 0
AH1_Liver_refined_nolncRNA_bed_5$minus[AH1_Liver_refined_nolncRNA_bed_5$minus < 0] <- 0
AH1_PC_Ribozero_refined_nolncRNA_bed_5$minus[AH1_PC_Ribozero_refined_nolncRNA_bed_5$minus < 0] <- 0
AH1_PC_refined_nolncRNA_bed_5$minus[AH1_PC_refined_nolncRNA_bed_5$minus < 0] <- 0
AH2_Liver_Ribozero_refined_nolncRNA_bed_5$minus[AH2_Liver_Ribozero_refined_nolncRNA_bed_5$minus < 0] <- 0
AH2_Liver_refined_nolncRNA_bed_5$minus[AH2_Liver_refined_nolncRNA_bed_5$minus < 0] <- 0
AH2_PC_Ribozero_refined_nolncRNA_bed_5$minus[AH2_PC_Ribozero_refined_nolncRNA_bed_5$minus < 0] <- 0
AH2_PC_refined_nolncRNA_bed_5$minus[AH2_PC_refined_nolncRNA_bed_5$minus < 0] <- 0

#order chr
AH1_Liver_Ribozero_refined_nolncRNA_bed_5 <- AH1_Liver_Ribozero_refined_nolncRNA_bed_5[with(
                                             AH1_Liver_Ribozero_refined_nolncRNA_bed_5,order(V1, minus)), ]
AH1_Liver_refined_nolncRNA_bed_5 <- AH1_Liver_refined_nolncRNA_bed_5[with(
                                    AH1_Liver_refined_nolncRNA_bed_5,order(V1, minus)), ]
AH1_PC_Ribozero_refined_nolncRNA_bed_5 <- AH1_PC_Ribozero_refined_nolncRNA_bed_5[with(
                                          AH1_PC_Ribozero_refined_nolncRNA_bed_5,order(V1, minus)), ]
AH1_PC_refined_nolncRNA_bed_5 <- AH1_PC_refined_nolncRNA_bed_5[with(
                                 AH1_PC_refined_nolncRNA_bed_5,order(V1, minus)), ]
AH2_Liver_Ribozero_refined_nolncRNA_bed_5 <- AH2_Liver_Ribozero_refined_nolncRNA_bed_5[with(
                                             AH2_Liver_Ribozero_refined_nolncRNA_bed_5,order(V1, minus)), ]
AH2_Liver_refined_nolncRNA_bed_5 <- AH2_Liver_refined_nolncRNA_bed_5[with(
                                    AH2_Liver_refined_nolncRNA_bed_5,order(V1, minus)), ]
AH2_PC_Ribozero_refined_nolncRNA_bed_5 <- AH2_PC_Ribozero_refined_nolncRNA_bed_5[with(
                                          AH2_PC_Ribozero_refined_nolncRNA_bed_5,order(V1, minus)), ]
AH2_PC_refined_nolncRNA_bed_5 <- AH2_PC_refined_nolncRNA_bed_5[with(
                                 AH2_PC_refined_nolncRNA_bed_5,order(V1, minus)), ]

write.table(AH1_Liver_Ribozero_refined_nolncRNA_bed_5, "AH1_Liver_Ribozero_refined_nolncRNA_5.bed", 
            row.names=F, col.names=F, sep = "\t")
write.table(AH1_Liver_refined_nolncRNA_bed_5, "AH1_Liver_refined_nolncRNA_5.bed", 
            row.names=F, col.names=F, sep = "\t")
write.table(AH1_PC_Ribozero_refined_nolncRNA_bed_5, "AH1_PC_Ribozero_refined_nolncRNA_5.bed",
            row.names=F, col.names=F, sep = "\t")
write.table(AH1_PC_refined_nolncRNA_bed_5, "AH1_PC_refined_nolncRNA_5.bed",
            row.names=F, col.names=F, sep = "\t")
write.table(AH2_Liver_Ribozero_refined_nolncRNA_bed_5, "AH2_Liver_Ribozero_refined_nolncRNA_5.bed",
            row.names=F, col.names=F, sep = "\t")
write.table(AH2_Liver_refined_nolncRNA_bed_5, "AH2_Liver_refined_nolncRNA_5.bed",
            row.names=F, col.names=F, sep = "\t")
write.table(AH2_PC_Ribozero_refined_nolncRNA_bed_5, "AH2_PC_Ribozero_refined_nolncRNA_5.bed",
            row.names=F, col.names=F, sep = "\t")
write.table(AH2_PC_refined_nolncRNA_bed_5, "AH2_PC_refined_nolncRNA_5.bed",
            row.names=F, col.names=F, sep = "\t")

#Isolate +1kb as 3' UTR
AH1_Liver_Ribozero_refined_nolncRNA_bed_3 <- AH1_Liver_Ribozero_refined_nolncRNA_bed[ ,c(
                                             "V1", "V3", "add", "V4", "V5","V6","V7","V8","V9","V10")]
AH1_Liver_refined_nolncRNA_bed_3 <- AH1_Liver_refined_nolncRNA_bed[ ,c(
                                    "V1", "V3", "add", "V4", "V5","V6","V7","V8","V9","V10")]
AH1_PC_Ribozero_refined_nolncRNA_bed_3 <- AH1_PC_Ribozero_refined_nolncRNA_bed[ ,c(
                                          "V1", "V3", "add", "V4", "V5","V6","V7","V8","V9","V10")]
AH1_PC_refined_nolncRNA_bed_3 <- AH1_PC_refined_nolncRNA_bed[ ,c(
                                 "V1", "V3", "add", "V4", "V5","V6","V7","V8","V9","V10")]
AH2_Liver_Ribozero_refined_nolncRNA_bed_3 <- AH2_Liver_Ribozero_refined_nolncRNA_bed[ ,c(
                                             "V1", "V3", "add", "V4", "V5","V6","V7","V8","V9","V10")]
AH2_Liver_refined_nolncRNA_bed_3 <- AH2_Liver_refined_nolncRNA_bed[ ,c(
                                    "V1", "V3", "add", "V4", "V5","V6","V7","V8","V9","V10")]
AH2_PC_Ribozero_refined_nolncRNA_bed_3 <- AH2_PC_Ribozero_refined_nolncRNA_bed[ ,c(
                                         "V1", "V3", "add", "V4", "V5","V6","V7","V8","V9","V10")]
AH2_PC_refined_nolncRNA_bed_3 <- AH2_PC_refined_nolncRNA_bed[ ,c(
                                 "V1", "V3", "add", "V4", "V5","V6","V7","V8","V9","V10")]
#order Chr
AH1_Liver_Ribozero_refined_nolncRNA_bed_3 <- AH1_Liver_Ribozero_refined_nolncRNA_bed_3[with(
  AH1_Liver_Ribozero_refined_nolncRNA_bed_3, order(V1, V3)), ]
AH1_Liver_refined_nolncRNA_bed_3 <- AH1_Liver_refined_nolncRNA_bed_3[with(
  AH1_Liver_refined_nolncRNA_bed_3, order(V1, V3)), ]
AH1_PC_Ribozero_refined_nolncRNA_bed_3 <- AH1_PC_Ribozero_refined_nolncRNA_bed_3[with(
  AH1_PC_Ribozero_refined_nolncRNA_bed_3, order(V1, V3)), ]
AH1_PC_refined_nolncRNA_bed_3 <- AH1_PC_refined_nolncRNA_bed_3[with(
  AH1_PC_refined_nolncRNA_bed_3, order(V1, V3)), ]
AH2_Liver_Ribozero_refined_nolncRNA_bed_3 <- AH2_Liver_Ribozero_refined_nolncRNA_bed_3[with(
  AH2_Liver_Ribozero_refined_nolncRNA_bed_3, order(V1, V3)), ]
AH2_Liver_refined_nolncRNA_bed_3 <- AH2_Liver_refined_nolncRNA_bed_3[with(
  AH2_Liver_refined_nolncRNA_bed_3, order(V1, V3)), ]
AH2_PC_Ribozero_refined_nolncRNA_bed_3 <- AH2_PC_Ribozero_refined_nolncRNA_bed_3[with(
  AH2_PC_Ribozero_refined_nolncRNA_bed_3, order(V1, V3)), ]
AH2_PC_refined_nolncRNA_bed_3 <- AH2_PC_refined_nolncRNA_bed_3[with(
  AH2_PC_refined_nolncRNA_bed_3, order(V1, V3)), ]

write.table(AH1_Liver_Ribozero_refined_nolncRNA_bed_3, "AH1_Liver_Ribozero_refined_nolncRNA_3.bed", 
            row.names = F, col.names = F, sep = "\t")
write.table(AH1_Liver_refined_nolncRNA_bed_3, "AH1_Liver_refined_nolncRNA_3.bed", 
            row.names = F, col.names = F, sep = "\t")
write.table(AH1_PC_Ribozero_refined_nolncRNA_bed_3, "AH1_PC_Ribozero_refined_nolncRNA_3.bed",
            row.names = F, col.names = F, sep = "\t")
write.table(AH1_PC_refined_nolncRNA_bed_3, "AH1_PC_refined_nolncRNA_3.bed",
            row.names = F, col.names = F, sep = "\t")
write.table(AH2_Liver_Ribozero_refined_nolncRNA_bed_3, "AH2_Liver_Ribozero_refined_nolncRNA_3.bed",
            row.names = F, col.names = F, sep = "\t")
write.table(AH2_Liver_refined_nolncRNA_bed_3, "AH2_Liver_refined_nolncRNA_3.bed",
            row.names = F, col.names = F, sep = "\t")
write.table(AH2_PC_Ribozero_refined_nolncRNA_bed_3, "AH2_PC_Ribozero_refined_nolncRNA_3.bed",
            row.names = F, col.names = F, sep = "\t")
write.table(AH2_PC_refined_nolncRNA_bed_3, "AH2_PC_refined_nolncRNA_3.bed",
            row.names = F, col.names = F, sep = "\t")

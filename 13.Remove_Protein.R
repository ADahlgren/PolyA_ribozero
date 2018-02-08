#R script to remove protein coding sequences using blastp, hmmersearch, and blastn files

setwd("~/research/polyA_ribozero/Filter3")

require(tidyr)
require(dplyr)
require(stringr)

#need to download blastp product files (end with sprot.outfmt6)

AH1_Liver_Ribozero_blastp <- read.table("683610_Liver_Ribozero_sprot.outfmt6", 
                                           header=FALSE, stringsAsFactors = FALSE, sep = "\t")
AH1_Liver_blastp <- read.table("683610_Liver_sprot.outfmt6", 
                                  header=FALSE, stringsAsFactors = FALSE, sep = "\t")
AH1_PC_Ribozero_blastp <- read.table("683610_Parietal_Cortex_Ribozero_sprot.outfmt6", 
                                        header=FALSE, stringsAsFactors = FALSE, sep = "\t")
AH1_PC_blastp <- read.table("683610_Parietal_Cortex_sprot.outfmt6", 
                               header=FALSE, stringsAsFactors = FALSE, sep = "\t")
AH2_Liver_Ribozero_blastp <- read.table("686521_Liver_Ribozero_sprot.outfmt6", 
                                           header=FALSE, stringsAsFactors = FALSE, sep = "\t")
AH2_Liver_blastp <- read.table("686521_Liver_sprot.outfmt6", 
                                  header=FALSE, stringsAsFactors = FALSE, sep = "\t")
AH2_PC_Ribozero_blastp <- read.table("686521_Parietal_Cortex_Ribozero_sprot.outfmt6", 
                                        header=FALSE, stringsAsFactors = FALSE, sep = "\t")
AH2_PC_blastp <- read.table("686521_Parietal_Cortex_sprot.outfmt6", 
                               header=FALSE, stringsAsFactors = FALSE, sep = "\t")

#parse out column V1 to get gene IDs
AH1_Liver_Ribozero_blastp_temp = data.frame(str_split_fixed(AH1_Liver_Ribozero_blastp$V1, ":", 4))
AH1_Liver_Ribozero_blastp_bed <- cbind(AH1_Liver_Ribozero_blastp_temp$X3, AH1_Liver_Ribozero_blastp[,c("V2","V11")])
names(AH1_Liver_Ribozero_blastp_bed) = c("V4", "target", "sig")

AH1_Liver_blastp_temp = data.frame(str_split_fixed(AH1_Liver_blastp$V1, ":", 4))
AH1_Liver_blastp_bed <- cbind(AH1_Liver_blastp_temp$X3, AH1_Liver_blastp[,c("V2","V11")])
names(AH1_Liver_blastp_bed) = c("V4", "target", "sig")

AH1_PC_Ribozero_blastp_temp = data.frame(str_split_fixed(AH1_PC_Ribozero_blastp$V1, ":", 4))
AH1_PC_Ribozero_blastp_bed <- cbind(AH1_PC_Ribozero_blastp_temp$X3, AH1_PC_Ribozero_blastp[,c("V2","V11")])
names(AH1_PC_Ribozero_blastp_bed) = c("V4", "target", "sig")

AH1_PC_blastp_temp = data.frame(str_split_fixed(AH1_PC_blastp$V1, ":", 4))
AH1_PC_blastp_bed <- cbind(AH1_PC_blastp_temp$X3, AH1_PC_blastp[,c("V2","V11")])
names(AH1_PC_blastp_bed) = c("V4", "target", "sig")

AH2_Liver_Ribozero_blastp_temp = data.frame(str_split_fixed(AH2_Liver_Ribozero_blastp$V1, ":", 4))
AH2_Liver_Ribozero_blastp_bed <- cbind(AH2_Liver_Ribozero_blastp_temp$X3, AH2_Liver_Ribozero_blastp[,c("V2","V11")])
names(AH2_Liver_Ribozero_blastp_bed) = c("V4", "target", "sig")

AH2_Liver_blastp_temp = data.frame(str_split_fixed(AH2_Liver_blastp$V1, ":", 4))
AH2_Liver_blastp_bed <- cbind(AH2_Liver_blastp_temp$X3, AH2_Liver_blastp[,c("V2","V11")])
names(AH2_Liver_blastp_bed) = c("V4", "target", "sig")

AH2_PC_Ribozero_blastp_temp = data.frame(str_split_fixed(AH2_PC_Ribozero_blastp$V1, ":", 4))
AH2_PC_Ribozero_blastp_bed <- cbind(AH2_PC_Ribozero_blastp_temp$X3, AH2_PC_Ribozero_blastp[,c("V2","V11")])
names(AH2_PC_Ribozero_blastp_bed) = c("V4", "target", "sig")

AH2_PC_blastp_temp = data.frame(str_split_fixed(AH2_PC_blastp$V1, ":", 4))
AH2_PC_blastp_bed <- cbind(AH2_PC_blastp_temp$X3, AH2_PC_blastp[,c("V2","V11")])
names(AH2_PC_blastp_bed) = c("V4", "target", "sig")

#write blastp tables into bed files (comparable format with hmmersearh results) 
write.table(AH1_Liver_Ribozero_blastp_bed, "AH1_Liver_Ribozero_blastp.bed", row.names = FALSE, 
            col.names = FALSE, quote = FALSE, sep = "\t")
write.table(AH1_Liver_blastp_bed, "AH1_Liver_blastp.bed", row.names = FALSE, 
            col.names = FALSE, quote = FALSE, sep = "\t")
write.table(AH1_PC_Ribozero_blastp_bed, "AH1_PC_Ribozero_blastp.bed", row.names = FALSE, 
            col.names = FALSE, quote = FALSE, sep = "\t")
write.table(AH1_PC_blastp_bed, "AH1_PC_blastp.bed", row.names = FALSE, 
            col.names = FALSE, quote = FALSE, sep = "\t")
write.table(AH2_Liver_Ribozero_blastp_bed, "AH2_Liver_Ribozero_blastp.bed", row.names = FALSE, 
            col.names = FALSE, quote = FALSE, sep = "\t")
write.table(AH2_Liver_blastp_bed, "AH2_Liver_blastp.bed", row.names = FALSE, 
            col.names = FALSE, quote = FALSE, sep = "\t")
write.table(AH2_PC_Ribozero_blastp_bed, "AH2_PC_Ribozero_blastp.bed", row.names = FALSE, 
            col.names = FALSE, quote = FALSE, sep = "\t")
write.table(AH2_PC_blastp_bed, "AH2_PC_blastp.bed", row.names = FALSE, 
            col.names = FALSE, quote = FALSE, sep = "\t")

##################################################################################################################

#upload hmmersearch files (end with pfam.tblout)
AH1_Liver_Ribozero_pfam <- read.table("683610_Liver_Ribozero_pfam.tblout", header = FALSE, stringsAsFactors = FALSE)
AH1_Liver_pfam <- read.table("683610_Liver_pfam.tblout", header = FALSE, stringsAsFactors = FALSE)
AH1_PC_Ribozero_pfam <- read.table("683610_Parietal_Cortex_Ribozero_pfam.tblout", header = FALSE, 
                                   stringsAsFactors = FALSE)
AH1_PC_pfam <- read.table("683610_Parietal_Cortex_pfam.tblout", header = FALSE, stringsAsFactors = FALSE)
AH2_Liver_Ribozero_pfam <- read.table("686521_Liver_Ribozero_pfam.tblout", header = FALSE, stringsAsFactors = FALSE)
AH2_Liver_pfam <- read.table("686521_Liver_pfam.tblout", header = FALSE, stringsAsFactors = FALSE)
AH2_PC_Ribozero_pfam <- read.table("686521_Parietal_Cortex_Ribozero_pfam.tblout", header = FALSE, 
                                   stringsAsFactors = FALSE)
AH2_PC_pfam <- read.table("686521_Parietal_Cortex_pfam.tblout", header = FALSE, stringsAsFactors = FALSE)

#manipulating .tblout to be compatible with blastp
AH1_Liver_Ribozero_pfam_temp=data.frame(str_split_fixed(AH1_Liver_Ribozero_pfam$V1, ":", 4))
AH1_Liver_Ribozero_pfam_bed <- cbind(AH1_Liver_Ribozero_pfam_temp$X3, AH1_Liver_Ribozero_pfam[,c("V3","V5")])
names(AH1_Liver_Ribozero_pfam_bed)=c("V4", "target", "sig")

AH1_Liver_pfam_temp=data.frame(str_split_fixed(AH1_Liver_pfam$V1, ":", 4))
AH1_Liver_pfam_bed <- cbind(AH1_Liver_pfam_temp$X3, AH1_Liver_pfam[,c("V3","V5")])
names(AH1_Liver_pfam_bed)=c("V4", "target", "sig")

AH1_PC_Ribozero_pfam_temp=data.frame(str_split_fixed(AH1_PC_Ribozero_pfam$V1, ":", 4))
AH1_PC_Ribozero_pfam_bed <- cbind(AH1_PC_Ribozero_pfam_temp$X3, AH1_PC_Ribozero_pfam[,c("V3","V5")])
names(AH1_PC_Ribozero_pfam_bed)=c("V4", "target", "sig")

AH1_PC_pfam_temp=data.frame(str_split_fixed(AH1_PC_pfam$V1, ":", 4))
AH1_PC_pfam_bed <- cbind(AH1_PC_pfam_temp$X3, AH1_PC_pfam[,c("V3","V5")])
names(AH1_PC_pfam_bed)=c("V4", "target", "sig")

AH2_Liver_Ribozero_pfam_temp=data.frame(str_split_fixed(AH2_Liver_Ribozero_pfam$V1, ":", 4))
AH2_Liver_Ribozero_pfam_bed <- cbind(AH2_Liver_Ribozero_pfam_temp$X3, AH2_Liver_Ribozero_pfam[,c("V3","V5")])
names(AH2_Liver_Ribozero_pfam_bed)=c("V4", "target", "sig")

AH2_Liver_pfam_temp=data.frame(str_split_fixed(AH2_Liver_pfam$V1, ":", 4))
AH2_Liver_pfam_bed <- cbind(AH2_Liver_pfam_temp$X3, AH2_Liver_pfam[,c("V3","V5")])
names(AH2_Liver_pfam_bed)=c("V4", "target", "sig")

AH2_PC_Ribozero_pfam_temp=data.frame(str_split_fixed(AH2_PC_Ribozero_pfam$V1, ":", 4))
AH2_PC_Ribozero_pfam_bed <- cbind(AH2_PC_Ribozero_pfam_temp$X3, AH2_PC_Ribozero_pfam[,c("V3","V5")])
names(AH2_PC_Ribozero_pfam_bed)=c("V4", "target", "sig")

AH2_PC_pfam_temp=data.frame(str_split_fixed(AH2_PC_pfam$V1, ":", 4))
AH2_PC_pfam_bed <- cbind(AH2_PC_pfam_temp$X3, AH2_PC_pfam[,c("V3","V5")])
names(AH2_PC_pfam_bed)=c("V4", "target", "sig")

#output rearranged table with tblout info
write.table(AH1_Liver_Ribozero_pfam_bed, "AH1_Liver_Ribozero_pfam.bed", 
            row.names = FALSE, col.names = FALSE, quote = FALSE, sep = "\t")
write.table(AH1_Liver_pfam_bed, "AH1_Liver_pfam.bed", 
            row.names = FALSE, col.names = FALSE, quote = FALSE, sep = "\t")
write.table(AH1_PC_Ribozero_pfam_bed, "AH1_PC_Ribozero_pfam.bed", 
            row.names = FALSE, col.names = FALSE, quote = FALSE, sep = "\t")
write.table(AH1_PC_pfam_bed, "AH1_PC_pfam.bed", 
            row.names = FALSE, col.names = FALSE, quote = FALSE, sep = "\t")
write.table(AH2_Liver_Ribozero_pfam_bed, "AH2_Liver_Ribozero_pfam.bed",
            row.names = FALSE, col.names = FALSE, quote = FALSE, sep = "\t")
write.table(AH2_Liver_pfam_bed, "AH2_Liver_pfam.bed",
            row.names = FALSE, col.names = FALSE, quote = FALSE, sep = "\t")
write.table(AH2_PC_Ribozero_pfam_bed, "AH2_PC_Ribozero_pfam.bed",
            row.names = FALSE, col.names = FALSE, quote = FALSE, sep = "\t")
write.table(AH2_PC_pfam_bed, "AH2_PC_pfam.bed",
            row.names = FALSE, col.names = FALSE, quote = FALSE, sep = "\t")

####################################################################################################

#download cdna.outfmt6 files from blastn
AH1_Liver_Ribozero_blastn <- read.table("683610_Liver_Ribozero_S51_f1_hg38_cdna.outfmt6", 
                                        header = FALSE, stringsAsFactors = FALSE, sep = "\t")
AH1_Liver_blastn <- read.table("683610_Liver_S61_f1_hg38_cdna.outfmt6", 
                                header = FALSE, stringsAsFactors = FALSE, sep = "\t")
AH1_PC_Ribozero_blastn <- read.table("683610_Parietal_Cortex_Ribozero_S53_f1_hg38_cdna.outfmt6", 
                                      header = FALSE, stringsAsFactors = FALSE, sep = "\t")
AH1_PC_blastn <- read.table("683610_Parietal_Cortex_S55_f1_hg38_cdna.outfmt6", 
                             header = FALSE, stringsAsFactors = FALSE, sep = "\t")
AH2_Liver_Ribozero_blastn <- read.table("686521_Liver_Ribozero_S50_f1_hg38_cdna.outfmt6",
                                         header = FALSE, stringsAsFactors = FALSE, sep = "\t")
AH2_Liver_blastn <- read.table("686521_Liver_S60_f1_hg38_cdna.outfmt6",
                                header = FALSE, stringsAsFactors = FALSE, sep = "\t")
AH2_PC_Ribozero_blastn <- read.table("686521_Parietal_Cortex_Ribozero_S52_f1_hg38_cdna.outfmt6",
                                      header = FALSE, stringsAsFactors = FALSE, sep = "\t")
AH2_PC_blastn <- read.table("686521_Parietal_Cortex_S54_f1_hg38_cdna.outfmt6",
                             header = FALSE, stringsAsFactors = FALSE, sep = "\t")

#select column V1 to get gene IDs
AH1_Liver_Ribozero_blastn_temp = data.frame(str_split_fixed(AH1_Liver_Ribozero_blastn$V1, ":", 4))
AH1_Liver_Ribozero_blastn_bed <- cbind(AH1_Liver_Ribozero_blastn_temp$X3, AH1_Liver_Ribozero_blastn[,c("V2", "V11")])
names(AH1_Liver_Ribozero_blastn_bed)=c("V4", "target", "sig")

AH1_Liver_blastn_temp = data.frame(str_split_fixed(AH1_Liver_blastn$V1, ":", 4))
AH1_Liver_blastn_bed <- cbind(AH1_Liver_blastn_temp$X3, AH1_Liver_blastn[,c("V2", "V11")])
names(AH1_Liver_blastn_bed)=c("V4", "target", "sig")

AH1_PC_Ribozero_blastn_temp = data.frame(str_split_fixed(AH1_PC_Ribozero_blastn$V1, ":", 4))
AH1_PC_Ribozero_blastn_bed <- cbind(AH1_PC_Ribozero_blastn_temp$X3, AH1_PC_Ribozero_blastn[,c("V2", "V11")])
names(AH1_PC_Ribozero_blastn_bed)=c("V4", "target", "sig")

AH1_PC_blastn_temp = data.frame(str_split_fixed(AH1_PC_blastn$V1, ":", 4))
AH1_PC_blastn_bed <- cbind(AH1_PC_blastn_temp$X3, AH1_PC_blastn[,c("V2", "V11")])
names(AH1_PC_blastn_bed)=c("V4", "target", "sig")

AH2_Liver_Ribozero_blastn_temp = data.frame(str_split_fixed(AH2_Liver_Ribozero_blastn$V1, ":", 4))
AH2_Liver_Ribozero_blastn_bed <- cbind(AH2_Liver_Ribozero_blastn_temp$X3, AH2_Liver_Ribozero_blastn[,c("V2", "V11")])
names(AH2_Liver_Ribozero_blastn_bed)=c("V4", "target", "sig")

AH2_Liver_blastn_temp = data.frame(str_split_fixed(AH2_Liver_blastn$V1, ":", 4))
AH2_Liver_blastn_bed <- cbind(AH2_Liver_blastn_temp$X3, AH2_Liver_blastn[,c("V2", "V11")])
names(AH2_Liver_blastn_bed)=c("V4", "target", "sig")

AH2_PC_Ribozero_blastn_temp = data.frame(str_split_fixed(AH2_PC_Ribozero_blastn$V1, ":", 4))
AH2_PC_Ribozero_blastn_bed <- cbind(AH2_PC_Ribozero_blastn_temp$X3, AH2_PC_Ribozero_blastn[,c("V2", "V11")])
names(AH2_PC_Ribozero_blastn_bed)=c("V4", "target", "sig")

AH2_PC_blastn_temp = data.frame(str_split_fixed(AH2_PC_blastn$V1, ":", 4))
AH2_PC_blastn_bed <- cbind(AH2_PC_blastn_temp$X3, AH2_PC_blastn[,c("V2", "V11")])
names(AH2_PC_blastn_bed)=c("V4", "target", "sig")

#output in bed format
write.table(AH1_Liver_Ribozero_blastn_bed, "AH1_Liver_Ribozero_blastn.bed", row.names = FALSE, 
            col.names = FALSE, quote = FALSE, sep = "\t" )
write.table(AH1_Liver_blastn_bed, "AH1_Liver_blastn.bed", row.names = FALSE, 
            col.names = FALSE, quote = FALSE, sep = "\t" )
write.table(AH1_PC_Ribozero_blastn_bed, "AH1_PC_Ribozero_blastn.bed", row.names = FALSE, 
            col.names = FALSE, quote = FALSE, sep = "\t" )
write.table(AH1_PC_blastn_bed, "AH1_PC_blastn.bed", row.names = FALSE, 
            col.names = FALSE, quote = FALSE, sep = "\t" )
write.table(AH2_Liver_Ribozero_blastn_bed, "AH2_Liver_Ribozero_blastn.bed", row.names = FALSE,
            col.names = FALSE, quote = FALSE, sep = "\t" )
write.table(AH2_Liver_blastn_bed, "AH2_Liver_blastn.bed", row.names = FALSE,
            col.names = FALSE, quote = FALSE, sep = "\t" )
write.table(AH2_PC_Ribozero_blastn_bed, "AH2_PC_Ribozero_blastn.bed", row.names = FALSE,
            col.names = FALSE, quote = FALSE, sep = "\t" )
write.table(AH2_PC_blastn_bed, "AH2_PC_blastn.bed", row.names = FALSE,
            col.names = FALSE, quote = FALSE, sep = "\t" )

##########################################################################################################

#use filter12.bed files next

AH1_Liver_Ribozero_bed_f12 = read.table("683610_Liver_Ribozero_S51_filter12.bed", header = F, 
                                    colClasses = "character", sep = "\t")
AH1_Liver_bed_f12 = read.table("683610_Liver_S61_filter12.bed", header = F, 
                           colClasses = "character", sep = "\t")
AH1_PC_Ribozero_bed_f12 = read.table("683610_Parietal_Cortex_Ribozero_S53_filter12.bed", header = F, 
                                 colClasses = "character", sep = "\t")
AH1_PC_bed_f12 = read.table("683610_Parietal_Cortex_S55_filter12.bed", header = F, 
                            colClasses = "character", sep = "\t")
AH2_Liver_Ribozero_bed_f12 = read.table("686521_Liver_Ribozero_S50_filter12.bed", header = F, 
                                        colClasses = "character", sep = "\t")
AH2_Liver_bed_f12 = read.table("686521_Liver_S60_filter12.bed", header = F, 
                               colClasses = "character", sep = "\t")
AH2_PC_Ribozero_bed_f12 = read.table("686521_Parietal_Cortex_Ribozero_S52_filter12.bed", header = F, 
                                     colClasses = "character", sep = "\t")
AH2_PC_bed_f12 = read.table("686521_Parietal_Cortex_S54_filter12.bed", header = F, 
                            colClasses = "character", sep = "\t")
                            
##############################################################################################################                            

#There is a combine step in Erica's script. However, this is a little different, so I am skipping this step

#merge protein findings (blastp + hmmer pfam + blastn)
AH1_Liver_Ribozero_protein <- rbind(AH1_Liver_Ribozero_blastp_bed, AH1_Liver_Ribozero_pfam_bed, 
                                    AH1_Liver_Ribozero_blastn_bed)
AH1_Liver_protein <- rbind(AH1_Liver_blastp_bed, AH1_Liver_pfam_bed, AH1_Liver_blastn_bed)
AH1_PC_Ribozero_protein <- rbind(AH1_PC_Ribozero_blastp_bed, AH1_PC_Ribozero_pfam_bed,
                                 AH1_PC_Ribozero_blastn_bed)
AH1_PC_protein <- rbind(AH1_PC_blastp_bed, AH1_PC_pfam_bed, AH1_PC_blastn_bed)
AH2_Liver_Ribozero_protein <- rbind(AH2_Liver_Ribozero_blastp_bed, AH2_Liver_Ribozero_pfam_bed,
                                    AH2_Liver_Ribozero_blastn_bed)
AH2_Liver_protein <- rbind(AH2_Liver_blastp_bed, AH2_Liver_pfam_bed, AH2_Liver_blastn_bed)
AH2_PC_Ribozero_protein <- rbind(AH2_PC_Ribozero_blastp_bed, AH2_PC_Ribozero_pfam_bed,
                                 AH2_PC_Ribozero_blastn_bed)
AH2_PC_protein <- rbind(AH2_PC_blastp_bed, AH2_PC_pfam_bed, AH2_PC_blastn_bed)

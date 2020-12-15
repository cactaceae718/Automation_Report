# Oncomine HEME create sample list on IR
# Oct. 27, 2020 (last modification)
# Yi Cai

library(readxl)
library(stringr)

args<-commandArgs(T)
run_id <- args[1]
target_dir <- args[2]

dir_vcf=paste0(target_dir, run_id,  sep ="/")

dir_sheet <- "/mnt/Z_drive/acc_pathology/molecular/MOLECULAR LAB ONLY/Validations/HEME Molecular/Heme NGS Panel/Oncomine Myeloid/Run Data/Worksheets-WETLAB/"

onc_1 <- as.data.frame(na.omit(read_excel(paste0(dir_sheet, run_id, ".xlsm"), sheet = "DNA" ,col_names = T,col_types = "text", range = cell_cols("A:D"))))
onc_2 <- as.data.frame(na.omit(read_excel(paste0(dir_sheet, run_id, ".xlsm"), sheet = "RNA" ,col_names = T,col_types = "text", range = cell_cols("A:D"))))

colnames(onc_1) <- onc_1[1,]
colnames(onc_2) <- onc_2[1,]

onc_1 <- onc_1[-1,]
onc_2 <- onc_2[-1,]

colnames(onc_1)[3] <- "DNA/RNA"
colnames(onc_2)[3] <- "DNA/RNA"

onc <- rbind(onc_1, onc_2)

onc$sample_ID <- paste0(onc[,2], "-", onc[,3])

write.table(onc[,c(5,4,1)], file = paste0(dir_vcf, run_id, "_samplelist_V2.txt"), row.names = F, col.names = F, quote = F)

cat(paste("\n finish samplelist generation  \n\n start calling APIs"))


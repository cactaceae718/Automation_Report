# Oncomine HEME reporting pipeline by extracting/parsing information from files called by API
# Nov. 11, 2020 (last modification)
# Yi Cai

library(readxl)
library(stringr)
library(stringi)
library(dplyr)
library(gtools)
library(writexl)

args<-commandArgs(T)
run_id <- args[1]
target_dir <- args[2]

out_dir <- paste0(target_dir, run_id , "/Report", sep ="/")
dir_zip <- paste0(target_dir, run_id , "/Variants")

##### Step 1. Parse to get information for filterset
path_files <- list.dirs(dir_zip)[-1]
files_maf <- list.files(path = path_files, pattern = "full.tsv")
files_onc <- list.files(path = path_files, pattern = "oncomine.tsv")[grepl("_Non-Filtered_", list.files(path = path_files, pattern = "oncomine.tsv"))]

combine_tsv <- data.frame()
report <- data.frame()
for (s in 1:length(path_files)) {
  sample_name <- unlist(str_split(path_files[s], "[/_]"))[14]
  
  tsv_maf <- NULL
  tsv <- NULL
  ## read in table with information of annotated MAF(ExAc Database)
  tsv_maf <- read.table(paste0(path_files[s], "/", files_maf[s]),  header = T, sep = "\t", stringsAsFactors = F, na.strings = c("", "NA", " "), skip =2, comment.char = "")
  tsv_maf_ <- tsv_maf[,c("X..locus", "gene", "type", "genotype", "maf", "ref", "length", "coverage", "allele_ratio", "X._frequency", "exac")]
  tsv_maf_$Sample <- rep(sample_name, nrow(tsv_maf))
  # to subset SNV from full tsv file containing maf information
  tsv_maf_ <- subset(tsv_maf_, grepl("SNV", tsv_maf_$type) | grepl("INDEL", tsv_maf_$type) | grepl("MNV", tsv_maf_$type) | grepl("FLT3ITD", tsv_maf_$type))
  tsv_maf_ <- tsv_maf_[order(tsv_maf_$X..locus),]
  
  ## parse maf number as single value in each cell
  if (dim(tsv_maf_)[1]!=0) {
    tsv_maf_final <- data.frame()
    
    for (l in 1:nrow(tsv_maf_)) {
      maf_single = as.numeric(unlist(str_split(tsv_maf_[l,"maf"],":")))
      alt <- unlist(str_split(tsv_maf_[l, "allele_ratio"], "[=,]"))
      alt_logi <- !alt[grepl("[ATCG]",alt)] %in% tsv_maf_$ref[l]
      alt_letter <- alt[grepl("[ATCG]",alt)][alt_logi]
      freq <- unlist(str_split(tsv_maf_[l,"X._frequency"], "[=,]"))
      freq_num <- as.numeric(freq[grepl("\\d+",freq)])
      
      exac <- unlist(str_split(tsv_maf_[l,"exac"],"[:;]"))
      exac_AF_OTH <- as.numeric(unlist(str_split(exac[grepl("OTH", exac)], "="))[2])
      exac_AF_Adj <- as.numeric(unlist(str_split(exac[grepl("Adj", exac)], "="))[2])
      exac_AF_AFR <- as.numeric(unlist(str_split(exac[grepl("AFR", exac)], "="))[2])
      exac_AF_EAS <- as.numeric(unlist(str_split(exac[grepl("EAS", exac)], "="))[2])
      exac_AF_FIN <- as.numeric(unlist(str_split(exac[grepl("FIN", exac)], "="))[2])
      exac_AF_SAS <- as.numeric(unlist(str_split(exac[grepl("SAS", exac)], "="))[2])
      exac_AF_NFE <- as.numeric(unlist(str_split(exac[grepl("NFE", exac)], "="))[2])
      exac_AF_AMR <- as.numeric(unlist(str_split(exac[grepl("AMR", exac)], "="))[2])
      
      na_convert <- function(str_char) {
        te <- ifelse(is.null(str_char), NA, str_char)
        return(te) }
      
      if (!is.na(maf_single) & length(maf_single)>1 & c(length(maf_single)!= length(freq_num))) {
        temp_maf <- data.frame(tsv_maf_[l,c(1:5)], maf_parse = NA, ref = tsv_maf_[l,6], ALT = alt_letter, frequency_percent = freq_num, 
                               ExAc_AF_OTH = na_convert(exac_AF_OTH), ExAc_AF_Adj = na_convert(exac_AF_Adj), ExAc_AF_AFR = na_convert(exac_AF_AFR), 
                               ExAc_AF_EAS = na_convert(exac_AF_EAS), ExAc_AF_FIN = na_convert(exac_AF_FIN), ExAc_AF_SAS = na_convert(exac_AF_SAS), 
                               ExAc_AF_NFE = na_convert(exac_AF_NFE), ExAc_AF_AMR = na_convert(exac_AF_AMR), tsv_maf_[l,c(7:8)],
                               Sample = tsv_maf_[l, 12], stringsAsFactors = F) } 
      else {
        temp_maf <- data.frame(tsv_maf_[l,c(1:5)], maf_parse = maf_single, ref = tsv_maf_[l,6], ALT = alt_letter, frequency_percent = freq_num, 
                               ExAc_AF_OTH = na_convert(exac_AF_OTH), ExAc_AF_Adj = na_convert(exac_AF_Adj), ExAc_AF_AFR = na_convert(exac_AF_AFR), 
                               ExAc_AF_EAS = na_convert(exac_AF_EAS), ExAc_AF_FIN = na_convert(exac_AF_FIN), ExAc_AF_SAS = na_convert(exac_AF_SAS), 
                               ExAc_AF_NFE = na_convert(exac_AF_NFE), ExAc_AF_AMR = na_convert(exac_AF_AMR), tsv_maf_[l,c(7:8)],
                               Sample = tsv_maf_[l, 12], stringsAsFactors = F) } 
        
      tsv_maf_final <- smartbind(tsv_maf_final, temp_maf) }
    } else {
      tsv_maf_final <- tsv_maf_ 
    }
    tsv_maf_final$maf_parse <- as.numeric(tsv_maf_final$maf_parse)
    
  ## read in tsv with annotations
  tsv <- read.table(paste0(path_files[s], "/", files_onc[s]),  header = T, sep = "\t", stringsAsFactors = F, na.strings = c("", "NA", " "), comment.char = "#")
  
  col_codings <- colnames(tsv)[grepl("coding",colnames(tsv))]
  col_funcs <- colnames(tsv)[grepl("function",colnames(tsv))]
  col_genes <- colnames(tsv)[grepl("gene",colnames(tsv))]
  col_trans <- colnames(tsv)[grepl("transcript",colnames(tsv))]
  col_pros <- colnames(tsv)[grepl("protein",colnames(tsv))]
  col_exon <- colnames(tsv)[grepl("exon",colnames(tsv))]
  col_gln <- colnames(tsv)[grepl("CLNSIG1", colnames(tsv))]
  col_oncoGene <- colnames(tsv)[grepl("oncomineGeneClass", colnames(tsv))]
  col_oncoVariant <- colnames(tsv)[grepl("oncomineVariantClass", colnames(tsv))]
  # select variables/columns which are useful for filter
  tsv_colnames_manifest <- c("vcf.rownum", "rowtype", "CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO...MISC", "INFO...READ_COUNT", "INFO...RPM", "INFO.0.HS", "INFO.0.IMPRECISE", "INFO.1.EXON_NUM",
                             "INFO.1.FAIL_REASON", "INFO.1.GENE_NAME", "INFO.1.PVAL", "INFO.1.SVTYPE", "INFO.1.DP", "INFO.1.FDP", "INFO.A.AO", "INFO.A.AF", "INFO.A.FAO", "INFO.A.FSAF", "INFO.A.FSAR",
                             "INFO.A.SAF", "INFO.A.SAR", col_codings, col_funcs, col_genes, col_trans, col_pros, col_exon, col_gln, col_oncoGene, col_oncoVariant)
  
  tsv_colnames <- tsv_colnames_manifest[tsv_colnames_manifest %in% colnames(tsv)]
  
  tsv_ <- tsv[, tsv_colnames]
  tsv_$Locus <- paste0(tsv_$CHROM, ":", tsv_$POS)
  tsv_$Sample <- rep(sample_name, nrow(tsv_))
  tsv_ <- tsv_[order(tsv_$Locus),]
  
  ####=========================================================================================================================
  ####=========================================================================================================================
  #### generate combined tsv with all samples all variables/columns
  pos <- setNames(data.frame(paste0(tsv$CHROM, ":", tsv$POS), stringsAsFactors = F), "Locus")
  s_id <- setNames(data.frame(rep(sample_name, nrow(tsv)),stringsAsFactors = F), "Sample")
  tsv <- cbind(s_id, pos, tsv)
  tsv <- tsv[order(tsv$Locus),]

  ## concatenate information coding, funcs, genes, trans, pros, exons
  con_info <- function(arg1, df) {
    info <- df[,grepl(arg1, colnames(df))]
    if (is.null(dim(info)[2])) {
      df[,arg1] <- info
      df[,arg1] <- str_replace_na(df[,arg1], "")
      return(df)
    } else {
      fn_info <- lapply(1:nrow(df), function(x) paste(info[x,][!is.na(info[x,])], collapse = ", "))
      df_info <- do.call(rbind, fn_info)
      df[,arg1] <- df_info
      return(df)
    }
  }

  tsv <- con_info(arg1 = "coding", tsv)
  tsv <- con_info(arg1 = "function", tsv)
  tsv <- con_info(arg1 = "gene", tsv)
  tsv <- con_info(arg1 = "transcript", tsv)
  tsv <- con_info(arg1 = "protein", tsv)
  tsv <- con_info(arg1 = "exon", tsv)
  tsv <- con_info(arg1 = "CLNSIG1", tsv)
  tsv <- con_info(arg1 = "oncomineGeneClass", tsv)
  tsv <- con_info(arg1 = "oncomineVariantClass", tsv)

  colnames(tsv)[colnames(tsv)=="function"] <- "Variant_Effect"
  colnames(tsv)[colnames(tsv)=="rowtype"] <- "Type"

  ## parse SNV variants
  combine_tsv_snv <- subset(tsv, subset = c(is.na(tsv$INFO.1.SVTYPE) | tsv$INFO.1.SVTYPE =="FLT3ITD"), select = c(1:5,12:ncol(tsv)))

  if (dim(combine_tsv_snv)[1]!=0) {
    merge_comb_snv <- merge.data.frame(tsv_maf_final, combine_tsv_snv, by.x = c("X..locus", "ALT"), by.y = c("Locus", "ALT"), all = T)

    snv_comb <- subset(merge_comb_snv, c(grepl("missense", merge_comb_snv$Variant_Effect) | grepl("frameshift", merge_comb_snv$Variant_Effect)
                                         | grepl("nonsense",  merge_comb_snv$Variant_Effect) | grepl("stoploss",  merge_comb_snv$Variant_Effect)
                                         | grepl("Function", merge_comb_snv$oncomineGeneClass)  | (merge_comb_snv$INFO.1.SVTYPE=="FLT3ITD" & merge_comb_snv$INFO.A.AF > 0)))

    snv_comb <- subset(snv_comb, c(snv_comb$INFO.A.AF >= 0.03 | grepl("Function", snv_comb$oncomineGeneClass)))


    snv_comb_ <- snv_comb[, c(which(colnames(snv_comb)=="Sample.y"), which(colnames(snv_comb)=="X..locus"), 7:8, 11:ncol(snv_comb))]
    colnames(snv_comb_)[colnames(snv_comb_)=="Sample.y"] <- "Sample"
    colnames(snv_comb_)[colnames(snv_comb_)=="X..locus"] <- "Locus"
    snv_comb_ <- subset(snv_comb_, select = !colnames(snv_comb_) %in% c(col_codings, col_funcs, col_genes, col_trans, col_pros, col_exon, col_oncoGene, col_oncoVariant))
  }
  else {
    snv_comb_ <- combine_tsv_snv
  }

  ## parse Fusion variants
  fusion_comb <- subset(tsv, subset = tsv$INFO.1.SVTYPE== "Fusion" | tsv$INFO.1.SVTYPE=="RNAExonVariant", select = c(1:5,12:ncol(tsv)))
  fusion_comb$INFO...READ_COUNT <- as.numeric(fusion_comb$INFO...READ_COUNT)
  fusion_comb$INFO...RPM <- as.numeric(fusion_comb$INFO...RPM)
  fusion_comb <- subset(fusion_comb, fusion_comb$INFO...READ_COUNT >=15 | grepl("Function", fusion_comb$oncomineGeneClass))

  if (dim(fusion_comb)[1]!=0) {
    for (l in 1:nrow(fusion_comb)) {
      fusion_comb$ID_parse_1[l] <-stri_sub(fusion_comb$ID[l], 1, nchar(fusion_comb$ID[l])-2)
      fusion_comb$ID_parse_2[l] <-stri_sub(fusion_comb$ID[l], -1, -1)
      fusion_comb$unique[l] <- paste0(fusion_comb$Sample[l], "_",fusion_comb$ID_parse_1[l])
    }
    fusion_comb_ls <- split(fusion_comb, f = fusion_comb$unique)

    fn_fusion_comb <- lapply(1:length(fusion_comb_ls),function(ls){
      Locus_m = paste(fusion_comb_ls[[ls]]$Locus[1], fusion_comb_ls[[ls]]$Locus[2], sep='-')
      Gene1 = paste(fusion_comb_ls[[ls]]$GENE_NAME[1], '(', fusion_comb_ls[[ls]]$exon[1], ')')
      Gene2 = paste(fusion_comb_ls[[ls]]$GENE_NAME[2], '(', fusion_comb_ls[[ls]]$exon[2], ')')
      Genes = paste(Gene1, Gene2, sep = '-')
      return(data.frame(Sample=fusion_comb_ls[[ls]]$Sample[1],
                        Locus=Locus_m, ID=Genes, fusion_comb_ls[[ls]][1,c(3:(ncol(fusion_comb_ls[[ls]])-3))], stringsAsFactors = F))
    })

    fusion_comb_ <- do.call(rbind, fn_fusion_comb)
  }
  else {
    fusion_comb_ <- fusion_comb
  }

  ## merge Fusion and SNV into one single file
  if (dim(fusion_comb)[1]+dim(snv_comb_)[1]==0) {
    med_Data <- c(sample_name, c("NONE", rep("0",ncol(snv_comb_)-2)))
    output <- structure(rbind(snv_comb_, med_Data), .Names = names(snv_comb_)) }
  else if (dim(fusion_comb)[1]==0 & dim(snv_comb_)[1]!=0) {
    output <- snv_comb_ }
  else {
    output <- smartbind(fusion_comb_, snv_comb_)
  }

  combine_tsv <- smartbind(combine_tsv, output)

  cat(paste("\n", sample_name, "combined tsv in ", out_dir, "\n\n"))

  ####=========================================================================================================================
  ####=========================================================================================================================
  ## to generate report for review  

  tsv_ <- con_info(arg1 = "coding", tsv_)
  tsv_ <- con_info(arg1 = "function", tsv_)
  tsv_ <- con_info(arg1 = "gene", tsv_)
  tsv_ <- con_info(arg1 = "transcript", tsv_)
  tsv_ <- con_info(arg1 = "protein", tsv_)
  tsv_ <- con_info(arg1 = "exon", tsv_)
  tsv_ <- con_info(arg1 = "CLNSIG1", tsv_)
  tsv_ <- con_info(arg1 = "oncomineGeneClass", tsv_)
  tsv_ <- con_info(arg1 = "oncomineVariantClass", tsv_)

  ## shorten colnames
  tsv_ <- subset(tsv_, select = !colnames(tsv_) %in% c(col_codings, col_funcs, col_genes, col_trans, col_pros, col_exon, col_gln, col_oncoGene, col_oncoVariant))

  for (col_names in 1:ncol(tsv_)) {
    new_col_names <- unlist(str_split(colnames(tsv_)[col_names], "[ .]"))[length(unlist(str_split(colnames(tsv_)[col_names], "[ .]")))]
    colnames(tsv_)[col_names] <- new_col_names
  }
  colnames(tsv_)[colnames(tsv_)=="rowtype"] <- "Type"
  colnames(tsv_)[colnames(tsv_)=="function"] <- "Variant_Effect"

  ##### Step 2.a separate into SNV
  info_snv <- c("rownum", "Type", "ID", "REF", "ALT", "QUAL", "FILTER", "HS", "IMPRECISE", "SVTYPE", "DP", "FDP", "AO", "AF", "FAO", "FSAF", "FSAR", "SAF", "SAR",
                "CLNSIG1", "coding", "codon", "exon", "Variant_Effect", "gene", "location", "protein", "transcript", "oncomineGeneClass", "oncomineVariantClass")
  tsv_snv <- subset(tsv_, subset = c(is.na(tsv_$SVTYPE) | tsv_$SVTYPE=="FLT3ITD"), select = c("Sample", "Locus", info_snv[info_snv %in% colnames(tsv_)]))

  tryCatch(
      {merge_snv <- merge.data.frame(tsv_maf_final, tsv_snv, by.x = c("X..locus", "ALT"), by.y = c("Locus", "ALT"), all = T)},
      error=function(msg)
      {cat("sample:", sample_name, "does not have SNV \n")}
    )
    
  ## apply filterset to SNV  
  if (nrow(tsv_snv)!=0){ 
    snv_out <- subset(merge_snv, c(grepl("missense", merge_snv$Variant_Effect) | grepl("frameshift", merge_snv$Variant_Effect)
                                   | grepl("nonsense", merge_snv$Variant_Effect) | grepl("stoploss", merge_snv$Variant_Effect)
                                   | grepl("Function", merge_snv$oncomineGeneClass) | (merge_snv$SVTYPE=="FLT3ITD" & merge_snv$AF > 0)))

    snv_out_ <- subset(snv_out, subset=c(snv_out$AF >= 0.03 & (snv_out$maf_parse <= 0.01 | is.na(snv_out$maf_parse)) | grepl("Function", snv_out$oncomineGeneClass)))

    snv_out_final <- snv_out_[, c(21, 1:2, 4:19, 25:ncol(snv_out_))]
    colnames(snv_out_final)[colnames(snv_out_final)=="Sample.y"] <- "Sample"
    colnames(snv_out_final)[colnames(snv_out_final)=="X..locus"] <- "Locus"
    colnames(snv_out_final)[colnames(snv_out_final)=="gene.y"] <- "ID"
  } else {
    snv_out_final <- tsv_snv
  }

  snv_out_final <- subset(snv_out_final, snv_out_final$ExAc_AF_Adj <= 0.04 | is.na(snv_out_final$ExAc_AF_Adj) | grepl("Function", snv_out_final$oncomineGeneClass))

  ##### Step 2.b separate into Fusion
  info_fusion <- c("type", "ID", "REF", "ALT", "QUAL", "FILTER", "HS", "IMPRECISE", "SVTYPE", "CDF_MAPD", "MISC", "READ_COUNT", "RPM", "EXON_NUM", "GENE_NAME", "CLNSIG1", "oncomineGeneClass", "oncomineVariantClass")
  tsv_fusion <- subset(tsv_, subset = tsv_$SVTYPE== "Fusion" | tsv_$SVTYPE=="RNAExonVariant",
                       select = c("Sample", "Locus", info_fusion[info_fusion %in% colnames(tsv_)]))
  
  ## apply filterset to Fusion
  tsv_fusion$READ_COUNT <- as.numeric(tsv_fusion$READ_COUNT)
  tsv_fusion$RPM <- as.numeric(tsv_fusion$RPM)

  colnames(tsv_fusion)[colnames(tsv_fusion)=="EXON_NUM"] <- "exon"

  fusion_out <- subset(tsv_fusion, tsv_fusion$READ_COUNT >=15 | grepl("Function", tsv_fusion$oncomineGeneClass))

  cat(paste("\n", sample_name, "was successfully processed in", out_dir, "\n\n"))

  fusion_pri <- fusion_out[(fusion_out$SVTYPE=="Fusion" | fusion_out$SVTYPE=="RNAExonVariant"),]

  ## parse fusion gene name and exon in one cell/row 
  if (dim(fusion_out)[1]+dim(snv_out_final)[1]==0) {
    med_Data <- c(sample_name, c("NONE", rep("0",ncol(snv_out_final)-2)))
    output_2 <- structure(rbind(snv_out_final, med_Data), .Names = names(snv_out_final))

    } else if (dim(fusion_pri)[1]==0 & dim(snv_out_final)[1]!=0) {
    output_2 <- snv_out_final

    } else if (dim(fusion_pri)[1]!=0) {
      for (l in 1:nrow(fusion_pri)) {
        fusion_pri$ID_parse_1[l] <-stri_sub(fusion_pri$ID[l], 1, nchar(fusion_pri$ID[l])-2)
        fusion_pri$ID_parse_2[l] <-stri_sub(fusion_pri$ID[l], -1, -1)
        fusion_pri$unique[l] <- paste0(fusion_pri$Sample[l], "_",fusion_pri$ID_parse_1[l])
      }
    re_group_fusion <- split(fusion_pri, f = fusion_pri$unique)

    fn_ref_fu <- lapply(1:length(re_group_fusion),function(ls){
      Locus_m = paste(re_group_fusion[[ls]]$Locus[1], re_group_fusion[[ls]]$Locus[2], sep=' - ')
      Gene1 = paste(re_group_fusion[[ls]]$GENE_NAME[1], '(', re_group_fusion[[ls]]$exon[1], ')', sep = "")
      Gene2 = paste(re_group_fusion[[ls]]$GENE_NAME[2], '(', re_group_fusion[[ls]]$exon[2], ')', sep = "")
      Genes = paste(Gene1, Gene2, sep = ' - ')
      return(data.frame(Sample=re_group_fusion[[ls]]$Sample[1],
                        Locus=Locus_m, ID=Genes, re_group_fusion[[ls]][1,c(3:(ncol(re_group_fusion[[ls]])-3))], stringsAsFactors = F))
    })
    final_regroup_fusion <- do.call(rbind, fn_ref_fu)
    colnames(final_regroup_fusion)[colnames(final_regroup_fusion)=="SVTYPE"] <- "type"
    output_2 <- smartbind(final_regroup_fusion, snv_out_final)
    }

    ## rename variables for final report
    colnames(output_2)[colnames(output_2)=="ID"] <- "Genes"
    colnames(output_2)[colnames(output_2)=="READ_COUNT"] <- "Read Counts"
    colnames(output_2)[colnames(output_2)=="RPM"] <- "Read Counts Per Million"
    colnames(output_2)[colnames(output_2)=="protein"] <- "Amino Acid Change"
    colnames(output_2)[colnames(output_2)=="coding"] <- "Coding"
    colnames(output_2)[colnames(output_2)=="exon"] <- "Exon"
    colnames(output_2)[colnames(output_2)=="transcript"] <- "Transcript"
    colnames(output_2)[colnames(output_2)=="genotype"] <- "Genotype"
    colnames(output_2)[colnames(output_2)=="length"] <- "Length"
    colnames(output_2)[colnames(output_2)=="type"] <- "Type"
    colnames(output_2)[colnames(output_2)=="coverage"] <- "Coverage"

    report <- smartbind(report, output_2)
}

if (length(grep("IMPRECISE", colnames(report)))==0) {
  report$HS <- ifelse(report$HS=="True", "HS", "0")
  colnames(report)[colnames(report)=="HS"] <- "Info"
} else {
  report$IMPRECISE <- ifelse(report$IMPRECISE=="False", "", "IMPRECISE")
  report$IMPRECISE <- ifelse(is.na(report$IMPRECISE), "", report$IMPRECISE)
  report$HS <- ifelse(report$HS=="True", "HS", "")
  report$Info <- paste0(report$HS, report$IMPRECISE)
}

report_colnames <- c("Sample", "Locus", "Genes", "Type", "Exon", "Transcript", "Coding", "Variant_Effect", "Genotype", "Info", "Length",
                    "frequency_percent", "Amino Acid Change", "Read Counts", "Read Counts Per Million", "Coverage")

## remove SC which should not in final report
report <- subset(report, subset=!grepl("SC", report$Sample), select = report_colnames[report_colnames %in% colnames(report)])
report <- replace(report, is.na(report), "NA")

##### read in tumor purity% and barcode# from samplelist generated from samplelist_create_V2.R
tumor <- read.table(paste0(target_dir, run_id,"/", run_id, "_samplelist_V2.txt"), stringsAsFactors = F)
colnames(tumor) <- c("Sample", "tumor%", "BC")
final_report <- merge.data.frame(report, tumor, by.x = "Sample", by.y = "Sample", all.x = T)
final_report <- final_report %>% relocate(BC, .after= Locus)

## add variables/columns for pathologists to fill in 
path_rev <- c("Validated", "Dx", "P/M", "Site", "Pri.Opi.", "Resident_Reivew",	"Reviewer",	"Mol.Consensus")
path_rev_df <- setNames(data.frame(matrix(ncol = length(path_rev), nrow = nrow(final_report))), path_rev)
path_rev_df <- replace(path_rev_df, is.na(path_rev_df), "")

final_report_PI <- cbind(final_report, path_rev_df)

cat(paste("\n is writing re-formated final report for pathologist review \n\n"))

write_xlsx(list("Final_report" = final_report_PI, "Combine_tsv" = combine_tsv), paste0(out_dir, run_id, "_report.xlsx"))

cat(paste("\n final report for pathologist review is generated in" ,out_dir, "\n\n"))

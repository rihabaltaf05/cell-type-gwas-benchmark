#!/usr/bin/env Rscript
# 05_cct.R — Cauchy Combination Test across enrichment methods
# Args: [--ewce FILE] [--fgsea FILE] [--cepo FILE] [--scdrs DIR] --outdir DIR
suppressPackageStartupMessages(library(data.table))
args_raw <- commandArgs(trailingOnly=TRUE)
get_arg  <- function(flag,default=NULL){ i<-which(args_raw==flag); if(!length(i)) return(default); args_raw[i+1] }
EWCE_FILE  <- get_arg("--ewce"); FGSEA_FILE <- get_arg("--fgsea")
CEPO_FILE  <- get_arg("--cepo"); SCDRS_DIR  <- get_arg("--scdrs"); OUTDIR <- get_arg("--outdir")
if(is.null(OUTDIR)) stop("--outdir required")
dir.create(OUTDIR,showWarnings=FALSE,recursive=TRUE)

cauchy_combine <- function(pvec){
    pvec <- pmax(pmin(pvec,1-1e-15),1e-15)
    T_cct <- sum(tan((0.5-pvec)*pi))/length(pvec)
    max(0.5-atan(T_cct)/pi, 1e-300)
}

all_pvals <- list()
if(!is.null(EWCE_FILE)&&file.exists(EWCE_FILE)){
    cat(sprintf("[05_cct] EWCE: %s\n",EWCE_FILE),flush=TRUE)
    df<-read.csv(EWCE_FILE,stringsAsFactors=FALSE)
    p_col<-intersect(c("p","p_adjust","P"),colnames(df))[1]; ct_col<-intersect(c("CellType","cell_type"),colnames(df))[1]
    if(!is.na(p_col)&&!is.na(ct_col))
        all_pvals[["ewce"]] <- data.frame(trait=df$trait,cell_type=df[[ct_col]],p_ewce=df[[p_col]],stringsAsFactors=FALSE)
}
if(!is.null(FGSEA_FILE)&&file.exists(FGSEA_FILE)){
    cat(sprintf("[05_cct] fGSEA: %s\n",FGSEA_FILE),flush=TRUE)
    df<-read.csv(FGSEA_FILE,stringsAsFactors=FALSE)
    p_col<-intersect(c("pval","pvalue","p.value"),colnames(df))[1]
    ct_col<-intersect(c("cell_type","CellType"),colnames(df))[1]; tr_col<-intersect(c("pathway","trait"),colnames(df))[1]
    if(!is.na(p_col)&&!is.na(ct_col)&&!is.na(tr_col))
        all_pvals[["fgsea"]] <- data.frame(trait=df[[tr_col]],cell_type=df[[ct_col]],p_fgsea=df[[p_col]],stringsAsFactors=FALSE)
}
if(!is.null(CEPO_FILE)&&file.exists(CEPO_FILE)){
    cat(sprintf("[05_cct] CEPO: %s\n",CEPO_FILE),flush=TRUE)
    df<-read.csv(CEPO_FILE,stringsAsFactors=FALSE)
    p_col<-intersect(c("P","p","pval"),colnames(df))[1]; ct_col<-intersect(c("VARIABLE","cell_type"),colnames(df))[1]
    if(!is.na(p_col)&&!is.na(ct_col)&&"trait"%in%colnames(df))
        all_pvals[["cepo"]] <- data.frame(trait=df$trait,cell_type=df[[ct_col]],p_cepo=df[[p_col]],stringsAsFactors=FALSE)
}
if(!is.null(SCDRS_DIR)&&dir.exists(SCDRS_DIR)){
    cat(sprintf("[05_cct] scDRS: %s\n",SCDRS_DIR),flush=TRUE)
    scdrs_rows <- list()
    for(sf in list.files(SCDRS_DIR,pattern="\\.tsv$",full.names=TRUE)){
        trait<-sub("\\.tsv$","",basename(sf))
        df<-tryCatch(read.table(sf,header=TRUE,sep="\t",stringsAsFactors=FALSE),error=function(e)NULL); if(is.null(df)) next
        p_col<-intersect(c("assoc_mcp","pval","p"),colnames(df))[1]; g_col<-intersect(c("group","cell_type"),colnames(df))[1]
        if(is.na(p_col)||is.na(g_col)) next
        scdrs_rows[[trait]] <- data.frame(trait=trait,cell_type=df[[g_col]],p_scdrs=df[[p_col]],stringsAsFactors=FALSE)
    }
    if(length(scdrs_rows)) all_pvals[["scdrs"]] <- do.call(rbind,scdrs_rows)
}
if(length(all_pvals)<2) stop(sprintf("Need >=2 methods; got %d",length(all_pvals)))
cat(sprintf("[05_cct] Combining: %s\n",paste(names(all_pvals),collapse=", ")),flush=TRUE)

merged <- Reduce(function(a,b) merge(a,b,by=c("trait","cell_type"),all=TRUE), all_pvals)
p_cols <- grep("^p_",colnames(merged),value=TRUE)
write.csv(merged,file.path(OUTDIR,"pval_matrix.csv"),row.names=FALSE)

merged$p_cct <- NA_real_; merged$n_methods <- NA_integer_
for(r in seq_len(nrow(merged))){
    pv<-as.numeric(merged[r,p_cols]); pv<-pv[!is.na(pv)]
    if(length(pv)<2) next
    merged$p_cct[r]<-cauchy_combine(pv); merged$n_methods[r]<-length(pv)
}
complete <- !is.na(merged$p_cct)
merged$fdr_cct <- NA_real_
merged$fdr_cct[complete] <- p.adjust(merged$p_cct[complete],method="BH")
merged <- merged[order(merged$p_cct,na.last=TRUE),]
write.csv(merged,file.path(OUTDIR,"cct_combined.csv"),row.names=FALSE)
sig <- sum(!is.na(merged$fdr_cct)&merged$fdr_cct<0.05,na.rm=TRUE)
cat(sprintf("[05_cct] %d significant (FDR<0.05) pairs\n",sig),flush=TRUE)
if(sig>0){ print(head(merged[!is.na(merged$p_cct),c("trait","cell_type","p_cct","fdr_cct","n_methods")],10),row.names=FALSE) }
cat("[05_cct] Done ✓\n",flush=TRUE)

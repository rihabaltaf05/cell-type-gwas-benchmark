#!/usr/bin/env Rscript
# 01_ewce.R — EWCE bootstrap cell-type enrichment
# Args: --datadir DIR --outdir DIR [--reps 1000] [--species mouse]
suppressPackageStartupMessages({ library(Matrix); library(EWCE); library(data.table) })
args_raw <- commandArgs(trailingOnly=TRUE)
get_arg <- function(flag, default=NULL) {
    i <- which(args_raw==flag); if(!length(i)) return(default); args_raw[i+1] }
DATADIR <- get_arg("--datadir"); OUTDIR <- get_arg("--outdir")
REPS    <- as.integer(get_arg("--reps","1000"))
SPECIES <- get_arg("--species","mouse")
if(is.null(DATADIR)) stop("--datadir required")
if(is.null(OUTDIR))  stop("--outdir required")
dir.create(OUTDIR, showWarnings=FALSE, recursive=TRUE)
cat(sprintf("[01_ewce] datadir=%s outdir=%s reps=%d species=%s\n",DATADIR,OUTDIR,REPS,SPECIES),flush=TRUE)

cat("[01_ewce] Loading matrix ...\n",flush=TRUE)
exp      <- readMM(file.path(DATADIR,"counts.mtx"))
genes    <- readLines(file.path(DATADIR,"genes.tsv"))
barcodes <- readLines(file.path(DATADIR,"barcodes.tsv"))
rownames(exp)<-genes; colnames(exp)<-barcodes
meta <- read.csv(file.path(DATADIR,"cell_metadata.csv"),row.names=1,stringsAsFactors=FALSE)
meta <- meta[colnames(exp),,drop=FALSE]
cat(sprintf("  %d genes x %d cells | %d cell types\n",nrow(exp),ncol(exp),length(unique(meta$cell_type))),flush=TRUE)

cat("[01_ewce] Building CTD ...\n",flush=TRUE)
annot <- list(level1class=setNames(meta$cell_type,rownames(meta)))
ctd   <- generate_celltype_data(exp=exp,annotLevels=annot,groupName="benchmark",
            outPath=OUTDIR,no.hits=FALSE,dendrograms=FALSE)
if(is.character(ctd)){env<-new.env();load(ctd[1],envir=env);ctd<-get(ls(env)[1],envir=env)}

spec_mat <- ctd[[1]]$specificity
spec_df  <- cbind(gene=rownames(spec_mat),as.data.frame(as.matrix(spec_mat)))
write.table(spec_df,file.path(OUTDIR,"ctd_specificity.tsv"),sep="\t",row.names=FALSE,quote=FALSE)
cat(sprintf("  Saved specificity: %d genes x %d cell types\n",nrow(spec_mat),ncol(spec_mat)),flush=TRUE)

cat("[01_ewce] Loading gene lists ...\n",flush=TRUE)
gl_lines   <- readLines(file.path(DATADIR,"gene_lists.tsv"))
gene_lists <- list()
for(line in gl_lines){
    parts <- strsplit(line,"\t")[[1]]
    if(length(parts)<2||nchar(trimws(parts[1]))==0) next
    gene_lists[[parts[1]]] <- strsplit(trimws(parts[2]),",")[[1]] }
bg <- readLines(file.path(DATADIR,"background_genes.tsv"))
cat(sprintf("  %d traits, %d background genes\n",length(gene_lists),length(bg)),flush=TRUE)

all_results <- list()
for(i in seq_along(gene_lists)){
    trait <- names(gene_lists)[i]; hits <- gene_lists[[trait]]
    cat(sprintf("[01_ewce] [%d/%d] %s — %d hits\n",i,length(gene_lists),trait,length(hits)),flush=TRUE)
    tryCatch({
        res <- bootstrap_enrichment_test(sct_data=ctd,hits=hits,bg=bg,reps=REPS,
                    annotLevel=1,genelistSpecies="human",sctSpecies=SPECIES)
        df <- res$results; df$trait <- trait; all_results[[trait]] <- df
    },error=function(e) cat(sprintf("  ERROR: %s\n",conditionMessage(e)),flush=TRUE))
}
if(length(all_results)>0){
    out_df<-do.call(rbind,all_results)
    write.csv(out_df,file.path(OUTDIR,"ewce_results.csv"),row.names=FALSE)
    cat(sprintf("[01_ewce] Saved %d rows\n",nrow(out_df)),flush=TRUE)
} else cat("[01_ewce] WARNING: no results\n",flush=TRUE)
cat("[01_ewce] Done ✓\n",flush=TRUE)

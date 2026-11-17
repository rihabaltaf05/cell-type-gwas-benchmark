#!/usr/bin/env Rscript
# 04_fgsea.R — fGSEA enrichment of GWAS gene sets against cell-type specificity
# Args: --specificity TSV --gs GS --outdir DIR [--min-size 5] [--max-size 500] [--nperm 1000]
suppressPackageStartupMessages({ library(fgsea); library(data.table) })
args_raw <- commandArgs(trailingOnly=TRUE)
get_arg  <- function(flag,default=NULL){ i<-which(args_raw==flag); if(!length(i)) return(default); args_raw[i+1] }
SPEC_FILE <- get_arg("--specificity"); GS_FILE <- get_arg("--gs"); OUTDIR <- get_arg("--outdir")
MIN_SIZE  <- as.integer(get_arg("--min-size","5"))
MAX_SIZE  <- as.integer(get_arg("--max-size","500"))
NPERM     <- as.integer(get_arg("--nperm","1000"))
if(is.null(SPEC_FILE)) stop("--specificity required")
if(is.null(GS_FILE))   stop("--gs required")
if(is.null(OUTDIR))    stop("--outdir required")
dir.create(OUTDIR,showWarnings=FALSE,recursive=TRUE)

cat("[04_fgsea] Loading specificity matrix ...\n",flush=TRUE)
spec_df   <- read.table(SPEC_FILE,header=TRUE,sep="\t",stringsAsFactors=FALSE,check.names=FALSE)
gene_col  <- colnames(spec_df)[1]
spec_mat  <- as.matrix(spec_df[,-1,drop=FALSE])
rownames(spec_mat) <- spec_df[[gene_col]]
cell_types <- colnames(spec_mat)
cat(sprintf("  %d genes x %d cell types\n",nrow(spec_mat),length(cell_types)),flush=TRUE)

cat("[04_fgsea] Parsing gene lists ...\n",flush=TRUE)
gene_sets <- list()
for(line in readLines(GS_FILE)){
    line <- trimws(line); if(!nchar(line)) next
    parts <- strsplit(line,"\t")[[1]]; if(length(parts)<2) next
    g <- sapply(strsplit(parts[2],",")[[1]],function(x)strsplit(x,":")[[1]][1])
    g <- g[nchar(g)>0]; if(length(g)) gene_sets[[parts[1]]] <- g
}
cat(sprintf("  %d traits\n",length(gene_sets)),flush=TRUE)

all_results <- list()
for(ct in cell_types){
    ct_safe <- gsub("[^A-Za-z0-9]","_",ct)
    ranked  <- spec_mat[,ct]; ranked <- ranked[!is.na(ranked)]; ranked <- sort(ranked,decreasing=TRUE)
    cat(sprintf("[04_fgsea] %s (%d genes) ...\n",ct,length(ranked)),flush=TRUE)
    tryCatch({
        res  <- fgsea(pathways=gene_sets,stats=ranked,minSize=MIN_SIZE,maxSize=MAX_SIZE,nPermSimple=NPERM)
        df   <- as.data.frame(res); df$cell_type <- ct
        df$leadingEdge <- sapply(df$leadingEdge,paste,collapse=",")
        all_results[[ct_safe]] <- df
    },error=function(e) cat(sprintf("  ERROR: %s\n",conditionMessage(e)),flush=TRUE))
}
if(length(all_results)>0){
    out_df<-do.call(rbind,all_results)
    write.csv(out_df,file.path(OUTDIR,"fgsea_results.csv"),row.names=FALSE)
    cat(sprintf("[04_fgsea] Saved %d rows\n",nrow(out_df)),flush=TRUE)
}
cat("[04_fgsea] Done ✓\n",flush=TRUE)

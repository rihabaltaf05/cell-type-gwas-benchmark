#!/usr/bin/env Rscript
# 02_cepo_magma.R — CEPO specificity + MAGMA linear gene-property test
# Args: --datadir DIR --genes-dir DIR --outdir DIR [--species mouse] [--magma-bin PATH]
suppressPackageStartupMessages({ library(Matrix); library(Cepo); library(org.Hs.eg.db); library(AnnotationDbi) })
args_raw <- commandArgs(trailingOnly=TRUE)
get_arg  <- function(flag,default=NULL){ i<-which(args_raw==flag); if(!length(i)) return(default); args_raw[i+1] }
DATADIR   <- get_arg("--datadir"); GENES_DIR <- get_arg("--genes-dir"); OUTDIR <- get_arg("--outdir")
SPECIES   <- get_arg("--species","mouse")
MAGMA_BIN <- get_arg("--magma-bin","/home/ayeshazubair/projects/scrna_selitti/magma/magma")
if(is.null(DATADIR))   stop("--datadir required")
if(is.null(GENES_DIR)) stop("--genes-dir required")
if(is.null(OUTDIR))    stop("--outdir required")
dir.create(OUTDIR,showWarnings=FALSE,recursive=TRUE)

cat("[02_cepo] Loading matrix ...\n",flush=TRUE)
exp<-readMM(file.path(DATADIR,"counts.mtx"))
genes<-readLines(file.path(DATADIR,"genes.tsv")); barcodes<-readLines(file.path(DATADIR,"barcodes.tsv"))
rownames(exp)<-genes; colnames(exp)<-barcodes
meta<-read.csv(file.path(DATADIR,"cell_metadata.csv"),row.names=1,stringsAsFactors=FALSE)
meta<-meta[colnames(exp),,drop=FALSE]
cat(sprintf("  %d genes x %d cells\n",nrow(exp),ncol(exp)),flush=TRUE)

cat("[02_cepo] Normalising ...\n",flush=TRUE)
cs<-Matrix::colSums(exp); cs[cs==0]<-1
exp_norm<-log1p(Matrix::t(Matrix::t(exp)/cs)*1e4)

cat("[02_cepo] Running CEPO ...\n",flush=TRUE)
cepo_res <- Cepo(exprsMat=exp_norm, cellTypes=meta$cell_type, minCells=5, exprsPct=0)
ds_mat   <- cepo_res$stats
cat(sprintf("  DS scores: %d genes x %d cell types\n",nrow(ds_mat),ncol(ds_mat)),flush=TRUE)

if(SPECIES=="mouse"){
    cat("[02_cepo] Mouse -> human ortholog mapping ...\n",flush=TRUE)
    if(!requireNamespace("homologene",quietly=TRUE)) install.packages("homologene",repos="https://cloud.r-project.org")
    library(homologene)
    hom <- homologene(rownames(ds_mat),inTax=10090,outTax=9606)
    colnames(hom) <- c("mouse","human","mouse_tax","human_tax")
    hom <- hom[!is.na(hom$human)&nchar(hom$human)>0,]
    hom <- hom[!duplicated(hom$mouse)&!duplicated(hom$human),]
    idx <- match(hom$mouse,rownames(ds_mat))
    ds_mat <- ds_mat[idx[!is.na(idx)],,drop=FALSE]; rownames(ds_mat) <- hom$human[!is.na(idx)]
    cat(sprintf("  After mapping: %d genes\n",nrow(ds_mat)),flush=TRUE)
}

spec_df  <- cbind(gene=rownames(ds_mat),as.data.frame(ds_mat))
write.table(spec_df,file.path(OUTDIR,"cepo_specificity.tsv"),sep="\t",row.names=FALSE,quote=FALSE)

cat("[02_cepo] Mapping to Entrez IDs ...\n",flush=TRUE)
entrez_map <- AnnotationDbi::mapIds(org.Hs.eg.db,keys=rownames(ds_mat),column="ENTREZID",keytype="SYMBOL",multiVals="first")
valid  <- !is.na(entrez_map)&nchar(entrez_map)>0
ds_ent <- ds_mat[valid,,drop=FALSE]; rownames(ds_ent) <- entrez_map[valid]
cat(sprintf("  %d genes with Entrez IDs\n",nrow(ds_ent)),flush=TRUE)
covar_df <- cbind(GENE=rownames(ds_ent),as.data.frame(ds_ent))
colnames(covar_df) <- gsub("[^A-Za-z0-9_]","_",colnames(covar_df))
covar_file <- file.path(OUTDIR,"cepo_gene_covar.txt")
write.table(covar_df,covar_file,sep=" ",row.names=FALSE,quote=FALSE)
cat(sprintf("[02_cepo] Gene covar saved: %s\n",covar_file),flush=TRUE)

genes_files <- list.files(GENES_DIR,pattern="\\.genes\\.out$",full.names=TRUE)
if(!length(genes_files)){
    cat(sprintf("[02_cepo] No .genes.out files in %s — CEPO only done\n",GENES_DIR),flush=TRUE)
    quit(save="no",status=0)
}
cat(sprintf("[02_cepo] Running MAGMA for %d traits ...\n",length(genes_files)),flush=TRUE)
magma_results <- list()
for(gf in genes_files){
    trait   <- sub("\\.genes\\.out$","",basename(gf))
    out_pfx <- file.path(OUTDIR,paste0("magma_cepo_",trait))
    cmd     <- sprintf("%s --gene-results %s --gene-covar %s --model linear --out %s",MAGMA_BIN,gf,covar_file,out_pfx)
    cat(sprintf("  [%s]\n",trait),flush=TRUE); system(cmd)
    gsa_file <- paste0(out_pfx,".gsa.out")
    if(!file.exists(gsa_file)) next
    df <- tryCatch(read.table(gsa_file,header=TRUE,comment.char="#",stringsAsFactors=FALSE),error=function(e)NULL)
    if(!is.null(df)){ df$trait<-trait; magma_results[[trait]]<-df }
}
if(length(magma_results)>0){
    out_df<-do.call(rbind,magma_results)
    write.csv(out_df,file.path(OUTDIR,"cepo_magma_results.csv"),row.names=FALSE)
    cat(sprintf("[02_cepo] Saved %d rows\n",nrow(out_df)),flush=TRUE)
}
cat("[02_cepo] Done ✓\n",flush=TRUE)

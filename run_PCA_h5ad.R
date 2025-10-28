
library(getopt)

# Define the specification for the options
spec = matrix(c(
  'h5ad', 'h', 1, 'character',
  'method', 'm', 1, 'character',
  'nthreads', 'n', 1, 'numeric',
  'ncells', 'c', 1, 'numeric',
  'nPC', 'k', 1, 'numeric',
  'out', 'o', 1, 'character'
), byrow=TRUE, ncol=4)

# Parse the arguments
opt = getopt(spec)

suppressPackageStartupMessages({
library(GenomicDataStream)
library(SingleCellExperiment)
library(zellkonverter)
library(BiocSingular)
})

# run_PCA_h5ad.R --h5ad ~/Downloads/4e6932db-5a78-40e4-b961-f87f66ba139a.h5ad

# Read data
sce = readH5AD(opt$h5ad, use_hdf5=TRUE, raw=TRUE, verbose=FALSE, uns=FALSE, obsp=FALSE, obsm=FALSE)

sce = sce[,seq(opt$ncells)]
counts(sce) = assay(sce, "X")

# Compute log CPM
sce$total_counts = colSums2(counts(sce))

lib.size = sce$total_counts
prior.count = 1
logcounts(sce) <- t(log2(t(counts(sce) + prior.count)) - log2(lib.size) + log2(1e6))

if( opt$method == "PCAstream" ){

  res_time = system.time({
    res <- PCAstream(sce, opt$nPC, assay="logcounts", threads=opt$nthreads)
  })

  res_time = data.frame(t(data.frame(res_time)), Method="PCAstream")
}

if( opt$method == "IRLBA" ){
  res_time = system.time({
    out <- runSVD( t(logcounts(sce)), k=opt$nPC, center=TRUE, scale=TRUE, BSPARAM=IrlbaParam())
  })
  res_time = data.frame(t(data.frame(res_time)), Method="IRLBA")
}

if( opt$method == "RSVD" ){
  res_time = system.time({
    out <- runSVD( t(logcounts(sce)), k=opt$nPC, center=TRUE, scale=TRUE, BSPARAM=RandomParam())
  })
  res_time = data.frame(t(data.frame(res_time)), Method="Randomized SVD")
}

if( opt$method == "exact" ){
  res_time = system.time({
    out <- runSVD( t(logcounts(sce)), k=opt$nPC, center=TRUE, scale=TRUE, BSPARAM=ExactParam())
  })
  res_time = data.frame(t(data.frame(res_time)), Method="Exact SVD")
}

res_time$ngenes = nrow(sce)
res_time$ncells = ncol(sce)
res_time$file = opt$h5ad

write.table(res_time, file=opt$out, quote=FALSE, sep="\t")







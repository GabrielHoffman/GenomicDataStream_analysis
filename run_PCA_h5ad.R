
library(getopt)

# Define the specification for the options
spec = matrix(c(
  'h5ad', 'h', 1, 'character',
  'method', 'm', 1, 'character',
  'nthreads', 'n', 1, 'numeric',
  'ncells', 'c', 1, 'numeric',
  'nPC', 'k', 1, 'numeric',
  'scaleAndCenter', 's', 0, "logical",
  'out', 'o', 1, 'character'
), byrow=TRUE, ncol=4)

# Parse the arguments
opt = getopt(spec)

if( is.null(opt$scaleAndCenter) ){
  opt$scaleAndCenter = FALSE
}

suppressPackageStartupMessages({
library(GenomicDataStream)
library(SingleCellExperiment)
library(BiocParallel)
library(BiocSingular)
library(HDF5Array)
})




# Read data
sce = readH5AD(opt$h5ad)

sce = sce[,seq(opt$ncells)]

# Compute log CPM
sce$total_counts = colSums2(counts(sce))

# 5000 genes
idx = as.integer(seq(1, nrow(sce), length.out=5000))
sce = sce[idx,]

lib.size = sce$total_counts
prior.count = 1
logcounts(sce) <- t(log2(t(counts(sce) + prior.count)) - log2(lib.size) + log2(1e6))

if( opt$method == "PCAstream" ){

  res_time = system.time({
    res <- PCAstream(sce, opt$nPC, assay="logcounts", threads=opt$nthreads, scaleAndCenter=opt$scaleAndCenter)
  })

  res_time = data.frame(t(data.frame(res_time)), Method="PCAstream")
}

if( opt$method == "IRLBA" ){
  res_time = system.time({
    out <- runSVD( t(logcounts(sce)), k=opt$nPC, center=opt$scaleAndCenter, scale=opt$scaleAndCenter, BSPARAM=IrlbaParam(), BPPARAM=SnowParam(opt$nthreads))
  })
  res_time = data.frame(t(data.frame(res_time)), Method="IRLBA")
}

if( opt$method == "RSVD" ){
  res_time = system.time({
    out <- runSVD( t(logcounts(sce)), k=opt$nPC, center=opt$scaleAndCenter, scale=opt$scaleAndCenter, BSPARAM=RandomParam(), BPPARAM=SnowParam(opt$nthreads))
  })
  res_time = data.frame(t(data.frame(res_time)), Method="RSVD")
}

if( opt$method == "Exact" ){
  res_time = system.time({
    out <- runSVD( t(logcounts(sce)), k=opt$nPC, center=opt$scaleAndCenter, scale=opt$scaleAndCenter, BSPARAM=ExactParam(), BPPARAM=SnowParam(opt$nthreads))
  })
  res_time = data.frame(t(data.frame(res_time)), Method="Exact")
}

if( opt$method == "Standard" ){
  library(RhpcBLASctl)
  omp_set_num_threads(opt$nthreads)
  blas_set_num_threads(opt$nthreads)
  res_time = system.time({
    X <- as.matrix(t(logcounts(sce)))
    if( opt$scaleAndCenter ){
      X <- scale(X)
    }
    out <- svd(X, nu=opt$nPC, nv=opt$nPC)
  })
  res_time = data.frame(t(data.frame(res_time)), Method="Standard")
}

res_time$ngenes = nrow(sce)
res_time$ncells = ncol(sce)
res_time$file = opt$h5ad
res_time$scaleAndCenter = opt$scaleAndCenter

write.table(res_time, file=opt$out, quote=FALSE, sep="\t", row.names=FALSE)







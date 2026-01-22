#!/usr/bin/env python3

import sys
import argparse
import anndata as ad
import numpy as np
import scipy.sparse as sp
from pathlib import Path


def main():

  # Define arguments
  #-------------------

  parser = argparse.ArgumentParser(
      description="Convert an AnnData .h5ad file so that X is stored in CSC sparse format."
  )

  parser.add_argument("--input",
      type=Path, 
      required=True,
      help="Input .h5ad file")

  parser.add_argument("--output",      
      type=Path, 
      required=True,
      help="Output .h5ad file")

  parser.add_argument("--ondisk", 
      action="store_true",
      help="Output .h5ad file")

  parser.add_argument("--sortBy", 
      default=None,
      help="Cols to sort by in _decreasing_ order of importance")

  parser.add_argument(
      "--compression",
      default=None,
      choices=[None, "None", "gzip", "lzf"],
      help="Optional compression for output file (gzip or lzf). Default: None")

  parser.add_argument(
      "--format",
      default="CSR",
      choices=["CSR", "CSC"],
      help="Store sparse count matrix in CSR or CSC format. Default: CSR")

  # parse args
  args = parser.parse_args()

  # if --ondisk set backed to "r"
  backed = None
  if args.ondisk:
    backed = "r"

  if args.compression == "None":
    args.compression = None

  # read file
  print("Read file...") 
  adata = ad.read_h5ad(args.input, backed=backed) 

  # Find counts entry
  if adata.X is not None:
      print("Using AnnData X matrix...")
  else:
    if 'X' in adata.layers: 
      print("Using AnnData layers/X matrix...")
      adata.X = adata.layers['X']
    elif 'counts' in adata.layers: 
      print("Using AnnData layers/counts matrix...")
      adata.X = adata.layers['counts']

    # sort cells by type  
  if args.sortBy != None:

    print("Sorting...") 
    fields = args.sortBy.split(",")

    # check if all fields are present
    missing = set(fields) - set(adata.obs.columns)

    if missing:
      print(f"Missing columns: {missing}")
      sys.exit(2)     

    # get sorted order
    # reverse since 1st sorted index is the last one for lexsort
    idx = np.lexsort( 
      keys = tuple(adata.obs[c].to_numpy() for c in fields[::-1]) )

    # apply reordering
    adata = adata[idx,:]

  if args.format == "CSC":
    if not sp.isspmatrix_csc(adata.X):
      print("Converting .X to CSC sparse format...")
      if backed is None:
        adata = adata.copy()
      else:
        adata = adata.to_memory()

      # convert matrix type
      adata.X = sp.csc_matrix(adata.X)

  if args.format == "CSR":
    if not sp.isspmatrix_csr(adata.X):
      print("Converting .X to CSR sparse format...")
      if backed is None:
        adata = adata.copy()
      else:
        adata = adata.to_memory()
      
      # convert matrix type
      adata.X = sp.csr_matrix(adata.X)

  if args.output.is_file():
    args.output.unlink(missing_ok=True)

  print("Writing H5AD...") 

  adata.write_h5ad( args.output, compression=args.compression )



if __name__ == "__main__":
  main()

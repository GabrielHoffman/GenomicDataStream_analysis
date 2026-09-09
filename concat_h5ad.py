#!/usr/bin/env python3

import sys, getopt, os
import anndata as ad
from numpy import loadtxt
import h5py
import numpy as np
from scipy import sparse
from anndata._core.sparse_dataset import sparse_dataset
from anndata.io import read_elem, write_elem

# From Donghoon Lee
def read_everything_but_X(pth) -> ad.AnnData:
    # read all keys but X and raw
    with h5py.File(pth) as f:
        attrs = list(f.keys())
        attrs.remove('X')
        if 'raw' in attrs:
            attrs.remove('raw')
        adata = ad.AnnData(**{k: read_elem(f[k]) for k in attrs})
        print(adata.shape)
    return adata

# From Donghoon Lee
def concat_on_disk(input_pths, output_pth, temp_pth='temp.h5ad'):
    """
    Params
    ------
    input_pths
        Paths to h5ad files which will be concatenated
    output_pth
        File to write as a result
    """
    annotations = ad.concat([read_everything_but_X(pth) for pth in input_pths])
    annotations.write_h5ad(output_pth)
    n_variables = annotations.shape[1]
    
    del annotations

    with h5py.File(output_pth, 'a') as target:
        
        # initiate empty X
        dummy_X = sparse.csr_matrix((0, n_variables), dtype=np.float32)
        dummy_X.indptr = dummy_X.indptr.astype(np.int64) # Guarding against overflow for very large datasets
        dummy_X.indices = dummy_X.indices.astype(np.int64) # Guarding against overflow for very large datasets
        write_elem(target, 'X', dummy_X)
        
        # append
        mtx = sparse_dataset(target['X'])
        for p in input_pths:
            with h5py.File(p, 'r') as src:
                
                # IF: src is in csc format, convert to csr and save to temp_pth
                if src['X'].attrs['encoding-type']=='csc_matrix':

                    # Convert to csr format
                    csc_mat = sparse.csc_matrix((src['X']['data'], src['X']['indices'], src['X']['indptr']))
                    csr_mat = csc_mat.tocsr()         
                    
                    # save to temp_pth
                    with h5py.File(temp_pth, 'w') as tmp:
                        write_elem(tmp, 'X', csr_mat)
                    
                    # read from temp_pth
                    with h5py.File(temp_pth, 'r') as tmp:
                        mtx.append(sparse_dataset(tmp['X']))
                        
                # ELSE: src is in csr format
                else:
                    mtx.append(sparse_dataset(src['X']))
                    
def main(argv):
   infile = ''
   outfile = ''
   ondisk = False
   try:
      opts, args = getopt.getopt(argv,"dhi:o:",["ifile=","ofile="])
   except getopt.GetoptError:
      print('concat_h5ad.py [-d] -i <file of h5ad files> -o <outfile h5ad>')
      sys.exit(2)
   for opt, arg in opts:
      if opt == '-h':
         print('concat_h5ad.py -i <file of h5ad files> --ondisk -o <outfile h5ad>')
         sys.exit()
      elif opt in ("-i", "--ifile"):
         infile = arg
      elif opt in ("-o", "--ofile"):
         outfile = arg
      elif opt in ("-d", "--ondisk"):
         ondisk = True

   if infile == None:
      print("Must specify infile")
      sys.exit(2)      

   if outfile == '':
      print("Must specify outfile")
      sys.exit(2)   

   # Run code data
   ################

   # read list of h5a files
   h5adfiles = open(infile, "r").read().split('\n')

   # remove blank lines
   h5adfiles = list(filter(None, h5adfiles))

   print( "Combine", len(h5adfiles), "H5AD files")

   if ondisk:
      concat_on_disk(h5adfiles, outfile)
   else:
      # Create anndata array
      print(" Reading files...")
      adatas = [ad.read_h5ad(file) for file in h5adfiles] 

      # Concatenate data
      print(" Concatenating data...")
      adata = ad.concat(adatas[:])

      # add .var data (i.e. rowData) to new object
      # adata.var = adatas[1].var

      # Write to disk
      print(" Writing to H5AD...")
      adata.write( outfile, compression="lzf" )

if __name__ == "__main__":
   main(sys.argv[1:])




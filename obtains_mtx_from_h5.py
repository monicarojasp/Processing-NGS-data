import h5py
import numpy as np
from scipy.sparse import csr_matrix
import pandas as pd
import gzip
import argparse
import os

def extract_h5_data(file_path):
    # Load data from the HDF5 file
    with h5py.File(file_path, 'r') as f:
        # Access the matrix group
        matrix_group = f['matrix']
        
        # Load barcodes
        barcodes = matrix_group['barcodes'][:].astype(str)

        # Access the features group
        features_group = matrix_group['features']
        
        # Load feature names and IDs
        feature_names = features_group['name'][:].astype(str)
        feature_ids = features_group['id'][:].astype(str)
        
        # Load sparse matrix components
        data = matrix_group['data'][:]
        indices = matrix_group['indices'][:]
        indptr = matrix_group['indptr'][:]
        shape = matrix_group['shape'][:]

      # Debugging output
    print(f"Data shape: {data.shape}")
    print(f"Indices shape: {indices.shape}")
    print(f"Indptr shape: {indptr.shape}")
    print(f"Shape: {shape}")
    print(f"Length of indptr: {len(indptr)}")
    print(f"Indptr first 10 values: {indptr[:10]}")
    print(f"Data first 10 values: {data[:10]}")
    print(f"Indices first 10 values: {indices[:10]}")


    # Create a sparse matrix
    sparse_matrix = csr_matrix((data, indices, indptr), shape=shape)

    # Save the matrix to a .mtx.gz file
    mtx_filename = 'matrix.mtx.gz'
    with gzip.open(mtx_filename, 'wb') as f:
        f.write(f"%%MatrixMarket matrix coordinate real general\n".encode())
        f.write(f"{sparse_matrix.shape[0]} {sparse_matrix.shape[1]} {sparse_matrix.nnz}\n".encode())
        for row, col in zip(*sparse_matrix.nonzero()):
            f.write(f"{row + 1} {col + 1} {sparse_matrix[row, col]}\n".encode())  # +1 for 1-based indexing

    # Save feature names and IDs to features.tsv.gz
    features_filename = 'features.tsv.gz'
    with gzip.open(features_filename, 'wt') as f:
        for feature_id, feature_name in zip(feature_ids, feature_names):
            f.write(f"{feature_id}\t{feature_name}\n")

    # Save barcodes to barcodes.tsv.gz
    barcodes_filename = 'barcodes.tsv.gz'
    with gzip.open(barcodes_filename, 'wt') as f:
        for barcode in barcodes:
            f.write(f"{barcode}\n")

    print(f"Files saved: {mtx_filename}, {features_filename}, {barcodes_filename}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Extract data from an HDF5 file.")
    parser.add_argument('file_path', type=str, help='Path to the HDF5 file (filtered_feature_bc_matrix.h5)')

    args = parser.parse_args()

    if not os.path.isfile(args.file_path):
        print(f"Error: The file '{args.file_path}' does not exist.")
    else:
        extract_h5_data(args.file_path)

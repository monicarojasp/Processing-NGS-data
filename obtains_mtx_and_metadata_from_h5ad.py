#! /usr/bin/env python

"""
This script extracts sample-level counts from a specified .h5ad file 
and writes the results into MTX format files. It processes each sample 
in the input file and generates separate output files for the matrix, 
barcodes, and features.

Input File Format:
    - The input file must be in the AnnData format (.h5ad).
    - It should contain a count matrix where rows correspond to genes and columns correspond to cells.
    - Sample information should be stored in `adata.obs` under the column 'Sample'.
    - The variable annotations (features) should be in `adata.var`.
    - Barcodes (cell identifiers) are expected to be the index of the `adata.obs` DataFrame, which should match the column names of the count matrix.

Usage:
    python process_mtx.py -i /path/to/input_file.h5ad -o /path/to/output_directory

Arguments:
    -i, --infile            Path to the input .h5ad file (required)
    -o, --outdir            Path to the directory where output files will be saved (required)
    -s, --split_by          Option to split the file into samples
    -l, --columns_metadata  Option to extract metadata specific columns
    -g, --gene_id           Option to choose second features gene name column
    -b, --barcode           Option to specify the name of the column containing barcode identifiers

Author: Monica L. Rojas-Pena
Email: monicarojasp@phenomic.ai
"""

import os
import sys
from pathlib import Path
import argparse
import scanpy as sc
import pandas as pd
from scipy.io import mmwrite
import shutil
import gzip
import numpy as np

# Add the data_curation module to the system path
home = str(Path.home())
sys.path.append(home + "/data_curation/src/data_curation")
from data_curation.data_curation_utils.python.utils.run_logs import RunLogs as run_logs

# Set up argument parser
usage = "usage: %prog [options]"
parser = argparse.ArgumentParser(description='Extracts sample-level counts from study.h5ad and writes MTX files')
parser.add_argument('-i', '--infile', type=str, required=True,
                    help='Path to the input .h5ad file')
parser.add_argument('-o', '--outdir', type=str, required=True,
                    help='Path to the directory where output files will be saved')
parser.add_argument("-s", "--split_by", type=str, required=False,
                    help="Indicates if MTX should be split by a specific column in adata.obs. Provide the column name as the value, e.g., 'sample_id'.")
parser.add_argument("-l", "--columns_metadata", type=str, required=False,
                    help="Path to a .txt file containing the list of columns to be extracted from adata.obs for output")
parser.add_argument("-g", "--gene_id", type=str, required=False,
                     help="Name of the column in adata.var containing gene identifiers.")
parser.add_argument("-b", "--barcode", type=str, required=False,
                     help="Name of the column in adata.obs containing barcode identifiers.")

# Parse the arguments
args = parser.parse_args()

# Get input file and output directory from arguments
input_file = args.infile
output_dir = args.outdir
split_by = args.split_by if args.split_by else 'Sample'  # Default to 'Sample' if not provided
gene_id = args.gene_id if args.gene_id else 'gene_id'

# Load your AnnData object
adata = sc.read_h5ad(input_file)

# Export metadata if the -l argument is provided
if args.columns_metadata:
    # Read the metadata columns from the provided file
    metadata_columns = pd.read_csv(args.columns_metadata, sep="\t", header=None).iloc[0].tolist()
    
    # Check if all requested columns exist in adata.obs
    missing_columns = [col for col in metadata_columns if col not in adata.obs.columns]
    if missing_columns:
        raise ValueError(f"Columns {missing_columns} not found in adata.obs.")
    
    # Export metadata to a compressed file
    metadata_outfile = os.path.join(os.path.dirname(input_file), "metadata.tsv.gz")  # Fix output path
    adata.obs[metadata_columns].drop_duplicates().to_csv(
        metadata_outfile, sep="\t", compression="gzip", index=False
    )
    print(f"Metadata exported to {metadata_outfile} successfully!")

# Write sample data to MTX files
def write_sample_data_to_mtx(adata, output_dir):
    # Iterate through each unique sample in the "sample" column
    for sample_name in adata.obs[split_by].unique():
        print(f"Writing {sample_name}")
        
        # Create the output directory for the sample
        outdir = os.path.join(output_dir, sample_name)
        os.makedirs(outdir, exist_ok=True)

        # Filter the data for the current sample
        sample_adata = adata[adata.obs[split_by] == sample_name,]

        # Extracting counts (ensure it's a sparse matrix)
        mtx = sample_adata.X.T.copy()  # Transpose to match MTX format

        # Prepare the barcodes and features data
        if args.barcode:
            # Use the barcode column name provided by the user
            barcode_column = args.barcode
            if barcode_column not in adata.obs.columns:
                raise ValueError(f"Column '{barcode_column}' not found in adata.obs.")
            barcodes = sample_adata.obs[barcode_column].copy()  # Use the barcode column specified by the user
        else:
            # Default to using the index if no --barcode argument is provided
            barcodes = sample_adata.obs.index.copy()

        # Convert barcodes (Index) to DataFrame
        barcodes_df = pd.DataFrame(barcodes)

        # Features (genes)
        features = pd.DataFrame(adata.var.index)
        features.columns = ["gene"]
        features["gene_id"] = adata.var[gene_id].values
        features["gene_exp"] = "Gene Expression"
        if gene_id not in adata.var.columns:
            raise ValueError(f"Column '{gene_id}' not found in adata.var.")

        # Define file paths for the matrix, barcodes, and features
        matrix_fp = os.path.join(outdir, "matrix.mtx")
        barcodes_fp = os.path.join(outdir, "barcodes.tsv.gz")
        features_fp = os.path.join(outdir, "features.tsv.gz")

        # Write the matrix (in MTX format) to a file
        mmwrite(target=matrix_fp, a=mtx, field="integer")
        print(f"Matrix shape: {mtx.shape}")

        # Write barcodes and features to their respective files
        barcodes_df.to_csv(barcodes_fp, sep="\t", header=None, index=False, compression="gzip")
        features.to_csv(features_fp, sep="\t", header=None, index=False, compression="gzip")

        # Compress the matrix file into .gz format
        with open(matrix_fp, "rb") as f_in:
            with gzip.open(matrix_fp + ".gz", "wb") as f_out:
                shutil.copyfileobj(f_in, f_out)

        # Remove the uncompressed matrix file after compression
        os.remove(matrix_fp)

        print(f"Files for sample {sample_name} written to {outdir}")

# Main script execution
if args.columns_metadata:
    write_sample_data_to_mtx(adata, output_dir)  # Run if metadata export occurred
else:
    write_sample_data_to_mtx(adata, output_dir)  # Run directly if no metadata export

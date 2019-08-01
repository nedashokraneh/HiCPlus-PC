import os
import sys
import numpy as np
import pandas as pd
import argparse
import cooler
import re
utils_path = "utils"
sys.path.insert(0, utils_path)
import utils

ap = argparse.ArgumentParser()
ap.add_argument("--input_type", help = "0 if it is cool and 1 if it is COO folder path", type = int)
ap.add_argument("--cool_file_path", help = "path of cool file to make frames from")
ap.add_argument("--resolution", help = "resolution of HiC matrix to be used", type = int)
ap.add_argument("--min_chrN", help = "minimum number of chromosome from cool file to be included in frames", type = int)
ap.add_argument("--max_chrN", help = "maximum number of chromosome from cool file to be included in frames", type = int)
ap.add_argument("--COO_folder_path", help = "path of folder containing COO files to make frames from")
ap.add_argument("--COO_format", help = "columns format in COO files(0 if bin number, 1 if genomic location)", type = int)
ap.add_argument("--genome_type", help = "hg19 or hg38")
ap.add_argument("--output_folder_path", help = "path of folder to save down sample data", required = True)
ap.add_argument("--output_file_name", help = "name of npy file to be saved", required = True)
args = vars(ap.parse_args())

if os.path.exists(os.path.join(args['output_folder_path'], args['output_file_name'] + ".npy")):
    print("data already exists")


if args['input_type'] == 0:
    cl = cooler.Cooler(args['cool_file_path'] + '::/resolutions/' + str(args['resolution']))
    frames_data = []
    index_data = []
    for i in range(args['min_chrN'], args['max_chrN'] + 1):
        chr_mat = cl.matrix(balance = False).fetch("chr" + str(i)).astype(float)
        chr_mat[np.isnan(chr_mat)] = 0
        chr_frames, chr_indices = utils.divide2(chr_mat,i)
        frames_data.extend(chr_frames)
        index_data.extend(chr_indices)

else:
    COO_folder_path = args['COO_folder_path']
    """
    # this works if format of file names in a COO_folder is not like chr[chr_num].txt but we should care about sequence of filenames to be same in
    # folder of high resolution files and folder of low resolution files
    chr_files_list = [f for f in os.listdir(COO_folder_path) if (not f.startswith('.')) & (not f.endswith('_npy_form_tmp.npy'))]
    chr_num_list = []
    for f in chr_files_list:
        m = re.search('chr(\d+|x)', f, re.IGNORECASE)
        chr_num_list.append(int(m.group(1)))
    """
    frames_data = []
    index_data = []
    for i in range(args['min_chrN'], args['max_chrN'] + 1):
        chr_COO_file_name = "chr" + str(i) + ".txt"
        chr_data_path = os.path.join(COO_folder_path, chr_COO_file_name)
        temp_frames, temp_index = utils.divide(chr_data_path, i, args['resolution'], args['genome_type'], args['COO_format'])
        frames_data.extend(temp_frames)
        index_data.extend(temp_index)
        print("chr" + str(i) + " is done!")

frames_data = np.stack(frames_data, axis = 0)
index_data = np.stack(index_data, axis = 0)
if not os.path.exists(args['output_folder_path']):
    os.makedirs(args['output_folder_path'])
np.save(os.path.join(args['output_folder_path'] , args['output_file_name'] + ".npy"), frames_data)
np.save(os.path.join(args['output_folder_path'] , args['output_file_name'] + "-index.npy"), index_data)

import os
import sys
import numpy as np
import pandas as pd
import argparse
import cooler


ap = argparse.ArgumentParser()
ap.add_argument("--input_type", help = "0 if it is cool and 1 if it is COO folder path", type = int)
ap.add_argument("--high_res_cool_path", help = "path of high resolution cool file")
ap.add_argument("--resolution", help = "resolution of HiC matrix to be used", type = int)
ap.add_argument("--COO_folder_path", help = "path of folder containing COO files to be down sampled")
ap.add_argument("--total_num_reads", help = "total number of reads in high resolution HiC data", type = int)
ap.add_argument("--num_sample_reads", help = "number of reads to be drawn as a sample data", type = int)
ap.add_argument("--low_res_cool_path", help = "path of low resolution cool file to get number of sample reads from")
ap.add_argument("--down_sample_ratio", help = "fraction of reads to be drawn as a sample data", type = int)
ap.add_argument("--output_folder_path", help = "path of folder to save down sample data", required = True)
args = vars(ap.parse_args())



if args['input_type'] == 0:
    high_res_cool = cooler.Cooler(args['high_res_cool_path'] + '::/resolutions/' + str(args['resolution']))
    # first it tries to find number of reads of Hi-C experiment (both intra and inter chromosomal) if it was not given
    try:
        if args['total_num_reads'] is not None:
            total_num_reads = args['total_num_reads']
        else:
            total_num_reads = sum(high_res_cool.pixels()[:,].iloc[:,2])
    except:
        print("high_res_cool_path is required when input_type is 0!")
    try:
        if args['num_sample_reads'] is not None:
            num_sample_reads = args['num_sample_reads']
        elif args['down_sample_ratio'] is not None:
            num_sample_reads = int(total_num_reads/args['down_sample_ratio'])
        else:
            low_res_cool = cooler.Cooler(args['low_res_cool_path'] + '::/resolutions/' + str(args['resolution']))
            num_sample_reads = sum(low_res_cool.pixels()[:,].iloc[:,2])
    except:
        print("At least one of num_sample_reads or low_res_cool_path or down_sample_ratio is required!")

    # vec_of_prob list is going to contain all non zero interactions from Hi-C matrices as a background distribution to sample from
    vec_of_prob = []
    for chrName in high_res_cool.chromnames:
        vec_of_prob.extend(high_res_cool.matrix(balance = False, as_pixels = True).fetch(chrName).iloc[:,2])
    num_inter_reads = total_num_reads - sum(vec_of_prob)
    vec_of_prob.append(num_inter_reads)
    vec_of_prob = [p/total_num_reads for p in vec_of_prob]
    print("start sampling!")
    # this function sample from calculated background distribution num_sample_reads we want
    down_sampled_counts = np.random.multinomial(num_sample_reads, vec_of_prob)
    print("sampling finished!")
    if not os.path.exists(args['output_folder_path']):
        os.makedirs(args['output_folder_path'])
    start_ind = 0
    # in this loop we are going to assign new values (sampled) of interactions to pixels
    for chrName in high_res_cool.chromnames:
        chr_pixel = high_res_cool.matrix(balance = False, as_pixels = True).fetch(chrName)
        pixel_size = chr_pixel.shape[0]
        offset = high_res_cool.extent(chrName)[0]
        new_pixel = np.column_stack((chr_pixel.iloc[:,0] - offset,
                                      chr_pixel.iloc[:,1] - offset,
                                      down_sampled_counts[start_ind:start_ind+pixel_size]))
        start_ind = start_ind + pixel_size
        np.savetxt(os.path.join(args['output_folder_path'], chrName + ".txt"), new_pixel, delimiter = "\t", fmt = "%i")
        print(chrName + " is done!")

else:

    try:
        total_num_reads = args['total_num_reads']
    except:
        print("total_num_reads is required when input type = 1!")
    try:
        if args['down_sample_ratio'] is not None:
            num_sample_reads = int(total_num_reads/args['down_sample_ratio'])
        else:
            num_sample_reads = args['num_sample_reads']
    except:
        print("At least one of num_sample_reads or down_sample_ratio is required!")
    try:
        COO_folder_path = args['COO_folder_path']
    except:
        print("Path of folder including high resolution COO files is required when input type = 1!")
    vec_of_prob = []
    chr_files_list = [f for f in os.listdir(COO_folder_path) if not f.startswith('.')]
    # it loops over all files in a folder to make down_sampled version of them
    for chr_file in chr_files_list:
        print("reading ", chr_file)
        chr_data = pd.read_csv(os.path.join(COO_folder_path, chr_file), delimiter = "\t", header = None)
        vec_of_prob.extend(chr_data.iloc[:,2])
    num_intra_reads = sum(vec_of_prob)
    num_inter_reads = total_num_reads - num_intra_reads
    vec_of_prob.append(num_inter_reads)
    vec_of_prob = [p/total_num_reads for p in vec_of_prob]
    print("start of sampling...")
    down_sampled_counts = np.random.multinomial(num_sample_reads, vec_of_prob)
    print("sampling finished!")
    if not os.path.exists(args['output_folder_path']):
        os.makedirs(args['output_folder_path'])
    start_ind = 0
    for chr_file in chr_files_list:
        chr_data = pd.read_csv(os.path.join(COO_folder_path, chr_file), delimiter = "\t", header = None)
        pixel_size = chr_data.shape[0]
        new_pixel = np.column_stack((chr_data.iloc[:,0],
                                      chr_data.iloc[:,1],
                                      down_sampled_counts[start_ind:start_ind+pixel_size]))
        start_ind = start_ind + pixel_size
        np.savetxt(os.path.join(args['output_folder_path'], chr_file), new_pixel, delimiter = "\t", fmt = "%i")
        print(chr_file + " is done!")






"""
Obs1: when we fetch a specific chromosome it means first columns belong to regions in that chromosome but
second column regions are through whole genome
Obs2: reads are not considered twice in files, for example when we fetch chr2 pixels, there are not interactions
between chr2 and chr1 any more.
Obs3: number of intra reads: 125015861, whole reads: 153752070 (in low resolution sample)
new_bins = high_res_cool.bins()
cooler.create_cooler(cool_uri = "/Users/neda/prostate-samples/PCa13266.down-sample.cool", bins = new_bins, pixels = new_pixel)
"""

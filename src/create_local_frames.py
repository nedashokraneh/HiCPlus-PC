### this script get cell type, resolution and list including chromosomes we want their data to include

import sys
import os
import numpy as np
import utils

dir_path = os.path.dirname(os.path.realpath(__file__))
root_path = dir_path + "/.."
data_path = root_path + "/data/"

resolution_dict = {
    5000: "5kb",
    10000: "10kb",
    25000: "25kb",
    50000: "50kb",
    100000: "100kb",
    250000: "250kb",
    500000: "500kb",
    1000000: "1mb"
}

def create_chr_frames(cell_type, chromosome_num, resolution_size): # for example GM12878_primary or GM12878_replicate could be cell_type there should be folder with this name in data folder

    input_file_name = "chr" + str(chromosome_num) + "_" + resolution_dict[resolution_size] + ".RAWobserved"
    input_file_path = data_path + cell_type + "/" + resolution_dict[resolution_size] + "_resolution_intrachromosomal/chr" + str(chromosome_num) + "/MAPQG0/" + input_file_name
    frames, index = utils.divide(input_file_path, chromosome_num, resolution_size)
    return (frames,index)

def create_frames(cell_type, chromosome_list, resolution_size):
    output_file_name = cell_type + "_chr" + str(min(chromosome_list)) + "-" + str(max(chromosome_list))
    output_file_path = data_path + "divided-data/" + output_file_name
    if os.path.exists(output_file_path):
        print("data already exists in ", output_file_path)
    else :
        frames_data = []
        index_data = []
        for chromosome_num in chromosome_list:
            temp_frames, temp_index = create_chr_frames(cell_type, chromosome_num, resolution_size)
            frames_data.append(temp_frames)
            index_data.append(index_data)
        frames_data = np.array(frames_data)
        index_data = np.array(index_data)
        np.save(output_file_path + ".npy", frames_data)
        np.save(output_file_path + "-index.npy", index_data)


def main():
    cell_type = "GM12878_primary"
    chromosome_list = [19,20,21,22]
    resolution_size = 10000
    create_frames(cell_type, chromosome_list, resolution_size)


if __name__ == "__main__":
    main()

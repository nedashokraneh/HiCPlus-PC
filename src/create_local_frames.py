### this script get cell type, resolution and list including chromosomes we want their data to include

import sys
import os
import numpy as np
import utils

dir_path = os.path.dirname(os.path.realpath(__file__))
root_path = dir_path + "/.."
data_path = root_path + "/data/"
is_down_sampled = 1
cell_type = "GM12878_primary"
training_chromosome_list = [i for i in range(1,18)]
test_chromosome_list = [i for i in range(18,23)]
resolution_size = 10000
down_sample_ratio = 16

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

def create_chr_frames(cell_type, chromosome_num): # for example GM12878_primary or GM12878_replicate could be cell_type there should be folder with this name in data folder
    if is_down_sampled == 1:
        input_file_name = "chr" + str(chromosome_num) + "_" + resolution_dict[resolution_size] + "_down.RAWobserved"
        input_file_path = data_path + cell_type + "/" + resolution_dict[resolution_size] + "_resolution_intrachromosomal_down" + str(down_sample_ratio) + "(rep2)/" + input_file_name
    else:
        input_file_name = "chr" + str(chromosome_num) + "_" + resolution_dict[resolution_size] + ".RAWobserved"
        input_file_path = data_path + cell_type + "/" + resolution_dict[resolution_size] + "_resolution_intrachromosomal/chr" + str(chromosome_num) + "/MAPQG0/" + input_file_name
    frames, index = utils.divide(input_file_path, chromosome_num, resolution_size)
    return (frames,index)

def create_frames(cell_type, chromosome_list):
    if is_down_sampled == 1:
        output_file_name = "chr" + str(min(chromosome_list)) + "-" + str(max(chromosome_list)) + "(down" + str(down_sample_ratio) + ")(rep2)"
    else:
        output_file_name = "chr" + str(min(chromosome_list)) + "-" + str(max(chromosome_list))
    output_file_path = data_path + "divided-data/" + cell_type + "/" + resolution_dict[resolution_size] + "_resolution/"
    if os.path.exists(output_file_path + output_file_name + ".npy"):
        print("data already exists in ", output_file_path + output_file_name)
    else :
        frames_data = []
        index_data = []
        for chromosome_num in chromosome_list:
            temp_frames, temp_index = create_chr_frames(cell_type, chromosome_num)
            frames_data.extend(temp_frames)
            index_data.extend(temp_index)
            print("chr", chromosome_num, " is done!")
        frames_data = np.stack(frames_data, axis = 0)
        index_data = np.stack(index_data, axis = 0)
        ### testing: cell type = "GM12878_primary", chromosome_list = [19,20,21,22], resolution_size = 10000 so following lines should result in 786
        ### and ['test' '21' '925' '925'] respectively.
        #print(frames_data[8719,16,16])
        #print(index_data[8719,:])
        if not os.path.exists(output_file_path):
            os.makedirs(output_file_path)
        np.save(output_file_path + output_file_name + ".npy", frames_data)
        np.save(output_file_path + output_file_name + "-index.npy", index_data)


def main():
    # create_frames(cell_type, training_chromosome_list)
    create_frames(cell_type, test_chromosome_list)

if __name__ == "__main__":
    main()

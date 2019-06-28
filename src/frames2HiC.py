import numpy as np
import os
from scipy.stats.stats import pearsonr
import matplotlib.pyplot as plt


chrN = 18
down_sample_ratio = 16
resolution_size = 10000
chrs_length = [249250621,243199373,198022430,191154276,180915260,171115067,159138663,146364022,141213431,135534747,135006516,133851895,115169878,107349540,102531392,90354753,81195210,78077248,59128983,63025520,48129895,51304566]
index_file_path = "../data/divided-data/GM12878_primary_down_chr18-22-index.npy"
enhanced_frames_path = "../data/enhanced-data/GM12878_primary_enhanced_chr18-22.npy"
low_res_HiC_file_path = "../data/GM12878_primary_down/10kb_resolution_intrachromosomal/chr" + str(chrN) + "_10kb_down.RAWobserved"
high_res_HiC_file_path = "../data/GM12878_primary/10kb_resolution_intrachromosomal/chr" + str(chrN) + "/MAPQG0/chr" + str(chrN) + "_10kb.RAWobserved"
low_res_HiC_matrix_file_path = "../data/GM12878_primary_down/10kb_resolution_intrachromosomal/chr" + str(chrN) + "_10kb_down.RAWobserved_npy_form_tmp.npy"
high_res_HiC_matrix_file_path = "../data/GM12878_primary/10kb_resolution_intrachromosomal/chr" + str(chrN) + "/MAPQG0/chr" + str(chrN) + "_10kb.RAWobserved_npy_form_tmp.npy"
total_length = int(chrs_length[chrN-1]/resolution_size) + 1
index = np.load(index_file_path)
enhanced_frames = np.load(enhanced_frames_path)
if os.path.exists(high_res_HiC_matrix_file_path):
    high_res_HiC_matrix = np.load(high_res_HiC_matrix_file_path)
else:
    high_res_HiC_matrix = utils.readSquareMatrix(high_res_HiC_file_path, total_length, resolution_size)
if os.path.exists(low_res_HiC_matrix_file_path):
    HiC_matrix = np.load(low_res_HiC_matrix_file_path)
else:
    HiC_matrix = utils.readSquareMatrix(low_res_HiC_file_path, total_length, resolution_size)
decoder = np.vectorize(lambda x: x.decode('UTF-8'))
index = decoder(index[:,1:]).astype(int)
chrN_index = np.where(index[:,0]==chrN)[0]
enhanced_HiC_matrix = HiC_matrix * down_sample_ratio
enhanced_HiC_matrix = enhanced_HiC_matrix.astype(float)
for i in chrN_index:
    x_pos = index[i,1]
    y_pos = index[i,2]
    enhanced_HiC_matrix[x_pos+6:x_pos+34,y_pos+6:y_pos+34] = enhanced_frames[i,:,:]

def vec_of_dist(matrix, x):
    return([matrix[i,i+x] for i in range(matrix.shape[1]-x)])

low_res_HiC_matrix = HiC_matrix * down_sample_ratio
highVSlow_corr_list = []
highVSenhanced_corr_list = []
for dist in range(100):
    low_res_vec = vec_of_dist(low_res_HiC_matrix, dist)
    high_res_vec = vec_of_dist(high_res_HiC_matrix, dist)
    enhanced_vec = vec_of_dist(enhanced_HiC_matrix, dist)
    highVSlow_corr_list.append(pearsonr(low_res_vec, high_res_vec)[0])
    highVSenhanced_corr_list.append(pearsonr(high_res_vec, enhanced_vec)[0])
plt.plot(highVSlow_corr_list, label = "highVSlow")
plt.plot(highVSenhanced_corr_list, label = "highVSenhanced")
plt.legend(loc='upper right', prop={'size': 5})
plt.show()


### creating N*3 array of coordinates list from enhanced matrix

output_file_path = "../data/enhanced-data/chr" + str(chrN) + "-enhanced.txt"
nonzero_indices = np.nonzero(enhanced_HiC_matrix)
source = nonzero_indices[0] * resolution_size
target = nonzero_indices[1] * resolution_size
weight = enhanced_HiC_matrix[nonzero_indices]
coordinate_list = np.transpose(np.array((source, target, weight)))
np.savetxt(output_file_path, coordinate_list, delimiter='\t')

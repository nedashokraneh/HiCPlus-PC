import os
import sys
#dir_path = os.path.dirname(os.path.realpath(__file__))
dir_path = "/Users/neda/HiCPlus_pytorch/src"
import numpy as np
import argparse
import cooler
import matplotlib.pyplot as plt
import torch
from torch.autograd import Variable
from scipy.stats.stats import pearsonr
model_path = dir_path + "/models"
utils_path = dir_path + "/utils"
sys.path.insert(0, model_path)
sys.path.insert(0, utils_path)
import model
import utils

ap = argparse.ArgumentParser()
ap.add_argument("--LowRes_matrix_type", help = "0 if it is a cool file and 1 if it is a COO", required = True, type = int)
ap.add_argument("--LowRes_matrix_path", help = "path of low resolution cool file or low resolution COO", required = True)
ap.add_argument("--HighRes_matrix_type", help = "0 if it is a cool file and 1 if it is a COO", required = True, type = int)
ap.add_argument("--HighRes_matrix_path", help = "path of high resolution cool file or high resolution COO", required = True)
ap.add_argument("--resolution", help = "resolution needed if input is cool file", type = int)
ap.add_argument("--chr_num", help = "chromosome number to be used for evaluation", type = int)
ap.add_argument("--genome_type", help = "hg19 or hg38")
ap.add_argument("--network_path", help = "path of ConvNet learned for enhancing data")
ap.add_argument("--result_path", help = "path of folder to save result in")
args = vars(ap.parse_args())

use_gpu = 0
Net = model.Net(40, 28)
Net.load_state_dict(torch.load(args['network_path']))
Net = Net.float()
if use_gpu:
    Net = Net.cuda()

if args['LowRes_matrix_type'] == 0:
    low_cool = cooler.Cooler(args['LowRes_matrix_path'] + '::/resolutions/' + str(args['resolution']))
    low_chr_mat = low_cool.matrix(balance = False).fetch("chr" + str(args['chr_num'])).astype(float)
    low_chr_mat[np.isnan(low_chr_mat)] = 0
    chr_frames, chr_indices = utils.divide2(low_chr_mat,args['chr_num'])
    enhanced_chr_mat = low_cool.matrix(balance = False).fetch("chr" + str(args['chr_num'])).astype(float)
    enhanced_chr_mat[np.isnan(enhanced_chr_mat)] = 0
    """
    average_chr_mat = low_cool.matrix(balance = False).fetch("chr" + str(args['chr_num'])).astype(float)
    average_chr_mat[np.isnan(average_chr_mat)] = 0
    """


else:
    chr_frames, chr_indices = utils.divide(args['LowRes_matrix_path'], args['chr_num'], args['resolution'], args['genome_type'])
    low_chr_mat = np.load(args['LowRes_matrix_path'] + '_npy_form_tmp.npy')
    enhanced_chr_mat = np.load(args['LowRes_matrix_path'] + '_npy_form_tmp.npy')
    # average_chr_mat = np.load(args['LowRes_matrix_path'] + '_npy_form_tmp.npy')
chr_frames = np.stack(chr_frames, axis = 0)
chr_indices = np.stack(chr_indices, axis = 0)
chr_frames = np.expand_dims(chr_frames, axis = 1)
lowres_set = torch.from_numpy(chr_frames).float()
enhanced_set = Net(Variable(lowres_set))
enhanced_set = enhanced_set.data.cpu().numpy()
enhanced_set = np.reshape(enhanced_set, (enhanced_set.shape[0], enhanced_set.shape[2], enhanced_set.shape[3]))

for i in range(chr_indices.shape[0]):
    x_pos = chr_indices[i,1]
    y_pos = chr_indices[i,2]
    enhanced_chr_mat[x_pos+6:x_pos+34,y_pos+6:y_pos+34] = enhanced_set[i,:,:]

iu = np.triu_indices(enhanced_chr_mat.shape[0],1)
il = (iu[1],iu[0])
enhanced_chr_mat[il]=enhanced_chr_mat[iu]
"""
chr_length = average_chr_mat.shape[0]
low_chr_mat2 = np.zeros((chr_length+6,chr_length+6))
low_chr_mat2[3:chr_length+3,3:chr_length+3] = low_chr_mat
for i1 in range(chr_length):
    for i2 in range(chr_length):
        average_chr_mat[i1,i2] = np.mean(low_chr_mat2[i1:i1+6,i2:i2+6])
"""
if args['HighRes_matrix_type'] == 0:
    high_cool = cooler.Cooler(args['HighRes_matrix_path'] + '::/resolutions/' + str(args['resolution']))
    high_chr_mat = high_cool.matrix(balance = False).fetch("chr" + str(args['chr_num'])).astype(float)
    high_chr_mat[np.isnan(high_chr_mat)] = 0
else:
    _, _ = utils.divide(args['HighRes_matrix_path'], args['chr_num'], args['resolution'], args['genome_type'])
    high_chr_mat = np.load(args['HighRes_matrix_path'] + '_npy_form_tmp.npy')


def vec_of_dist(matrix, x):
    return([matrix[i,i+x] for i in range(matrix.shape[1]-x)])

highVSlow_corr_list = []
highVSenhanced_corr_list = []
#highVSaverage_corr_list = []
for dist in range(100):
    low_res_vec = vec_of_dist(low_chr_mat, dist)
    high_res_vec = vec_of_dist(high_chr_mat, dist)
    enhanced_vec = vec_of_dist(enhanced_chr_mat, dist)
    #average_vec = vec_of_dist(average_chr_mat, dist)
    highVSlow_corr_list.append(pearsonr(low_res_vec, high_res_vec)[0])
    highVSenhanced_corr_list.append(pearsonr(high_res_vec, enhanced_vec)[0])
    #highVSaverage_corr_list.append(pearsonr(high_res_vec, average_vec)[0])
fig = plt.figure()
plt.plot(highVSlow_corr_list, label = "highVSlow")
plt.plot(highVSenhanced_corr_list, label = "highVSenhanced")
#plt.plot(highVSaverage_corr_list, label = "highVSaverage")
plt.legend(loc='upper right', prop={'size': 5})
plt.show()
fig.savefig(os.path.join(args['result_path'],'correlation-plot.png'))

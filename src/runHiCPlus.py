# Author: Yan Zhang
# Email: zhangyan.cse (@) gmail.com

import os
import sys
dir_path = os.path.dirname(os.path.realpath(__file__))

import numpy as np
import argparse
import cooler

model_path = dir_path + "/models"
utils_path = dir_path + "/utils"
sys.path.insert(0, model_path)
sys.path.insert(0, utils_path)
import model
import utils


ap = argparse.ArgumentParser()
ap.add_argument("-i", "--input_path", required = True, help = "path of cool file need to be enhanced")
ap.add_argument("-o", "--output_path", required = True, help = "path of output cool file to be saved") #default: where the input file was
ap.add_argument("-m", "--model_path", required = True, help = "path of learned model for corresponding down sample ratio")
ap.add_argument("-r", "--resolution", required = True, help = "resolution size of HiC contact matrix", default = "10000")
args = vars(ap.parse_args())


print("check: ", model.conv2d1_filters_numbers)
use_gpu = 0

epochs = 10
HiC_max_value = 100
for i in range(23):

cooler_path = args['input_path'] + "::/resolutions/" + args['resolution']
c = cooler.Cooler(cooler_path)
chr_mat = c.matrix(balance = False).fetch("chr" + arg['chr_num'])
chr_mat[np.isnan(chr18_mat)] = 0
low_resolution_samples, indices = utils.divide2(chr_mat)
low_resolution_samples = np.stack(low_resolution_samples, axis = 0)
indices = np.stack(indices, axis = 0)
low_resolution_samples = np.expand_dims(low_resolution_samples, axis=1)
batch_size = low_resolution_samples.shape[0]

# Reshape the high-quality Hi-C sample as the target value of the training.
sample_size = low_resolution_samples.shape[-1]
padding = conv2d1_filters_size + conv2d2_filters_size + conv2d3_filters_size - 3
half_padding = padding / 2
output_length = sample_size - padding


print(low_resolution_samples.shape)

lowres_set = data.TensorDataset(torch.from_numpy(low_resolution_samples), torch.from_numpy(np.zeros(low_resolution_samples.shape[0])))
lowres_loader = torch.utils.data.DataLoader(lowres_set, batch_size=batch_size, shuffle=False)



model = model.Net(40, 28)
model.load_state_dict(torch.load(args['model_path']))

if use_gpu:
    model = model.cuda()

model = model.float()
frames_prediction = model(torch.from_numpy(low_resolution_samples).float())
frames_prediction = frames_prediction.data.cpu().numpy()
frames_prediction = np.reshape(frames_prediction, (frames_prediction.shape[0], frames_prediction.shape[2], frames_prediction.shape[3]))
enhanced_chr_mat = chr_mat.astype(float)
for i in range(indices.shape[0]):
    x_pos = indices[i,0]
    y_pos = indices[i,1]
    enhanced_chr18_mat[x_pos+6:x_pos+34,y_pos+6:y_pos+34] = y_prediction[i,:,:]


np.save(output_file_path, y_predict)

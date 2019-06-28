import os
import sys
dir_path = os.path.dirname(os.path.realpath(__file__))

model_path = dir_path + "/models"
sys.path.insert(0, model_path)
import model

import numpy as np
import matplotlib.pyplot as plt
import pickle
import gzip
from torch.utils import data
import torch
import torch.optim as optim
from torch.autograd import Variable
from time import gmtime, strftime
import torch.nn as nn
from scipy.stats.stats import pearsonr
import argparse

ap = argparse.ArgumentParser()
ap.add_argument("-d", "--high_res_path", required = True, help = "path of high resolution npy file")
ap.add_argument("-l", "--low_res_path", required = True, help = "path of corresponding low resolution npy file") #default: where the input file was
ap.add_argument("-m", "--model_path", required = True, help = "location where network should be saved in")
ap.add_argument("-n", "--model_name", required = True, help = "name of model to be saved")
args = vars(ap.parse_args())

use_gpu = 0
down_sample_ratio = 16
epochs = 3
HiC_max_value = 100 # ?????
batch_size = 256

def corr_highVSlow(index,data1,data2):
    return pearsonr(data1[index,0,:,:].flatten(),data2[index,0,:,:].flatten())[0]

low_resolution_samples = np.load(args['low_res_path'], "r").astype(np.float32) * down_sample_ratio
low_resolution_samples = np.expand_dims(low_resolution_samples, axis=1)
high_resolution_samples = np.load(args['high_res_path'], "r").astype(np.float32)
high_resolution_samples = np.expand_dims(high_resolution_samples, axis=1)
#low_resolution_samples = np.minimum(HiC_max_value, low_resolution_samples)
#high_resolution_samples = np.minimum(HiC_max_value, high_resolution_samples)
sample_size = low_resolution_samples.shape[-1]
high_resolution_samples = high_resolution_samples[:,:,model.half_padding:(sample_size-model.half_padding),model.half_padding:(sample_size-model.half_padding)]
lowres_set = data.TensorDataset(torch.from_numpy(low_resolution_samples))
lowres_loader = torch.utils.data.DataLoader(lowres_set, batch_size=batch_size, shuffle=False)
hires_set = data.TensorDataset(torch.from_numpy(high_resolution_samples))
hires_loader = torch.utils.data.DataLoader(hires_set, batch_size=batch_size, shuffle=False)

Net = model.Net(40, 28)
if use_gpu:
    Net = Net.cuda()

optimizer = optim.SGD(Net.parameters(), lr = 0.00001)
_loss = nn.MSELoss()
Net.train()

running_loss = 0.0
losslist = []

# write the log file to record the training process
log = open('HindIII_train.txt', 'w')
for epoch in range(0, epochs):
    for i, (v1, v2) in enumerate(zip(lowres_loader, hires_loader)):
        if (i == len(lowres_loader) - 1):
            continue
        _lowRes = Variable(v1[0])
        _highRes = Variable(v2[0])
        if use_gpu:
            _lowRes = _lowRes.cuda()
            _highRes = _highRes.cuda()
        optimizer.zero_grad()
        y_prediction = Net(_lowRes)
        loss = _loss(y_prediction, _highRes)
        loss.backward()
        optimizer.step()
        running_loss += loss.item()

    print '-------', i, epoch, running_loss/i, strftime("%Y-%m-%d %H:%M:%S", gmtime())
    losslist.append(running_loss/i)
    log.write(str(epoch) + ', ' + str(running_loss/i,) + '\n')
    running_loss = 0.0

if not os.path.exists(args['model_path']):
    os.makedirs(args['model_path'])
torch.save(Net.state_dict(), args['model_path'] + "/" + args['model_name'])

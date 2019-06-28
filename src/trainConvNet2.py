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

use_gpu = 0
conv2d1_filters_numbers = 8
conv2d1_filters_size = 9
conv2d2_filters_numbers = 8
conv2d2_filters_size = 1
conv2d3_filters_numbers = 1
conv2d3_filters_size = 5
down_sample_ratio = 16
epochs = 101
HiC_max_value = 100 # ?????
batch_size = 256
resolution_size = 10000
cell_type = "GM12878_primary"
chr_num_min = 1
chr_num_max = 17
model_path = "../models/19-06-17(res10kb,down180)/"
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

low_resolution_path = "../data/divided-data/" + cell_type + "/" + resolution_dict[resolution_size] + "_resolution/chr" + str(chr_num_min) + "-" + str(chr_num_max) + "(down" + str(down_sample_ratio) + ").npy"
high_resolution_path = "../data/divided-data/" + cell_type + "/" + resolution_dict[resolution_size] + "_resolution/chr" + str(chr_num_min) + "-" + str(chr_num_max) + ".npy"
low_resolution_samples = np.load(low_resolution_path).astype(np.float32) * down_sample_ratio
low_resolution_samples = np.expand_dims(low_resolution_samples, axis=1)
high_resolution_samples = np.load(high_resolution_path).astype(np.float32)
high_resolution_samples = np.expand_dims(high_resolution_samples, axis=1)
high_resolution_samples = np.minimum(high_resolution_samples, HiC_max_value)
low_resolution_samples = np.minimum(low_resolution_samples, HiC_max_value)

sample_size = low_resolution_samples.shape[-1]
padding = conv2d1_filters_size + conv2d2_filters_size + conv2d3_filters_size - 3
half_padding = padding / 2
high_resolution_samples = high_resolution_samples[:,:,half_padding:(sample_size-half_padding),half_padding:(sample_size-half_padding)]

lowres_set = data.TensorDataset(torch.from_numpy(low_resolution_samples), torch.from_numpy(np.zeros(low_resolution_samples.shape[0])))
lowres_loader = torch.utils.data.DataLoader(lowres_set, batch_size=batch_size, shuffle=False)

hires_set = data.TensorDataset(torch.from_numpy(high_resolution_samples), torch.from_numpy(np.zeros(high_resolution_samples.shape[0])))
hires_loader = torch.utils.data.DataLoader(hires_set, batch_size=batch_size, shuffle=False)



Net = model.Net(40, 28)
if use_gpu:
    Net = Net.cuda()

optimizer = optim.SGD(Net.parameters(), lr = 0.00001)
_loss = nn.MSELoss()
Net.train()
running_loss = 0.0
losslist = []
log = open('HindIII_train.txt', 'w')

for epoch in range(0, epochs):
    # iterate over two lists and their indices using enumerate together with zip
    # lowres_loader is list of batches
    for i, (v1, v2) in enumerate(zip(lowres_loader, hires_loader)):
        # probably it is for skipping last incomplete batch
        if (i == len(lowres_loader) - 1):
            continue


        # v1 is list with length = 2. v1[0] is data tensor so with shape 256*1*40*40. v1[1] is vector of 256 zeros because pf line 85 but what's the reason?
        _lowRes, _ = v1
        _highRes, _ = v2
        # print "_lowres:", _lowRes, "\n shape: ", _lowRes.shape

        _lowRes = Variable(_lowRes)
        _highRes = Variable(_highRes)


        if use_gpu:
            _lowRes = _lowRes.cuda()
            _highRes = _highRes.cuda()
        optimizer.zero_grad()
        y_prediction = Net(_lowRes)
        loss = _loss(y_prediction, _highRes)
        loss.backward()
        optimizer.step()
        #print(loss.item())
        running_loss += loss.item()

    print '-------', i, epoch, running_loss/i, strftime("%Y-%m-%d %H:%M:%S", gmtime())
    losslist.append(running_loss/i)
    log.write(str(epoch) + ', ' + str(running_loss/i,) + '\n')
    running_loss = 0.0
    # save the model every 10 epoches

if not os.path.exists(model_path):
    os.makedirs(model_path)
torch.save(Net.state_dict(), model_path + str(epoch))

plt.plot(losslist)
#plt.axis([0, 6, 0, 20])
plt.show()


"""

test_y_prediction = Net(torch.from_numpy(test_low_resolution_samples)).detach().numpy()
    train_y_prediction = Net(torch.from_numpy(train_low_resolution_samples)).detach().numpy()
    ### evaluating part ###
    test_cor = np.nanmean([corr_highVSlow(i,test_y_prediction,test_high_resolution_samples) for i in range(test_low_resolution_samples.shape[0])])
    train_cor = np.nanmean([corr_highVSlow(i,train_y_prediction,train_high_resolution_samples) for i in range(train_low_resolution_samples.shape[0])])
    print("test data correlation: ", test_cor)
    print("train data correlation: ", train_cor)
    test_cor_list.append(test_cor)
    train_cor_list.append(train_cor)

"""

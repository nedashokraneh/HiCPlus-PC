# Author: Yan Zhang
# Email: zhangyan.cse (@) gmail.com



import sys
import numpy as np
import matplotlib.pyplot as plt
import pickle
import os
import gzip
import model
from torch.utils import data
import torch
import torch.optim as optim
from torch.autograd import Variable
from time import gmtime, strftime
import sys
import torch.nn as nn
from scipy.stats.stats import pearsonr

use_gpu = 0

conv2d1_filters_numbers = 8
conv2d1_filters_size = 9
conv2d2_filters_numbers = 8
conv2d2_filters_size = 1
conv2d3_filters_numbers = 1
conv2d3_filters_size = 5


down_sample_ratio = 16
epochs = 10
HiC_max_value = 100 # ?????
batch_size = 256

def corr_highVSlow(index,data1,data2):
    return pearsonr(data1[index,0,:,:].flatten(),data2[index,0,:,:].flatten())[0]
# This block is the actual training data used in the training. The training data is too large to put on Github, so only toy data is used.
# cell = "GM12878_replicate"
# chrN_range1 = '1_8'
# chrN_range = '1_8'

# low_resolution_samples = np.load(gzip.GzipFile('/home/zhangyan/SRHiC_samples/'+cell+'down16_chr'+chrN_range+'.npy.gz', "r")).astype(np.float32) * down_sample_ratio
# high_resolution_samples = np.load(gzip.GzipFile('/home/zhangyan/SRHiC_samples/original10k/'+cell+'_original_chr'+chrN_range+'.npy.gz', "r")).astype(np.float32)

# low_resolution_samples = np.minimum(HiC_max_value, low_resolution_samples)
# high_resolution_samples = np.minimum(HiC_max_value, high_resolution_samples)


#low_resolution_samples = np.load(gzip.GzipFile('../../data/GM12878_replicate_down16_chr19_22.npy.gz', "r")).astype(np.float32) * down_sample_ratio
#high_resolution_samples = np.load(gzip.GzipFile('../../data/GM12878_replicate_original_chr19_22.npy.gz', "r")).astype(np.float32)

#low_resolution_samples = np.load(gzip.GzipFile('/home/zhangyan/SRHiC_samples/IMR90_down_HINDIII16_chr1_8.npy.gz', "r")).astype(np.float32) * down_sample_ratio
#high_resolution_samples = np.load(gzip.GzipFile('/home/zhangyan/SRHiC_samples/original10k/_IMR90_HindIII_original_chr1_8.npy.gz', "r")).astype(np.float32)

low_resolution_samples = np.load('../data/GM12878_replicatedown16_chr19_22.npy', "r").astype(np.float32) * down_sample_ratio
high_resolution_samples = np.load('../data/GM12878_replicate_original_chr19_22.npy', "r").astype(np.float32)
num_samples = low_resolution_samples.shape[0]
print(low_resolution_samples.shape[0])

# most of the high resolution sample entries are bigger than HiC_max_value. What it is supposed to do?
#low_resolution_samples = np.minimum(HiC_max_value, low_resolution_samples)
#high_resolution_samples = np.minimum(HiC_max_value, high_resolution_samples)


# Reshape the high-quality Hi-C sample as the target value of the training.
sample_size = low_resolution_samples.shape[-1]
padding = conv2d1_filters_size + conv2d2_filters_size + conv2d3_filters_size - 3
half_padding = padding / 2
output_length = sample_size - padding

"""
why like that? takes time and reduce dimension which cause warning
Y = []
for i in range(high_resolution_samples.shape[0]):
    no_padding_sample = high_resolution_samples[i][0][half_padding:(sample_size-half_padding) , half_padding:(sample_size - half_padding)]
    Y.append(no_padding_sample)
Y = np.array(Y).astype(np.float32)
"""
high_resolution_samples = high_resolution_samples[:,:,half_padding:(sample_size-half_padding),half_padding:(sample_size-half_padding)]
#print low_resolution_samples.shape, Y.shape

# what's the reason of np.zeros?
# length of loader is length of lower_set divided by 256 which is the batch size
test_indices = np.arange(0,num_samples, 5)
train_indices = np.setdiff1d(np.arange(num_samples), np.arange(0,num_samples, 5))
train_low_resolution_samples = low_resolution_samples[train_indices,:,:,:]
test_low_resolution_samples = low_resolution_samples[test_indices,:,:,:]
train_high_resolution_samples = high_resolution_samples[train_indices,:,:,:]
test_high_resolution_samples = high_resolution_samples[test_indices,:,:,:]

lowres_set = data.TensorDataset(torch.from_numpy(train_low_resolution_samples), torch.from_numpy(np.zeros(train_low_resolution_samples.shape[0])))
lowres_loader = torch.utils.data.DataLoader(lowres_set, batch_size=batch_size, shuffle=False)

hires_set = data.TensorDataset(torch.from_numpy(train_high_resolution_samples), torch.from_numpy(np.zeros(train_high_resolution_samples.shape[0])))
hires_loader = torch.utils.data.DataLoader(hires_set, batch_size=batch_size, shuffle=False)

Net = model.Net(40, 28)
if use_gpu:
    Net = Net.cuda()

optimizer = optim.SGD(Net.parameters(), lr = 0.00001)
_loss = nn.MSELoss()
Net.train()

running_loss = 0.0
running_loss_validate = 0.0
reg_loss = 0.0
losslist = []
train_cor_list = []
test_cor_list = []
# write the log file to record the training process
log = open('HindIII_train.txt', 'w')

for epoch in range(0, 100):
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

    print '-------', i, epoch, running_loss/batch_size, strftime("%Y-%m-%d %H:%M:%S", gmtime())
    losslist.append(running_loss/batch_size)
    log.write(str(epoch) + ', ' + str(running_loss/batch_size,) + '\n')
    running_loss = 0.0
    running_loss_validate = 0.0
    # save the model every 100 epoches
    if (epoch % 100 == 0):
        if not os.path.exists('../models/19-05-29/'):
            os.makedirs('../models/19-05-29/')
        torch.save(Net.state_dict(), '../models/19-05-29/' + str(epoch))
        pass
    test_y_prediction = Net(torch.from_numpy(test_low_resolution_samples)).detach().numpy()
    train_y_prediction = Net(torch.from_numpy(train_low_resolution_samples)).detach().numpy()
    ### evaluating part ###
    test_cor = np.nanmean([corr_highVSlow(i,test_y_prediction,test_high_resolution_samples) for i in range(test_low_resolution_samples.shape[0])])
    train_cor = np.nanmean([corr_highVSlow(i,train_y_prediction,train_high_resolution_samples) for i in range(train_low_resolution_samples.shape[0])])
    print("test data correlation: ", test_cor)
    print("train data correlation: ", train_cor)
    test_cor_list.append(test_cor)
    train_cor_list.append(train_cor)


plt.plot(test_cor_list, label = "test_cor")
plt.plot(train_cor_list, label = "train_cor")
#plt.axis([0, 6, 0, 20])
plt.show()

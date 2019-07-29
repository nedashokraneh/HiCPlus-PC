# input of this script is npy files including frames which are samples


import os
import sys
dir_path = os.path.dirname(os.path.realpath(__file__))

model_path = dir_path + "/models"
sys.path.insert(0, model_path)
import model
import model2
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
ap.add_argument("--LowRes_frames_path", help = "path of npy file including LowRes frames", required = True)
ap.add_argument("--HighRes_frames_path", help = "path of npy file including HighRes frames", required = True)
ap.add_argument("--batch_size", help = "batch size in learning process", required = True, type = int)
ap.add_argument("--epochs", help = "number of epochs", required = True, type = int)
ap.add_argument("--output_network_path", help = "path of folder to save learned network", required = True)
ap.add_argument("--network_name", help = "name of network for saving", required = True)
args = vars(ap.parse_args())

use_gpu = 0
# down_sample_ratio = 16 learn with scaling or not????
epochs = args['epochs']
#HiC_max_value = 100 # ?????
batch_size = args['batch_size']

low_resolution_samples = np.load(args['LowRes_frames_path']).astype(np.float32)
low_resolution_samples = np.expand_dims(low_resolution_samples, axis=1)
high_resolution_samples = np.load(args['HighRes_frames_path']).astype(np.float32)
high_resolution_samples = np.expand_dims(high_resolution_samples, axis=1)
#high_resolution_samples = np.minimum(high_resolution_samples, HiC_max_value)
#low_resolution_samples = np.minimum(low_resolution_samples, HiC_max_value)
print(low_resolution_samples.shape)
print(high_resolution_samples.shape)
sample_size = low_resolution_samples.shape[-1]
half_padding = model.half_padding
lb = int(half_padding)
ub = int(sample_size - half_padding)
high_resolution_samples = high_resolution_samples[:,:,lb:ub,lb:ub]
lowres_loader = torch.utils.data.DataLoader(torch.from_numpy(low_resolution_samples), batch_size=batch_size, shuffle=False)
hires_loader = torch.utils.data.DataLoader(torch.from_numpy(high_resolution_samples), batch_size=batch_size, shuffle=False)



Net = model.Net(40, 28)
if use_gpu:
    Net = Net.cuda()
optimizer = optim.SGD(Net.parameters(), lr = 0.00001)
_loss = nn.MSELoss()
Net.train()
running_loss = 0.0
losslist = []
for epoch in range(0, epochs):
    # iterate over two lists and their indices using enumerate together with zip
    # lowres_loader is list of batches
    for i, (v1, v2) in enumerate(zip(lowres_loader, hires_loader)):
        # probably it is for skipping last incomplete batch
        if (i == len(lowres_loader) - 1):
            continue
        _lowRes = Variable(v1)
        _highRes = Variable(v2)
        if use_gpu:
            _lowRes = _lowRes.cuda()
            _highRes = _highRes.cuda()
        optimizer.zero_grad()
        y_prediction = Net(_lowRes)
        loss = _loss(y_prediction, _highRes)
        loss.backward()
        optimizer.step()
        running_loss += loss.item()
    print ('-------', i, epoch, running_loss/i, strftime("%Y-%m-%d %H:%M:%S", gmtime()))
    losslist.append(running_loss/i)
    running_loss = 0.0

if not os.path.exists(args['output_network_path']):
    os.makedirs(args['output_network_path'])
torch.save(Net.state_dict(), os.path.join(args['output_network_path'],args['network_name']))
fig = plt.figure()
plt.plot(losslist)
fig.savefig(os.path.join(args['output_network_path'] ,args['network_name'] + "-loss.png"))

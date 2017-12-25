import numpy as np
import pandas as pd 
import os
from functions import *
import scipy.io 
import neuroseries as nts
import matplotlib.pyplot as plt
data_directory = 'data/'
files = os.listdir(data_directory)
generalinfo = scipy.io.loadmat(data_directory+'Mouse12-120806_GeneralInfo.mat')
shankstructure = loadShankStructure(generalinfo)
spikes,shanks = loadSpikeData(data_directory+'Mouse12-120806_SpikeData.mat', shankstructure['thalamus'])
thalamus_index = list(spikes.keys())

#compute firing rates for each neuron for 0.1 second time bins during the first rem epoch
HD_index = loadHDCellInfo(data_directory+'Mouse12-120806_HDCells.mat', thalamus_index)
HD_spikes = {}
for i in HD_index:
	HD_spikes[i] = spikes[i]
rem_ep = loadEpoch(data_directory, 'rem')
first_rem = rem_ep.loc[0]
first_rem = nts.IntervalSet(start = first_rem['start'], end = first_rem['end'])
binsize = 0.1
rem_start = first_rem.as_units('s')['start'].values[0]
rem_end = first_rem.as_units('s')['end']
duration = rem_end-rem_start
duration = duration + 0.1
num_points = duration/binsize
num_points = int(num_points)
for j in HD_index:
	my_neuron = HD_spikes[j]
	my_neuron = my_neuron.restrict(first_rem)
	first_spike = my_neuron.as_units('s').index[0]
	last_spike = my_neuron.as_units('s').index[-1]
	firing_rate = np.zeros(num_points)
	for i in range(num_points):
		start = first_rem.as_units('s')['start'] + i*binsize
		end = start + binsize
		if ((start.values[0]>= first_spike) and (last_spike>=end.values[0])):
			spikes_in_interval = my_neuron.as_units('s').loc[start.values[0]:end.values[0]]
			firing_rate[i]= len(spikes_in_interval)
	firing_rate = nts.Tsd(t = np.arange(rem_start, rem_end, binsize), d = firing_rate, time_units = 's')
	if j == HD_index[0]:
		sleep_spikes = firing_rate
	else:
		sleep_spikes = np.vstack([sleep_spikes, firing_rate])
sleep_spikes = sleep_spikes.transpose()
time_bins = np.arange(rem_start, rem_end, binsize)
sleep_spikes = nts.TsdFrame(t = time_bins, d = sleep_spikes, time_units = 's')
sleep_spikes.columns = HD_index

#Compute tuning curves for each neuron using the wake epoch
mouse_position = np.genfromtxt(data_directory+'Mouse12-120806_PosHD.txt')
mouse_HD = nts.TsdFrame(t = mouse_position[:,0], d = mouse_position[:,3], time_units = 's')
wake_ep = loadEpoch(data_directory, 'wake')
wake_start = wake_ep.as_units('s')['start'].values[0]
wake_end = wake_ep.as_units('s')['end'].values[0]
duration = wake_end - wake_start
duration = duration + 0.1
num_points = duration/binsize
num_points = int(num_points)
head_direction = np.zeros(num_points)
#head direction at each time bin
for i in range(num_points):
	start = wake_start + i*binsize
	end = start + binsize
	HD_in_interval = mouse_HD.as_units('s').loc[start:end]
	avg_HD = np.mean(HD_in_interval)
	head_direction[i] = avg_HD
#firing rate for each time bin
for j in HD_index:
	my_neuron = HD_spikes[j]
	my_neuron = my_neuron.restrict(wake_ep)
	first_spike = my_neuron.as_units('s').index[0]
	last_spike = my_neuron.as_units('s').index[-1]
	firing_rate = np.zeros(num_points)
	for i in range(num_points):
		start = wake_start + i*binsize
		end = start + binsize
		if ((start>= first_spike) and (last_spike>=end)):
			spikes_in_interval = my_neuron.as_units('s').loc[start:end]
			firing_rate[i] = len(spikes_in_interval)
	firing_rate = nts.Tsd(t = np.arange(wake_start, wake_end, binsize), d = firing_rate, time_units = 's')
	if j == HD_index[0]:
		wake_spikes = firing_rate
	else:
		wake_spikes = np.vstack([wake_spikes, firing_rate])
data = nts.TsdFrame(t = np.arange(wake_start, wake_end, binsize), d = np.vstack([head_direction, wake_spikes]).transpose(), time_units = 's')
data = data[~data.isnull()[0].values]
columns_array = np.zeros(HD_index.size + 1)
columns_array = columns_array.astype(str)
columns_array[0] = 'Head Direction'
for i in range(HD_index.size):
	columns_array[i+1] = HD_index[i]
data.columns = columns_array
angular_bins = np.linspace(0, 2*np.pi+0.001, 61)
tuning_curves = np.zeros((HD_index.size, 60))
for i in range(HD_index.size):
	for j in range(60):
		left_border = angular_bins[j]
		right_border = angular_bins[j+1]
		index = np.logical_and(data['Head Direction']>left_border, data['Head Direction']<=right_border).values
		my_bin = data[index][HD_index[i].astype(str)]
		tuning_curves[i,j] = np.mean(my_bin)
tuning_curves = pd.DataFrame(data = tuning_curves)
angles = np.zeros(60)
for i in range(angles.size):
	angles[i] = np.mean([angular_bins[i], angular_bins[i+1]])

tuning_curves.columns = angles
tuning_curves.index = HD_index
angle_probabilities = np.ones(60)
for i in range(60):
	for j in HD_index:
		num_spikes = sleep_spikes[j][3905618400]
		angle_probabilities[i] = angle_probabilities[i]*(tuning_curves[angles[i]].loc[j]**num_spikes)
max_probability = np.max(angle_probabilities)

	

#tuning curves is a data frame with the neurons as the rows and the angular bins as the colums 
#showing probability of each neuron firing at each head direction





		












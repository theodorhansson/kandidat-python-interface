#%%

import os, sys
if os.path.dirname(os.path.dirname(os.path.realpath(__file__))) not in sys.path:
    sys.path.append(os.path.dirname(os.path.dirname(os.path.realpath(__file__))))
    
    
import matplotlib.pyplot as plt
import numpy as np

from scipy.signal import butter, filtfilt
from spr_project.plot_measurement_data.spr_plot_functions import load_measurement_data, plot_single_trace, plot_all_traces, savgol

def butterworth(data):
    sampling_freq=1/5
    cutoff_freq=1/300
    order=1
    nyquist_freq = 0.5 * sampling_freq
    normalized_cutoff_freq = cutoff_freq/nyquist_freq
    b, a = butter(order, normalized_cutoff_freq, btype='lowpass')
    data = filtfilt(b, a, data)
    return data


#%% Attempt 5 at bio-measurements

## Load all data
spr_data_folder = 'spr_measurements_240123'
vcsels = ['VCSEL_0', 'VCSEL_1', 'VCSEL_2']


files_to_plot = np.array([4], dtype=int)

mov_avg_window = 6
sav_gol_window = 15
sav_gol_order  = 5
plot_raw, plot_moving_avg, plot_savgol, frame_time = load_measurement_data(spr_data_folder, files_to_plot, vcsels,
                                                               mov_avg_window=mov_avg_window,
                                                               sav_gol_window=sav_gol_window, sav_gol_order=sav_gol_order)

time_scale = 60
savgol_window = 100
savgol_order  = 9

remove_air = 20

## Time trace
time_trace = plot_raw[vcsels[0]][0, :]
time_trace = time_trace[0:-remove_air]

## Li
reaction_start_time = 3600
reaction_stopp_time = 5300
reaction_start = np.argmin(np.abs(time_trace - reaction_start_time))
reaction_stopp = np.argmin(np.abs(time_trace - reaction_stopp_time))

time_trace_reaction = time_trace[reaction_start:reaction_stopp]

## Set start time to zero
time_trace_reaction_zeroed = (time_trace_reaction - time_trace_reaction[0])/time_scale


# High concentration
spr_high   = plot_raw[vcsels[0]][1, :]
spr_high   = spr_high[0:-remove_air]
spr_high_reaction = spr_high[reaction_start:reaction_stopp]
spr_high_reaction = spr_high_reaction - spr_high_reaction[0]
spr_high_reaction_savgol = savgol(spr_high_reaction, savgol_window, savgol_order)
spr_high_reaction_butter = butterworth(spr_high_reaction)


# Low concentration
spr_low    = plot_raw[vcsels[1]][1, :]
spr_low    = spr_low[0:-remove_air]
spr_low_reaction = spr_low[reaction_start:reaction_stopp]
spr_low_reaction = spr_low_reaction - spr_low_reaction[0]
spr_low_reaction_savgol = savgol(spr_low_reaction, savgol_window, savgol_order)
spr_low_reaction_butter = butterworth(spr_low_reaction_savgol)


## Reference
spr_ref    = plot_raw[vcsels[2]][1, :]
spr_ref    = spr_ref[0:-remove_air]
spr_ref_reaction = spr_ref[reaction_start:reaction_stopp]
spr_ref_reaction = spr_ref_reaction - spr_ref_reaction[0]
spr_ref_reaction_savgol = savgol(spr_ref_reaction, savgol_window, savgol_order)
spr_ref_reaction_butter = butterworth(spr_ref_reaction)


## Figure object
fig = plt.figure(1, figsize=(8,5))
ax = fig.add_subplot(111)


## Plot raw data
marker_size = 0.1
marker_type = '-x'
marker_color  = 'black'
linewidth_raw = 0.1
plt.plot(time_trace_reaction_zeroed, spr_high_reaction, marker_type, color=marker_color, ms=marker_size, linewidth=linewidth_raw)
plt.plot(time_trace_reaction_zeroed, spr_low_reaction, marker_type, color=marker_color, ms=marker_size, linewidth=linewidth_raw)
plt.plot(time_trace_reaction_zeroed, spr_ref_reaction, marker_type, color=marker_color, ms=marker_size, linewidth=linewidth_raw)


## Plot filtered data
linewidth_filtered = 0.8
plt.plot(time_trace_reaction_zeroed, spr_high_reaction_butter, color='blue', label=r'1.0 nM miR-122', linewidth=linewidth_filtered)
plt.plot(time_trace_reaction_zeroed, spr_low_reaction_butter, color='red', label=r'0.5 nM miR-122', linewidth=linewidth_filtered)
plt.plot(time_trace_reaction_zeroed, spr_ref_reaction_butter, color='green', label=r'DNA reference', linewidth=linewidth_filtered)


plt.grid(True)
plt.legend(fontsize=12, framealpha=0.3, edgecolor='black', loc='upper left')
plt.xlabel(r'Time [min]')
plt.ylabel(r'SPR shift [$\mu$m]')

ax.grid(linewidth=1, alpha=0.3)
x_tick_dx = 3
ax.set_xticks(np.arange(time_trace_reaction_zeroed[0], time_trace_reaction_zeroed[-1] + x_tick_dx, x_tick_dx))
ax.set_xticklabels(np.array(ax.get_xticks()).astype(int), rotation=0)

image_name = 'miRNA.svg'
image_format = 'svg'
plt.savefig(image_name, format=image_format, dpi=600)


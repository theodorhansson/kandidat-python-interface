#%%
import os, sys
if os.path.dirname(os.path.dirname(os.path.realpath(__file__))) not in sys.path:
    sys.path.append(os.path.dirname(os.path.dirname(os.path.realpath(__file__))))
    
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

from scipy.optimize import curve_fit

from spr_project.plot_measurement_data.spr_plot_functions import load_measurement_data, plot_single_trace, plot_all_traces, butterworth
from spr_project.spr_calc.spr_calculations import ResonantN, SPR_ang, ResonantAngle, SPR_loc, ref_idx

from scipy.signal import butter, filtfilt

def butterworth(data):
    sampling_freq=1/5
    cutoff_freq=1/300
    order=1
    nyquist_freq = 0.5 * sampling_freq
    normalized_cutoff_freq = cutoff_freq/nyquist_freq
    b, a = butter(order, normalized_cutoff_freq, btype='lowpass')
    data = filtfilt(b, a, data)
    return data

## Load all data
spr_data_folder = 'spr_measurements_240222'
parent_path = Path(__file__).resolve().parents[1]
measurement_path = Path(parent_path, 'spr_measurement_data', spr_data_folder)


names = os.listdir(measurement_path)
IPV_path = Path(measurement_path, names[3])

VCSEL_0 = Path(IPV_path, 'VCSEL_0')
VCSEL_1 = Path(IPV_path, 'VCSEL_1')

data = np.genfromtxt(str(Path(VCSEL_0, 'data.txt')), delimiter=',')
time_trace = data[:,0]

Ib_0 = np.loadtxt(str(Path(VCSEL_0, 'Ib.txt'))) - 1
Ib_1 = np.loadtxt(str(Path(VCSEL_1, 'Ib.txt')))

Vb_0 = np.loadtxt(str(Path(VCSEL_0, 'Vb.txt'))) - 0.05
Vb_1 = np.loadtxt(str(Path(VCSEL_1, 'Vb.txt')))



# plt.plot(time_trace, Ib_0[:-1])
# # plt.plot(time_trace, Vb_0[:-1])
# plt.grid(linewidth=1, alpha=0.3)

# plt.xlabel(r'Time [s]')
# plt.ylabel(r'$I_b$ [mA]')
# plt.ylim([0.999995, 1.001])

# plt.ylabel(r'$V_b$ [V]')
# plt.ylim([2.41, 2.45])


plt.tight_layout()
# plt.ylim([0, 3])


time_trace_dt = time_trace[2] - time_trace[1]
time_trace_long = np.arange(0, 6*59*60, time_trace_dt)/2

repeats = len(time_trace_long)//len(Ib_0)
Ib_long = np.repeat(Ib_0, repeats)
Vb_long = np.repeat(Vb_0, repeats)
np.random.shuffle(Ib_long)
np.random.shuffle(Vb_long)
Ib_long = butterworth(Ib_long)
Vb_long = butterworth(Vb_long)

Ib_breakdown = np.ones(shape=50)*0
Vb_breakdown = np.ones(shape=50)*5.00545


Ib_long = np.concatenate((Ib_long, Ib_breakdown))
Vb_long = np.concatenate((Vb_long, Vb_breakdown))

remove = len(time_trace_long) - len(Ib_long)

plt.plot(time_trace_long[:-remove], Ib_long)
plt.xlabel(r'Time [min]')
plt.ylabel(r'$I_b$ [mA]')
# plt.ylim([0.999995, 1.001])

# plt.plot(time_trace_long[:-remove], Vb_long)
# plt.ylabel(r'$V_b$ [V]')
# plt.xlabel(r'Time [min]')

plt.grid(linewidth=1, alpha=0.3)








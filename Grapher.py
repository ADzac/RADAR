import matplotlib.pyplot as plt
from scipy.interpolate import make_interp_spline, BSpline
import numpy as np

freqs = np.array([2.1218, 2.1538, 2.1619, 2.165, 2.171, 2.17256, 2.175,  2.17796, 2.18, 2.183,  2.184,2.1841, 2.1853, 2.185, 2.18610])
Volts = np.array([0, 3, 5 , 6, 9, 10, 12, 15, 18, 20, 21, 24, 25, 27, 30])
V2 = np.array([0, 3, 6, 9, 12, 15, 18, 21, 24, 27, 30])
Vd = np.array([.01, 1.18, 2.3, 3.5, 4.7, 5.9, 7.1, 8.3, 9.5, 10, 11])
f2 = np.array([2.1218, 2.1538,  2.165, 2.171, 2.175,  2.17796, 2.18, 2.184,2.1841, 2.185, 2.18610])



Vnew = np.linspace(Volts.min(), Volts.max(), 200) 
spl = make_interp_spline(Volts, freqs, k=3)
y_smooth = spl(Vnew)

Vdnew = np.linspace(Vd.min(), Vd.max(), 200) 
spl2 = make_interp_spline(Vd, f2, k=3)
y_smooth2 = spl2(Vdnew)

plt.figure()
plt.plot(Vnew, y_smooth)
plt.xlabel("Volts(V)")
plt.ylabel("Frequences (MHz)")
plt.title("Courbe les frequences en fonction du VE")
plt.show()

# Plotting the second graph and highlighting VD = 1.99
highlight_VD = 1.99

# Interpolating to find the corresponding frequency
highlight_freq = spl2(highlight_VD)

plt.figure()
plt.plot(Vdnew, y_smooth2, label='Frequency vs VD')
plt.scatter(highlight_VD, highlight_freq, color='red', zorder=5, label=f'VD = {highlight_VD}')
plt.text(highlight_VD, highlight_freq, f'({highlight_VD},{highlight_freq:.4f})', fontsize=12, ha='right')
plt.xlabel("VD(V)")
plt.ylabel("Frequence (MHz)")
plt.title("Courbe frequence en fonction du VD")
plt.show()
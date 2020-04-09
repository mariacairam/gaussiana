import matplotlib.ticker as ticker
import numpy as np
from matplotlib import gridspec
from matplotlib.ticker import AutoMinorLocator
import matplotlib.pyplot as plt
import scipy as scipy
from scipy.optimize import curve_fit


y_array = open('y_array.txt', 'w')
imq = np.empty(30, dtype=object)
istart = np.empty(30, dtype=object)
istop = np.empty(30, dtype=object)
icounts = np.empty(5000, dtype=object)
# aqui estou apenas abrindo arrays (vetores) vazios e definindo os tamanhos

IT = "Intervalo_TOF.txt"
imq, istart, istop = np.loadtxt(IT, delimiter=None, unpack=True)

TT = "Table_TOF.txt"
icounts = np.loadtxt(TT)


istart = istart.astype(int)
istop = istop.astype(int)
i = istart[0]
while (i <= istop[0]):
    y_array.write('' + str(icounts[i]) + '\n')
    i += 1
y_array.close()
xp = np.linspace(istart[0], istop[0], (istop[0] - istart[0] + 1))
yp = np.loadtxt("y_array.txt")
print(yp)


fig = plt.figure(figsize=(4,3))
gs = gridspec.GridSpec(1,1)

ax1 = fig.add_subplot(gs[0])
ax1.plot(xp, yp, "ro")
fig.savefig("gau.png", format="png", dpi=1000)
cen = (istop[0]+istart[0])/2
print (cen)
center = int(cen)
amp= icounts[center]
sigma= 2

def gaussiana1(x, amp, cen, sigma):
    return amp * (1 / (sigma * (np.sqrt(2 * np.pi)))) * (np.exp((-1.0 / 2.0) * (((xp - cen) / sigma) ** 2)))
popt, pcov = scipy.optimize.curve_fit(gaussiana1, xp, yp, p0=[amp, cen, sigma])
perr = np.sqrt(np.diag(pcov))
print("amplitude = %0.2f (+/-) %0.2f" % (popt[0], perr[0]))
print("center = %0.2f (+/-) %0.2f" % (popt[1], perr[1]))
print("sigma = %0.2f (+/-) %0.2f\n\n" % (popt[2], perr[2]))


fig = plt.figure(figsize=(4,3))
gs = gridspec.GridSpec(1,1)
ax1 = fig.add_subplot(gs[0])
ax1.plot(xp, yp, "ro")
ax1.plot(xp, gaussiana1(xp, *popt), 'k--')
fig.savefig("fitgau.png", format="png", dpi=1000)
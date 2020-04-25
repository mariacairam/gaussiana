import matplotlib.ticker as ticker
import numpy as np
from matplotlib import gridspec
from matplotlib.ticker import AutoMinorLocator
import matplotlib.pyplot as plt
import scipy as scipy
import scipy.signal
from scipy.optimize import curve_fit

dados = open('dadosfit.log', 'w')
y_array = open('y_array.txt', 'w')

MQ = "mq1.txt"
mq = np.loadtxt(MQ)
print(mq)
i = 0
with open(MQ, 'r') as f:
    for line in f:
        i += 1
        nLinhas = i
print(nLinhas)

TT = "Table_TOF.txt"
icounts = np.loadtxt(TT)
print(icounts[12])

i = 12
while i < nLinhas + 12:
    y_array.write('' + str(icounts[i]) + '\n')
    i += 1
y_array.close()
yp = np.loadtxt("y_array.txt")
dados.write('xp e yp pico \n')
dados.write('' + str(mq) + '\n' + str(yp) + '\n')


def gaussiana1(x, amp, cen, sigma):
    return amp * (1 / (sigma * (np.sqrt(2 * np.pi)))) * (np.exp((-1.0 / 2.0) * (((mq - cen) / sigma) ** 2)))


picos, _ = scipy.signal.find_peaks(yp, height=0, threshold=None, distance=200)
print(picos)
print(yp[picos])
fig = plt.figure(figsize=(4, 3))
gs = gridspec.GridSpec(1, 1)
ax1 = fig.add_subplot(gs[0])
ax1.plot(mq, yp, "ro")
ax1.plot(mq[picos], yp[picos], "x")
fig.savefig("gau1.png", format="png", dpi=1000)

cen = mq[picos[0]]
amp = yp[picos[0]]
sigma = 2
print(cen)
print(amp)
popt, pcov = scipy.optimize.curve_fit(gaussiana1, mq, yp, p0=[amp, cen, sigma])
perr = np.sqrt(np.diag(pcov))
print("amplitude = %0.2f (+/-) %0.2f" % (popt[0] / (popt[2] * (np.sqrt(2 * np.pi))), perr[0]))
print("center = %0.2f (+/-) %0.8f" % (popt[1], perr[1]))
print("sigma = %0.2f (+/-) %0.8f\n\n" % (popt[2], perr[2]))
dados.write('\nAjustes pico')
dados.write('\nAmplitude pico : ' + str(popt[0]) + ' (+/-) ' + str(perr[0]) +
            '\nCentro pico :   ' + str(popt[1]) + ' (+/-) ' + str(perr[1]) +
            '\nSigma pico :    ' + str(popt[2]) + ' (+/-) ' + str(perr[2]) + '\n\n')

xpfit = np.linspace(mq[0], mq[nLinhas-1], 500)

def gaussiana1fit(x, amp, cen, sigma):
    return amp * (1 / (sigma * (np.sqrt(2 * np.pi)))) * (np.exp((-1.0 / 2.0) * (((xpfit - cen) / sigma) ** 2)))

fig = plt.figure(figsize=(4, 3))
gs = gridspec.GridSpec(1, 1)
ax1 = fig.add_subplot(gs[0])
ax1.plot(mq, yp, "ro")
ax1.plot(xpfit, gaussiana1fit(xpfit, *popt), 'k--')
ax1.set_xlabel("x: m/q", family="serif", fontsize=12)
ax1.set_ylabel("y: icounts", family="serif", fontsize=12)
fig.savefig("fitgau1.png", format="png", dpi=1000)

dados.close()

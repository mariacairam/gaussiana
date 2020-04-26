import matplotlib.pyplot as plt
import numpy as np
import scipy.signal
from matplotlib import gridspec
from scipy.optimize import curve_fit

dados = open('dadosfit1.log', 'w')
y_array = open('y_array.txt', 'w')

MQ = "mq2.txt"
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

def gaussiana2(x, amp1, cen1, sigma1, amp2, cen2, sigma2):
    return amp1 * (1 / (sigma1 * (np.sqrt(2 * np.pi)))) * (np.exp((-1.0 / 2.0) * (((mq - cen1) / sigma1) ** 2))) + \
           amp2 * (1 / (sigma2 * (np.sqrt(2 * np.pi)))) * (np.exp((-1.0 / 2.0) * (((mq - cen2) / sigma2) ** 2)))

picos, _ = scipy.signal.find_peaks(yp, height=100, threshold=None)
print(picos)
print(yp[picos])
fig = plt.figure(figsize=(4, 3))
gs = gridspec.GridSpec(1, 1)
ax1 = fig.add_subplot(gs[0])
ax1.plot(mq, yp, "ro")
ax1.plot(mq[picos], yp[picos], "x")
fig.savefig("gau2.png", format="png", dpi=1000)



sigma1 = 0.02
sigma2 = 0.04
cen1 = mq[picos[0]]
amp1 = yp[picos[0]]
cen2 = mq[picos[1]]
amp2 = yp[picos[1]]

print(cen1)
print(amp1)
print(cen2)
print(amp2)
popt2, pcov2 = scipy.optimize.curve_fit(gaussiana2, mq, yp, p0=[amp1, cen1, sigma1, amp2, cen2, sigma2])
perr = np.sqrt(np.diag(pcov2))
print(popt2)
print("amplitude1 = %0.2f (+/-) %0.2f" % (popt2[0] / (popt2[2] * (np.sqrt(2 * np.pi))), perr[0]))
print("center1 = %0.2f (+/-) %0.8f" % (popt2[1], perr[1]))
print("sigma1 = %0.2f (+/-) %0.8f\n\n" % (popt2[2], perr[2]))
print("amplitude2 = %0.2f (+/-) %0.2f" % (popt2[3] / (popt2[5] * (np.sqrt(2 * np.pi))), perr[3]))
print("center2 = %0.2f (+/-) %0.8f" % (popt2[4], perr[4]))
print("sigma2 = %0.2f (+/-) %0.8f\n\n" % (popt2[5], perr[5]))
dados.write('\nAjustes pico')
dados.write('\nAmplitude pico 1: ' + str(popt2[0]) + ' (+/-) ' + str(perr[0]) +
            '\nCentro pico 1:   ' + str(popt2[1]) + ' (+/-) ' + str(perr[1]) +
            '\nSigma pico 1:    ' + str(popt2[2]) + ' (+/-) ' + str(perr[2]) + '\n\n'
            '\nAmplitude pico 2: ' + str(popt2[3]) + ' (+/-) ' + str(perr[3]) +
            '\nCentro pico 2:   ' + str(popt2[4]) + ' (+/-) ' + str(perr[4]) +
            '\nSigma pico 2:    ' + str(popt2[5]) + ' (+/-) ' + str(perr[5]) + '\n\n'
            )

fig = plt.figure(figsize=(4, 3))
gs = gridspec.GridSpec(1, 1)
ax1 = fig.add_subplot(gs[0])
ax1.plot(mq, yp, "ro")
ax1.plot(mq, gaussiana2(mq, *popt2), 'k--')
ax1.set_xlabel("x: m/q", family="serif", fontsize=12)
ax1.set_ylabel("y: icounts", family="serif", fontsize=12)
fig.savefig("fitgau2.png", format="png", dpi=1000)

dados.close()

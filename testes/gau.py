import matplotlib.pyplot as plt
import numpy as np
import scipy.signal
from scipy.optimize import curve_fit
from matplotlib import gridspec

dados = open('dadosfit.log', 'w')
y_array = open('y_array.txt', 'w')

MQ = "MassWave.txt"
mq = np.loadtxt(MQ)
# print(mq)
i = 0
with open(MQ, 'r') as f:
    for line in f:
        i += 1
        nLinhas = i
print(nLinhas)

TT = "Table_TOF.txt"
icounts = np.loadtxt(TT)
print(icounts[12])

PA = "parametros.txt"
par = open(PA, "r")
parametros = [(line.strip()).split() for line in par]
par.close()

for i in range(12, nLinhas + 12, 1):
    y_array.write('' + str(icounts[i]) + '\n')
y_array.close()
yp = np.loadtxt("y_array.txt")


def linear(x, a, b):
    return a * x + b


mqlin = np.zeros(82)
yplin = np.zeros(82)
for i in range(0, 41, 1):
    mqlin[i] = mq[109 + i]
    mqlin[81 - i] = mq[nLinhas - 109 - i]
    yplin[i] = yp[109 + i]
    yplin[81 - i] = yp[nLinhas - 109 - i]
print(yplin)

popt_linear, pcov_linear = scipy.optimize.curve_fit(linear, mqlin, yplin)

fig = plt.figure(figsize=(4, 3))
plt.plot(mq, yp)
plt.plot(mq, linear(mq, *popt_linear), 'k--')
plt.show()

base = linear(mq, *popt_linear)
yp = yp - base
f = scipy.interpolate.interp1d(mq, yp)
mq = np.arange(mq[0], mq[nLinhas - 1] - 0.1, 0.01)
yp = f(mq)
base = linear(mq, *popt_linear)

picos, _ = scipy.signal.find_peaks(yp, height=100 - base, threshold=None, distance=50, width=1.5)

picos = np.delete(picos, [10])
print(picos)
print(mq[picos])
print(yp[picos])

plt.figure(figsize=(4, 3))
plt.title('Picos Espectro')
plt.plot(mq, yp)
plt.scatter(mq[picos], yp[picos], color='red', marker='+', label='Picos encontrados')
plt.legend()
plt.savefig("gau.png", format="png", dpi=1000)
plt.show()


def gaussiana1(x, amp, cen, sigma):
    return amp * (1 / (sigma * (np.sqrt(2 * np.pi)))) * (np.exp((-1.0 / 2.0) * (((mq - cen) / sigma) ** 2)))


def gaussiana2(x, *parametros):
    y = np.zeros_like(x)
    for c in range(0, (3 * len(picos)), 3):
        amp = parametros[c]
        cen = parametros[c + 1]
        sigma = parametros[c + 2]
        y = y + amp * (1 / (sigma * (np.sqrt(2 * np.pi)))) * (np.exp((-1.0 / 2.0) * (((mq - cen) / sigma) ** 2)))
    return y


fwhm_rel, h_eval, left_ips, right_ips = scipy.signal.peak_widths(yp, picos, rel_height=0.5)
esq = left_ips.astype(int)
dire = right_ips.astype(int)
sigma = (mq[dire] - mq[esq]) / 2 * np.sqrt(2 * np.log(2))
cen = mq[picos]
amp = yp[picos]

chute = np.zeros(3 * len(picos))
boundinf = np.zeros(3 * len(picos))
boundsup = np.zeros(3 * len(picos))
c = 0
i = 0
while c < (3 * len(picos)) and i < len(picos):
    chute[c] = amp[i]
    chute[c + 1] = cen[i]
    chute[c + 2] = sigma[i]
    boundinf[c] = -np.inf
    boundinf[c + 1] = cen[i] - 0.2
    boundinf[c + 2] = -np.inf
    boundsup[c] = np.inf
    boundsup[c + 1] = cen[i] + 0.2
    boundsup[c + 2] = np.inf
    c += 3
    i += 1

popt2, pcov2 = scipy.optimize.curve_fit(gaussiana2, mq, yp, p0=chute,
                                        bounds=(boundinf, boundsup))
perr = np.sqrt(np.diag(pcov2))

dados.write('Ajustes picos\n')
dados.write('Numero de picos: ' + str(len(picos)) + '\n')

i = 0
c = 0
area = np.zeros(len(picos))
areatot = 0
while c < (3 * len(picos)) and i < len(picos):
    pico_gauss = gaussiana1(mq, popt2[c], popt2[c + 1], popt2[c + 2])
    area[i] = np.trapz(pico_gauss, x=mq, dx=0.001, axis=0)
    areatot = areatot + area[i]
    dados.write('\nPico de m/q experimental = ' + str(mq[picos[i]]) +
                '\nArea:      ' + str(area[i]) +
                '\nAmplitude: ' + str((popt2[c] / (popt2[c + 2] * (np.sqrt(2 * np.pi))))) + ' (+/-) ' + str(perr[c]) +
                '\nCentro:    ' + str(popt2[c + 1]) + ' (+/-) ' + str(perr[c + 1]) +
                '\nSigma:     ' + str(popt2[c + 2]) + ' (+/-) ' + str(perr[c + 2]) + '\n\n'
                )
    i += 1
    c += 3
dados.write('Area total: ' + str(areatot) + '\n\n')

for i in range(0, len(picos), 1):
    dados.write('Pico de m/q experimental = ' + str(mq[picos[i]]) +
                '\nArea relativa: ' + str((area[i] / areatot * 100)) + '\n\n')

fig = plt.figure(figsize=(4, 3))
gs = gridspec.GridSpec(1, 1)
ax1 = fig.add_subplot(gs[0])
plt.plot(mq, yp)
plt.plot(mq, gaussiana2(mq, *popt2), 'k--')
pico_gauss = gaussiana1(mq, popt2[39], popt2[40], popt2[41])
ax1.plot(mq, pico_gauss, "b")
ax1.fill_between(mq, pico_gauss.min(), pico_gauss, facecolor="blue", alpha=0.5)
plt.xlabel("x: m/q", family="serif", fontsize=12)
plt.ylabel("y: icounts", family="serif", fontsize=12)
plt.tick_params(axis='both', which='major', direction="out", top="on", right="on", bottom="on", length=8, labelsize=8)
plt.tick_params(axis='both', which='minor', direction="out", top="on", right="on", bottom="on", length=5, labelsize=8)
fig.tight_layout()
fig.savefig("fitgau.png", format="png", dpi=1000)
plt.show()

dados.close()

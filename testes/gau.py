import matplotlib.pyplot as plt
import numpy as np
import scipy.signal
from scipy.optimize import curve_fit

dados = open('dadosfit.log', 'w')
y_array = open('y_array.txt', 'w')

MQ = "MassWave.txt"
mq = np.loadtxt(MQ)
#print(mq)
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
#dados.write('xp e yp pico \n')
#dados.write('' + str(mq) + '\n' + str(yp) + '\n')

def linear(x, a, b):
    return a * x + b


aa = 0
bb = 0
ix1 = 109
ix2 = nLinhas - 109
for i in range(1, 42, 1):
    aa = aa + yp[109 + i]
    bb = bb + yp[nLinhas - 109 - i]

aa = aa / 41
bb = bb / 41

bcoef = (bb - aa) / (ix2 - ix1)
acoef = aa - bcoef * ix1

popt_linear, pcov_linear = scipy.optimize.curve_fit(linear, mq, yp, p0=[acoef, bcoef])

fig = plt.figure(figsize=(4, 3))
plt.plot(mq, yp)
plt.plot(mq, linear(mq, *popt_linear), 'k--')
plt.show()


picos, _ = scipy.signal.find_peaks(yp, height=100, threshold=None, distance=10, width=1.5)


#print(picos)
#print(mq[picos])
#print(yp[picos])
picos = np.delete(picos, [10])
print(picos)
print(mq[picos])
print(yp[picos])

fig = plt.figure(figsize=(4, 3))
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

chute = np.column_stack((amp, cen, sigma))
print(chute)

popt2, pcov2 = scipy.optimize.curve_fit(gaussiana2, mq, yp, p0=chute)
perr = np.sqrt(np.diag(pcov2))

dados.write('Ajustes picos\n')
dados.write('Numero de picos: ' + str(len(picos)) + '\n')

i = 0
c = 0
area = np.zeros(len(picos))
areatot = 0
while c < (3 * len(picos)) and i < len(picos):
    pico_gauss = gaussiana1(mq, popt2[c], popt2[c + 1], popt2[c + 2])
    area[i] = np.trapz(pico_gauss, x=mq, dx=0.001)
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
plt.plot(mq, yp)
plt.plot(mq, gaussiana2(mq, *popt2), 'k--')
plt.xlabel("x: m/q", family="serif", fontsize=12)
plt.ylabel("y: icounts", family="serif", fontsize=12)
plt.tick_params(axis='both', which='major', direction="out", top="on", right="on", bottom="on", length=8, labelsize=8)
plt.tick_params(axis='both', which='minor', direction="out", top="on", right="on", bottom="on", length=5, labelsize=8)
fig.tight_layout()
fig.savefig("fitgau.png", format="png", dpi=1000)

dados.close()

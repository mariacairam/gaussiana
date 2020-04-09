import matplotlib.ticker as ticker
import numpy as np
from matplotlib import gridspec
from matplotlib.ticker import AutoMinorLocator
import matplotlib.pyplot as plt
import scipy as scipy
from scipy.optimize import curve_fit

dados = open('dadosajuste.log', 'w')

y_array = open('y_array.txt', 'w')
imq = np.empty(30, dtype=object)
istart = np.empty(30, dtype=object)
istop = np.empty(30, dtype=object)
icounts = np.empty(5000, dtype=object)
# aqui estou apenas abrindo arrays (vetores) vazios e definindo os tamanhos

IT = "Intervalo_TOF.txt"
imq, istart, istop = np.loadtxt(IT, delimiter=None, unpack=True)
# np.loadtxt é uma função que lê o txt e escreve as informações no array, onde cada linha do array é uma linha do txt,
# definindo delimiter=None ele vai separar os elementos por espaço e unpack=True vai colocar cada elemento da linha
# em um array diferente.

TT = "Table_TOF.txt"
icounts = np.loadtxt(TT)
# a mesma coisa mas só tem 1 elemento por linha então as opções que a função oferece está no default

i=0
with open(IT, 'r') as f:
    for line in f:
        i += 1
        nFrag = i
print (nFrag)


c= 0
while(c<nFrag):
    y_array = open('y_array.txt', 'w')
    y_array.truncate(0)
    istart = istart.astype(int)
    istop = istop.astype(int)
    # estou transformando os elementos dos arrays e inteiros para poder utilizar no icounts[i]
    i = istart[c]
    while (i <= istop[c]):
        y_array.write('' + str(icounts[i]) + '\n')
        i += 1
    y_array.close()
    xp = np.linspace(istart[c], istop[c], (istop[c] - istart[c] + 1))
    yp = np.loadtxt("y_array.txt")
    open('y_array.txt', 'w')
    dados.write('xp e yp pico '+str(c+1)+'\n')
    dados.write(''+str(xp)+'\n'+str(yp)+'\n')

    fig = plt.figure(figsize=(4,3))
    gs = gridspec.GridSpec(1,1)
    ax1 = fig.add_subplot(gs[0])
    ax1.plot(xp, yp, "ro")
    fig.savefig("gau.png", format="png", dpi=1000)

    cen = (istop[c]+istart[c])/2
    print (cen)
    center = int(cen)
    amp= icounts[center]
    sigma= 2
    # estou supondo que o meu centro é a média entre istart e istop, que minha amplitude é o meu icounts no meu centro
    # suposto e que o sigma é dois.
    #dados.write(''+str(amp)+'      '+str(cen)+'     '+str(sigma)+'\n')
    def gaussiana1(x, amp, cen, sigma):
        return amp*(1/(sigma*(np.sqrt(2*np.pi))))*(np.exp((-1.0/2.0)*(((xp-cen)/sigma)**2)))

    popt, pcov = scipy.optimize.curve_fit(gaussiana1, xp, yp, p0=[amp, cen, sigma])
    perr = np.sqrt(np.diag(pcov))
    print("amplitude = %0.2f (+/-) %0.2f" % (popt[0], perr[0]))
    print("center = %0.2f (+/-) %0.2f" % (popt[1], perr[1]))
    print("sigma = %0.2f (+/-) %0.2f\n\n" % (popt[2], perr[2]))
    dados.write('\nAjustes pico'+str(c+1))
    dados.write('\nAplitude pico '+str(c+1) +': ' +str(popt[0]) +' (+/-) ' +str(perr[0]) +
                '\nCentro pico ' +str(c+1) +':   ' + str(popt[1]) +' (+/-) ' + str(perr[1])+
                '\nSigma pico ' +str(c+1) +':    ' +str(popt[2]) +' (+/-) ' + str(perr[2]) + '\n\n')

    xpfit= np.linspace(istart[c], istop[c], 50)
    def gaussiana1fit(x, amp, cen, sigma):
        return amp*(1/(sigma*(np.sqrt(2*np.pi))))*(np.exp((-1.0/2.0)*(((xpfit-cen)/sigma)**2)))

    fig = plt.figure(figsize=(4, 3))
    gs = gridspec.GridSpec(1, 1)
    ax1 = fig.add_subplot(gs[0])
    ax1.plot(xp, yp, "ro")
    ax1.plot(xpfit, gaussiana1fit(xpfit, *popt), 'k--')
    ax1.set_xlabel("x: tvoo", family="serif", fontsize=12)
    ax1.set_ylabel("y: icounts", family="serif", fontsize=12)
    fig.savefig("fitgaunovo.png", format="png", dpi=1000)

    c +=1

dados.close()
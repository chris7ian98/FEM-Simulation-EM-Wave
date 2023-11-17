import numpy as np
from settings_domain1D import *
from matplotlib import pyplot as plt
import shapeFunctions
import main_FEM
import plotResult
import integration

# Zeitbereich [s]
t_0 = 1e-10
c_0 = 299792458
# Kreisfrequenz [Hz]
omega = 4 * np.pi * c_0 / b_  #genau zwei Wellenlaengen im Problemgebiet

#Damit Anfangswert und Anfangsgeschwindigkeit Null sind
#bessere Ergebnisse für Zeitschrittverfahren
t_phi= 2/omega
a= 1/np.sqrt(2)/t_phi**2
phi= np.pi/4-omega*t_phi

# Anregung
def w(t):
    if(t<t_phi):
        return a*t**2
    return np.sin(omega *t+phi)

# Ableitung der Anregung
def w_der(t):
    if (t < t_phi):
        return 2*a * t
    return np.cos(omega * t+phi) * omega
    
'''
Funktion eps(x) und sigma(x) darstellen
'''
if (1):
    x = [b_ / 1000 * i for i in range(0, 1001)]
    y_eps = [eps(xi) / eps0 for xi in x]
    y_sigma = [sigma(xi) for xi in x]
    
    fig = plt.figure()
    ax1 = fig.add_subplot(111)
    ax1.plot(x, y_eps, 'r')
    ax1.set_ylabel('$\epsilon_r(x)$', fontsize= 15)
    ax1.set_xlabel('$x [m]$', fontsize=15)

    ax2 = ax1.twinx()
    ax2.plot(x, y_sigma, 'b--')
    ax2.set_ylabel('$\sigma(x) \,\, [(\Omega m)^{-1}]$', color='b', fontsize= 15)
    for tl in ax2.get_yticklabels():
        tl.set_color('b')
    fig.suptitle('Materialparameter $\epsilon_r(x), \sigma(x)$', fontsize=15)
    ax1.tick_params(axis="x", labelsize=15)
    ax1.tick_params(axis="y", labelsize=15)
    ax2.tick_params(axis="y", labelsize=15) 
    plt.show()

'''
Anregung w(t) darstellen
'''

if (1):
    x = [t_0 / 1000 * i for i in range(0, 1001)]
    y = [w(xi) for xi in x]
    plt.plot(x, y)
    plt.xlabel('$t [s]$', fontsize= 15)
    plt.ylabel('$w(t) [V/m]$', fontsize= 15)
    plt.xticks(fontsize= 15)
    plt.yticks(fontsize=15)
    plt.title('Anregung w(t) (linker Rand)', fontsize=15)
    plt.grid()
    plt.show()



#Basisfunktionen(Shape Functions)
SF =  1     #0->Lagrange, 1->Legendre
m =   1000 #Anzahl der Zeitschritte
nFE=  100   #Anzahl der finiten Elemente
nSF = 2    #Ordnung der Basisfunktionen
integrationOrderM = nSF + 1
integrationOrderC = nSF + 1
n = nFE * nSF + 1
h_t = t_0 / m
h_x = (b_ - a_) / nFE


'''
Animation
'''
if (1):
    U = np.zeros(shape=(n, m))
    U = main_FEM.getU_main(a_, b_, nFE, n, m, nSF, h_t, h_x, integrationOrderM, integrationOrderC, w, w_der, SF)
    plotResult.animation2D(U, m, nSF, h_x, nFE, SF)

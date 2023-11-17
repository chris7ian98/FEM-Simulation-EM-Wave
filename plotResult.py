from matplotlib import animation
from matplotlib import pyplot as plt
from settings_domain1D import *
import numpy as np
import shapeFunctions

resolution_x= 200
resolution_m= 1000

'''
Da die Spaltenvektoren der Naeherungsloesung U die Koeffizienten entalten
muss um den Loesungswert bei x zu erhalten die relevanten Basisfunktionen
bei x ausgewertet und mit den jewiligen Koeffizienten gewichtet werden
'''
def getSolutionValueFromU(U, timestep, x, nSF, h_x, nFE, SF):
    coeffs= np.zeros(nSF+1)
    k= int(x/h_x)
    if(k>=nFE):
        k= k-1

    for i in range(0, nSF+1):
        coeffs[i]= U[k*nSF+i][timestep]
    X= (x-k*h_x)*2/h_x -1

    sum= 0
    if(SF==0):
        for i in range(0, nSF+1):
            sum+= coeffs[i]*shapeFunctions.Lagrange(nSF, i, X, h_x)
    else:
        for i in range(0, nSF+1):
            sum+= coeffs[i]*shapeFunctions.Legendre(nSF, i, X, h_x)
    return sum



'''
Animation um die Loesung U darzustellen
Spaltenvektoren von U enthalten Koeffizienten der Basisfunktionen
Python Animations Bsp.:
(vgl. https://jakevdp.github.io/blog/2012/08/18/matplotlib-animation-tutorial/)
'''
def animation2D(U, m, nSF, h_x, nFE, SF):
    if m>=resolution_m:
        m_new= resolution_m
    else:
        m_new= m
    U_nodes = np.zeros(shape=(resolution_x+1, m_new))
    X = np.zeros(shape=(resolution_x+1, m_new))

    #Problemgebiet [a,b]
    for i in range(0, resolution_x+1):
        for j in range(0, m_new):
            X[i][j]= a_+ (b_-a_)/resolution_x*i

    for i in range(0, resolution_x+1):
        for j in range(0, m_new):
            U_nodes[i][j]= getSolutionValueFromU(U, int(j*m/m_new), X[i][j], nSF, h_x, nFE, SF)

    fig = plt.figure()
    ax  = plt.axes(xlim=(a_ - (b_ - a_) * 0.05, b_ + (b_ - a_) * 0.05), ylim=(-2, 2)) 
    data, = ax.plot([], [], lw=2)
    text_for_timestep = ax.text(0.045, 0.90, '', transform=ax.transAxes, fontsize= 15)

    def initialisation():
        data.set_data([], [])
        text_for_timestep.set_text('')
        return data, text_for_timestep

    def anima(j):
        data.set_data(X[:, j], U_nodes[:, j])
        text_for_timestep.set_text('$t/t_0=$'+ str((round(j/m_new * 10)) / 10))
        return data, text_for_timestep

    anm= animation.FuncAnimation(fig, anima, init_func=initialisation,
                                   frames=m_new, interval= 20, blit=True)

    ax.set_xlabel('$x [m]$', fontsize= 15)
    ax.set_ylabel('$u(x,t) [V/m]$', fontsize=15)
    ax.tick_params(axis="x", labelsize=14)
    ax.tick_params(axis="y", labelsize=14)

    plt.grid(True)
    plt.show()

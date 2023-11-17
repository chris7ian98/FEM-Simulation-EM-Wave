import numpy as np
import integration
from settings_domain1D import *

'''
Loesen eines tridiagonalen Gleichungssystems 
(nur bei linearen Basisfunktionen verwendbar)
'''
def TDSolver(A, x, d):
    n= len(x)
    array_ci= np.zeros(n)
    array_di= np.zeros(n)

    array_ci[0]= A[0][1]/A[0][0] #c1/b1
    array_di[0]= d[0]/A[0][0]   #c1/b1
    for i in range(1, len(x)-1):
        array_ci[i]= A[i][i+1]/(A[i][i]- A[i][i-1]*array_ci[i-1])
    for i in range(1, len(x)):
        array_di[i]= (d[i]-A[i][i-1]*array_di[i-1])/(A[i][i]-A[i][i-1]*array_ci[i-1])
    x[n-1]= array_di[n-1]
    for i in range(n-2, -1, -1):
        x[i]= array_di[i]-array_ci[i]*x[i+1]


'''
NEWMARK-BETA-VERFAHREN:
Loesungsmatrix U, Geschwindigkeitsmatrix V, Beschleunigungsmatrix A
Systemmatrizen: Steifigkeitsmatrix K, Massenmatrix M, Daempfungsmatrix C
Loesung fuer einen Zeitschritt  wird ineinem Spaltenvektor der Matrix U, V, A gespeichert.
'''
def newmarkBetaMethod(U, V, A, K, M, C, w, w_der, h_t, m, nSF):
    gamma= 1/2
    beta= 1/4

    b = np.zeros(np.shape(K)[0])

    U[0][0] = w(0)      #Anfangszustand
    V[0][0] = w_der(0)  #Anfangsgeschwindigkeit

    MatrixAi = M / beta / h_t ** 2 + gamma / beta / h_t * C + K

    for i in range(1, m):
        b[0]= 2*(np.sqrt(integration.mue_eps(a_)))*w_der(i*h_t)
        VektorBi= b+ M.dot(U[:,i-1]/beta/h_t**2+ V[:,i-1]/beta/h_t+(1/2/beta-1)*A[:,i-1]) +C.dot(gamma*U[:,i-1]/beta/h_t- V[:,i-1]*(1-gamma/beta)-h_t*A[:,i-1]*(1-gamma/2/beta))
        if (nSF > 1):  # Ordnung der Formfunktionen
            U[:, i] = np.linalg.solve(MatrixAi, VektorBi)
        else:
            TDSolver(MatrixAi, U[:, i], VektorBi)
        A[:,i]= (U[:,i]-U[:,i-1])/beta/h_t**2- V[:,i-1]/beta/h_t-(1/2/beta-1)*A[:,i-1]
        V[:,i]= gamma *(U[:, i]-U[:, i-1])/beta/h_t+ V[:, i-1]*(1-gamma/beta) + h_t*A[:, i-1]*(1-gamma/2/beta)
    return U


import time
import numpy as np
import timeIntegrationMethod
import integration


'''
FEM + NEWMARK-BETA-VERFAHREN
Rueckgabe: U
Spaltenvektoren von U enthalten die Koeffizienten der Basisfunktionen fuer einen Zeitschritt
'''
def getU_main(a_, b_, nFE, n, m, nSF, h_t, h_x, integrationOrderM, integrationOrderC, w, w_der, SF):
    '''LAUFZEIT'''
    startTime = time.time()
    '''
    DEFINITION DER MATRIZEN UND VEKTOREN
    nFE...Anzahl der Finiten Elemente
    nSF...Ordnung der Basisfunktionen (Shape functions)
    n= nFE*nSF+1
    (n)... Anzahl aller Basisfunktionen
    '''

    #System-Matrizen
    K = np.zeros(shape=(n, n))  #Steifigkeits-Matrix
    M = np.zeros(shape=(n, n))  #Massenmatrix
    C = np.zeros(shape=(n, n))  #Daempfungsmatrix

    #Loesungsvektoren
    U = np.zeros(shape=(n, m))  #elektrische Feldstaerke
    V = np.zeros(shape=(n, m))  #Geschwindigkeit
    A = np.zeros(shape=(n, m))  #Beschleunigung

    #Lokale Matrix
    A_local = np.zeros(shape=(nSF+1, nSF+1)) #lokale Matrix

    '''
    FINITE ELEMENTE METHODE:
    Eintraege der Matrizen K,M und C werden berechnet:
    Dabei werden die lokalen Matritzen A fuer jedes Element k der nFE 
    Finiten Elemente berechnet und anschliessend in 
    die globale Matrix (K, M, C) eingetragen. 
    Parameter:
    -k-tes Element(wird fuer Transformation benoetigt, eps(x), sigma(x))
    -{0,1} es wird ueber die normalen (abgeleiteten) Basisfunktionen integriert  
    -Gewichte der Gauss-Legendre Quadraturregel
    -Gewichte der Gauss-Legendre Quadraturregel
    -Funktion (1, mue_esp(x), sig_eps(x)) welche im Integraden vorkommt
    '''

    '''
    STEIFIGKEITSMATRIX K
    '''
    gaussQuadK= integration.GaussLegendreQuadratureRule(nSF+1)

    for k in range(0, nFE):  #Finite Elemente
        integration.calculateLocalMatrix(k, 1, gaussQuadK[0], gaussQuadK[1], A_local,  integration.Eins, h_x, nSF, SF)
        integration.local_global_table(k, A_local, K, nSF)

    '''
    MASSENMATRIX M
    '''
    gaussQuadM = integration.GaussLegendreQuadratureRule(integrationOrderM)

    for k in range(0, nFE):#Finite Elemente
        integration.calculateLocalMatrix(k, 0, gaussQuadM[0], gaussQuadM[1], A_local, integration.mue_eps, h_x, nSF, SF)
        integration.local_global_table(k, A_local, M, nSF)

    '''
    Daempfungsmatrix
    '''
    
    gaussQuadC = integration.GaussLegendreQuadratureRule(integrationOrderC)

    for k in range(0, nFE):#Finite Elemente
        integration.calculateLocalMatrix(k, 0, gaussQuadC[0], gaussQuadC[1], A_local, integration.mue_sigma, h_x, nSF, SF)
        integration.local_global_table(k, A_local, C, nSF)


    #Rechter Rand, linker Rand (nichtreflektierende Randbedingungen)
    C[0][0]= C[0][0]+ np.sqrt(integration.mue_eps(a_)) 
    C[n-1][n-1]= C[n-1][n-1]+ np.sqrt(integration.mue_eps(b_))
    
    print('Laufzeit (K,M,C) berechnet:', time.time()-startTime, 's')
    print('K', K, '\nM', M, '\nC', C)


    '''
    ZEITSCHRITTVERFAHREN:
    Newmark-Beta-Verfahren
    '''
    U= timeIntegrationMethod.newmarkBetaMethod(U, V, A, K, M, C, w, w_der, h_t, m, nSF)

    print('Laufzeit nach FEM und Zeitschrittverfahren:', time.time()-startTime, 's')

    return U

import numpy as np
import shapeFunctions
from settings_domain1D import *

'''
Legendre Polynome ueber 3-Term-REkursion berechnet
wird in Pderived benoetigt
Matrix T= [L0, L1, L2, ..., Ln] 
'''
def LegndrePolynomials(T):
    n= T.shape[0]

    T[0][0]= 1 #erster Spaltenvektor (1)
    T[1][1]= 1 #zweiter Spaltenvektor (x)
    for j in range(1, n-1):# Spalte
        for i in range(0, n): #Zeile
            T[i][j+1]= (2*j+1)/(j+1)*T[i-1][j]- (j)/(j+1)*T[i][j-1]
    return T


'''
n-te Legendre Polynom abgeleitet und an der Stelle x ausgewertet
um die Gewichte der Gauss Legendre Quadraturformel zu bestimmen
'''
def Pderived(n, x):
    S= 0
    T= np.zeros(shape=(n+1, n+1))#+1  (0. bis n.-Legendre Polynom)
    LegndrePolynomials(T)
    for i in range(1, n+1):# wegen Ableitung 0->1 (shift)
        S+= T[i][n]* (i)* x**(i-1)
    return S


'''
GAUSS-LEGENDRE-QUADRATURE_RULE
I=[-1, 1]
(1)Stuetzstellen ->Nullstellen des n-ten Legendre-Polynoms.
3-Term-Rekursion-> Eigenwertproblem->Stuetzstellen
(2): aus Nullstellen Gewichte berechnen ueber Formel (benoetigt Pderived)
Ruekgabe: 2 Listen -> [Liste der Stuetzstellem, Liste der Gewichte] 
'''
def GaussLegendreQuadratureRule(n): #order
    T = np.zeros(shape=(n, n))

    T[0][0] = 0
    T[0][1] = 1
    for i in range(1, n-1):
        #
        ai = ((i+1) -1)/(2*(i+1) - 1)
        bi = ((i+1))/ (2*(i+1) - 1)
        #
        T[i][i - 1] = ai
        T[i][i] = 0
        T[i][i+1] = bi
    #
    an = (n -1)/(2*n - 1)
    T[n - 1][n - 2] = an


    points = np.linalg.eigvals(T)
    points.sort()

    #--weights--
    weights = np.zeros(len(points))
    for i in range(0, len(points)):
        weights[i] = (2 / (1 - points[i] ** 2) / (Pderived(n, points[i])) ** 2)
    return [points, weights]


'''
Transformation von X auf x des k-ten Elements
[-1, 1] auf [xk, xk+1]
'''
def transform_X_to_x(k, X, h_x):
    return k*h_x +(h_x)*(X+1)/2


'''
Einsfunktion
wird fuer die Berechnung der Eintraege der Steifigkeitsmatrix benoetigt
'''
def Eins(x):
    return 1


'''
Funktion
wird fuer die Berechnung der Eintraege der Massematrix benoetigt
'''
def mue_eps(x):
    return mue0*eps(x)


'''
Funktion
wird fuer die Berechnung der Eintraege der Daempfungsmatrix benoetigt
'''
def mue_sigma(x):
    return mue0*sigma(x)

##########################################

'''
Berechnung der lokalen Elementmatrix des k-ten Elements
-Berechnung der k-ten Elementmatrix: 
Jede lokale Basisfunktion p_i mit jeder lokalen Basisfunktion p_j multipliziert 
und ueber das k-te Finite Element integriert.
'''
def calculateLocalMatrix(element_k, derived, xi, w, A, f, h_x, nSF, SF):
    for i in range(0, np.shape(A)[0]):
        for j in range(0, np.shape(A)[0]):
            A[i][j]= quadrature(element_k, derived, i, j, f, xi, w, nSF, h_x, SF) #Gauss-quad
    return A


'''
Zusammenhang zwischen Eintraegen der lokalen Element-Matrizen und Eintraegen der globalen Matrix. 
k... k-tes finite Element
'''
def local_global_table(k, A_local, A_global, nSF):
    for i in range(0, nSF+1):
        for j in range(0, nSF+1):
            A_global[nSF*k +i][nSF*k +j]+= A_local[i][j]
    return 1


'''
Gauss-Legendre Quadratur: 
Uebergabeparameter:
-element_k...k-te finte Element
-derived...Unterscheidung ob Basisfunktion od. abgeleitete Basisfuktion
-(i,j)... Eintrag in der lokalen Matrix der berechnet wird->(zugriff auf jeweilige Basisfunktionen)
-(xi, w)... Listen an Stuetzstellen und Gewichten
Ruekgabe: Naeherungslloesung des Integrals des Eintrags(i,j) der k-ten lokalen Elementmatrix A_k
'''
def quadrature(element_k, derived, i, j, f, xi, w, nSF, h_x, SF):
    Q= 0
    for k in range(0, len(w)):
        if(derived== 0):
            if(SF==0):
                Q+=  w[k] * f(transform_X_to_x(element_k, xi[k], h_x)) *h_x/2* shapeFunctions.Lagrange(nSF, i, xi[k], h_x) * shapeFunctions.Lagrange(nSF, j, xi[k], h_x)
            else:
                Q+=  w[k] * f(transform_X_to_x(element_k, xi[k], h_x)) *h_x/2* shapeFunctions.Legendre(nSF, i, xi[k], h_x) * shapeFunctions.Legendre(nSF, j, xi[k], h_x)
        else:
            if(SF==0):
                Q+= w[k] * f(transform_X_to_x(element_k, xi[k], h_x))  *h_x/2*shapeFunctions.Lagrange_der(nSF, i, xi[k], h_x) *shapeFunctions.Lagrange_der(nSF, j, xi[k], h_x)
            else:
                Q+= w[k] * f(transform_X_to_x(element_k, xi[k], h_x))  *h_x/2*shapeFunctions.Legendre_der(nSF, i, xi[k], h_x) *shapeFunctions.Legendre_der(nSF, j, xi[k], h_x)
    return Q
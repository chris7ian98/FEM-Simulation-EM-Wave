#Problemgebiet [m]
a_= 0
b_= 0.01

mue0= 1.25663706212*1e-6  #magn. Feldkonstante
eps0= 8.8541878128*1e-12  #elektr. Feldkonstante

#Permittivitaet
def eps(x):
    if x> b_*2/5: #and x<b_*3/5:
        return 3*eps0
    return 1*eps0

#Leitfaehigkeit
def sigma(x):
    if x> b_*4/8 and x<b_*5/8:
        return 0#1e16
    else:
        return 0

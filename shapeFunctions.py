'''
Lagrange Formfunktionen beliebiger Ordnung n 
[-1,1]
'''
def Lagrange(n, j, X, h_x):
    product = 1
    for k in range(0, n+1):
        if (k != j):
            product *= (X - (-1 + k * 2 / n)) / ((-1 + j * 2 / n) - (-1 + k * 2 / n))
    return product

'''
abgeleitete Lagrange Formfunktionen beliebiger Ordnung n 
[-1,1]
innere Ableitung (2/h_x), bei Integration benoetigt
'''
def Lagrange_der(n, j, X, h_x):
    sum = 0
    for k in range(0, n+1):
        product = 1
        if (k != j):
            for m in range(0, n + 1):
                if (m != j and m != k):
                    product *= (X - (-1 + m * 2 / n)) / ((-1 + j * 2 / n) - (-1 + m * 2 / n))

            sum += product / ((-1 + j * 2 / n) - (-1 + k * 2 / n))
    return sum*(2/h_x)


'''
Legendre Formfunktionen beliebiger Ordnung n 
[-1,1]
'''
def Legendre(order, n, X, h_x):
    l1 = (1 - X) / 2
    l2 = (X + 1) / 2
    if (n == 0):
        return l1

    if(n== order):
        return l2

    li_minus_1= -1
    li= X

    for i in range(2, n + 2):
        l = li
        li = (2 * i - 3) / (i) * X * li - (i - 3) / (i) * li_minus_1
        li_minus_1 = l
    return li


'''
abgeleitete Legendre Formfunktionen beliebiger Ordnung n
[-1,1]
innere Ableitung (2/h_x), bei Integration benoetigt
'''
def Legendre_der(order, n, X, h_x):
    l1 = (-1) /2
    l2 = (1) /2
    if (n == 0):
        return l1 *(2/h_x)
    if (n == order):
        return l2 *(2/h_x)

    # n>1
    ln_minus1=1
    ln = X
    if(n == 1):
        return ln *(2/h_x)

    for i in range(2, n+1):
        l = ln
        ln = (2*i-1) / (i) * X * ln - (i - 1) / i * ln_minus1
        ln_minus1 = l
    return ln *(2/h_x)
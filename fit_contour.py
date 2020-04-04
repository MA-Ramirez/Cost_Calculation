import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

def fun(x,y,a,b):
    #return (a*x)/(float(b*y))
    return a*y-b*x

def fun1(data,a,b):
    return a*data[1] - b*data[0]

#3x3
#KA
KKA = np.linspace(1,3,3)

#KB
KKB = np.linspace(1,3,3)

#Generar grid - KKA es y KKB es x
#Ahi hay confusion de ejes
X, Y = np.meshgrid(KKA, KKB)

Ans = np.zeros((3,3))

take = [-5.858933808835822382e-05,-2.112965386516043112e-02,-3.903618175338666174e-02,2.083671057859595374e-02,-1.891626158611038784e-03,-2.965555453410631537e-02,3.772909134794374258e-02,2.811379580390168872e-02,-9.342721789747449523e-04]

z =0
for i in range(3):
    for j in range(3):
       Ans[i,j] = take[z]
       z+=1

print Ans

ink = np.array([[1.0,2.0,3.0,1.0,2.0,3.0,1.0,2.0,3.0],[1,1,1,2,2,2,3,3,3]])
print ink[1]


    #Fit lineal A

    #Retorna pendiente de recta y varianza (same for A and B)
    #LA y LB tienen mismo len
mAA, covAA = curve_fit(fun1, ink, take)
#mA = mAA[0]
#cA = mAA[1]
#covA = covAA[0][0]

print "Param"
print mAA
print covAA
avg = (mAA[0]+mAA[1])/2.0
print "Avg"
print avg

Fi = np.ones((len(KKA),len(KKB)))

for i in range(len(KKA)):
    for j in range(len(KKB)):
        #print inI
        Fi[i,j] = fun(KKB[j], KKA[i],0.02,0.02)
        #Fi[i,j] = fun(KKB[j], KKA[i],1,1)
        print Fi[i,j], KKB[j], KKA[i]


plt.contourf(X, Y, Fi, 16, cmap="RdBu_r")
plt.colorbar()
plt.savefig("Ajuste.png")
plt.show()
plt.clf()


plt.contourf(X, Y, Ans, 16, cmap="RdBu_r")
plt.colorbar()
#plt.savefig("y_x2.png")
#plt.show()
plt.clf()

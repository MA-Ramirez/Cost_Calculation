import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

#math.log(number, base)
#natural logarithm
#np.log(number)

def fun(x,y,a,b):
    #return (a*y**0.45) - (b*x**0.45)
    return (a*np.log(y)) - (b*np.log(x))


def fun1(data,a,b):
    return a*np.log(data[0]) - b*np.log(data[1])

#Normal
#kh2 = np.array([2,6,10,14,18])
#kh3 = np.array([4,10,16,22,28,34])
#kh4 = np.array([5,13,21,29,37,45])

#Long
kh2 = np.array([2,4,6,8,10,12,14,16,18])
kh3 = np.array([4,7,10,13,16,19,22,25,28,31,34])
kh4 = np.array([5,9,13,17,21,25,29,33,37,41,45])

#Short
#kh3 = np.array([1,2,3])


def fitCC(h):
    #First column of cost matrix
    ColIni = np.genfromtxt("Results/DataCost_"+str(h)+"_W.txt", delimiter=",",usecols=0)

    #The cost matrix in nxn, so #columns = #rows
    #Lenght of first colum
    n = np.size(ColIni)

    #Define size of matrix
    if h==2:
        #KA
        KKA = kh2
        #KB
        KKB = kh2
        kh = kh2
    elif h==23:
        #KA
        KKA = kh23
        #KB
        KKB = kh23
        kh = kh23
    elif h==27:
        #KA
        KKA = kh27
        #KB
        KKB = kh27
        kh = kh27
    elif h==3:
        #KA
        KKA = kh3
        #KB
        KKB = kh3
        kh = kh3
    elif h==4:
        #KA
        KKA = kh4
        #KB
        KKB = kh4
        kh = kh4


    #Generate grid - KKA is Y and KKB is X
    X, Y = np.meshgrid(KKA, KKB)

    for i in range(1,n):
        Col = np.genfromtxt("Results/DataCost_"+str(h)+"_W.txt", delimiter=",",usecols=i)
        ColIni = np.append(ColIni, Col)

    Ans = np.zeros((n,n))

    z=0
    for j in range(n):
        for i in range(n):
           Ans[i,j] = ColIni[z]
           z+=1

    #print Ans

    fir = []
    sec = []

    for i in kh:
        for j in kh:
            fir.append(j)
            sec.append(i)

    ink = np.array([fir,sec])

    #Fit Contour

    #Retorna pendiente de recta y varianza (same for A and B)
    #LA y LB tienen mismo len
    mAA, covAA = curve_fit(fun1, ink, ColIni)
    p1 = mAA[0]
    p2 = mAA[1]

    print "Param"
    print mAA
    print covAA
    Cavg = (p1+p2)/2.0
    print "CAvg"
    print Cavg

    #Print in Cost cte file
    text_file = open("Results/Fit/CCte.csv", "a+")
    n = text_file.write(str(h)+","+ str(Cavg)+"\n")
    text_file.close()

    #Print Fit Var in file
    text_file = open("Results/Fit/Var_Fit.csv", "a+")
    n = text_file.write(str(h)+","+str(covAA[0][0])+","+str(covAA[0][1])+"\n")
    text_file.close()


    Fi = np.ones((len(KKA),len(KKB)))

    for i in range(len(KKA)):
        for j in range(len(KKB)):
            Fi[i,j] = fun(KKB[j], KKA[i],p1,p2)

    if h ==2:
        #Med and long
        L =[-0.056,-0.048,-0.04,-0.032,-0.024,-0.016,-0.008,0.0,0.008,0.016,0.024,0.032,0.04,0.048,0.056]
    elif h ==3:
        #Med
        #L = 16
        #Long
        L = [-0.105,-0.09,-0.075,-0.06,-0.045,-0.03,-0.015,0.0,0.015,0.03,0.045,0.06,0.075,0.09,0.105]
        #Short
        #L = [-0.04,-0.035,-0.03,-0.025,-0.02,-0.015,-0.01,-0.005,0.0,0.005,0.01,0.015,0.02,0.025,0.03,0.035,0.04]
    elif h ==4:
        #Med
        #L = [-0.175,-0.15,-0.125,-0.1,-0.075,-0.05,-0.025,0.0,0.025,0.05,0.075,0.1,0.125,0.15,0.175]
        #Long
        L=[-0.18,-0.16,-0.14,-0.12,-0.1,-0.08,-0.06,-0.04,-0.02,0.0,0.02,0.04,0.06,0.08,0.1,0.12,0.14,0.16,0.18]

    """
    plt.title("Fit Cost for plasmids type A ($h$ = "+str(h)+")")
    plt.contourf(X, Y, Fi, L, cmap="RdBu_r")
    plt.colorbar()
    plt.xlabel("$K_{A}$")
    plt.ylabel("$K_{B}$")
    plt.xticks(kh)
    plt.yticks(kh)
    plt.savefig("Results/Fit/Fit_"+str(h)+"_W.png")
    #plt.show()
    plt.clf()

    plt.title("Check Cost for plasmids type A ($h$ = "+str(h)+")")
    plt.contourf(X, Y, Ans, 16, cmap="RdBu_r")
    plt.colorbar()
    plt.xlabel("$K_{A}$")
    plt.ylabel("$K_{B}$")
    plt.xticks(kh)
    plt.yticks(kh)
    plt.savefig("Results/Fit/Check_"+str(h)+"_W.png")
    #plt.show()
    plt.clf()
    """

fitCC(2)
fitCC(3)
fitCC(4)


#Fit of ctes of fit function
def Fit_Cte():
    HS = np.genfromtxt("Results/Fit/CCte.csv", delimiter=",",usecols=0)
    CV = np.genfromtxt("Results/Fit/CCte.csv", delimiter=",",usecols=1)

    plt.title("Fit Cost Constant ($A$) vs Hill Coefficient ($h$)")
    plt.plot(HS,CV,c="r")
    plt.xlabel("Hill Coefficient ($h$)")
    plt.ylabel("Fit cost constant ($A$)")
    plt.xticks(HS)
    plt.savefig("Results/Fit/Cte_Fit.png")
    plt.clf()

#Fit_Cte()

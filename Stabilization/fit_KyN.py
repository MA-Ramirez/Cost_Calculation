import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

#Linear function
#m: slope
#c: y-intercept
def fun(x,m,c):
    return m*x + c

#Graph and fit of stabilization data
#Data is reported by pairs
def graphSta(h):
    #k: repression coefficient
    #nn: average number of plasmids
    K = np.genfromtxt(str(h)+"Stable.csv", delimiter=",", usecols=0)
    N = np.genfromtxt(str(h)+"Stable.csv", delimiter=",", usecols=1)

    #Ph: parameters obtained of fit
    #Ph[0]: slope, Ph[1]: y-intercept
    #coh: covariance of parameters
    Ph, coh = curve_fit(fun, K, N)
    print Ph, coh

    #Visualization of data and fit
    #Given a K obtain NN
    plt.plot(K,N,c="r",label="Data")
    plt.plot(K,fun(K,Ph[0],Ph[1]),c="g",label="Fit")
    plt.title("Stabilization Fit $h=$"+str(h))
    plt.xlabel("Repression Coefficient ($K$)")
    plt.ylabel("Average amount of plasmids ($\overline{n}$)")
    plt.plot(0,0,c="white",label="$\overline{n}$="+str(round(Ph[0],3))+"$K$+"+str(round(Ph[1],3)))
    plt.legend(loc=7,fontsize="small")
    #plt.text(29, 150, "$\overline{n}$="+str(round(Ph[0],3))+"$K$+"+str(round(Ph[1],3)), bbox=dict(facecolor='g', alpha=0.2))
    plt.savefig("StableGraphs/Fit_Stabi_h"+str(h)+".png")
    plt.clf()

    #Returns slope and y-intercept fit parameters
    return Ph

#Execute code
Ph2=graphSta(2)
Ph23=graphSta(2.3)
Ph27=graphSta(2.7)
Ph3=graphSta(3)
Ph33=graphSta(3.3)
Ph37=graphSta(3.7)
Ph4=graphSta(4)
Ph43=graphSta(4.3)
Ph47=graphSta(4.7)
Ph5=graphSta(5)

#------General Fit----------
#m: slope
#c: y-intercept
def funG(x,a,c):
    return a*np.exp(-x) +c

#Fit of slope for each Hill coefficient array (y)
m = np.array([Ph2[0],Ph23[0],Ph27[0],Ph3[0],Ph33[0],Ph37[0],Ph4[0],Ph43[0],Ph47[0],Ph5[0]])

#Hill coefficient array (x)
HH = np.array([2.0,2.3,2.7,3.0,3.3,3.7,4.0,4.3,4.7,5.0])

#Phm: parameters obtained of fit
#com: covariance of parameters
Phm, cohm = curve_fit(funG, HH, m)
print Phm,cohm

#plot
plt.plot(HH,m,c="indigo", label="Slopes")
plt.plot(HH,funG(HH,Phm[0],Phm[1]),c="cyan",label="Fit")
plt.title("Slope ($m_{nK}$) vs Hill Coefficient ($h$)")
plt.xlabel("Hill Coefficient ($h$)")
plt.ylabel("Slope of graph $\overline{n}$ vs $K$ ($m_{nK}$)")
plt.legend(loc=1)
plt.text(4, 9, "$m_{nK}$="+str(round(Phm[0],3))+"$e^{(-h)}$+"+str(round(Phm[1],3)), bbox=dict(facecolor='cyan', alpha=0.2))
plt.savefig("StableGraphs/GeneralFit_m_vs_h.png")
plt.clf()

cc = np.array([Ph2[1],Ph23[1],Ph27[1],Ph3[1],Ph33[1],Ph37[1],Ph4[1],Ph43[1],Ph47[1],Ph5[1]])
plt.plot(HH,cc,c="Orange")
plt.title("y-intercept ($c_{nK}$) vs Hill Coefficient ($h$)")
plt.xlabel("Hill Coefficient ($h$)")
plt.ylabel("y-intercept of graph $\overline{n}$ vs $K$ ($c_{nK}$)")
#plt.legend(loc=1)
plt.savefig("StableGraphs/GeneralFit_c_vs_h.png")

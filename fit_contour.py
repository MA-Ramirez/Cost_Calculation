import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

#Code by Maria Alejandra Ramirez
#Script: fit_contour.py

#---------------CONTOUR_GRAPH_FIT--------------

#Performs the fit of countour graphs that show the cost of each
# one-on-one competition
#The fit allows to find the relationship between c,KA and KB

#Proposed function for the fit
#Function to graph results (2)
def fun(x,y,a,b):
    return (a*np.log(y)) - (b*np.log(x))

#Function to obtain the fit (1)
def fun1(data,a,b):
    return a*np.log(data[0]) - b*np.log(data[1])

#RANGE OF KS USED FOR EACH HILL COEFFICIENT (H)

#Short
#kh3 = np.array([1,2,3])

#Medium
#kh2 = np.array([2,6,10,14,18])
#kh3 = np.array([4,10,16,22,28,34])
#kh4 = np.array([5,13,21,29,37,45])

#Long Sample
kh2 = np.array([2,4,6,8,10,12,14,16,18])
kh3 = np.array([4,7,10,13,16,19,22,25,28,31,34])
kh4 = np.array([5,9,13,17,21,25,29,33,37,41,45])
kh23 = np.array([3,6,9,12,15,18,21,24])
kh27 = np.array([3,6,9,12,15,18,21,24,27,30])
kh33 = np.array([4,8,12,16,20,24,28,32,36,40])
kh37 = np.array([4,10,16,22,28,34,40,46])
kh43 = np.array([5,10,15,20,25,30,35,40,45,50])
kh47 = np.array([5,12,19,26,33,40,47,54])
kh5 = np.array([6,12,18,24,30,36,42,48,54])

#Function to graph the fit result contour graph
# with the same color bar range than in the results, in order to
# make a clear visual comparisons of the results with the fit
#Param: st Negative start of the color bar range
#Param: hop Interval to the next color value limit
#Return: L Array with the values of the color bar
def barRange(st,hop):
    s = st
    ini = s
    L = []
    br = True
    while br == True:
        if s == (-1*ini):
            br = False
        L.append(s)
        s = round(s+hop,3)
    return L

#Performs the fit of the cost results given a Hill coefficient value
#Param: h Hill coefficient
#Return: Graph of the fit results
# CSV file with param of fit (CCte.csv)
# and with variance of the fit (Var_Fit.csv)
def fitCC(h):
    #Obtains the results from the txt files
    #First column of cost matrix txt file
    ColIni = np.genfromtxt("Results/DataCost_"+str(h)+".txt", delimiter=",",usecols=0)

    #The cost matrix in nxn, so #columns = #rows
    #Lenght of first colum
    n = np.size(ColIni)

    #Define size of matrix
    #The size is different for each Hill coefficient that depends
    # on the Ks range stated previously
    if h==2:
        #KA
        KKA = kh2
        #KB
        KKB = kh2
        kh = kh2
    elif h==2.3:
        #KA
        KKA = kh23
        #KB
        KKB = kh23
        kh = kh23
    elif h==2.7:
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
    elif h==3.3:
        #KA
        KKA = kh33
        #KB
        KKB = kh33
        kh = kh33
    elif h==3.7:
        #KA
        KKA = kh37
        #KB
        KKB = kh37
        kh = kh37
    elif h==4:
        #KA
        KKA = kh4
        #KB
        KKB = kh4
        kh = kh4
    elif h==4.3:
        #KA
        KKA = kh43
        #KB
        KKB = kh43
        kh = kh43
    elif h==4.7:
        #KA
        KKA = kh47
        #KB
        KKB = kh47
        kh = kh47
    elif h==5:
        #KA
        KKA = kh5
        #KB
        KKB = kh5
        kh = kh5


    #Generate grid - KKA is Y and KKB is X
    X, Y = np.meshgrid(KKA, KKB)

    #Obtains the rest of cost results data
    for i in range(1,n):
        Col = np.genfromtxt("Results/DataCost_"+str(h)+".txt", delimiter=",",usecols=i)
        ColIni = np.append(ColIni, Col)

    #Array where the cost results from the txt file are arranged.
    # This array helps to establish how to properly array the txt data in
    # numpy arrays
    Ans = np.zeros((n,n))

    #Loop to organize the data from the txt files into a numpy array
    z=0
    for j in range(n):
        for i in range(n):
           Ans[i,j] = ColIni[z]
           z+=1

    #K range for the x and y axis obtained from the arrays previously stated
    # This range is different for each Hill coefficient
    fir = []
    sec = []

    #Organizes the K ranges appropiately
    for i in kh:
        for j in kh:
            fir.append(j)
            sec.append(i)

    ink = np.array([fir,sec])

    #-----FIT----------
    #Fit Contour
    #Performs the fit of the data from the txt files of cost results
    # curve_fit(function, xdata, ydata)
    # xdata refers to the independent variable and
    # ydata to the depende variable
    # xdata--> Ks range    ydata--> cost results
    mAA, covAA = curve_fit(fun1, ink, ColIni)
    p1 = mAA[0]
    p2 = mAA[1]

    #Print the fit parameters and their variance
    print "Param"
    print mAA
    print covAA

    #Calculate and print the average of the parameters
    Cavg = (p1+p2)/2.0
    print "CAvg"
    print Cavg

    #Print Cavg in Cost cte file
    text_file = open("Results/Fit/CCte.csv", "a+")
    n = text_file.write(str(h)+","+ str(Cavg)+"\n")
    text_file.close()

    #Print Fit Var in file
    text_file = open("Results/Fit/Var_Fit.csv", "a+")
    n = text_file.write(str(h)+","+str(covAA[0][0])+","+str(covAA[0][1])+"\n")
    text_file.close()


    #Graphs the fit results
    Fi = np.ones((len(KKA),len(KKB)))

    for i in range(len(KKA)):
        for j in range(len(KKB)):
            Fi[i,j] = fun(KKB[j], KKA[i],p1,p2)

    #Assigns the correct color bar for each contour graph
    # This allows to make a direct visual comparison between the results
    # and their fit
    if h ==2:
        #Med and long
        L = barRange(-0.056,0.008)
    elif h ==2.3:
        L = barRange(-0.064,0.008)
    elif h ==2.7:
        L = barRange(-0.09,0.01)
    elif h ==3:
        #Short
        #L = barRange(-0.04,0.005)

        #Med and long
        L = barRange(-0.105,0.015)
    elif h ==3.3:
        L = barRange(-0.12,0.015)
    elif h ==3.7:
        L = barRange(-0.16,0.02)
    elif h ==4:
        #Med
        #L = barRange(-0.175,0.025)

        #Long sample
        L = barRange(-0.18,0.02)
    elif h ==4.3:
        L = barRange(-0.2,0.025)
    elif h ==4.7:
        L = barRange(-0.24,0.03)
    elif h ==5:
        L = barRange(-0.27,0.03)

    #Graph of the fit
    plt.title("Fit Cost for plasmids type A ($h$ = "+str(h)+")")
    plt.contourf(X, Y, Fi, L, cmap="RdBu_r")
    plt.colorbar()
    plt.xlabel("$K_{A}$")
    plt.ylabel("$K_{B}$")
    plt.xticks(kh)
    plt.yticks(kh)
    plt.savefig("Results/Fit/Fit_"+str(h)+".png")
    #plt.show()
    plt.clf()

    #Graph of the Results
    # This graph is used to check that the values are being shown in the
    # correct way
    """
    plt.title("Check Cost for plasmids type A ($h$ = "+str(h)+")")
    plt.contourf(X, Y, Ans, 16, cmap="RdBu_r")
    plt.colorbar()
    plt.xlabel("$K_{A}$")
    plt.ylabel("$K_{B}$")
    plt.xticks(kh)
    plt.yticks(kh)
    plt.savefig("Results/Fit/Check_"+str(h)+".png")
    #plt.show()
    plt.clf()
    """

#Executes the code for each set of results
# Each set is characterized by its Hill function
fitCC(2)
fitCC(2.3)
fitCC(2.7)
fitCC(3)
fitCC(3.3)
fitCC(3.7)
fitCC(4)
fitCC(4.3)
fitCC(4.7)
fitCC(5)


#--------------------FIT_OF_CTE-----------------------------

#Performs the fit needed to find the relationship between the Constant
# previously found with the Hill coefficient
#Return: Graphs Cte vs h
def Fit_Cte():
    #Obtains the fit cte recorded in a CSV file, along with its
    # corresponding Hill coefficient
    HS = np.genfromtxt("Results/Fit/CCte.csv", delimiter=",",usecols=0)
    CV = np.genfromtxt("Results/Fit/CCte.csv", delimiter=",",usecols=1)

    #Proposed fit Function
    #Linear function
    #m: slope
    #c: y-intercept
    def funA(x,m,c):
        return m*x + c

    #Ph: parameters obtained of fit
    #Ph[0]: slope, Ph[1]: y-intercept
    #coh: covariance of parameters
    #curve_fit(function, xdata, ydata)
    Ph1, coh1 = curve_fit(funA, HS, CV)
    print Ph1, coh1

    #Graphs Cte vs h
    plt.title("Fit Cost Constant ($A$) vs Hill Coefficient ($h$)")
    plt.plot(HS,CV,c="r",label="Data")

    #Shows the fit results in the graph
    plt.plot(HS,funA(HS,Ph1[0],Ph1[1]),c="g",label="Fit")
    plt.plot(0,0,c="white",label="$A$="+str(round(Ph1[0],3))+"$h$"+str(round(Ph1[1],3)))

    plt.xlabel("Hill Coefficient ($h$)")
    plt.ylabel("Fit cost constant ($A$)")
    plt.xticks(HS)
    plt.xlim((2,5))
    plt.legend(loc=2,fontsize="small")
    plt.savefig("Results/Fit/Cte_Fit.png")
    plt.clf()

#Execute the code
Fit_Cte()

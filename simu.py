import numpy as np
import matplotlib.pyplot as plt
import time
from scipy.optimize import curve_fit

#------WHOLE SIMULATION-------
#Whole simulation to obtain the contour graph of the cost

#HILL FUNCTION FOR A REPRESSOR
#Param: h Hill coefficient
#Param: k repression coefficient
#Param: x TOTAL amount of plasmids
#Return: Bmax, Y
def repression(h,k,x):
    Bmax = 10.0
    a = float(1 + (x/k)**h)
    y = Bmax / a
    return Bmax,y

#------STABILIZATION-------------------
#Stabilization data is done in another script for optimization

#-----MORAN_FIT----------------
#BIRTH
#Generates a reproduction decision A,B or None
#Param: h Hill coefficient
#Param: kA K of plasmid type A
#Param: kB K of plasmid type B
#Param: A amount of plasmids type A
#Param: B amount of plasmids type B
#Return: amount of A, B after reproduction decision
def Birth(h,kA,kB,A,B):

    #Total amount of plasmids
    tot = float(A+B)

    #Type A (kA)
    #nor: Bmax to normalize y, so that it has same range as random.random()
    norA, yA = repression(h,kA,tot)

    #Nomalized reproduction probability
    probA = float(yA/norA)

    #Type B (kB)
    norB, yB =repression(h,kB,tot)
    probB = float(yB/norB)

    #Threshold A
    TA = (probA*A) / tot

    #Threshold B
    TB = (probB*B) / tot

    #Random number between 0-1
    Q1 = np.random.random()
    Q2 = np.random.random()

    #Birth A
    if Q1 <= TA:
        A+= 1
    #Birth B
    if Q2 <= TB:
        B += 1
    #"else" ninguno

    return A,B,probA,probB

#Death
#Generates a death decision A,B
#Param: A amount of plasmids type A
#Param: B amount of plasmids type B
#Return: amount of A, B after death decision
def Death(A,B):

    #Normalized
    #Threshold A
    TDA = (0.5*A) / (0.5*A+0.5*B)

    #Threshold B
    TDB = (0.5*B) / (0.5*A+0.5*B)

    #random number between 0-1
    QQ = np.random.random()

    #Either a death of A or B should happen
    #Death A
    if QQ <= TDA:
        A = A-1

    #Death B
    else:
        B = B-1

    return A,B

#Normalize
#Normalize the amount of plasmids to 1
#Param: AA amount of plasmids type A
#Param: BB amount of plasmids type B
#Return: cA, cB normalized amounts to 1
def nor(AA,BB):
    cA = AA / float(AA+BB)
    cB = BB / float(AA+BB)
    return cA, cB

#Process
#The whole process is performed
#Param: inA initial amount of plasmids type A
#Param: inB initial amount of plasmids type B
#Param: h Hill coefficient
#Param: kA K of plasmid type A
#Param: kB K of plasmid type B
#Return: pA, pB arrays with data of each event for the plasmids type A and B
#Return: npA, npB normalized arrays with data of each event for the plasmids type A and B
def Go(inA, inB, h, kA, kB):
    #Stop marker to avoid infinite loops
    S = 0

    #Arrays with data of each event for the plasmids
    pA = np.array([inA])
    pB = np.array([inB])
    #Normalization
    cinA, cinB = nor(inA,inB)
    npA = np.array([cinA])
    npB = np.array([cinB])

    #Upper and lower cuts for the normalized amount of plasmids
    #  This avoids the code to run for amounts where the reproduction
    #  probability is determined by 1/N (irrelevant for this project)
    #tUP = 1 - (1/float(inA+inB))
    #tDw = (1/float(inA+inB))
    tUP = 0.75
    tDw = 0.25

    #Limit of number of events for the simulation
    # this limit is called LEN
    LEN = 400

    #Markers to break loop
    C1 = True
    C2 = True

    #Code will run while the length of the array is lower than the LEN limit
    while len(npA) < LEN and S == 0:
        #Condition 1 (Upper cuts both types)
        if cinA >=tUP or cinB >= tUP:
            C1 = False

        #Condition 2 (Lower cuts both types)
        if cinA <= tDw or cinB <= tDw:
            C2 = False

        #Condition 3
        #Code will run as long as no plasmid takes completely over the population
        if C1 == True and C2 == True and S ==0:
            #Run whole process
            #Primer for condition below
            A,B,probA,probB = Birth(h,kA,kB,inA,inB)

            #--------------BIRTH----------------

            #There will be reproduction while reproduction probability
            # is above 2%
            while (probA >= 0.02 or probB >= 0.02) and S == 0:

                A,B,probA,probB = Birth(h,kA,kB,inA,inB)
                inA = A
                inB = B
                pA = np.append(pA,inA)
                pB = np.append(pB,inB)

                #Normalization
                cinA, cinB = nor(inA,inB)
                npA = np.append(npA,cinA)
                npB = np.append(npB,cinB)

                #Code will run while the length of the array is lower
                # than the LEN limit
                if len(npA) >= LEN:
                    S = 1
                    break

                #Condition 1 (Upper cuts both types)
                if cinA >= tUP or cinB >= tUP:
                    C1 = False
                #Condition 2 (Lower cuts both types)
                if cinA <= tDw or cinB <= tDw:
                    C2 = False
                #If markers to break loop have been modified, break process
                if (C1 != True or C2 != True):
                    S = 1

            #TOTAL amount of plasmids after REPRODUCTIVE PHASE
            # this amount is called fixed
            fixed = inA + inB


            #--------------DEATH----------------

            #Random death happens until half of fixed is reached
            while (inA+inB > (fixed/2.0)) and S == 0:

                A,B = Death(inA,inB)
                inA = A
                inB = B
                pA = np.append(pA,inA)
                pB = np.append(pB,inB)

                #Normalization
                cinA, cinB = nor(inA,inB)
                npA = np.append(npA,cinA)
                npB = np.append(npB,cinB)

                #Code will run while the length of the array is lower
                # than the LEN limit
                if len(npA) >= LEN:
                    S = 1
                    break

                #Condition 1 (Upper cuts both types)
                if cinA >= tUP or cinB >= tUP:
                    C1 = False
                #Condition 2 (Lower cuts both types)
                if cinA <= tDw or cinB <= tDw:
                    C2 = False

                #If markers to break loop have been modified, break process
                if (C1 == False or C2 == False):
                    S = 1


        #Creation of normalized arrays
        #Limits avoid the code to run for amounts where the reproduction
        #  probability is determined by 1/N (irrelevant for this project)
        #Upper normalized limit of plasmids
        Unl = tUP
        #Lower normalized limit of plasmids
        Lnl = tDw

        #Markers and conditions to break the loop
        if C1 == False and C2 == False:
            S = 0
            if cinA <= Lnl:
                 npA = np.append(npA,tDw)
                 npB = np.append(npB,tUP)
            elif cinA >= Unl:
                 npA = np.append(npA,tUP)
                 npB = np.append(npB,tDw)
            elif cinB <= Lnl:
                npA = np.append(npA,tUP)
                npB = np.append(npB,tDw)
            elif cinB >= Unl:
                npA = np.append(npA,tDw)
                npB = np.append(npB,tUP)

            #General process
            #Code will run while the length of the array is lower
            # than the LEN limit
            if len(npA) >= LEN:
                S = 1
                break

    #Results of the whole process
    #pA, pB arrays with amount of plasmids for each event (type A and B)
    #npA, npB arrays of pA and pB normalized to 1
    return pA, pB, npA, npB


#--------REPETITION OF THE PROCESS & GRAPHS------------
#Repetition of the process and corresponding
#Param: rounds times that the process is repeated +1
#Param: rep repetitions of the same general simulation
#Param: cc iteration of the current repetition
#Param: inA initial amount of plasmids type A
#Param: inB initial amount of plasmids type B
#Param: h Hill coefficient
#Param: kA K of plasmid type A
#Param: kB K of plasmid type B
#rep and cc are used to only graph one general simulation (optimization)
def repetitionHist(rounds, rep, cc, inA, inB, h, kA, kB):

    #Start measuring simulation time
    t0 = time.time()

    #-------Limits for the average and the corresponding fit--------
    #Data of average of total plasmids is calculated in another script
    # this data is recorded on csv files for optimization
    NKb = np.genfromtxt("Total/"+str(h)+"_"+str(kA)+"Total.csv", delimiter=",",usecols=1)
    Navg = np.genfromtxt("Total/"+str(h)+"_"+str(kA)+"Total.csv", delimiter=",",usecols=2)

    #Average of total plasmids for this competition
    inNKb = np.where(NKb==kB)[0]
    NNavg = Navg[inNKb]

    #Cuts of fit
    CutIniFit = 0
    #print CutIniFit
    FinalFit = int(CutIniFit+50)
    #print FinalFit

    #Cut final for average line
    LENF = 400

    #Return of Go function (check Go comments)
    pA, pB, npA, npB = Go(inA, inB, h, kA, kB)


    #-----First round-----
    #Array of average per event point
    # initialized with data of first round
    primerA = npA[:LENF]
    primerB = npB[:LENF]

    #Graph for first round simulation
    # type A = blue   type B = green
    #Condition to graph only 1 general simulation
    if (rep == cc) == True:
        plt.plot(np.linspace(0,len(npA), num = len(npA)), npA, c = "b")
        plt.plot(np.linspace(0,len(npB), num = len(npB)), npB, c = "g")

    #Subsequent rounds
    for i in range(rounds):

        #Return of Go function (check Go comments)
        pA, pB, npA, npB = Go(inA, inB, h, kA, kB)

        #Graph for simulations
        # type A = blue   type B = green
        #Condition to graph only 1 general simulation
        if (rep == cc) == True:
            plt.plot(np.linspace(0,len(npA), num = len(npA)), npA, c = "b")
            plt.plot(np.linspace(0,len(npB), num = len(npB)), npB, c = "g")

        #Sum of plasmid amounts per event point
        primerA += npA[:LENF]
        primerB += npB[:LENF]

        #Done simulation round


    #------Calculation of average simulation line------
    primerA = primerA/float(rounds+1)
    primerB = primerB/float(rounds+1)
    meanA = primerA
    meanB = primerB

    #--------Fit-------
    #Cuts of fit
    #CutIniFit = 0
    #FinalFit = CutIniFit+50

    #Cut final for average line
    #LENF = 250

    #Fit is done in the average simulation line
    #Array with delimitation that it is important for the fit
    LA = meanA[CutIniFit:FinalFit]
    LB = meanB[CutIniFit:FinalFit]

    #Proposed linear fit function
    #Param: m slope
    #Param: c y-intercept
    #  should be around 0.5 since there are the same number of plasmids
    #  in the beginning of the simulation
    def fun(x,m,c):
        return m*x + c

    #Linear fit type A
    #Returns the parameters of the fit

    #LA y LB are of the same length
    mAA, covAA = curve_fit(fun, np.linspace(CutIniFit, FinalFit, num = len(LA)), LA)

    #The slope is the variable of interest
    #Slope parameter
    mA = mAA[0]
    #y-intercept parameter
    cA = mAA[1]
    #Variance of the slope as a fit paramater
    covA = covAA[0][0]

    #Print in Log txt ***
    text_file = open("Results/Log_"+str(h)+".txt", "a+")
    n = text_file.write("Fit Params A: "+str(mAA)+", "+str(covA)+"\n")
    text_file.close()
    #print mAA, covA

    #Linear fit type B
    #Same but for type B
    mBB, covBB = curve_fit(fun, np.linspace(CutIniFit, FinalFit, num = len(LB)), LB)
    #Slope parameter
    mB = mBB[0]
    #y-intercept parameter
    cB = mBB[1]
    #Variance of the slope as a fit paramater
    covB = covBB[0][0]

    #Print in Log txt ***
    text_file = open("Results/Log_"+str(h)+".txt", "a+")
    n = text_file.write("Fit Params B: "+str(mBB)+", "+str(covB)+"\n")
    text_file.close()
    #print mBB, covB


    #----------COST_CALCULATION---------
    #Average of total plasmids for this competition
    #NNavg = Navg[inNKb]

    #Cost of type A
    FA = 2*mA*(NNavg)

    #Cost of type B
    FB = 2*mB*(NNavg)

    #Print in Log txt ***
    text_file = open("Results/Log_"+str(h)+".txt", "a+")
    n = text_file.write("Fitness/Cost A: "+str(FA)+"\n")
    text_file.close()

    #-----------------GRAPHS---------------------
    #Graph of average simulation line and corresponding fit

    #Condition to graph only 1 general simulation
    if (rep == cc) == True:
        #Labels of the simulation
        # type A = blue   type B = green
        plt.plot(0,0, c = "b", label = "$K_{A}$ = " + str(kA))
        plt.plot(0,0, c = "g", label = "$K_{B}$ = " + str(kB))

        #----Plot average simulation line---
        plt.plot(np.linspace(0,len(meanA), num = len(meanA)), meanA, c = "deepskyblue", label="Avg $K_{A}$", linewidth = 3)
        plt.plot(np.linspace(0,len(meanB), num = len(meanB)), meanB, c = "lime", label="Avg $K_{B}$", linewidth = 3)

        #----Plot fit----
        #x axis info for plotting fit
        xxxa = np.linspace(CutIniFit, FinalFit, num = len(LA))
        xxxb = np.linspace(CutIniFit, FinalFit, num = len(LB))

        #Plot fit of average line simulation
        plt.plot(xxxa, fun(xxxa, mA*np.ones(len(LA)) , cA*np.ones(len(LA))) , '--',c = "deeppink", label="Fit $K_{A}$", linewidth = 4)
        plt.plot(xxxb, fun(xxxb, mB*np.ones(len(LB)), cB*np.ones(len(LB))) ,'--', c = "y", label="Fit $K_{B}$", linewidth = 4)

        #Info of the general graph
        # info of cost for type A
        plt.plot(0,0, c = "white", label = "$u_{A}$ = " + str(round(mA,5)))
        plt.plot(0,0, c = "white", label = "$u_{B}$ = " + str(round(mB,5)))
        #plt.plot(0,0, c = "k", label = "Rounds = "+str(rounds+1))

        #Upper limit should match with Go LEN
        #General simulation
        plt.xlim((0,400))
        plt.title("General simulation $K_{A}$ = "+str(kA)+" vs $K_{B}$ = "+str(kB)+" ($h$="+str(h)+")")
        plt.xlabel("Events")
        plt.ylabel("Normalized amount of plasmids ($\eta$)")
        plt.legend(loc = 4, fontsize = "x-small")
        plt.savefig("SimulationGraphs/Simu_h"+str(h)+"/Graph_" + str(kA) + "_" + str(kB) + "_" + str(int(h)) + ".png")
        plt.clf()

    #Return
    #FA: the calculated cost of type A
    #covA: variance of the slope fit parameter of type A
    return FA, covA


#Start measuring simulation time
t0 = time.time()

#Generates repetitions of the whole simulation to find an average of the Cost
# for the same general simulation. It also gives the average of the Variance
# of the fit of the slope
#Param: h Hill coefficient
#Param: kA K of plasmid type A
#Param: kB K of plasmid type B
#Param: inI Initial amount of plasmids (same for type A and B)
#Param: rep Repetitions of the whole simulation to find the average
#Return: avgFA Average of cost for the repetitions of the same simulation
#Return: avgvA Average of variance of the fit of the slope for the repetitions
# of the same simulation
def full(h,kA,kB,inI,rep):

     #Arrays of fitness y variance
     lFA = np.array([])
     lva = np.array([])

     #Print in Log txt ***
     text_file = open("Results/Log_"+str(h)+".txt", "a+")
     n = text_file.write("Competition: (kA,kB) "+str(kA)+", "+str(kB)+"\n")
     text_file.close()

     #General rounds are rounds +1
     rounds = 599

     #Marker that helps to graph only one general simulation (optimization)
     cc = 0

     #Perform the general simulation rep times and find the average of their
     # results
     for i in range(rep):
        cc += 1
        #Perform the general simulation rep times
        FA, va = repetitionHist(rounds,rep,cc,inI,inI,h,kA,kB)
        lFA = np.append(lFA,FA)
        lva = np.append(lva,va)

     #Average of cost for type A
     avgFA = np.mean(lFA)
     #Average of variance of the fit slope for type A
     avgva = np.mean(lva)

     #Returns the average of cost and variance of the fit slope
     # for type A
     return avgFA, avgva


#----------------COST_CONTOUR_FIGURE----------------

#Creates the contour figure for the cost of type A (nxn graph)
#The cost for each competition is the average of rep repetitions
# for the whole simulation of the given competition
#Param: start Initial K
#Param: stop Final K
#Param: hop Amount of values between start and stop
#Param: Hill coefficient
#Param: rep Repetitions for the whole simulation of a given competition
#Return: Graph of Cost, Graph of variance of slope parameter
#Return: txt with info of cost and variance of slope parameter
def contourG(start,stop,hop,h,rep):
    #Arrays for kA y kB that are gonna be used (nxn)
    # Example: np.linspace(1,3,3) = 1,2,3
    # Example: np.linspace(2,4,3) = 2,3,4
    KKA = np.linspace(start,stop,hop)
    KKB = np.linspace(start,stop,hop)

    #Print in Log txt ***
    text_file = open("Results/Log_"+str(h)+".txt", "a+")
    n = text_file.write("RANGE OF CALCULATION: "+str(KKA)+"\n")
    text_file.close()

    #Generate Grid - KKA is x KKB is y
    X, Y = np.meshgrid(KKA, KKB)

    #Results of the cost of A = FA
    Fi = np.ones((len(KKA),len(KKB)))

    #Print in Log txt ***
    text_file = open("Results/Log_"+str(h)+".txt", "a+")
    n = text_file.write("LEN OF RANGE: (KKA,KKB) "+str(len(KKA))+","+str(len(KKB))+"\n \n")
    text_file.close()

    #Variance of the slope for type A
    vi = np.ones((len(KKA),len(KKB)))

    #Data of stabilization is calculated in another script
    # this data is recorded on csv files for optimization
    Sta = np.genfromtxt("Stabilization/"+str(h)+"Stable.csv", delimiter=",",usecols=1)
    kIndex = np.genfromtxt("Stabilization/"+str(h)+"Stable.csv", delimiter=",",usecols=0)

    #Initial inA y inB (same) is the average of both stabilizations, half half
    #Simulation every competition
    #Counter for location
    xF = 0
    for i in KKA:
        #Stabilization value for kA
        inKA = np.where(kIndex==i)[0]
        Stable1 = Sta[inKA]

        #Print in Log txt ***
        text_file = open("Results/Log_"+str(h)+".txt", "a+")
        n = text_file.write("STABLE1: "+ str(Stable1)+"   i = "+str(i)+"\n")
        text_file.close()

        #Counter for location
        yF = 0

        for j in KKB:
            #Stabilization value for kB
            inKB = np.where(kIndex==j)[0]
            Stable2 = Sta[inKB]

            #Print in Log txt ***
            text_file = open("Results/Log_"+str(h)+".txt", "a+")
            n = text_file.write("Stable2: "+ str(Stable2)+"   j = "+str(j)+"\n")
            text_file.close()

            #Average of stabilization
            AvgStable = (Stable1 + Stable2)/2.0
            #Turns it into int, no decimal amount of plasmids
            inI = int(AvgStable/2.0)

            #Print in Log txt ***
            text_file = open("Results/Log_"+str(h)+".txt", "a+")
            n = text_file.write("Average Stable: "+ str(AvgStable)+"   InI: "+str(inI)+"\n")
            text_file.close()

            #Assigns calculated cost and variance
            #Check comments for full - to only calculate cost use full
            Fi[xF,yF], vi[xF,yF] = full(h,i,j,inI,rep)

            #Print in Log txt ***
            text_file = open("Results/Log_"+str(h)+".txt", "a+")
            n = text_file.write("---------//---------\n")
            text_file.close()

            #Counter for location
            yF+=1

        #Print in Log txt ***
        text_file = open("Results/Log_"+str(h)+".txt", "a+")
        n = text_file.write("\n")
        text_file.close()
        #Counter for location
        xF+=1

    #txt with the data of the cost and of the variance of slope parameter
    np.savetxt("Results/DataCost_"+str(h)+".txt", Fi, delimiter=',')
    np.savetxt("Results/DataVariance_"+str(h)+".txt", vi, delimiter=',')

    #COST_CONTOUR_FIGURE
    #x,y,z,levels = number of divisions, cmap = color
    plt.contourf(X, Y, Fi, 16, cmap="RdBu_r")
    plt.colorbar()
    plt.xlabel("$K_{A}$")
    plt.ylabel("$K_{B}$")
    plt.xticks(KKA)
    plt.yticks(KKB)
    plt.title("Cost for plasmids type A ($h$ = "+str(h)+")")
    plt.savefig("Results/CostA_"+str(h)+".png")
    plt.clf()

    #Variance contour figure
    plt.contourf(X, Y, vi, cmap="YlOrRd")
    plt.colorbar()
    plt.xlabel("$K_{A}$")
    plt.ylabel("$K_{B}$")
    plt.xticks(KKA)
    plt.yticks(KKB)
    plt.title("Variance of the slope fit parameter type A ($h$ = "+str(h)+")",fontsize="small")
    plt.savefig("Results/VarianceA_"+str(h)+".png")
    plt.clf()

    #Time of all the total simulation
    tsim = round(time.time()-t0,3)

    #Print in Log txt ***
    text_file = open("Results/Log_"+str(h)+".txt", "a+")
    n = text_file.write("Total Simu Time "+ str(tsim))
    text_file.close()

#contourG(start,stop,hop,h,rep)
#h=2
#contourG(2,18,9,2,3)
#h=3
#contourG(4,34,11,3,3)
#contourG(1,3,3,3,1)
#h=4
#contourG(5,45,11,4,3)

#h2.3
contourG(3,24,8,2.3,3)
#h=2.7
contourG(3,30,10,2.7,3)
#h=3.3
contourG(4,40,10,3.3,3)
#h=3.7
contourG(4,46,8,3.7,3)
#h=4.3
contourG(5,50,10,4.3,3)
#h=4.7
contourG(5,54,8,4.7,3)
#h=5
contourG(6,54,9,5,3)

import numpy as np
import matplotlib.pyplot as plt
import time
from scipy.optimize import curve_fit


#------SINGLE COMPETITION-------
#Code to obtain the graph of a single competition
# with rep repetitions

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

    #Threshold none
    #TN = 1.0 - TA - TB

    #Random number between 0-1
    Q = np.random.random()

    #Birth A
    if Q <= TA:
        A+= 1
    #Birth B
    elif Q > TA and Q <= TA + TB:
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
def Go(inA, inB, h, kA, kB):
    #Stop marker to avoid infinite loops
    S = 0
    #Win count
    winA = 0
    winB = 0

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
    tUP = 1 - (1/float(inA+inB))
    tDw = (1/float(inA+inB))

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
            #Primer Birth
            A,B,probA,probB = Birth(h,kA,kB,inA,inB)
            inA = A
            inB = B
            pA = np.append(pA,inA)
            pB = np.append(pB,inB)

            #Normalization
            cinA, cinB = nor(inA,inB)
            npA = np.append(npA,cinA)
            npB = np.append(npB,cinB)


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
                if (C1 != True and C2 != True):
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
                if (C1 == False and C2 == False):
                    S = 1

                #If one plasmid has already took completely over the population
                # the process is stopped
                elif inA == 0 or cinA <= tDw:
                    winB += 1
                    break
                elif inB == 0 or cinB >= tUP:
                    winA += 1
                    break

        #Creation of normalized arrays
        #Limits avoid the code to run for amounts where the reproduction
        #  probability is determined by 1/N (irrelevant for this project)
        #Upper normalized limit of plasmids
        Unl = 0.75
        #Lower normalized limit of plasmids
        Lnl = 0.25

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
    #winA, winB win counts = amount of times each plasmids
    # takes over the population (partially)
    return pA, pB, npA, npB, winA, winB


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
    #Counts win A (partial dominations)
    GwinA = 0
    #Counts win B (partial dominations)
    GwinB = 0

    #Start measuring simulation time
    t0 = time.time()

    #Array for all events of type A, in terms of plasmid amount, for all rounds
    AA = np.array([])

    #Array for all events of type B, in terms of plasmid amount, for all rounds
    BB = np.array([])

    #Top limit for the average and the corresponding fit
    LENF = 250

    #Return of Go function (check Go comments)
    pA, pB, npA, npB, winA, winB = Go(inA, inB, h, kA, kB)


    #-----First round-----
    #Array of average per event point
    # initialized with data of first round
    primerA = npA[:LENF]
    primerB = npB[:LENF]

    #Graph for first round simulation
    # type A = blue   type B = green
    #Condition to graph only 1 general simulation
    if (rep == cc) == True:
        plt.plot(np.linspace(0,len(npA[:(inA+inB)]), num = len(npA)), npA, c = "b")
        plt.plot(np.linspace(0,len(npB[:(inA+inB)]), num = len(npB)), npB, c = "g")

    #Subsequent rounds
    for i in range(rounds):

        #Return of Go function (check Go comments)
        pA, pB, npA, npB, winA, winB = Go(inA, inB, h, kA, kB)

        #Sum of wins for each type of plasmid
        GwinA += winA
        GwinB += winB

        #Graph for simulations
        # type A = blue   type B = green
        #Condition to graph only 1 general simulation
        if (rep == cc) == True:
            plt.plot(np.linspace(0,len(npA), num = len(npA)), npA, c = "b")
            plt.plot(np.linspace(0,len(npB), num = len(npB)), npB, c = "g")

        #Sum of plasmid amounts per event point
        primerA += npA[:LENF]
        primerB += npB[:LENF]

        #Length of arrays npA and npB are the same
        sizea = npA.size

        #Addition of events of each round to the general array,
        # of all events in terms of plasmid amount, for both types
        AA = np.append(AA,npA)
        BB = np.append(BB,npB)

        #Done simulation round


    #------Calculation of average simulation line------
    primerA = primerA/float(rounds+1)
    primerB = primerB/float(rounds+1)
    meanA = primerA
    meanB = primerB

    #----Fit----
    #Limits for performing the fit
    #Top cut is LENF
    #print "Fitcut: " + str(LENF)
    #Initial limit is set so that 1 generation is allowed
    # (mainly done for competitions with equal k)
    # the system to globally stabilize
    cutIN = (inA)
    print cutIN,LENF

    #Fit is done in the average simulation line
    #Array with delimitation that it is important for the fit
    LA = meanA[(cutIN):(LENF)]
    LB = meanB[(cutIN):(LENF)]

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
    mAA, covAA = curve_fit(fun, np.linspace(cutIN, LENF, num = len(LA)), LA)

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
    mBB, covBB = curve_fit(fun, np.linspace(cutIN, LENF, num = len(LA)), LB)
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


    #----------COST_CALCULATION---------
    #Cost of type A
    FA = 2*mA*(inA+inB)
    #Cost of type B
    FB = 2*mB*(inA+inB)

    #Print in Log txt ***
    text_file = open("Results/Log_"+str(h)+".txt", "a+")
    n = text_file.write("Fitness/Cost A: "+str(FA)+"\n")
    text_file.close()
    #print "Fitness/Cost A: " +str(FA)


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
        xxxa = np.linspace(cutIN, LENF, num = len(LA))
        xxxb = np.linspace(cutIN, LENF, num = len(LB))

        #Plot fit of average line simulation
        plt.plot(xxxa, fun(xxxa, mA*np.ones(len(LA)) , cA*np.ones(len(LA))) , '--',c = "deeppink", label="Fit $K_{A}$", linewidth = 4)
        plt.plot(xxxb, fun(xxxb, mB*np.ones(len(LB)), cB*np.ones(len(LB))) ,'--', c = "y", label="Fit $K_{B}$", linewidth = 4)

        #Info of the general graph
        # info of cost for type A
        plt.plot(0,0, c = "white", label = "$u_{A}$ = " + str(round(mA,5)))
        plt.plot(0,0, c = "white", label = "$u_{B}$ = " + str(round(mB,5)))
        #plt.plot(0,0, c = "k", label = "Rounds = "+str(rounds+1))
        #plt.plot(0,0, c = "k", label = "Win " + str(kA) + " = " + str(GwinA) + "   Win " + str(kB) + " = " + str(GwinB))

        #Upper limit should match with Go LEN
        #No se muestra parte de estabilizacion
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
     #print kA,kB

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


full(3,6.0,22.0,33.0,1)

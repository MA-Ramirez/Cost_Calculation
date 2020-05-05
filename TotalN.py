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
    #return A,B,TA,TB

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
    #tUP = 1 - (1/float(inA+inB))
    tUP = 0.75
    tDw = 0.25
    #tDw = (1/float(inA+inB))
    #Lnl = 0.25
    #print tUP
    #print tDw

    #Limit of number of events for the simulation
    # this limit is called LEN
    LEN = 400

    #Markers to break loop
    C1 = True
    C2 = True

    #Marker to cut into it reaches SS
    CAVG = 0
    Round1 = True
    Birth1 = False

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
            #inA = A
            #inB = B
            #pA = np.append(pA,inA)
            #pB = np.append(pB,inB)

            #Normalization
            #cinA, cinB = nor(inA,inB)
            #npA = np.append(npA,cinA)
            #npB = np.append(npB,cinB)

            #--------------BIRTH----------------

            #There will be reproduction while reproduction probability
            # is above 2%
            while (probA >= 0.02 and probB >= 0.02) and S == 0:
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

                #If we are in 1st birth-death round
                if Round1 == True:
                    #If Birth phase occurs first than a Death phase
                    # then True just by going inside of Birth round
                    Birth1 = True


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

                #If we are in 1st birth-death round
                # and if birth didnt happen first
                if Round1 == True and Birth1 == False:
                    #Count how many death events happen of stabilization
                    CAVG +=1

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

            #End of first Birth-Death round
            Round1 = False

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
    #CAVG cut of fuchsia average of total amount of plasmid to avoid
    # initial death phase of stabilization
    return pA, pB, npA, npB, winA, winB, CAVG


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

    #Top limit for the average and the corresponding fit
    LENF = 250

    #Return of Go function (check Go comments)
    pA, pB, npA, npB, winA, winB, CAVG = Go(inA, inB, h, kA, kB)

    aCAVG = np.array([CAVG])

    #-----First round-----
    #Array of average per event point
    # initialized with normalized data of first round
    #primerA = npA[:LENF]
    #primerB = npB[:LENF]


    #pA: number of plasmids type A (no normalization)
    #pB: number of plasmids type B (no normalization)
    #TOTALG: total number of plamids - typeA + type B
    TOTALG = pA[:400]+pB[:400]

    #Graphs total number of plasmids per event
    plt.plot(np.linspace(0,len(TOTALG), num = len(TOTALG)), TOTALG, c  = "darkorange")


    #Graph for first round simulation
    # type A = blue   type B = green
    #Condition to graph only 1 general simulation
    if (rep == cc) == True:
        plt.plot(np.linspace(0,len(pA), num = len(pA)), pA, c = "b")
        plt.plot(np.linspace(0,len(pB), num = len(pB)), pB, c = "g")

    #Subsequent rounds
    for i in range(rounds):

        #Return of Go function (check Go comments)
        pA, pB, npA, npB, winA, winB, CAVG = Go(inA, inB, h, kA, kB)
        aCAVG = np.append(aCAVG,CAVG)

        #Sum of wins for each type of plasmid
        GwinA += winA
        GwinB += winB

        #Graph for simulations
        # type A = blue   type B = green
        #Condition to graph only 1 general simulation
        if (rep == cc) == True:
            plt.plot(np.linspace(0,len(pA), num = len(pA)), pA, c = "b")
            plt.plot(np.linspace(0,len(pB), num = len(pB)), pB, c = "g")

        #Sum of normalized plasmid amounts per event point
        #primerA += npA[:LENF]
        #primerB += npB[:LENF]

        #TOTALG: total number of plamids - typeA + type B
        TOTALG2 = pA[:400]+pB[:400]

        #Graphs total number of plasmids per event
        plt.plot(np.linspace(0,len(TOTALG2), num = len(TOTALG2)), TOTALG2, c = "darkorange")

        #Sum of total plasmids amounts per event point
        #TOTALG += TOTALG2

        if TOTALG2.size < TOTALG.size:
            c = np.copy(TOTALG)
            part = c[:TOTALG2.size]
            insert = np.mean([part, TOTALG2], axis=0)
            TOTALOK = insert
            TOTALOK = np.append(TOTALOK,c[TOTALG2.size:])
            TOTALG = TOTALOK
        else:
            c = np.copy(TOTALG2)
            part = c[:TOTALG.size]
            insert = np.mean([part, TOTALG], axis=0)
            TOTALOK = insert
            TOTALOK = np.append(TOTALOK,c[TOTALG.size:])
            TOTALG = TOTALOK


        #***Done simulation round***


    #------Calculation of average simulation line------
    #primerA = primerA/float(rounds+1)
    #primerB = primerB/float(rounds+1)
    #meanA = primerA
    #meanB = primerB

    #------Average total plasmids amount------
    #Average of total plasmids amount for all rounds
    #Graph of previous quantity
    plt.plot(np.linspace(0,len(TOTALG), num = len(TOTALG)), TOTALG, c = "orangered",linewidth = 3, label = "Average total amount of plasmids")



    CCAVG = round(np.mean(aCAVG))
    print CCAVG
    print "--------"
    #Average of Global Average
    #Uses cut found previously to not take into account stabilization to SS
    TTN = round(np.mean(TOTALG[CCAVG:]),1)
    plt.plot(np.linspace(CCAVG,len(TOTALG[CCAVG:]), num = len(TOTALG[CCAVG:])), np.ones(len(TOTALG[CCAVG:]))*TTN, c = "fuchsia",linewidth = 2.5, label = "Total average = "+str(TTN))

    #First column kA, second column kB, average of Global Average
    text_file = open("Total/"+str(h)+"_"+str(kA)+"Total.csv", "a+")
    n = text_file.write(str(kA)+","+str(kB)+","+str(TTN)+"\n")
    text_file.close()


    #-----------------GRAPHS---------------------
    #Graph of average simulation line and corresponding fit

    #Condition to graph only 1 general simulation
    if (rep == cc) == True:
        #Labels of the simulation
        # type A = blue   type B = green
        plt.plot(0,0, c = "b", label = "$K_{A}$ = " + str(kA))
        plt.plot(0,0, c = "g", label = "$K_{B}$ = " + str(kB))

        #General label of total amount of plasmids
        plt.plot(0,0, c="darkorange",label="Total amount of plasmids")

        #Upper limit should match with Go LEN
        plt.xlim((0,400))
        plt.title("Total Amount plasmids $K_{A}$ = "+str(kA)+" vs $K_{B}$ = "+str(kB)+" ($h$="+str(h)+")")
        plt.xlabel("Events")
        plt.ylabel("Amount of plasmids")
        plt.legend(loc = 4, fontsize = "x-small")
        plt.savefig("Total/Total_h"+str(h)+"/Total_" + str(kA) + "_" + str(kB) + "_" + str(int(h)) + ".png")
        plt.clf()


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
     text_file = open("Total/Log_"+str(h)+".txt", "a+")
     n = text_file.write("Competition: (kA,kB) "+str(kA)+", "+str(kB)+"\n")
     text_file.close()
     #print kA,kB

     #General rounds are rounds +1
     rounds = 299

     #Marker that helps to graph only one general simulation (optimization)
     cc = 0

     #Perform the general simulation rep times and find the average of their
     # results
     for i in range(rep):
        cc += 1
        repetitionHist(rounds,rep,cc,inI,inI,h,kA,kB)

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
    text_file = open("Total/Log_"+str(h)+".txt", "a+")
    n = text_file.write("RANGE OF CALCULATION: "+str(KKA)+"\n")
    text_file.close()

    #print "RANGE OF CALCULATION:"
    #print KKA

    #Generate Grid - KKA is x KKB is y
    X, Y = np.meshgrid(KKA, KKB)

    #Results of the cost of A = FA
    Fi = np.ones((len(KKA),len(KKB)))

    #Print in Log txt ***
    text_file = open("Total/Log_"+str(h)+".txt", "a+")
    n = text_file.write("LEN OF RANGE: (KKA,KKB) "+str(len(KKA))+","+str(len(KKB))+"\n \n")
    text_file.close()
    #print "LEN OF RANGE: "
    #print len(KKA), len(KKB)

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
        text_file = open("Total/Log_"+str(h)+".txt", "a+")
        n = text_file.write("STABLE1: "+ str(Stable1)+"   i = "+str(i)+"\n")
        text_file.close()

        #print "Stable1: " + str(Stable1)
        #print i

        #Counter for location
        yF = 0

        for j in KKB:
            #Stabilization value for kB
            inKB = np.where(kIndex==j)[0]
            Stable2 = Sta[inKB]

            #Print in Log txt ***
            text_file = open("Total/Log_"+str(h)+".txt", "a+")
            n = text_file.write("Stable2: "+ str(Stable2)+"   j = "+str(j)+"\n")
            text_file.close()

            #print "Stable2: " + str(Stable2)
            #print j

            #Average of stabilization
            print i,Stable1
            print j,Stable2
            AvgStable = (Stable1 + Stable2)/2.0
            #Turns it into int, no decimal amount of plasmids
            #int cuz it rounds down .5, more accurate for low amount of P
            inI = int(AvgStable/2.0)

            #Print in Log txt ***
            text_file = open("Total/Log_"+str(h)+".txt", "a+")
            n = text_file.write("Average Stable: "+ str(AvgStable)+"   InI: "+str(inI)+"\n")
            text_file.close()

            #print "Average Stable: " + str(AvgStable)
            #print "InI: " + str(inI)

            full(h,i,j,inI,rep)

            #Counter for location
            yF+=1

        #Print in Log txt ***
        text_file = open("Total/Log_"+str(h)+".txt", "a+")
        n = text_file.write("\n")
        text_file.close()
        #Counter for location
        xF+=1

    #Time of all the total simulation
    tsim = round(time.time()-t0,3)

    #Print in Log txt ***
    text_file = open("Total/Log_"+str(h)+".txt", "a+")
    n = text_file.write("Total Simu Time "+ str(tsim))
    text_file.close()

#contourG(start,stop,hop,h,rep)
contourG(2,22,6,3,1)

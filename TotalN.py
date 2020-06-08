import numpy as np
import matplotlib.pyplot as plt
import time

#Code by Maria Alejandra Ramirez
#Script: TotalN.py

#------SINGLE COMPETITION (KA vs KB)-------
#Non-relative population growth to calculate NT
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

#--------------------------------------
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

    #"else" none (implicit)

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
#Normalizes the amount of plasmids to 1
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
            # is above 2%, for either one type
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
            #Code will run while the length of the array is less
            # than the LEN limit
            if len(npA) >= LEN:
                S = 1
                break

    #Results of the whole process
    #pA, pB arrays with amount of plasmids for each event (type A and B)
    #npA, npB arrays of pA and pB normalized to 1
    return pA, pB, npA, npB


#--------REPETITION OF THE PROCESS & GRAPHS------------
#Repetition of the process and corresponding graphs
#Param: rounds times that the process is repeated +1
#Param: rep repetitions of the same general simulation
#Param: cc iteration of the current repetition
#Param: inA initial amount of plasmids type A
#Param: inB initial amount of plasmids type B
#Param: h Hill coefficient
#Param: kA K of plasmid type A
#Param: kB K of plasmid type B
#rep and cc are used to only graph one general simulation (optimization)
#Return: graphs of the simualtion
#Return: CVS file for the corresponding quantity NT of each competition
def repetitionHist(rounds, rep, cc, inA, inB, h, kA, kB):

    #Return of Go function (check Go comments)
    pA, pB, npA, npB = Go(inA, inB, h, kA, kB)

    #-----First round-----

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
        pA, pB, npA, npB = Go(inA, inB, h, kA, kB)

        #Graph for simulations
        # type A = blue   type B = green
        #Condition to graph only 1 general simulation
        if (rep == cc) == True:
            plt.plot(np.linspace(0,len(pA), num = len(pA)), pA, c = "b")
            plt.plot(np.linspace(0,len(pB), num = len(pB)), pB, c = "g")

        #TOTALG: total number of plamids - typeA + type B
        TOTALG2 = pA[:400]+pB[:400]

        #Graphs total number of plasmids per event
        plt.plot(np.linspace(0,len(TOTALG2), num = len(TOTALG2)), TOTALG2, c = "darkorange")

        #Average of total N
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

    #------Average total plasmids amount------
    #Average of total plasmids amount for all rounds
    #Graph of previous quantity
    plt.plot(np.linspace(0,len(TOTALG), num = len(TOTALG)), TOTALG, c = "orangered",linewidth = 3, label = "Average total amount of plasmids")

    #Average of Global Average
    #Uses cut found previously to not take into account stabilization to SS
    TTN = round(np.mean(TOTALG),1)
    plt.plot(np.linspace(0,len(TOTALG), num = len(TOTALG)), np.ones(len(TOTALG))*TTN, c = "fuchsia",linewidth = 2.5, label = "Total average = "+str(TTN))

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


#Generates repetitions of the whole simulation to find an average of the results
# for the same general simulation.
#Param: h Hill coefficient
#Param: kA K of plasmid type A
#Param: kB K of plasmid type B
#Param: inI Initial amount of plasmids (same for type A and B)
#Param: rep Repetitions of the whole simulation to find the average
def full(h,kA,kB,inI,rep):

     #Print in Log txt ***
     text_file = open("Total/Log_"+str(h)+".txt", "a+")
     n = text_file.write("Competition: (kA,kB) "+str(kA)+", "+str(kB)+"\n")
     text_file.close()

     #General rounds are rounds +1
     rounds = 499

     #Marker that helps to graph only one general simulation (optimization)
     cc = 0

     #Perform the general simulation rep times and find the average of their
     # results
     for i in range(rep):
        cc += 1
        repetitionHist(rounds,rep,cc,inI,inI,h,kA,kB)

#----------GENERATES COUPLES FOR COMPETITIONS---------------------------

#Creates the one-on-one competitions
#The quantity NT for each competition can be the average of rep repetitions
# for the whole simulation of the given competition
#Param: start Initial K
#Param: stop Final K
#Param: hop Amount of values between start and stop
#Param: Hill coefficient
#Param: rep Repetitions for the whole simulation of a given competition
def contourG(start,stop,hop,h,rep):

    #Start measuring simulation time
    t0 = time.time()

    #Arrays for kA y kB that are gonna be used (nxn)
    # Example: np.linspace(1,3,3) = 1,2,3
    # Example: np.linspace(2,4,3) = 2,3,4
    KKA = np.linspace(start,stop,hop)
    KKB = np.linspace(start,stop,hop)

    #Print in Log txt ***
    text_file = open("Total/Log_"+str(h)+".txt", "a+")
    n = text_file.write("RANGE OF CALCULATION: "+str(KKA)+"\n")
    text_file.close()

    #Print in Log txt ***
    text_file = open("Total/Log_"+str(h)+".txt", "a+")
    n = text_file.write("LEN OF RANGE: (KKA,KKB) "+str(len(KKA))+","+str(len(KKB))+"\n \n")
    text_file.close()

    #Data of stabilization is calculated in another script
    # this data is recorded on csv files for optimization

    #Note: to perform a simulation, the Stabilization data for the Ks involved
    # should be available **********

    Sta = np.genfromtxt("Stabilization/"+str(h)+"Stable.csv", delimiter=",",usecols=1)
    kIndex = np.genfromtxt("Stabilization/"+str(h)+"Stable.csv", delimiter=",",usecols=0)

    #Initial inA y inB (same) is the average of both stabilizations, half half
    #Simulation  ofevery competition

    for i in KKA:
        #Stabilization value for kA
        inKA = np.where(kIndex==i)[0]
        Stable1 = Sta[inKA]

        #Print in Log txt ***
        text_file = open("Total/Log_"+str(h)+".txt", "a+")
        n = text_file.write("STABLE1: "+ str(Stable1)+"   i = "+str(i)+"\n")
        text_file.close()

        for j in KKB:
            #Stabilization value for kB
            inKB = np.where(kIndex==j)[0]
            Stable2 = Sta[inKB]

            #Print in Log txt ***
            text_file = open("Total/Log_"+str(h)+".txt", "a+")
            n = text_file.write("Stable2: "+ str(Stable2)+"   j = "+str(j)+"\n")
            text_file.close()

            #Average of stabilization
            AvgStable = (Stable1 + Stable2)/2.0
            #Turns it into int, no decimal amount of plasmids
            #int cuz it rounds down .5, more accurate for low amount of P
            inI = int(AvgStable/2.0)

            #Print in Log txt ***
            text_file = open("Total/Log_"+str(h)+".txt", "a+")
            n = text_file.write("Average Stable: "+ str(AvgStable)+"   InI: "+str(inI)+"\n")
            text_file.close()

            #Executes shortcut for Simulation
            #if (i-j) != 0:
                #full(h,i,j,inI,rep)

            #Executes Simulation
            full(h,i,j,inI,rep)


        #Print in Log txt ***
        text_file = open("Total/Log_"+str(h)+".txt", "a+")
        n = text_file.write("\n")
        text_file.close()

    #Time of all the total simulation
    tsim = round(time.time()-t0,3)

    #Print in Log txt ***
    text_file = open("Total/Log_"+str(h)+".txt", "a+")
    n = text_file.write("Total Simu Time "+ str(tsim)+"\n")
    text_file.close()

#Example
#contourG(start,stop,hop,h,rep)
contourG(41,45,2,4,1)

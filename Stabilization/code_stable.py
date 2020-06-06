import numpy as np
import matplotlib.pyplot as plt
import time

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

#BIRTH 1 TYPE
#Generates a decision of reproduction A,B or None
#Param: h Hill coefficient
#Param: kA of the plasmid
#Param: A amount of plasmids
#Return: amount of A, B after the reproduction decision
def Birth(h,kA,A):

    #Total amount of plasmids
    tot = float(A)

    #Type A (kA)
    #nor: Bmax to normalize y, should fit with random.random()
    norA, yA = repression(h,kA,tot)

    #Probability of reproduction normalized
    probA = float(yA/norA)

    #Threshold A
    TA = (probA*A) / tot

    #Random number between 0-1
    Q = np.random.random()

    #Birth A
    if Q <= probA:
        A+= 1

    return A,probA

#Death
#Generates a decision of death A,B
#Param: A amount of plasmids type A
#Return: amount of A after the decision of death
def Death(A):

    #Normalized
    #Threshold A
    TDA = (0.5*A) / (1*A)

    #random number between 0-1
    QQ = np.random.random()

    #Could be optimize, since there is only one type of plasmids
    #Death A
    if QQ <= TDA:
        A = A-1

    return A

#Process
#All the process is done together
#Reproduction-Death cycle
#Param: inA initial amount of plasmids type A
#Param: h Hill coefficient
#Param: kA of plasmid type A
#Return: pA array with data of each event of the plasmids type A
def Go(inA, h, kA):

    #Arrays with data of each event of the plasmids type A
    pA = np.array([inA])

    #Amount of events
    rounds = 90000
    #It stops with a determined amount of events
    while len(pA) < rounds:

        #Primer Birth
        A,probA = Birth(h,kA,inA)

        #While the probability of reproduction is >= 2% it is performed:
        while probA >= 0.02:

            A,probA = Birth(h,kA,inA)
            inA = A
            pA = np.append(pA,inA)

            #break to avoid an infinite loop, based on number of events
            if len(pA) >= rounds:
                break

        #Total amount of plasmids after the reproductive phase = fixed
        fixed = inA

        #Random death until half of the fixed amount is reached
        while inA > (fixed/2.0):
            A = Death(inA)
            inA = A
            pA = np.append(pA,inA)

            #break to avoid an infinite loop, based on number of events
            if len(pA) >= rounds:
                break

    return pA

#Execute
#It executes the process, graphs and reports the average stabilization number
#Param: inA initial amount of plasmids type A
#Param: h Hill coefficient
#Param: kA of plasmid type A
def exe(inA, h, kA):
    #Start measuring the time
    t0 = time.time()

    #Array with amount of plasmids per event
    pA = Go(inA, h, kA)

    #Time of simulation
    tsim = round(time.time()-t0,3)

    #average stabilization number
    avgpA = round(np.mean(pA))

    #Prints h, k and average of stabilization
    print h,kA,avgpA

    #Records in .csv the data
    #First column K, second column average stabilization number
    text_file = open(str(h)+"Stable.csv", "a+")
    n = text_file.write(str(kA)+","+str(avgpA)+"\n")
    text_file.close()

    #Graph
    plt.plot(np.linspace(0,len(pA), num = len(pA)), pA, label = "$K$ = " + str(kA) +", $h$ = " + str(h))
    #Legends
    plt.plot(np.linspace(0,len(pA), num = len(pA)),np.ones(len(pA))*avgpA, c = "c", label = "$\overline{n}$ = " + str(avgpA))
    plt.plot(0,0, c = "y", label = "Time Simu= " + str(tsim) +"s")
    plt.legend(loc = 4, fontsize = "small")
    plt.title("Average Stabilization Number ($\overline{n}$)")
    plt.xlabel("Events")
    plt.ylabel("Number of plasmids")
    plt.savefig("h"+str(h)+"/graph_" + str(kA)+ "_" + str(h) + ".png")
    plt.clf()

#Example range of K
#np.linspace(start,stop,hop)
SS = np.linspace(3,50,48)

#Execute code for same h for a given range of k
for i in SS:
    Always start with 1 plasmids
    #exe(inA,h,kA)
    exe(1,2.7,i)

#Example 1 single K
exe(1,5,70)

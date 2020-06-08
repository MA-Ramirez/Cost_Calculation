# Cost_Calculation
- Code and results to calculate the cost of an altruistic act for the system of plasmids in bacteria

###### The code for this simulation consists on 4 main parts.

1. **Stabilization**

On the folder *Stabilization*, the script *code_stable.py* can be found. This script runs the simulation for the evolution of a single plasmid
type, in order to obtain its steady-state. The results are recorded on CSV files and are graphed, when indicated, in folders like *h3*.

The script *fit_KyN.py* can be also found in this folder. This script graphs the steady-states in function of the repression coefficient
(K) for each Hill coefficient (h) considered. Afterwards, it performs a fit to obtain the relationship between the steady-state and K. 
These results are located in the folder *StableGraphs*. Similarly, the relationship between the parameters of the previous fit and 
the Hill coefficient was obtained via *fit_KyN.py* and its results can be also found in the folder *StableGraphs*.

In particular, the steady-state is the main quantity that determines the initial amount of plasmids in each competition between plasmid types.

2. **Total**

The quantity NT used to calculate the cost of an altruistic act is obtained with the script *TotalN.py*. The main goal of this script is to
run one-on-one competitions of different plasmid types to record the general total population size (NT). The corresponding results are recorded
in the folder *Total* via CSV files and graphs, such graphs are organized in different folders according to the Hill coefficient of the plasmid types.

3. **Simulation**

The whole simulation of the one-one-competitions of different plasmid types is performed by the script *simu.py*. In this case, the results
of the simulation are relative population growth graphs from which the cost of an altruistic act can be graphically calculated. The demonstration
of the mathematical expression used to find the cost is explained on the document related to this bachelor's thesis.

In general, the graphical results of the population growth graphs are located in the folder *SimulationGraphs* and are organized according to the Hill 
coefficient of the plasmid types.

Furthermore, in the script *simu.py* the calculation of the cost is performed. The results of the cost can be found in the folder *Results*. 
That folder contains the following results: filled contour plots representing the cost related to each one-on-one competition, txt files with the 
cost results, filled contour plots illustrating the variance associated with the linear fits performed on the previous graphs and txt files with such 
information.

4. **Fit**

The fit performed on the cost results was done through the script *fit_contour.py*. Therefore, this script takes the results from the previous subsection
and carries out the corresponding fit to find the relationship between the cost, KA and KB. The results of the fits are located in the folder
*Results/Fit*. Particularly, the visual representations of the fits are saved as PNG files and the results of the fit parameters are recorded
in a CSV files, along with the variance related to those fits (in separate files).

Finally, through the script *fit_contour.py* the relationship between the fit parameters previously found and the Hill coefficient was determined.
Such relation is documented in the file *Cte_Fit.png*.

Note: the files *Log_#.txt* are used to record information about the progress of the simulations.

- Overall, to run a competition between plasmid types the scripts should be run in the following order to calculate the cost:
1. *code_stable.py*
2. *TotalN.py*
3. *simu.py*

This order should be followed, because *code_stable.py* establishes the initial amount of plasmids in the competitions by calculating the steady-states.
Subsequently, *TotalN.py* computes NT, a quantity that is necessary to determine the cost via the calculation perfomed in *simu.py*.

It should be noted that the competitions take place between pairs of plasmid types, which means that this fact should be taken into account while 
running the code. 

The scripts *fit_KyN.py* and *fit_contour.py* are not needed to run a competition, since they are mainly used to establish mathematical relationships among variables.

>For further information about this research, please contact:
###### Author: Maria Alejandra Ramirez - B.S. in Physics (ma.ramirez@uniandes.edu.co)      

###### Supervisor: Juan Manuel Pedraza - Professor at Universidad de los Andes

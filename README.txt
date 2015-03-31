MyAMI Specific Ion Interaction Model (Version 1.0):

MyAMI is Python code to calculate thermodynamic pK's and conditional pK's

Author: Mathis P. Hain -- m.p.hain@soton.ac.uk

Reference:
Hain, M.P., Sigman, D.M., Higgins, J.A., and Haug, G.H. (2015) The effects of secular calcium and magnesium concentration changes on the thermodynamics of seawater acid/base chemistry: Implications for Eocene and Cretaceous ocean carbon chemistry and buffering, Global Biogeochemical Cycles, 29, doi:10.1002/2014GB004986

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

HOW TO:
Calculate the parameters that define temperature and salinity dependences of the conditional equilibrium constants ? at arbitrary [Ca2+] and [Mg2+].

(1) open command line
(2) navigate to folder (cd ~/path/to/folder/with/MyAMI/code)
(3) run MyAMI (python MyAMI_V1.py)
(4) follow instructions
(5) output is saved to file [OUTPUT_INTERACTIVE.txt]

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

HOW TO:
Call Python code from MATLAB to calculate conditional and thermodynamic equilibrium constants at a given [Ca2+], [Mg2+], temperature and salinity.

(1) open Matlab
(2) navigate to folder that contains [PyMyAMI.m]
(3) >>PITZERpath = '~/path/to/folder/of/MyAMI/code/Pitzer.py';
(4) >>T = 20;
(5) >>S = 35;
(6) >>Ca = 0.0102821 * S/35;
(7) >>Mg = 0.0528171* S/35;
(8) >>B = xB*S / 1.80655 * 0.0000219;   
(9) >>[K] = PyMyAMI(PITZERpath,num2str(T),num2str(S),num2str(Ca),num2str(Mg));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

REQUIREMENTS: (1) you need Python (2.7.5 verified); (2) you need the standard Python modules ?NumPy? and ?SciPy?; (3) have fun.

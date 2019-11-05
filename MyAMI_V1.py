# MyAMI Specific Ion Interaction Model (Version 1.0):
# This is a Python script to calculate thermodynamic pK's and conditional pK's
# Author: Mathis P. Hain -- m.p.hain@soton.ac.uk
#
# Reference:
# Hain, M.P., Sigman, D.M., Higgins, J.A., and Haug, G.H. (2015) The effects of secular calcium and magnesium concentration changes on the thermodynamics of seawater acid/base chemistry: Implications for Eocene and Cretaceous ocean carbon chemistry and buffering, Global Biogeochemical Cycles, 29, doi:10.1002/2014GB004986
#
# For general context on the calculations see Millero, 2007 (Chemical Reviews) and Millero and Pierrot, 1998 (Aquatic Geochemistry)

# NOTE: in Python log is natural-logarithm & log10 is log-base-10

# header information needed to run
# [[ REQUIRED: NUMPY, SCIPY, MATPLOTLIB ]]
import numpy
from scipy.optimize import curve_fit
from multiprocessing import Pool
import sys
import pKs
import K_thermo_conditional

## for plotting
import matplotlib
import matplotlib.pyplot as plt

# Determine if code is run under INTERACTIVE USER or NONINTERACTIVE MODEL mode
INTERACTIVE = 1
if (len(sys.argv)>1):
    try:
        XmCa = float(sys.argv[1])
        INTERACTIVE = 0
    except Exception:
        INTERACTIVE = 1
        print("Error: The first argument ([Ca2+]) provided is not a number.")

if (len(sys.argv)>2) and (INTERACTIVE == 0) :
    try:
        XmMg = float(sys.argv[2])
    except Exception:
        INTERACTIVE = 1
        print("Error: The first argument ([Mg2+]) provided is not a number.")

if (INTERACTIVE == 0):
    if (XmCa>0.1) or (XmCa<0) or (XmMg>0.1) or (XmMg<0):
        INTERACTIVE = -1

# Some documentation printed to screen when code is being run
if (INTERACTIVE == 1) or (INTERACTIVE == -1):
    print()
    print()
    print
    print("--------------------------------------------------------")
    print("--- MyAMI specific ion interaction model Version 1.0 ---")
    print("--- Author: Dr. Mathis P. Hain (m.p.hain@soton.ac.uk) --")
    print("--------------------------------------------------------")
    print()
    print("PURPOSE:")
    print("This code calculates conditional equilibrium constants in seawater under variable [Ca2+] and [Mg2+] concentrations. The equilibrium constants are presented in functional form including their temperature and salinity dependence.")
    print()
    print("BASIS:")
    print("The 'MyAMI' code used here is based on the 'MIAMI' model, as described in detail by Millero and Pierrot (Aquatic Geochemistry 4: 153-199, 1998). However, some of the PitZre parameters used here differ from those tabulated by these authors.")
    print()
    print("NOTE:")
    print("(#1) For modern seawater composition (e.g., modern [Ca2+] and [Mg2+] the model output is consistent with the modern conditional equilibrium constants summarized by Dickson, Sabine & Christian (Guide to best practises for ocean CO2 measurements, PICES Special Publication, 2007). This is by design because the Pitzer-type model is used to estimate the sensitivity of the equilibrium constants to [Ca2+] and [Mg2+] change, which is then applied to the modern reference equilibrium constants.")
    print("(#2): The code assumes a simplified seawater composition, with Na, K, H, Ca, Mg, Sr, Cl, OH, HCO3, CO3, HSO4, SO4, B(OH)4, CO2* and B(OH)3. Concentrations of these species are according to 'standard seawater' of Millero et al. (Deep-Sea Research I: 55: 50-72, 2008)")
    print("(#3): All equilibrium constants are on the 'total pH scale' with pH = -log10(H+HSO4).")
    print("(#4): Salinity is taken to be a master variable, determining (1) ionic strength according to I = 19.924*Sal/(1000-1.005*Sal); and (2) determining the dilution/concentration of the dissolved species according to Conc(Sal=X) = Conc(Sal=35) * X / 35.")
    print()
    print("USAGE:")
    print("(#1): If you see this text you are running the code in INTERACTIVE USER mode (e.g., python MyAMI_V1.py). The code will (first) ask you to input [Ca2+] and [Mg2+] manually and (second) output the Temp&Sal dependent formulas of the equilibrium constants in human readable form on screen.")
    print("(#2): If you call the code including [Ca2+] and [Mg2+] (in mol/kg) as command line arguments (e.g., python MyAMI_V1.py 0.0102821 0.0528171) the code will use these values and output only the coefficients for the Temp&Sal dependent formulas of the equilibrium constants. This NONINTERACTIVE MODEL mode is useful if the coefficients are inteded to be used directly by some other computer code (e.g. carbon cycle model).")
    print()
    print()
    print()

# Experiment (X) concentration (m) of Ca and Mg
if (INTERACTIVE == 1) or (INTERACTIVE == -1):
    print()
    print("------------------INTERACTIVE_MODE----------------------")
    print("--- MyAMI specific ion interaction model Version 0.9 ---")
    print("------------------INTERACTIVE_MODE----------------------")
    print()

if (INTERACTIVE == -1):
    print("Your choice of [Ca2+] and/or [Mg2+] falls outside the reasonabel range of 0 to 100 mM. Please reconsider.")
    print()

if (INTERACTIVE == 1) or (INTERACTIVE == -1):
    XmCa = input("Enter [Ca2+] in mol/kg for nominal Sal=35 (modern 0.0102821): ")
    XmMg = input("Enter [Mg2+] in mol/kg for nominal Sal=35 (modern 0.0528171): ")
    XmCa = float(XmCa)
    XmMg = float(XmMg)
    print()
    print("CODE IS RUNNING ... be patient. Sometimes it takes a while to fit the Temp/Sal function to the Model output.")

# Modern (M) concentration (m) of Ca and Mg case T=25C, I=0.7, seawatercomposition
MmMg = 0.0528171 # Mg Millero et al., 2008; Dickson OA-guide
MmCa = 0.0102821 # Ca Millero et al., 2008; Dickson OA-guide

# number of Temp and Sal steps used as the basis dataset for the fitting of the pK's
n =21 # number Temp and Sal levels
N = n*n # number of distinct combinations of Temp/Sal

# create list of Temp's and Sal's defining the grid for fitting pK's
TempC = numpy.linspace(0,40,n) # 0-40degC in N steps
Sal = numpy.linspace(30,40,n) # 30-40 Sal
TempC_M,Sal_M = numpy.meshgrid(TempC,Sal) # generate grid in matrix form
TempC = numpy.ravel(TempC_M) # "unroll" matrix into vectors
Sal = numpy.ravel(Sal_M) # "unroll" matrix into vectors
TempK = TempC + 273.15
TempK_M = TempC_M + 273.15
lnTempK = numpy.log(TempK)
lnTempK_M = numpy.log(TempK_M)

# Allocate memory for modern conditional K's
KspC_mod = numpy.zeros(N)
K1_mod = numpy.zeros(N)
K2_mod = numpy.zeros(N)
Kw_mod = numpy.zeros(N)
Kb_mod = numpy.zeros(N)
KspA_mod = numpy.zeros(N)
K0_mod = numpy.zeros(N)
KSO4_mod = numpy.zeros(N)

# Calculate K's for modern seawater composition
for idx in range(0,N):
    KspC_mod[idx], K1_mod[idx], K2_mod[idx], Kw_mod[idx], Kb_mod[idx], KspA_mod[idx], K0_mod[idx], KSO4_mod[idx] = K_thermo_conditional.CalculateKcond(TempC[idx],Sal[idx])

# Allocate memory for the aggregated gammas that go with the various K's (gK's)
gKspC_mod = numpy.zeros(N)
gK1_mod = numpy.zeros(N)
gK2_mod = numpy.zeros(N)
gKw_mod = numpy.zeros(N)
gKb_mod = numpy.zeros(N)
gKspA_mod = numpy.zeros(N)
gK0_mod = numpy.zeros(N)
gKSO4_mod = numpy.zeros(N)
############################
gKspC_X = numpy.zeros(N)
gK1_X = numpy.zeros(N)
gK2_X = numpy.zeros(N)
gKw_X = numpy.zeros(N)
gKb_X = numpy.zeros(N)
gKspA_X = numpy.zeros(N)
gK0_X = numpy.zeros(N)
gKSO4_X = numpy.zeros(N)

# Calculate gK's for modern (mod) and experimental (x) seawater composition
for idx in range(0,N):
    [gKspC_mod[idx], gK1_mod[idx], gK2_mod[idx], gKw_mod[idx], gKb_mod[idx], gKspA_mod[idx], gK0_mod[idx], gKSO4_mod[idx]] = pKs.calculate_gKs(TempC[idx],Sal[idx], MmCa, MmMg)

for idx in range(0,N):
    [gKspC_X[idx], gK1_X[idx], gK2_X[idx], gKw_X[idx], gKb_X[idx], gKspA_X[idx], gK0_X[idx], gKSO4_X[idx]] = pKs.calculate_gKs(TempC[idx],Sal[idx], XmCa, XmMg)

# Calculate conditional K's predicted for seawater composition X
KspC_X = KspC_mod * gKspC_X / gKspC_mod
K1_X = K1_mod * gK1_X / gK1_mod
K2_X = K2_mod * gK2_X / gK2_mod
Kw_X = Kw_mod * gKw_X / gKw_mod
Kb_X = Kb_mod * gKb_X / gKb_mod
KspA_X = KspA_mod * gKspA_X / gKspA_mod
K0_X = K0_mod * gK0_X / gK0_mod
KSO4_X = KSO4_mod * gKSO4_X / gKSO4_mod


# ABOVE: the PITZER model was applied to predict N-times (across Temp&Sal) conditional constants applicable to some arbitrary seawater composition X
# BELOW: the above dataset is fit to the identical functional form of the constants for modern SW by using a non-linear least-square method
maxfevN=2000000  # number of optimiztion timesteps allowed to reach convergence


# functional form used for K0
guess_K0 = (1, 1, 1, 1, 1, 1) # initial guess for fitting algorithm taken to be unity when modern params (i.e., normalized for nummeric reasons)
def func_K0(TkS, a, b, c, d, e, f):
    TempK = TkS[0]
    Sal = TkS[1]
    return numpy.exp( (-60.2409*a)+(93.4517*b)*100/TempK+(23.3585*c)*numpy.log(TempK/100) + Sal*((0.023517*d)+(-0.023656*e)*TempK/100+(0.0047036*f)*TempK/100*TempK/100) )
# fitting the params for K0
popt, pcov = curve_fit(func_K0,(TempK,Sal), K0_X, guess_K0,maxfev=maxfevN)
params_K0 = [[-60.2409*popt[0], 93.4517*popt[1], 23.3585*popt[2], 0.023517*popt[3], -0.023656*popt[4],0.0047036*popt[5]],[popt[0], popt[1], popt[2], popt[3], popt[4], popt[5]],[-60.2409, 93.4517, 23.3585, 0.023517, -0.023656, 0.0047036]]
test = func_K0((TempK,Sal),params_K0[1][0],params_K0[1][1],params_K0[1][2],params_K0[1][3],params_K0[1][4],params_K0[1][5])
if (INTERACTIVE != 0):
    print()
    print("lnK0 = {0:.4} (+) {1:.4}(100/T) (+) {2:.4}*ln(T/100) (+) Sal*[ {3:.6} (+) {4:.6}*(T/100) (+) {5:.7}*(T/100)*(T/100)]".format(params_K0[0][0],params_K0[0][1],params_K0[0][2],params_K0[0][3],params_K0[0][4],params_K0[0][5]))
    print("Model fitted parameters at user defined [Ca2+] and [Mg2+] relative to parameters fit to experimental data for modern seawater:", popt)
    print("Average accuracy of fitting the function to the model output:", numpy.sum(numpy.sqrt((test/K0_X)*(test/K0_X)))/N)
    print()


# functional form used for K1
guess_K1 = (1, 1, 1, 1, 1) # initial guess for fitting algorithm taken to be unity when modern params (i.e., normalized for nummeric reasons)
def func_K1(TkS, a, b, c, d, e):
    TempK = TkS[0]
    Sal = TkS[1]
    return numpy.power(10,((61.2172*a) + (-3633.86*b)/TempK + (-9.6777*c)*numpy.log(TempK) + (0.011555*d)*Sal + (-0.0001152*e)*Sal*Sal))
# fitting the params for K1
popt, pcov = curve_fit(func_K1,(TempK,Sal), K1_X, guess_K1,maxfev=maxfevN)
params_K1 = [[61.2172*popt[0], -3633.86*popt[1], -9.6777*popt[2], 0.011555*popt[3], -0.0001152*popt[4]],[popt[0], popt[1], popt[2], popt[3], popt[4]],[61.2172, -3633.86, -9.6777, 0.011555, -0.0001152]]
test = func_K1((TempK,Sal),params_K1[1][0],params_K1[1][1],params_K1[1][2],params_K1[1][3],params_K1[1][4])
if (INTERACTIVE != 0):
    print("log10K1 = {0:.4f} (+) {1:.2f}/T (+) {2:.4f}*lnT (+) {3:.6f}*Sal (+) {4:.7f}*Sal*Sal".format(params_K1[0][0],params_K1[0][1],params_K1[0][2],params_K1[0][3],params_K1[0][4]))
    print("Model fitted parameters at user defined [Ca2+] and [Mg2+] relative to parameters fit to experimental data for modern seawater:", popt)
    print("Average accuracy of fitting the function to the model output:", numpy.sum(numpy.sqrt((test/K1_X)*(test/K1_X)))/N)
    print()


# functional form used for K2
guess_K2 = (1, 1, 1, 1, 1) # initial guess for fitting algorithm taken to be unity when modern params (i.e., normalized for nummeric reasons)
def func_K2(TkS, a, b, c, d, e):
    TempK = TkS[0]
    Sal = TkS[1]
    return numpy.power(10,((9) + (-25.9290*a) + (-471.78*b)/TempK + (3.16967*c)*numpy.log(TempK) + (0.01781*d)*Sal + (-0.0001122*e)*Sal*Sal))
# fitting the params for K2
popt, pcov = curve_fit(func_K2,(TempK,Sal), K2_X*1e9, guess_K2,maxfev=maxfevN)
params_K2 = [[-25.9290*popt[0], -471.78*popt[1], 3.16967*popt[2], 0.01781*popt[3], -0.0001122*popt[4]],[popt[0], popt[1], popt[2], popt[3], popt[4]],[-25.9290, -471.78, 3.16967, 0.01781, -0.0001122]]
test = func_K2((TempK,Sal),params_K2[1][0],params_K2[1][1],params_K2[1][2],params_K2[1][3],params_K2[1][4])
if (INTERACTIVE != 0):
    print("log10K2 = {0:.4f} (+) {1:.2f}/T (+) {2:.5f}*lnT (+) {3:.5f}*Sal (+) {4:.7f}*Sal*Sal".format(params_K2[0][0],params_K2[0][1],params_K2[0][2],params_K2[0][3],params_K2[0][4]))
    print("Model fitted parameters at user defined [Ca2+] and [Mg2+] relative to parameters fit to experimental data for modern seawater:", popt)
    print("Average accuracy of fitting the function to the model output:", numpy.sum(numpy.sqrt((test/(K2_X*1e9))*(test/(K2_X*1e9))))/N)
    print()


# functional form used for Kb
guess_Kb = (1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1) # initial guess for fitting algorithm taken to be unity when modern params (i.e., normalized for nummeric reasons)
def func_Kb(TkS, a, b, c, d, e, f, g, h, i, j, k, l):
    TempK = TkS[0]
    Sal = TkS[1]
    sqrtSal = numpy.sqrt(Sal)
    lnTempK = numpy.log(TempK)
    return numpy.exp( (148.0248*a) + (137.1942*b)*sqrtSal + (1.62142*c)*Sal + (1/TempK)*((-8966.90*d)+(-2890.53*e)*sqrtSal+(-77.942*f)*Sal+(1.728*g)*Sal*sqrtSal+(-0.0996*h)*Sal*Sal) + lnTempK*((-24.4344*i)+(-25.085*j)*sqrtSal+(-0.2474*k)*Sal) + (0.053105*l)*sqrtSal*TempK )
# fitting the params for Kb
popt, pcov = curve_fit(func_Kb,(TempK,Sal), Kb_X, guess_Kb,maxfev=maxfevN)
params_Kb = [[148.0248*popt[0], 137.1942*popt[1], 1.62142*popt[2], -8966.90*popt[3], -2890.53*popt[4], -77.942*popt[5], 1.728*popt[6], -0.0996*popt[7], -24.4344*popt[8], -25.085*popt[9], -0.2474*popt[10], 0.053105*popt[11]],[popt[0], popt[1], popt[2], popt[3], popt[4], popt[5], popt[6], popt[7], popt[8], popt[9], popt[10], popt[11]],[148.0248, 137.1942, 1.62142, -8966.90, -2890.53, -77.942, 1.728, -0.0996, -24.4344, -25.085, -0.2474, 0.053105]]
test = func_Kb((TempK,Sal),params_Kb[1][0],params_Kb[1][1],params_Kb[1][2],params_Kb[1][3],params_Kb[1][4],params_Kb[1][5],params_Kb[1][6],params_Kb[1][7],params_Kb[1][8],params_Kb[1][9],params_Kb[1][10],params_Kb[1][11])
if (INTERACTIVE != 0):
    print("lnKb = {0:.4f} (+) {1:.4f}*Sal^0.5 (+) {2:.5f}*Sal (+) (1/T*[{3:.2f} (+) {4:.2f}*sqrtSal (+) {5:.3f}*Sal (+) {6:.3f}*Sal^1.5 (+) {7:.4f}*Sal^2] (+) lnT*[{8:.4f} (+) {9:.3f}*Sal^0.5 (+) {10:.4f}*Sal] (+) {11:.6f}*T*Sal^0.5".format(params_Kb[0][0],params_Kb[0][1],params_Kb[0][2],params_Kb[0][3],params_Kb[0][4],params_Kb[0][5],params_Kb[0][6],params_Kb[0][7],params_Kb[0][8],params_Kb[0][9],params_Kb[0][10],params_Kb[0][11]))
    print("Model fitted parameters at user defined [Ca2+] and [Mg2+] relative to parameters fit to experimental data for modern seawater:", popt)
    print("Average accuracy of fitting the function to the model output:", numpy.sum(numpy.sqrt((test/Kb_X)*(test/Kb_X)))/N)
    print()


# functional form used for Kw
guess_Kw = (1, 1, 1, 1, 1, 1, 1) # initial guess for fitting algorithm taken to be unity when modern params (i.e., normalized for nummeric reasons)
def func_Kw(TkS, a, b, c, d, e, f, g):
    TempK = TkS[0]
    Sal = TkS[1]
    sqrtSal = numpy.sqrt(Sal)
    lnTempK = numpy.log(TempK)
    return numpy.exp( (148.9652*a) + (-13847.26*b)/TempK + (-23.6521*c)*lnTempK + sqrtSal*( (118.67*d)/TempK + (-5.977*e) + (1.0495*f)*lnTempK ) + (-0.01615*g)*Sal  )
# fitting the params for Kw
popt, pcov = curve_fit(func_Kw,(TempK,Sal), Kw_X, guess_Kw,maxfev=maxfevN)
params_Kw = [[148.9652*popt[0], -13847.26*popt[1], -23.6521*popt[2], 118.67*popt[3], -5.977*popt[4], 1.0495*popt[5], -0.01615*popt[6]],[popt[0], popt[1], popt[2], popt[3], popt[4], popt[5], popt[6]],
             [148.9652, -13847.26, -23.6521, 118.67, -5.977, 1.0495, -0.01615]]
test = func_Kw((TempK,Sal),params_Kw[1][0],params_Kw[1][1],params_Kw[1][2],params_Kw[1][3],params_Kw[1][4],params_Kw[1][5],params_Kw[1][6])
if (INTERACTIVE != 0):
    print("lnKw = {0:.4f} (+) {1:.2f}/T (+) {2:.4f}*lnT (+) Sal^0.5*[{3:.2f}/T (+) {4:.3f} (+) {5:.4f}*lnT] (+) {6:.5f}*Sal".format(params_Kw[0][0],params_Kw[0][1],params_Kw[0][2],params_Kw[0][3],params_Kw[0][4],params_Kw[0][5],params_Kw[0][6]))
    print("Model fitted parameters at user defined [Ca2+] and [Mg2+] relative to parameters fit to experimental data for modern seawater:", popt)
    print("Average accuracy of fitting the function to the model output:", numpy.sum(numpy.sqrt((test/Kw_X)*(test/Kw_X)))/N)
    print()


# functional form used for KspC
guess_KspC = (1, 1, 1, 1, 1, 1, 1, 1, 1) # initial guess for fitting algorithm taken to be unity when modern params (i.e., normalized for nummeric reasons)
def func_KspC(TkS, a, b, c, d, e, f, g, h, i):
    TempK = TkS[0]
    Sal = TkS[1]
    sqrtSal = numpy.sqrt(Sal)
    log10TempK = numpy.log10(TempK)
    return numpy.power(10, (-171.9065*a) + (-0.077993*b)*TempK + (2839.319*c)/TempK + (71.595*d)*log10TempK + sqrtSal*((-0.77712*e)+(0.0028426*f)*TempK+(178.34*g)/TempK) + (-0.07711*h)*Sal+(0.0041249*i)*Sal*sqrtSal )
# fitting the params for KspC
popt, pcov = curve_fit(func_KspC,(TempK,Sal), KspC_X, guess_KspC,maxfev=maxfevN)
params_KspC = [[-171.9065*popt[0], -0.077993*popt[1], 2839.319*popt[2], 71.595*popt[3], -0.77712*popt[4], 0.0028426*popt[5], 178.34*popt[6], -0.07711*popt[7], 0.0041249*popt[8]],[popt[0], popt[1], popt[2], popt[3], popt[4], popt[5], popt[6], popt[7], popt[8]],[-171.9065, -0.077993, 2839.319, 71.595, -0.77712, 0.0028426, 178.34, -0.07711, 0.0041249]]
test = func_KspC((TempK,Sal),params_KspC[1][0],params_KspC[1][1],params_KspC[1][2],params_KspC[1][3],params_KspC[1][4],params_KspC[1][5],params_KspC[1][6],params_KspC[1][7],params_KspC[1][8])
if (INTERACTIVE != 0):
    print("log10KspC = {0:.4f} (+) {1:.6f}*T (+) {2:.3f}/T (+) {3:.3f}*log10T (+) Sal^0.5*({4:.6f} (+) {5:.7f}*T (+) {6:.3f}/T) (+) {7:.5f}*Sal (+) {8:.7f}*Sal^1.5".format(params_KspC[0][0],params_KspC[0][1],params_KspC[0][2],params_KspC[0][3],params_KspC[0][4],params_KspC[0][5],params_KspC[0][6],params_KspC[0][7],params_KspC[0][8]))
    print("Model fitted parameters at user defined [Ca2+] and [Mg2+] relative to parameters fit to experimental data for modern seawater:", popt)
    print("Average accuracy of fitting the function to the model output:", numpy.sum(numpy.sqrt((test/KspC_X)*(test/KspC_X)))/N)
    print()


# functional form used for KspA
guess_KspA = (1, 1, 1, 1, 1, 1, 1, 1, 1) # initial guess for fitting algorithm taken to be unity when modern params (i.e., normalized for nummeric reasons)
def func_KspA(TkS, a, b, c, d, e, f, g, h, i):
    TempK = TkS[0]
    Sal = TkS[1]
    sqrtSal = numpy.sqrt(Sal)
    log10TempK = numpy.log10(TempK)
    return numpy.power(10, (-171.945*a) + (-0.077993*b)*TempK + (2903.293*c)/TempK + (71.595*d)*log10TempK + sqrtSal*((-0.068393*e)+(0.0017276*f)*TempK+(88.135*g)/TempK) + (-0.10018*h)*Sal+(0.0059415*i)*Sal*sqrtSal )
# fitting the params for KspA
popt, pcov = curve_fit(func_KspA,(TempK,Sal), KspA_X, guess_KspA,maxfev=maxfevN)
params_KspA = [[-171.945*popt[0], -0.077993*popt[1], 2903.293*popt[2], 71.595*popt[3], -0.068393*popt[4], 0.0017276*popt[5], 88.135*popt[6], -0.10018*popt[7], 0.0059415*popt[8]],[popt[0], popt[1], popt[2], popt[3], popt[4], popt[5], popt[6], popt[7], popt[8]],[-171.945, -0.077993, 2903.293, 71.595, -0.068393, 0.0017276, 88.135, -0.10018, 0.0059415]]
test = func_KspA((TempK,Sal),params_KspA[1][0],params_KspA[1][1],params_KspA[1][2],params_KspA[1][3],params_KspA[1][4],params_KspA[1][5],params_KspA[1][6],params_KspA[1][7],params_KspA[1][8])
if (INTERACTIVE != 0):
    print("log10KspA = {0:.4f} (+) {1:.6f}*T (+) {2:.3f}/T (+) {3:.3f}*log10T (+) Sal^0.5*({4:.6f} (+) {5:.7f}*T (+) {6:.3f}/T) (+) {7:.5f}*Sal (+) {8:.7f}*Sal^1.5".format(params_KspA[0][0],params_KspA[0][1],params_KspA[0][2],params_KspA[0][3],params_KspA[0][4],params_KspA[0][5],params_KspA[0][6],params_KspA[0][7],params_KspA[0][8]))
    print("Model fitted parameters at user defined [Ca2+] and [Mg2+] relative to parameters fit to experimental data for modern seawater:", popt)
    print("Average accuracy of fitting the function to the model output:", numpy.sum(numpy.sqrt((test/KspA_X)*(test/KspA_X)))/N)
    print()


# functional form used for KSO4
guess_KSO4 = (1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1) # initial guess for fitting algorithm taken to be unity when modern params (i.e., normalized for nummeric reasons)
def func_KSO4(TkS, a, b, c, d, e, f, g, h, i, j, k):
    TempK = TkS[0]
    Sal = TkS[1]
    I = 19.924*Sal/(1000-1.005*Sal)
    sqrtI = numpy.sqrt(I)
    lnTempK = numpy.log(TempK)
    return numpy.exp( (141.328*a) + (-4276.1*b)/TempK + (-23.093*c)*lnTempK + sqrtI*( (-13856*d)/TempK  + (324.57*e) + (-47.986*f)*lnTempK) + I*( (35474*g)/TempK  + (-771.54*h) + (114.723*i)*lnTempK) + (-2698*j)/TempK *I*sqrtI + (1776*k)/TempK *I*I + numpy.log(1-0.001005*Sal) )
# fitting the params for KSO4
popt, pcov = curve_fit(func_KSO4,(TempK,Sal), KSO4_X, guess_KSO4,maxfev=maxfevN)
params_KSO4 = [[141.328*popt[0], -4276.1*popt[1], -23.093*popt[2], -13856*popt[3], 324.57*popt[4], -47.986*popt[5], 35474*popt[6], -771.54*popt[7], 114.723*popt[8], -2698*popt[9], 1776*popt[10]],[popt[0], popt[1], popt[2], popt[3], popt[4], popt[5], popt[6], popt[7], popt[8], popt[9], popt[10]],[141.328, -4276.1, -23.093, -13856, 324.57, -47.986, 35474, -771.54, 114.723, -2698, 1776]]
test = func_KSO4((TempK,Sal),params_KSO4[1][0],params_KSO4[1][1],params_KSO4[1][2],params_KSO4[1][3],params_KSO4[1][4],params_KSO4[1][5],params_KSO4[1][6],params_KSO4[1][7],params_KSO4[1][8],params_KSO4[1][9],params_KSO4[1][10])
if (INTERACTIVE != 0):
    print("lnKSO4 = {0:.3f} (+) {1:.1f}/T (+) {2:.3f}*lnT (+) (I^0.5*[{3:.0f}/T (+) {4:.2f} (+) {5:.3f}lnT] (+) (I*[{6:.0f}/T (+) {7:.2f} (+) {8:.3f}*lnT] (+) {9:.0f}/T*I^1.5 (+) {10:.0f}/T*I^2 + ln(1-0.001005*Sal)".format(params_KSO4[0][0],params_KSO4[0][1],params_KSO4[0][2],params_KSO4[0][3],params_KSO4[0][4],params_KSO4[0][5],params_KSO4[0][6],params_KSO4[0][7],params_KSO4[0][8],params_KSO4[0][9],params_KSO4[0][10]))
    print("Model fitted parameters at user defined [Ca2+] and [Mg2+] relative to parameters fit to experimental data for modern seawater:", popt)
    print("Average accuracy of fitting the function to the model output:", numpy.sum(numpy.sqrt((test/KSO4_X)*(test/KSO4_X)))/N) #,  test, KSO4_X
    print()
    print()
    print("FINISHED !!")
    print("The parameters for all pK's have been written to the file [OUTPUT_INTERACTIVE.txt]")
    print()
    print("----------------------FINISHED--------------------------")
    print("--- MyAMI specific ion interaction model Version 0.9 ---")
    print("----------------------FINISHED--------------------------")
    print()
    print()


if (INTERACTIVE == 0):
    OUTPUTFILE = open('OUTPUT_NONINTERACTIVE.txt','w')
    

    OUTPUTFILE.write(str(params_K0[0][0])+' '+str(params_K0[0][1])+' '+str(params_K0[0][2])+' '+str(params_K0[0][3])+' '+str(params_K0[0][4])+' '+str(params_K0[0][5])+'\n')

    OUTPUTFILE.write(str(params_K1[0][0])+' '+str(params_K1[0][1])+' '+str(params_K1[0][2])+' '+str(params_K1[0][3])+' '+str(params_K1[0][4])+'\n')

    OUTPUTFILE.write(str(params_K2[0][0])+' '+str(params_K2[0][1])+' '+str(params_K2[0][2])+' '+str(params_K2[0][3])+' '+str(params_K2[0][4])+'\n')
    
    OUTPUTFILE.write(str(params_Kb[0][0])+' '+str(params_Kb[0][1])+' '+str(params_Kb[0][2])+' '+str(params_Kb[0][3])+' '+str(params_Kb[0][4])+' '+str(params_Kb[0][5])+' '+str(params_Kb[0][6])+' '+str(params_Kb[0][7])+' '+str(params_Kb[0][8])+' '+str(params_Kb[0][9])+' '+str(params_Kb[0][10])+' '+str(params_Kb[0][11])+'\n')
    
    OUTPUTFILE.write(str(params_Kw[0][0])+' '+str(params_Kw[0][1])+' '+str(params_Kw[0][2])+' '+str(params_Kw[0][3])+' '+str(params_Kw[0][4])+' '+str(params_Kw[0][5])+' '+str(params_Kw[0][6])+'\n')

    OUTPUTFILE.write(str(params_KspC[0][0])+' '+str(params_KspC[0][1])+' '+str(params_KspC[0][2])+' '+str(params_KspC[0][3])+' '+str(params_KspC[0][4])+' '+str(params_KspC[0][5])+' '+str(params_KspC[0][6])+' '+str(params_KspC[0][7])+' '+str(params_KspC[0][8])+'\n')
    
    OUTPUTFILE.write(str(params_KspA[0][0])+' '+str(params_KspA[0][1])+' '+str(params_KspA[0][2])+' '+str(params_KspA[0][3])+' '+str(params_KspA[0][4])+' '+str(params_KspA[0][5])+' '+str(params_KspA[0][6])+' '+str(params_KspA[0][7])+' '+str(params_KspA[0][8])+'\n')
    
    OUTPUTFILE.write(str(params_KSO4[0][0])+' '+str(params_KSO4[0][1])+' '+str(params_KSO4[0][2])+' '+str(params_KSO4[0][3])+' '+str(params_KSO4[0][4])+' '+str(params_KSO4[0][5])+' '+str(params_KSO4[0][6])+' '+str(params_KSO4[0][7])+' '+str(params_KSO4[0][8])+' '+str(params_KSO4[0][9])+' '+str(params_KSO4[0][10])+'\n')

#    OUTPUTFILE.write(str(params_K0[0])+'\n'+'\n')
#    OUTPUTFILE.write(str(params_K1[0])+'\n'+'\n')
#    OUTPUTFILE.write(str(params_K2[0])+'\n'+'\n')
#    OUTPUTFILE.write(str(params_Kb[0])+'\n'+'\n')
#    OUTPUTFILE.write(str(params_Kw[0])+'\n'+'\n')
#    OUTPUTFILE.write(str(params_KspC[0])+'\n'+'\n')
#    OUTPUTFILE.write(str(params_KspA[0])+'\n'+'\n')
#    OUTPUTFILE.write(str(params_KSO4[0])+'\n'+'\n')
    
    OUTPUTFILE.close() # you can omit in most cases as the destructor will call if


if (INTERACTIVE != 0):
    OUTPUTFILE = open('OUTPUT_INTERACTIVE.txt','w')


    OUTPUTFILE.write('------------------------OUTPUT--------------------------\n')
    OUTPUTFILE.write('--- MyAMI specific ion interaction model Version 0.9 ---\n')
    OUTPUTFILE.write('--- Author: Dr. Mathis P. Hain (m.p.hain@soton.ac.uk) --\n')
    OUTPUTFILE.write('------------------------OUTPUT--------------------------\n')
    OUTPUTFILE.write('[Ca2+] @ Sal=35 was set to: ' + str(XmCa) +' \n')
    OUTPUTFILE.write('[Mg2+] @ Sal=35 was set to: ' + str(XmMg) +' \n'+' \n')

    
    OUTPUTFILE.write('pK0 parameters:\n')
    OUTPUTFILE.write(str(params_K0[0])+'\n'+'\n')
    OUTPUTFILE.write('pK1 parameters:\n')
    OUTPUTFILE.write(str(params_K1[0])+'\n'+'\n')
    OUTPUTFILE.write('pK2 parameters:\n')
    OUTPUTFILE.write(str(params_K2[0])+'\n'+'\n')
    OUTPUTFILE.write('pKb parameters:\n')
    OUTPUTFILE.write(str(params_Kb[0])+'\n'+'\n')
    OUTPUTFILE.write('pKw parameters:\n')
    OUTPUTFILE.write(str(params_Kw[0])+'\n'+'\n')
    OUTPUTFILE.write('pKspC parameters:\n')
    OUTPUTFILE.write(str(params_KspC[0])+'\n'+'\n')
    OUTPUTFILE.write('pKspA parameters:\n')
    OUTPUTFILE.write(str(params_KspA[0])+'\n'+'\n')
    OUTPUTFILE.write('pKSO4 parameters:\n')
    OUTPUTFILE.write(str(params_KSO4[0])+'\n'+'\n')
    
    OUTPUTFILE.close() # you can omit in most cases as the destructor will call if



# SOME PLOTTING

K1_X2_M = numpy.zeros((n,n))
K1_X_M = numpy.zeros((n,n))

for idx in range(0,n):
    for idy in range (0,n):
        K1_X2_M[idx,idy] = func_K1((TempK_M[idx,idy],Sal_M[idx,idy]),popt[0], popt[1], popt[2], popt[3], popt[4])
        K1_X_M[idx,idy] = numpy.reshape(K1_X,(n,n))[idx,idy]
#print Kout

#K1_X_M = numpy.reshape(K1_X,(n,n))

#print numpy.sum(numpy.sqrt(numpy.power(K1_X2_M/K1_X_M,2)))/N



plt.close('all')
# Four axes, returned as a 2-d array
f, axarr = plt.subplots(3, 1)
CS = axarr[0, ].contour(Sal_M, TempC_M, -numpy.log10(K1_X_M))
plt.clabel(CS, inline=1, fontsize=10)
axarr[0, ].set_title('pK1_X')

CS = axarr[1, ].contour(Sal_M, TempC_M, -numpy.log10(K1_X2_M))
plt.clabel(CS, inline=1, fontsize=10)
axarr[1, ].set_title('pK1_X2')

CS = axarr[2, ].contour(Sal_M, TempC_M, 10000*K1_X2_M/K1_X_M-10000)
plt.clabel(CS, inline=1, fontsize=10)
axarr[2, ].set_title('K1_X2/K1_X-1')

#plt.show()


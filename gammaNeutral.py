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
import math
import numpy

def gammaCO2(Tc, m_an, m_cat):
    T = Tc + 273.15
    lnT = math.log(T)
    
    m_ion = numpy.array([m_cat[0],m_cat[1],m_cat[2],m_cat[3],m_cat[4],m_an[1],m_an[6]])

    param_lamdaCO2 = numpy.zeros([7,5])
    param_lamdaCO2[0,:]=[0,0,0,0,0] #H
    param_lamdaCO2[1,:]=[-5496.38465, -3.326566, 0.0017532, 109399.341, 1047.021567] #Na
    param_lamdaCO2[2,:]=[2856.528099, 1.7670079, -0.0009487, -55954.1929, -546.074467] #K
    param_lamdaCO2[3,:]=[-479.362533, -0.541843, 0.00038812, 3589.474052, 104.3452732] #Mg
    #param_lamdaCO2[3,:]=[9.03662673e+03, 5.08294701e+00, -2.51623005e-03, -1.88589243e+05, -1.70171838e+03] #Mg refitted
    param_lamdaCO2[4,:]=[-12774.6472, -8.101555, 0.00442472, 245541.5435, 2452.50972] #Ca
    #param_lamdaCO2[4,:]=[-8.78153999e+03, -5.67606538e+00, 3.14744317e-03, 1.66634223e+05, 1.69112982e+03] #Ca refitted
    param_lamdaCO2[5,:]=[1659.944942, 0.9964326, -0.00052122, -33159.6177, -315.827883] #Cl
    param_lamdaCO2[6,:]=[2274.656591, 1.8270948, -0.00114272, -33927.7625, -457.015738] #SO4
    
    param_zetaCO2 = numpy.zeros([2,6,5])
    param_zetaCO2[0,0,:] = [-804.121738,-0.470474,0.000240526,16334.38917,152.3838752] #Cl & H
    param_zetaCO2[0,1,:] = [-379.459185,-0.258005,0.000147823,6879.030871,73.74511574] #Cl & Na
    param_zetaCO2[0,2,:] = [-379.686097,-0.257891,0.000147333,6853.264129,73.79977116] #Cl & K
    param_zetaCO2[0,3,:] = [-1342.60256,-0.772286,0.000391603,27726.80974,253.62319406] #Cl & Mg
    param_zetaCO2[0,4,:] = [-166.06529,-0.018002,-0.0000247349,5256.844332,27.377452415] #Cl & Ca
    param_zetaCO2[1,1,:] = [67030.02482,37.930519,-0.0189473,-1399082.37,-12630.27457] #SO4 & Na
    param_zetaCO2[1,2,:] = [-2907.03326,-2.860763,0.001951086,30756.86749,611.37560512] #SO4 & K
    param_zetaCO2[1,3,:] = [-7374.24392,-4.608331,0.002489207,143162.6076,1412.302898] #SO4 & Mg
    
    lamdaCO2 = numpy.zeros(7)
    for ion in range(0,7):
        lamdaCO2[ion] = param_lamdaCO2[ion,0] + param_lamdaCO2[ion,1]*T +param_lamdaCO2[ion,2]*T*T + param_lamdaCO2[ion,3]/T + param_lamdaCO2[ion,4]*lnT

    zetaCO2 = numpy.zeros([2,5])
    for ion in range(0,5):
        zetaCO2[0,ion] = param_zetaCO2[0,ion,0] + param_zetaCO2[0,ion,1]*T +param_zetaCO2[0,ion,2]*T*T + param_zetaCO2[0,ion,3]/T + param_zetaCO2[0,ion,4]*lnT
    for ion in range(1,4):
        zetaCO2[1,ion] = param_zetaCO2[1,ion,0] + param_zetaCO2[1,ion,1]*T +param_zetaCO2[1,ion,2]*T*T + param_zetaCO2[1,ion,3]/T + param_zetaCO2[1,ion,4]*lnT

    ln_gammaCO2 = 0
    for ion in range(0,7):
        ln_gammaCO2 = ln_gammaCO2 + m_ion[ion]*2*lamdaCO2[ion]
#for cat in range(0,5):
#ln_gammaCO2 = ln_gammaCO2 + m_ion[5]*m_ion[cat]*zetaCO2[0,cat] + m_ion[6]*m_ion[cat]*zetaCO2[1,cat]

    gammaCO2 = math.exp(ln_gammaCO2) # as according to He and Morse 1993
            #gammaCO2 = math.pow(10,ln_gammaCO2) # pK1 is "correct if log-base 10 is assumed


    gammaCO2gas = math.exp(1/(8.314462175*T*(0.10476 - 61.0102/T -660000/T/T/T -2.47E27/math.pow(T,12))))


    ##########################
    # CALCULATION OF gammaB
    
    lamdaB = numpy.array([0, -0.097, -0.14, 0, 0, 0.091, 0.018]) #Felmy and Wear 1986
    #lamdaB = numpy.array([0.109, 0.028, -0.026, 0.191, 0.165, 0, -0.205]) #Chanson and Millero 2006
    ln_gammaB = m_ion[1]*m_ion[6]*0.046 #tripple ion interaction Na-SO4
    for ion in range(0,7):
        ln_gammaB = ln_gammaB + m_ion[ion]*2*lamdaB[ion]

    gammaB = math.exp(ln_gammaB) # as according to Felmy and Wear 1986
#print gammaB


    return gammaCO2, gammaCO2gas, gammaB
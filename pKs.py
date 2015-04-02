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

import gammaANDalpha
import K_thermo_conditional
import K_HSO4_thermo
import gammaNeutral

# "calculate_gKs" calculates the activity coefficients only (i.e., gamma)
def calculate_gKs(Tc,Sal, mCa, mMg):
    I = 19.924*Sal/(1000-1.005*Sal)
    
    m_cation = numpy.zeros(6)
    m_cation[0] = 0.00000001*Sal/35 # H ion; pH of about 8
    m_cation[1] = 0.4689674*Sal/35 # Na Millero et al., 2008; Dickson OA-guide
    m_cation[2] = 0.0102077*Sal/35 # K Millero et al., 2008; Dickson OA-guide
    m_cation[3] = mMg*Sal/35 # Mg Millero et al., 2008; Dickson OA-guide
    m_cation[4] = mCa*Sal/35 # Ca Millero et al., 2008; Dickson OA-guide
    m_cation[5] = 0.0000907*Sal/35 # Sr Millero et al., 2008; Dickson OA-guide
    
    m_anion = numpy.zeros(7)
    m_anion[0] = 0.0000010*Sal/35 # OH ion; pH of about 8
    m_anion[1] = 0.5458696*Sal/35 # Cl Millero et al., 2008; Dickson OA-guide
    m_anion[2] = 0.0001008*Sal/35 # BOH4 Millero et al., 2008; Dickson OA-guide; pH of about 8 -- borate, not Btotal
    m_anion[3] = 0.0017177*Sal/35 # HCO3 Millero et al., 2008; Dickson OA-guide
    m_anion[4] = 0.0282352*1e-6*Sal/35 # HSO4 Millero et al., 2008; Dickson OA-guide
    m_anion[5] = 0.0002390*Sal/35 # CO3 Millero et al., 2008; Dickson OA-guide
    m_anion[6] = 0.0282352*Sal/35 # SO4 Millero et al., 2008; Dickson OA-guide
    
    
    [gamma_cation, gamma_anion, alpha_Hsws, alpha_Ht, alpha_OH, alpha_CO3] = gammaANDalpha.CalculateGammaAndAlphas(Tc,Sal,I,m_cation,m_anion)
    
    gammaT_OH = gamma_anion[0]*alpha_OH
    gammaT_BOH4 = gamma_anion[2]
    gammaT_HCO3 = gamma_anion[3]
    gammaT_CO3 = gamma_anion[5]*alpha_CO3
    
    gammaT_Hsws = gamma_cation[0]*alpha_Hsws
    gammaT_Ht = gamma_cation[0]*alpha_Ht
    gammaT_Ca = gamma_cation[4]
    
    [gammaCO2, gammaCO2gas, gammaB] = gammaNeutral.gammaCO2(Tc, m_anion, m_cation)
    
    gKspC = 1/gammaT_Ca/gammaT_CO3
    gKspA = 1/gammaT_Ca/gammaT_CO3
    gK1 = 1/gammaT_Ht/gammaT_HCO3*gammaCO2
    gK2 = 1/gammaT_Ht/gammaT_CO3*gammaT_HCO3
    gKw = 1/gammaT_Ht/gammaT_OH
    gKb = 1/gammaT_BOH4/gammaT_Ht*(gammaB)
    gK0 = 1/gammaCO2*gammaCO2gas
    gKHSO4 = 1/gamma_anion[6]/gammaT_Ht*gamma_anion[4]

    
    return gKspC, gK1, gK2, gKw, gKb, gKspA, gK0, gKHSO4

# "calculate_Ks" first calculates the activity coefficients (gammma), second the thermodynamic equilibrium constants, and third the predicted conditional equilibrium constants
def calculate_Ks(Tc,Sal, mCa, mMg):
    I = 19.924*Sal/(1000-1.005*Sal)

    m_cation = numpy.zeros(6)
    m_cation[0] = 0.00000001*Sal/35 # H ion; pH of about 8
    m_cation[1] = 0.4689674*Sal/35 # Na Millero et al., 2008; Dickson OA-guide
    m_cation[2] = 0.0102077*Sal/35 # K Millero et al., 2008; Dickson OA-guide
    m_cation[3] = mMg*Sal/35 # Mg Millero et al., 2008; Dickson OA-guide
    m_cation[4] = mCa*Sal/35 # Ca Millero et al., 2008; Dickson OA-guide
    m_cation[5] = 0.0000907*Sal/35 # Sr Millero et al., 2008; Dickson OA-guide


    m_anion = numpy.zeros(7)
    m_anion[0] = 0.0000010*Sal/35 # OH ion; pH of about 8
    m_anion[1] = 0.5458696*Sal/35 # Cl Millero et al., 2008; Dickson OA-guide
    m_anion[2] = 0.0001008*Sal/35 # BOH4 Millero et al., 2008; Dickson OA-guide; pH of about 8 -- borate, not Btotal
    m_anion[3] = 0.0017177*Sal/35 # HCO3 Millero et al., 2008; Dickson OA-guide
    m_anion[4] = 0.0282352*1e-6*Sal/35 # HSO4 Millero et al., 2008; Dickson OA-guide
    m_anion[5] = 0.0002390*Sal/35 # CO3 Millero et al., 2008; Dickson OA-guide
    m_anion[6] = 0.0282352*Sal/35 # SO4 Millero et al., 2008; Dickson OA-guide


    [gamma_cation, gamma_anion, alpha_Hsws, alpha_Ht, alpha_OH, alpha_CO3] = gammaANDalpha.CalculateGammaAndAlphas(Tc,Sal,I,m_cation,m_anion)

    gammaT_OH = gamma_anion[0]*alpha_OH
    gammaT_BOH4 = gamma_anion[2]
    gammaT_HCO3 = gamma_anion[3]
    gammaT_CO3 = gamma_anion[5]*alpha_CO3

    gammaT_Hsws = gamma_cation[0]*alpha_Hsws
    gammaT_Ht = gamma_cation[0]*alpha_Ht
    gammaT_Ca = gamma_cation[4]

    [KspCcond, K1cond, K2cond, Kwcond, Kbcond, KspAcond, K0cond, KHSO4cond] = K_thermo_conditional.CalculateKcond(Tc,Sal)
    [KspC,K1, K2, Kw, Kb, KspA, K0, KHSO4] = K_thermo_conditional.CalculateKthermo(Tc,Sal)
    
    # moved to "K_thermo_cond"[KHSO4cond, KHSO4] = K_HSO4_thermo.supplyK(Tc+273.15,I,Sal)#*(gamma_anion[4]/gamma_anion[6]/gamma_cation[0])
    [gammaCO2, gammaCO2gas, gammaB] = gammaNeutral.gammaCO2(Tc, m_anion, m_cation)
    #print gammaT_OH, gamma_cation[0], gammaT_Ht, gammaT_Hsws, gammaCO2,gammaCO2gas, gammaT_HCO3, gammaT_CO3, gamma_cation, gamma_anion
    
    KspCcondcalc = KspC/gammaT_Ca/gammaT_CO3
    KspAcondcalc = KspA/gammaT_Ca/gammaT_CO3
    K1condcalc = K1/gammaT_Ht/gammaT_HCO3*gammaCO2
    K2condcalc = K2/gammaT_Ht/gammaT_CO3*gammaT_HCO3
    Kwcondcalc = Kw/gammaT_Ht/gammaT_OH
    Kbcondcalc = Kb/gammaT_BOH4/gammaT_Ht*(gammaB)
    K0condcalc = K0/gammaCO2*gammaCO2gas
    KHSO4condcalc = KHSO4/gamma_anion[6]/gammaT_Ht*gamma_anion[4]
    #print KHSO4condcalc/KHSO4cond, "HSO4",KHSO4cond, KHSO4,KHSO4condcalc,gamma_anion[6],gammaT_Ht,gamma_anion[4]

    Kall = numpy.zeros([9,3])
    Kall[0,:] = numpy.array([0.0, Tc, Sal])
    Kall[1,:] = numpy.array([Kw, Kwcond, Kwcondcalc])
    Kall[2,:] = numpy.array([K1, K1cond, K1condcalc])
    Kall[3,:] = numpy.array([K2, K2cond, K2condcalc])
    Kall[4,:] = numpy.array([KspC, KspCcond, KspCcondcalc])
    Kall[5,:] = numpy.array([Kb, Kbcond, Kbcondcalc])
    Kall[6,:] = numpy.array([KspA, KspAcond, KspAcondcalc])
    Kall[7,:] = numpy.array([K0, K0cond, K0condcalc])
    Kall[8,:] = numpy.array([KHSO4, KHSO4cond, KHSO4condcalc])

    #print Tc, Sal, (Kwcondcalc)/Kwcond, (K1condcalc)/K1cond, (K2condcalc)/K2cond, (KspCcondcalc)/KspCcond, Kbcondcalc/Kbcond, gammaB

    return Kall






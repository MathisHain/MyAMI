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


import PitzerParams
import K_HSO4_thermo
import K_HF_cond

def CalculateGammaAndAlphas(Tc,S,I,m_cation,m_anion):
    #Testbed case T=25C, I=0.7, seawatercomposition
    T = Tc + 273.15
    sqrtI = math.sqrt(I)
    
    Z_cation = numpy.zeros(6)
    Z_cation[0] = 1
    Z_cation[1] = 1
    Z_cation[2] = 1
    Z_cation[3] = 2
    Z_cation[4] = 2
    Z_cation[5] = 2

    Z_anion = numpy.zeros(7)
    Z_anion[0] = -1
    Z_anion[1] = -1
    Z_anion[2] = -1
    Z_anion[3] = -1
    Z_anion[4] = -1
    Z_anion[5] = -2
    Z_anion[6] = -2

    #######################################################################################
    [beta_0, beta_1, beta_2, C_phi, Theta_negative, Theta_positive, Phi_NNP, Phi_PPN, C1_HSO4] = PitzerParams.SupplyParams(T)

    
    A_phi = 3.36901532E-01 - 6.32100430E-04*T + 9.14252359/T - 1.35143986E-02*math.log(T) + 2.26089488E-03/(T-263) + 1.92118597E-6*T*T + 4.52586464E+01/(680-T) # note correction of last parameter, E+1 instead of E-1
    #A_phi = 8.66836498e1 + 8.48795942e-2*T - 8.88785150e-5*T*T + 4.88096393e-8*T*T*T -1.32731477e3/T - 1.76460172e1*math.log(T) #Spencer et al 1990

    f_gamma = -A_phi * (sqrtI/(1+1.2*sqrtI) +(2/1.2)*math.log(1+1.2*sqrtI))

    #E_cat = sum(m_cation*Z_cation)
    E_an = -sum(m_anion*Z_anion)
    E_cat = -E_an

    #BMX_phi
    BMX_phi = numpy.zeros((6,7))
    BMX = numpy.zeros((6,7))
    BMX_apostroph = numpy.zeros((6,7))
    CMX = numpy.zeros((6,7))

    for cat in range(0,6):
        for an in range(0,7):
            BMX_phi[cat,an] = beta_0[cat,an] + beta_1[cat,an] * math.exp(-2*sqrtI)
            BMX[cat,an] = beta_0[cat,an] + (beta_1[cat,an]/(2*I)) * (1 - (1+2*sqrtI)*math.exp(-2*sqrtI))
            BMX_apostroph[cat,an] = (beta_1[cat,an]/(2*I*I)) * (-1 + (1+(2*sqrtI)+(2*sqrtI))*math.exp(-2*sqrtI))
            CMX[cat,an] = C_phi[cat,an]/(2*math.sqrt(-Z_anion[an]*Z_cation[cat]))

    # BMX* and CMX are calculated differently for 2:2 ion pairs, corrections below # alpha2= 6 for borates ... see Simonson et al 1988
    cat = 3; an = 2 # MgBOH42
    BMX_phi[cat,an] = beta_0[cat,an] + beta_1[cat,an] * math.exp(-1.4*sqrtI) + beta_2[cat,an] * math.exp(-6*sqrtI)
    BMX[cat,an] = beta_0[cat,an] + (beta_1[cat,an]/(0.98*I)) * (1 - (1+1.4*sqrtI)*math.exp(-1.4*sqrtI))+ (beta_2[cat,an]/(18*I)) * (1 - (1+6*sqrtI)*math.exp(-6*sqrtI))
    BMX_apostroph[cat,an] = (beta_1[cat,an]/(0.98*I*I)) * (-1+(1+1.4*sqrtI+0.98*I)*math.exp(-1.4*sqrtI)) + (beta_2[cat,an]/(18*I)) * (-1-(1+6*sqrtI+18*I)*math.exp(-6*sqrtI))
    cat = 3; an = 6 # MgSO4
    BMX_phi[cat,an] = beta_0[cat,an] + beta_1[cat,an] * math.exp(-1.4*sqrtI) + beta_2[cat,an] * math.exp(-12*sqrtI)
    BMX[cat,an] = beta_0[cat,an] + (beta_1[cat,an]/(0.98*I)) * (1 - (1+1.4*sqrtI)*math.exp(-1.4*sqrtI))+ (beta_2[cat,an]/(72*I)) * (1 - (1+12*sqrtI)*math.exp(-12*sqrtI))
    BMX_apostroph[cat,an] = (beta_1[cat,an]/(0.98*I*I)) * (-1+(1+1.4*sqrtI+0.98*I)*math.exp(-1.4*sqrtI)) + (beta_2[cat,an]/(72*I*I)) * (-1-(1+12*sqrtI+72*I)*math.exp(-12*sqrtI))
    #BMX_apostroph[cat,an] = (beta_1[cat,an]/(0.98*I)) * (-1+(1+1.4*sqrtI+0.98*I)*math.exp(-1.4*sqrtI)) + (beta_2[cat,an]/(72*I)) * (-1-(1+12*sqrtI+72*I)*math.exp(-12*sqrtI)) #not 1/(0.98*I*I) ... compare M&P98 equation A17 with Pabalan and Pitzer 1987 equation 15c/16b
    cat = 4; an = 2 # CaBOH42
    BMX_phi[cat,an] = beta_0[cat,an] + beta_1[cat,an] * math.exp(-1.4*sqrtI) + beta_2[cat,an] * math.exp(-6*sqrtI)
    BMX[cat,an] = beta_0[cat,an] + (beta_1[cat,an]/(0.98*I)) * (1 - (1+1.4*sqrtI)*math.exp(-1.4*sqrtI))+ (beta_2[cat,an]/(18*I)) * (1 - (1+6*sqrtI)*math.exp(-6*sqrtI))
    BMX_apostroph[cat,an] = (beta_1[cat,an]/(0.98*I*I)) * (-1+(1+1.4*sqrtI+0.98*I)*math.exp(-1.4*sqrtI)) + (beta_2[cat,an]/(18*I)) * (-1-(1+6*sqrtI+18*I)*math.exp(-6*sqrtI))
    cat = 4; an = 6 # CaSO4
    BMX_phi[cat,an] = beta_0[cat,an] + beta_1[cat,an] * math.exp(-1.4*sqrtI) + beta_2[cat,an] * math.exp(-12*sqrtI)
    BMX[cat,an] = beta_0[cat,an] + (beta_1[cat,an]/(0.98*I)) * (1 - (1+1.4*sqrtI)*math.exp(-1.4*sqrtI))+ (beta_2[cat,an]/(72*I)) * (1 - (1+12*sqrtI)*math.exp(-12*sqrtI))
    BMX_apostroph[cat,an] = (beta_1[cat,an]/(0.98*I*I)) * (-1+(1+1.4*sqrtI+0.98*I)*math.exp(-1.4*sqrtI)) + (beta_2[cat,an]/(72*I)) * (-1-(1+12*sqrtI+72*I)*math.exp(-12*sqrtI))


    cat = 5; an = 2 # SrBOH42
    BMX_phi[cat,an] = beta_0[cat,an] + beta_1[cat,an] * math.exp(-1.4*sqrtI) + beta_2[cat,an] * math.exp(-6*sqrtI)
    BMX[cat,an] = beta_0[cat,an] + (beta_1[cat,an]/(0.98*I)) * (1 - (1+1.4*sqrtI)*math.exp(-1.4*sqrtI))+ (beta_2[cat,an]/(18*I)) * (1 - (1+6*sqrtI)*math.exp(-6*sqrtI))
    BMX_apostroph[cat,an] = (beta_1[cat,an]/(0.98*I*I)) * (-1+(1+1.4*sqrtI+0.98*I)*math.exp(-1.4*sqrtI)) + (beta_2[cat,an]/(18*I)) * (-1-(1+6*sqrtI+18*I)*math.exp(-6*sqrtI))
    
    # BMX* is calculated with T-dependent alpha for H-SO4; see Clegg et al., 1994 --- Millero and Pierrot are completly off for this ion pair
    xClegg = (2-1842.843*(1/T - 1/298.15))*sqrtI
    #xClegg = (2)*sqrtI
    gClegg = 2*(1-(1+xClegg)*math.exp(-xClegg))/(xClegg*xClegg)
    BMX[0,6] = beta_0[0,6] + beta_1[0,6] * gClegg # alpha = (2-1842.843*(1/T - 1/298.15)) see Table 6 in Clegg et al 1994
    BMX_apostroph[0,6] = beta_1[0,6]/I * (math.exp(-xClegg) -gClegg)

    CMX[0,6] = C_phi[0,6] + 4*C1_HSO4* (6-(6+2.5*sqrtI*(6+3*2.5*sqrtI+2.5*sqrtI*2.5*sqrtI))*math.exp(-2.5*sqrtI))/(2.5*sqrtI*2.5*sqrtI*2.5*sqrtI*2.5*sqrtI)  # w = 2.5 ... see Clegg et al., 1994

    # unusual alpha=1.7 for Na2SO4
    #BMX[1,6] = beta_0[1,6] + (beta_1[1,6]/(2.89*I)) * 2*(1 - (1+1.7*sqrtI)*math.exp(-1.7*sqrtI))
    #BMX[1,6] = beta_0[1,6] + (beta_1[1,6]/(1.7*I)) * (1 - (1+1.7*sqrtI)*math.exp(-1.7*sqrtI))

    #BMX[4,6] =BMX[4,6]*0 #knock out Ca-SO4
    
    R = 0
    S = 0
    for cat in range(0,6):
        for an in range(0,7):
            R = R + m_anion[an]*m_cation[cat]*BMX_apostroph[cat,an]
            S = S + m_anion[an]*m_cation[cat]*CMX[cat,an]

#print BMX_apostroph[4,:], "BMX'"

#print CMX[4,:], "CMX"

    gamma_anion = numpy.zeros(7)
    ln_gamma_anion = numpy.zeros(7)
    #ln_gammaCl = Z_anion[1]*Z_anion[1]*f_gamma + R - S
    #print math.exp(ln_gammaCl), ln_gammaCl
    XX = 99
    for an in range(0,7):
        ln_gamma_anion[an] = Z_anion[an]*Z_anion[an]*(f_gamma + R) + Z_anion[an]*S
        if an==XX:
            print(ln_gamma_anion[an], "init")
        for cat in range(0,6):
            ln_gamma_anion[an] = ln_gamma_anion[an] + 2*m_cation[cat]*(BMX[cat,an]+ E_cat*CMX[cat,an])
            if an==XX:
                print(ln_gamma_anion[an], cat)
        for an2 in range(0,7):
            ln_gamma_anion[an] = ln_gamma_anion[an] + m_anion[an2]*(2*Theta_negative[an,an2])
            if an==XX:
                print(ln_gamma_anion[an], an2)
        for an2 in range(0,7):
            for cat in range(0,6):
                ln_gamma_anion[an] = ln_gamma_anion[an] + m_anion[an2] * m_cation[cat] * Phi_NNP[an,an2,cat]
                if an==XX:
                    print(ln_gamma_anion[an], an2, cat)
        for cat in range(0,6):
            for cat2 in range(cat+1,6):
                ln_gamma_anion[an] = ln_gamma_anion[an] + m_cation[cat]*m_cation[cat2]*Phi_PPN[cat,cat2,an]
                if an==XX:
                    print(ln_gamma_anion[an], cat, cat2)

    gamma_cation = numpy.zeros(6)
    ln_gamma_cation = numpy.zeros(6)
    #ln_gammaCl = Z_anion[1]*Z_anion[1]*f_gamma + R - S
    #print math.exp(ln_gammaCl), ln_gammaCl
    XX = 99
    for cat in range(0,6):
        ln_gamma_cation[cat] = Z_cation[cat]*Z_cation[cat]*(f_gamma + R) + Z_cation[cat]*S
        if cat==XX:
            print(ln_gamma_cation[cat], "init")
        for an in range(0,7):
            ln_gamma_cation[cat] = ln_gamma_cation[cat] + 2*m_anion[an]*(BMX[cat,an]+ E_cat*CMX[cat,an])
            if cat==XX:
                print(ln_gamma_cation[cat], an , BMX[cat,an], E_cat*CMX[cat,an])
        for cat2 in range(0,6):
            ln_gamma_cation[cat] = ln_gamma_cation[cat] + m_cation[cat2]*(2*Theta_positive[cat,cat2])
            if cat==XX:
                print(ln_gamma_cation[cat], cat2)
        for cat2 in range(0,6):
            for an in range(0,7):
                ln_gamma_cation[cat] = ln_gamma_cation[cat] +  m_cation[cat2] * m_anion[an] * Phi_PPN[cat,cat2,an]
                if cat==XX:
                    print(ln_gamma_cation[cat], cat2, an)
        for an in range(0,7):
            for an2 in range(an+1,7):
                ln_gamma_cation[cat] = ln_gamma_cation[cat] + m_anion[an]*m_anion[an2]*Phi_NNP[an,an2,cat]
                if cat==XX:
                    print(ln_gamma_cation[cat], an, an2)
    for an in range(0,7):
        gamma_anion[an] = math.exp(ln_gamma_anion[an])

    for cat in range(0,6):
        gamma_cation[cat] = math.exp(ln_gamma_cation[cat])

    # choice of pH-scale = total pH-scale [H]T = [H]F + [HSO4]
    # so far gamma_H is the [H]F activity coefficient (= free-H pH-scale)
    # thus, conversion is required
    [K_HSO4_conditional, K_HSO4] = K_HSO4_thermo.supplyK(T,I,S)#*(gamma_anion[4]/gamma_anion[6]/gamma_cation[0])
    #print K_HSO4_conditional
    #print gamma_anion[4], gamma_anion[6], gamma_cation[0]
    #alpha_H = 1/(1+ m_anion[6]/K_HSO4_conditional + 0.0000683/(7.7896E-4*1.1/0.3/gamma_cation[0]))
    alpha_Hsws = 1/(1+ m_anion[6]/K_HSO4_conditional + 0.0000683/(K_HF_cond.supplyK(T,sqrtI)))
    alpha_Ht = 1/(1+ m_anion[6]/K_HSO4_conditional)
    #alpha_H = 1/(1+ m_anion[6]/K_HSO4_conditional)
    


    # A number of ion pairs are calculated explicitly: MgOH, CaCO3, MgCO3, SrCO3
    # since OH and CO3 are rare compared to the anions the anion alpha (free/total) are assumed to be unity
    gamma_MgCO3 = 1
    gamma_CaCO3 = gamma_MgCO3
    gamma_SrCO3 = gamma_MgCO3

    b0b1CPhi_MgOH = numpy.array([-0.1, 1.658, 0, 0.028])
    BMX_MgOH = b0b1CPhi_MgOH[0] + (b0b1CPhi_MgOH[1]/(2*I)) * (1 - (1+2*sqrtI)*math.exp(-2*sqrtI))
    ln_gamma_MgOH = 1*(f_gamma + R) + (1)*S
    ln_gamma_MgOH = ln_gamma_MgOH + 2*m_anion[1]*(BMX_MgOH+ E_cat*b0b1CPhi_MgOH[2]) # interaction between MgOH-Cl affects MgOH gamma
    ln_gamma_MgOH = ln_gamma_MgOH  +  m_cation[3] * m_anion[1] * b0b1CPhi_MgOH[3] # interaction between MgOH-Mg-OH affects MgOH gamma
    gamma_MgOH = math.exp(ln_gamma_MgOH)

    K_MgOH = math.pow(10,-(3.87-501.6/T))/(gamma_cation[3]*gamma_anion[0]/gamma_MgOH)
    K_MgCO3 = math.pow(10,-(1.028+0.0066154*T))/(gamma_cation[3]*gamma_anion[5]/gamma_MgCO3)
    K_CaCO3 = math.pow(10,-(1.178+0.0066154*T))/(gamma_cation[4]*gamma_anion[5]/gamma_CaCO3)
#K_CaCO3 = math.pow(10,(-1228.732 - 0.299444*T + 35512.75/T +485.818*math.log10(T)))/(gamma_cation[4]*gamma_anion[5]/gamma_CaCO3) #Plummer and Busenberg82
#K_MgCO3 = math.pow(10,(-1228.732 +(0.15) - 0.299444*T + 35512.75/T +485.818*math.log10(T)))/(gamma_cation[4]*gamma_anion[5]/gamma_CaCO3)#Plummer and Busenberg82
    K_SrCO3 = math.pow(10,-(1.028+0.0066154*T))/(gamma_cation[5]*gamma_anion[5]/gamma_SrCO3)

    alpha_OH = 1/(1+ (m_cation[3]/K_MgOH))
    alpha_CO3 = 1/(1+ (m_cation[3]/K_MgCO3) + (m_cation[4]/K_CaCO3) + (m_cation[5]/K_SrCO3))
    
    
    return gamma_cation, gamma_anion, alpha_Hsws, alpha_Ht, alpha_OH, alpha_CO3




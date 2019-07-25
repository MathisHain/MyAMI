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


def SupplyParams(T): # assumes T [K] -- not T [degC]
    
    Tinv = 1/T
    lnT = math.log(T)
    ln_of_Tdiv29815 = math.log(T/298.15)
    Tpower2 = T*T
    Tpower3 = Tpower2*T
    Tpower4 = Tpower2*Tpower2


    # PART 1 -- calculate thermodynamic pK's for acids, gases and complexes

    #paramerters [A, B, C, D] according to Millero (2007) Table 11
    param_HF = [-12.641, 1590.2, 0, 0]
    param_H2S = [225.8375, -13275.324, -34.64354, 0]
    param_H2O = [148.9802, -13847.26, -23.6521, 0]
    param_BOH3 = [148.0248, -8966.901, -24.4344, 0]
    param_HSO4 = [141.411, -4340.704, -23.4825, 0.016637]
    param_NH4 = [-0.25444, -6285.33, 0, 0.0001635]
    param_H2CO3 = [290.9097, -14554.21, -45.0575, 0]
    param_HCO3 = [207.6548, -11843.79, -33.6485, 0]
    param_H2SO3 = [554.963, -16700.1, -93.67, 0.1022]
    param_HSO3 = [-358.57, 5477.1, 65.31, -0.1624]
    param_H3PO4 = [115.54, -4576.7518, -18.453, 0]
    param_H2PO4 = [172.1033, -8814.715, -27.927, 0]
    param_HPO4 = [-18.126, -3070.75, 0, 0]
    param_CO2 = [-60.2409, 9345.17, 18.7533, 0]
    param_SO2 = [-142.679, 8988.76, 19.8967, -0.0021]
    param_Aragonite = [303.5363, -13348.09, -48.7537, 0]
    param_Calcite = [303.1308, -13348.09, -48.7537, 0]


    # definition of the function that takes (Temp, param) as input and returns the lnK at that temp
    def Eq_lnK_calcABCD(T,paramABCD):
        return paramABCD[0] + paramABCD[1]/T + paramABCD[2]*math.log(T) + paramABCD[3]*T
    # How to use:  ln_of_K_HCO3_at_18degC = lnK_calcABCD(18,param_HCO3)


    #paramerters [A, B, C] according to Millero (2007) Table 12
    param_MgOH = [3.87, -501.6, 0]
    param_MgF = [3.504, -501.6, 0]
    param_CaF = [3.014, -501.6, 0]
    param_MgCO3 = [1.028, 0, 0.0066154]
    param_CaCO3 = [1.178, 0, 0.0066154]
    param_SrCO3 = [1.028, 0, 0.0066154]
    param_MgH2PO4 = [1.13, 0, 0]
    param_CaH2PO4 = [1, 0, 0]
    param_MgHPO4 = [2.7, 0, 0]
    param_CaHPO4 = [2.74, 0, 0]
    param_MgPO4 = [5.63, 0, 0]
    param_CaPO4 = [7.1, 0, 0]

    # definition of the function that takes (Temp, param) as input and returns the lnK at that temp
    def lnK_calcABC(T,paramABC):
        return paramABC[0] + paramABC[1]/T  + paramABC[2]*T
    # How to use:  ln_of_K_CaHPO4_at_18degC = lnK_calcABC(18,param_CaHPO4)



    ################################################################################
    # PART 2 -- Pitzer equation (based on Millero and Pierrot (1998))

    # Table A1 (Millero and Pierrot, 1998; after Moller, 1988 & Greenberg and Moller, 1989) valid 0 to 250degC
    param_NaCl = numpy.array([(1.43783204E01, 5.6076740E-3, -4.22185236E2, -2.51226677E0, 0.0, -2.61718135E-6, 4.43854508, -1.70502337),(-4.83060685E-1, 1.40677470E-3, 1.19311989E2, 0.0, 0.0, 0.0, 0.0, -4.23433299),(-1.00588714E-1, -1.80529413E-5, 8.61185543E0, 1.2488095E-2, 0.0, 3.41172108E-8, 6.83040995E-2, 2.93922611E-1)]) # note that second value is changed to original ref (e-3 instead e01)
    param_KCl = numpy.array([[2.67375563E1, 1.00721050E-2, -7.58485453E2, -4.70624175, 0.0, -3.75994338E-6, 0.0, 0.0],[-7.41559626, 0.0, 3.22892989E2, 1.16438557, 0.0, 0.0, 0.0, -5.94578140],[-3.30531334, -1.29807848E-3, 9.12712100E1, 5.864450181E-1, 0.0, 4.95713573E-7, 0.0, 0.0]])
    param_K2SO4 = numpy.array([[4.07908797E1, 8.26906675E-3, -1.418242998E3, -6.74728848, 0.0, 0.0, 0.0, 0.0],[-1.31669651E1, 2.35793239E-2, 2.06712592E3, 0.0, 0.0, 0.0, 0.0, 0.0],[-1.88E-2, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]])
    param_CaCl2 = numpy.array([[-9.41895832E1, -4.04750026E-2, 2.34550368E3, 1.70912300E1, -9.22885841E-1, 1.51488122E-5, -1.39082000E0, 0.0],[3.4787, -1.5417E-2, 0.0, 0.0, 0.0, 3.1791E-5, 0.0, 0.0],[1.93056024E1, 9.77090932E-3, -4.28383748E2, -3.57996343, 8.82068538E-2, -4.62270238E-6, 9.91113465, 0.0]])#[-3.03578731e1, 1.36264728e-2, 7.64582238e2, 5.50458061e0, -3.27377782e-1, 5.69405869e-6, -5.36231106e-1, 0]])
    param_CaCl2_Spencer = numpy.array([[-5.62764702e1, -3.00771997e-2, 1.05630400e-5, 3.3331626e-9, 1.11730349e3, 1.06664743e1],[3.4787e0, -1.5417e-2, 3.1791e-5,0,0,0],[2.64231655e1, 2.46922993e-2, -2.48298510e-5, 1.22421864e-8, -4.18098427e2, -5.35350322e0]])
    param_CaSO4 = numpy.array([[0.015, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],[3.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],[0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]]) #corrected after Greenberg and Moller 1989 0.015 instead of 0.15
    param_SrSO4 = param_CaSO4
    param_CaSO4_Spencer = numpy.array([[0.795e-1, -0.122e-3, 0.5001e-5, 0.6704e-8, -0.15228e3, -0.6885e-2], [0.28945e1, 0.7434e-2, 0.5287e-5, -0.101513e-6, -0.208505e4, 0.1345e1]])
    
    def Equation_TabA1(T, Tinv, lnT, a):
        return a[:,0] + a[:,1]*T + a[:,2]*Tinv + a[:,3]*lnT + a[:,4]/(T-263) + a[:,5]*T*T + a[:,6]/(680-T) + a[:,7]/(T-227)

    def EquationSpencer(T, lnT, q):
        return q[:,0] + q[:,1]*T + q[:,2]*T*T + q[:,3]*T*T*T + q[:,4]/T + q[:,5]*lnT
    

    
    
    # Table A2 (Millero and Pierrot, 1998; after Pabalan and Pitzer, 1987) valid 25 to 200degC
    param_MgCl2 = numpy.array([[0.576066, -9.31654E-04, 5.93915E-07],[2.60135, -0.0109438, 2.60169E-05],[0.059532, -2.49949E-04, 2.41831E-07]])
    param_MgSO4 = numpy.array([[-1.0282, 8.4790E-03, -2.33667E-05, 2.1575E-08, 6.8402E-04, 0.21499],[-2.9596E-01, 9.4564E-04, 0.0, 0.0, 1.1028E-02, 3.3646],[4.2164E-01, -3.5726E-03, 1.0040E-05, -9.3744E-09, -3.5160E-04, 2.7972E-02]])
    #param_MgSO4 = numpy.array([[-1.0282, 8.4790E-03, -2.33667E-05, 2.1575E-08, 6.8402E-04, 0.21499],[-2.9596E-01, 9.4564E-04, 0.0, 0.0, 1.1028E-02, 3.3646],[1.0541E-01, -8.9316E-04, 2.51E-06, -2.3436E-09, -8.7899E-05, 0.006993]]) #Cparams corrected after Pabalan and Pitzer ... but note that column lists Cmx not Cphi(=4xCmx) ... MP98 is correct
    
    def Equation1_TabA2(T, q):
        return q[:,0] + q[:,1]*T + q[:,2]*T*T

    def Equation2_TabA2(T,Tpower2, Tpower3,Tpower4, q):
        return q[:,0]*((T/2)+(88804)/(2*T)-298) + q[:,1]*((Tpower2/6)+(26463592)/(3*T)-(88804/2)) + q[:,2]*(Tpower3/12 + 88804*88804/(4*T)-26463592/3) +q[:,3]*((Tpower4/20)+88804*26463592/(5*T)-88804*88804/4) + q[:,4]*(298-(88804/T)) + q[:,5]

    
    # Table A3 (Millero and Pierrot, 1998; after mutiple studies, at least valid 0 to 50degC)
    #param_NaHSO4 = numpy.array([[0.030101, -0.362E-3, 0.0],[0.818686, -0.019671, 0.0],[0.0, 0.0, 0.0]]) # corrected after Pierrot et al., 1997
    param_NaHSO4 = numpy.array([[0.0544, -1.8478e-3, 5.3937e-5],[0.3826401, -1.8431e-2, 0.0],[0.003905, 0.0, 0.0]]) # corrected after Pierrot and Millero, 1997
    param_NaHCO3 = numpy.array([[0.028, 1.0E-3, -2.6E-5/2],[0.044, 1.1E-3, -4.3E-5/2],[0.0, 0.0, 0.0]]) # corrected after Peiper and Pitzer 1982
    param_Na2SO4 = numpy.array([[6.536438E-3, -30.197349, -0.20084955],[0.8742642, -70.014123, 0.2962095],[7.693706E-3, 4.5879201, 0.019471746]])# corrected according to Hovey et al 1993; note also that alpha = 1.7, not 2
    param_Na2SO4_Moller = numpy.array([[81.6920027 + 0.0301104957*T - 2321.93726/T - 14.3780207*lnT - 0.666496111/(T-263) - 1.03923656e-05*T*T],[1004.63018 + 0.577453682*T - 21843.4467/T - 189.110656*lnT - 0.2035505488/(T-263) - 0.000323949532*T*T + 1467.72243/(680-T)],[-80.7816886 - 0.0354521126*T + 2024.3883/T + 14.619773 * lnT - 0.091697474/(T-263) + 1.43946005e-05*T*T - 2.42272049/(680-T)]]) # Moller 1988 parameters as used in Excel MIAMI code !!!!!! careful this formula assumes alpha1=2 as opposed to alpha1=1.7 for the Hovey parameters
    # XXXXX --> need to go to the calculation of beta's (to switch Hovey/Moller) and of B et al (to switch alpha1
    
    
    #param_Na2CO3 = numpy.array([[0.0362, 1.79E-3, 1.694E-21],[1.51, 2.05E-3, 1.626E-19],[0.0052, 0.0, 0.0]]) # Millero and Pierrot referenced to Peiper and Pitzer
    param_Na2CO3 = numpy.array([[0.0362, 1.79E-3, -4.22E-5/2],[1.51, 2.05E-3, -16.8E-5/2],[0.0052, 0.0, 0.0]]) # Peiper and Pitzer 1982
    # XXXX check below if Haynes 2003 is being used.

    param_NaBOH4 = numpy.array([[-0.051, 5.264E-3, 0.0],[0.0961, -1.068E-2, 0.0],[0.01498, -15.7E-4, 0.0]]) #corrected after Simonson et al 1987 5th param should be e-2

    def Equation_TabA3andTabA4andTabA5(T, a):
        return a[:,0] + a[:,1]*(T-298.15) + a[:,2]*(T-298.15)*(T-298.15)

    def Equation_Na2SO4_TabA3(T, ln_of_Tdiv29815, a):
        return a[:,0] + a[:,1]*((1/T)-(1/298.15)) + a[:,2]*ln_of_Tdiv29815


    # Table A4 (Millero and Pierrot, 1998; after mutiple studies, at least valid 5 to 45degC)
    param_KHCO3 = numpy.array([[-0.0107, 0.001, 0.0],[0.0478, 0.0011, 6.776E-21],[0.0, 0.0, 0.0]])
    param_K2CO3 = numpy.array([[0.1288, 1.1E-3, -5.1E-6],[1.433, 4.36E-3, 2.07E-5],[0.0005, 0.0, 0.0]])
    param_KBOH4 = numpy.array([[0.1469, 2.881E-3, 0.0],[-0.0989, -6.876E-3, 0.0],[-56.43/1000, -9.56E-3, 0.0]]) #corrected after Simonson et al 1988
    #same function as TabA3 "Equation_TabA3andTabA4andTabA5(T,a)"


    # Table A5 (Millero and Pierrot, 1998; after Simonson et al, 1987b; valid 5-55degC
    param_MgBOH42 = numpy.array([[-0.623, 6.496E-3, 0.0],[0.2515, -0.01713, 0.0],[0.0, 0.0, 0.0]]) #corrected after Simonson et al 1988 first param is negative
    param_CaBOH42 = numpy.array([[-0.4462, 5.393E-3, 0.0],[-0.868, -0.0182, 0.0],[0.0, 0.0, 0.0]])
    param_SrBOH42 = param_CaBOH42 # see Table A6 
    def Equation_TabA3andTabA4andTabA5_Simonson(T, a):
        return a[:,0] + a[:,1]*(T-298.15) + a[:,2]*(T-303.15)*(T-303.15)


    # Table A7 (Millero and Pierrot, 1998; after multiple studies; valid 0-50degC
    param_KOH = numpy.array([[0.1298, -0.946E-5, 9.914E-4],[0.32, -2.59E-5, 11.86E-4],[0.0041, 0.0638E-5, -0.944E-4]])
    param_SrCl2 = numpy.array([[0.28575, -0.18367E-5, 7.1E-4],[1.66725, 0.0E-5, 28.425E-4],[-0.0013, 0.0E-5, 0.0E-4]])
    def Equation_TabA7(T, P):
        return P[:,0] + P[:,1]*(8834524.639 - 88893.4225*P[:,2])*(1/T-(1/298.15)) + P[:,1]/6*(T*T-88893.4225)


    # Table A8 --- Pitzer parameters unknown; beta's known for 25degC
    Equation_KHSO4 = numpy.array([-0.0003, 0.1735, 0.0]) 
    #Equation_MgHSO42 = numpy.array([0.4746, 1.729, 0.0]) #  XX no Cphi #from Harvie et al 1984 as referenced in MP98
    Equation_MgHSO42 = numpy.array([-0.61656 - 0.00075174 * (T - 298.15), 7.716066 - 0.0164302 * (T - 298.15), 0.43026 + 0.00199601 * (T - 298.15)]) #from Pierrot and Millero 1997 as used in the Excel file
    
    #Equation_MgHCO32 = numpy.array([0.329, 0.6072, 0.0]) # Harvie et al 1984
    Equation_MgHCO32 = numpy.array([0.03, 0.8, 0.0]) #Millero and Pierrot redetermined after Thurmond and Millero 1982
    Equation_CaHSO42 = numpy.array([0.2145, 2.53, 0.0]) 
    
    ### Correction of parameters:
    Equation_CaHCO32 = numpy.array([0.2, 0.3,0])# XXX after He and Morse 1993 ### original version of MyAMI was using suspect paramerters numpy.array([0.4, 2.977, 0.0]) from Harvie et al 1984 (same as in Pitzeretal85) numpy.array([0.4, 2.977, 0.0])
    ### See Comment Zeebe and Tyrrell (2018) and Response by Hain et al (2018) -- both published in GBC
    
    Equation_CaOH2 = numpy.array([-0.1747, -0.2303, -5.72]) # according to Harvie84, the -5.72 should be for beta2, not Cphi (which is zero) -- but likely typo in original ref since 2:1 electrolytes don't usually have a beta2
    Equation_SrHSO42 = Equation_CaHSO42 
    Equation_SrHCO32 = Equation_CaHCO32 
    Equation_SrOH2 = Equation_CaOH2
    Equation_MgOHCl = numpy.array([-0.1, 1.658, 0.0])
    Equation_NaOH = numpy.array([0.0864, 0.253, 0.0044]) # Rai et al 2002 ref to Pitzer91(CRC Press)
    Equation_CaSO4_PnM74 = numpy.array([0.2, 2.65, 0]) #Pitzer and Mayorga74
                                 
    # Table A9 --- (Millero and Pierrot, 1998; after multiple studies; valid 0-50degC
    param_HCl = numpy.array([[1.2859, -2.1197e-3, -142.58770],[-4.4474, 8.425698E-3, 665.7882],[-0.305156, 5.16E-4, 45.521540]]) # beta1 first param corrected to negative according to original reference (Campbell et al)
    param_HSO4 = numpy.array([[0.065, 0.134945, 0.022374, 7.2E-5],[-15.009, -2.405945, 0.335839, -0.004379],[0.008073, -0.113106, -0.003553, 3.57E-5]]) #XXXXX two equations for C
    param_HSO4_Clegg94 = numpy.array([[0.0348925351, 4.97207803, 0.317555182, 0.00822580341],[-1.06641231, -74.6840429, -2.26268944, -0.0352968547],[0.00764778951, -0.314698817, -0.0211926525, 0.000586708222],[0.0, -0.176776695, -0.731035345, 0.0]])
    def Equation_HCl(T, a):
        return a[:,0] + a[:,1]*T +a[:,2]/T
    def Equation_HSO4(T, a):
        return a[:,0] + (T-328.15)*1E-3 * (a[:,1] + (T-328.15)*((a[:,2]/2)+(T-328.15)*(a[:,3]/6)))

    def Equation_HSO4_Clegg94(T, a):
        return a[:,0] + (T-328.15) * (1E-3 *a[:,1]+(T-328.15) * ((1e-3*a[:,2]/2)+(T-328.15)*1e-3*a[:,3]/6))

    ############################################################
    # beta_0, beta_1 and C_phi values arranged into arrays
    N_cations = 6;  # H+=0; Na+=1; K+=2; Mg2+=3; Ca2+=4; Sr2+=5
    N_anions = 7;   # OH-=0; Cl-=1; B(OH)4-=2; HCO3-=3; HSO4-=4; CO3-=5; SO4-=6;

    beta_0 = numpy.zeros((N_cations,N_anions)) # creates empty array
    beta_1 = numpy.zeros((N_cations,N_anions)) # creates empty array
    C_phi = numpy.zeros((N_cations,N_anions)) # creates empty array

    # H = cation
    #[beta_0[0,0],beta_1[0,0],C_phi[0,0]] = n/a
    [beta_0[0,1],beta_1[0,1],C_phi[0,1]] = Equation_HCl(T,param_HCl)
    #[beta_0[0,2],beta_1[0,2],C_phi[0,2]] = n/a
    #[beta_0[0,3],beta_1[0,3],C_phi[0,3]] = n/a
    #[beta_0[0,4],beta_1[0,4],C_phi[0,4]] = n/a
    #[beta_0[0,5],beta_1[0,5],C_phi[0,5]] = n/a
    #[beta_0[0,6],beta_1[0,6],C_phi[0,6]] = Equation_HSO4(T,param_HSO4)
    #[beta_0[0,6],beta_1[0,6],C_phi[0,6], C1_HSO4] = Equation_HSO4_Clegg94(T,param_HSO4_Clegg94)
    C1_HSO4 =0
    #print beta_0[0,:],beta_1[0,:]#,beta_2[0,:]
    
    # Na = cation
    [beta_0[1,0],beta_1[1,0],C_phi[1,0]] = Equation_NaOH
    [beta_0[1,1],beta_1[1,1],C_phi[1,1]] = Equation_TabA1(T, Tinv, lnT, param_NaCl)
    [beta_0[1,2],beta_1[1,2],C_phi[1,2]] = Equation_TabA3andTabA4andTabA5(T,param_NaBOH4)
    [beta_0[1,3],beta_1[1,3],C_phi[1,3]] = Equation_TabA3andTabA4andTabA5(T,param_NaHCO3)
    [beta_0[1,4],beta_1[1,4],C_phi[1,4]] = Equation_TabA3andTabA4andTabA5(T,param_NaHSO4)
    [beta_0[1,5],beta_1[1,5],C_phi[1,5]] = Equation_TabA3andTabA4andTabA5(T,param_Na2CO3)
    [beta_0[1,6],beta_1[1,6],C_phi[1,6]] = param_Na2SO4_Moller #Equation_Na2SO4_TabA3(T,ln_of_Tdiv29815, param_Na2SO4)

    # K = cation
    [beta_0[2,0],beta_1[2,0],C_phi[2,0]] = Equation_TabA7(T,param_KOH)
    [beta_0[2,1],beta_1[2,1],C_phi[2,1]] = Equation_TabA1(T, Tinv, lnT, param_KCl)
    [beta_0[2,2],beta_1[2,2],C_phi[2,2]] = Equation_TabA3andTabA4andTabA5(T,param_KBOH4)
    [beta_0[2,3],beta_1[2,3],C_phi[2,3]] = Equation_TabA3andTabA4andTabA5(T,param_KHCO3)
    [beta_0[2,4],beta_1[2,4],C_phi[2,4]] = Equation_KHSO4
    [beta_0[2,5],beta_1[2,5],C_phi[2,5]] = Equation_TabA3andTabA4andTabA5(T,param_K2CO3)
    [beta_0[2,6],beta_1[2,6],C_phi[2,6]] = Equation_TabA1(T, Tinv, lnT, param_K2SO4)

    # Mg = cation
    #[beta_0[3,0],beta_1[3,0],C_phi[3,0]] = n/a
    [beta_0[3,1],beta_1[3,1],C_phi[3,1]] = Equation1_TabA2(T,param_MgCl2)
    [beta_0[3,2],beta_1[3,2],C_phi[3,2]] = Equation_TabA3andTabA4andTabA5_Simonson(T,param_MgBOH42)
    [beta_0[3,3],beta_1[3,3],C_phi[3,3]] = Equation_MgHCO32
    [beta_0[3,4],beta_1[3,4],C_phi[3,4]] = Equation_MgHSO42
    #[beta_0[3,5],beta_1[3,5],C_phi[3,5]] = n/a
    [beta_0[3,6],beta_1[3,6],C_phi[3,6]] = Equation2_TabA2(T,Tpower2, Tpower3, Tpower4, param_MgSO4)
    #print beta_0[3,6],beta_1[3,6],C_phi[3,6]

    # Ca = cation
    [beta_0[4,0],beta_1[4,0],C_phi[4,0]] = Equation_CaOH2
    [beta_0[4,1],beta_1[4,1],C_phi[4,1]] = Equation_TabA1(T, Tinv, lnT, param_CaCl2)
    [beta_0[4,2],beta_1[4,2],C_phi[4,2]] = Equation_TabA3andTabA4andTabA5_Simonson(T,param_CaBOH42)
    [beta_0[4,3],beta_1[4,3],C_phi[4,3]] = Equation_CaHCO32
    [beta_0[4,4],beta_1[4,4],C_phi[4,4]] = Equation_CaHSO42
    #[beta_0[4,5],beta_1[4,5],C_phi[4,5]] = n/a
    [beta_0[4,6],beta_1[4,6],C_phi[4,6]] = Equation_CaSO4_PnM74 # Equation_TabA1(T, Tinv, lnT, param_CaSO4) #

    # Sr = cation
    [beta_0[5,0],beta_1[5,0],C_phi[5,0]] = Equation_SrOH2
    [beta_0[5,1],beta_1[5,1],C_phi[5,1]] = Equation_TabA7(T,param_SrCl2)
    [beta_0[5,2],beta_1[5,2],C_phi[5,2]] = Equation_TabA3andTabA4andTabA5_Simonson(T,param_SrBOH42)
    [beta_0[5,3],beta_1[5,3],C_phi[5,3]] = Equation_SrHCO32
    [beta_0[5,4],beta_1[5,4],C_phi[5,4]] = Equation_SrHSO42
    #[beta_0[5,5],beta_1[5,5],C_phi[5,5]] = n/a
    [beta_0[5,6],beta_1[5,6],C_phi[5,6]] = Equation_CaSO4_PnM74 #Equation_TabA1(T, Tinv, lnT, param_SrSO4)
    
    # for 2:2 ion pairs beta_2 is needed
    beta_2 = numpy.zeros((N_cations,N_anions))
    b2_param_MgSO4 = numpy.array([-13.764, 0.12121, -2.7642e-4, 0, -0.21515, -32.743])
    def Eq_b2_MgSO4(T,Tpower2, Tpower3,Tpower4, q):
        return q[0]*((T/2)+(88804)/(2*T)-298) + q[1]*((Tpower2/6)+(26463592)/(3*T)-(88804/2)) + q[2]*(Tpower3/12 + 88804*88804/(4*T)-26463592/3) +q[3]*((Tpower4/20)+88804*26463592/(5*T)-88804*88804/4) + q[4]*(298-(88804/T)) + q[5]
    b2_param_MgBOH42 = numpy.array([-11.47, 0.0, -3.24e-3])
    b2_param_CaBOH42 = numpy.array([-15.88, 0.0, -2.858e-3])
    def Eq_b2_MgANDCaBOH42(T, a):
        return a[0] + a[1]*(T-298.15) + a[2]*(T-303.15)*(T-303.15)
    b2_param_CaSO4 = numpy.array([-55.7,0]) # Pitzer and Mayorga74 # [-1.29399287e2, 4.00431027e-1]) Moller88
    def Eq_b2_CaSO4(T, a):
        return a[0] + a[1]*T


    beta_2[3,6] = Eq_b2_MgSO4(T,Tpower2,Tpower3,Tpower4,b2_param_MgSO4)
    beta_2[3,2] = Eq_b2_MgANDCaBOH42(T,b2_param_MgBOH42)
    beta_2[4,2] = Eq_b2_MgANDCaBOH42(T,b2_param_CaBOH42)
    beta_2[4,6] = Eq_b2_CaSO4(T,b2_param_CaSO4)
    beta_2[5,2] = beta_2[4,2]


    #############################################################################
    #############################################################################
    # Data and T-based calculations to create arrays holding Theta and Phi values
    # based on Table A10 and A11

    # Theta of positive ions H+=0; Na+=1; K+=2; Mg2+=3; Ca2+=4; Sr2+=5
    Theta_positive = numpy.zeros((6,6)) # Array to hold Theta values between ion two ions (for numbering see list above)

    # H-Sr
    Theta_positive[0,5] = 0.0591 + 4.5*1E-4*(T-298.15)
    Theta_positive[5,0] = Theta_positive[0,5]

    # H-Na
    Theta_positive[0,1] = 0.03416 - 2.09*1E-4*(T-298.15)
    Theta_positive[1,0] = Theta_positive[0,1]

    # H-K
    Theta_positive[0,2] = 0.005 - 2.275*1E-4*(T-298.15)
    Theta_positive[2,0] = Theta_positive[0,2]

    # H-Mg
    Theta_positive[0,3] = 0.062 + 3.275*1E-4*(T-298.15)
    Theta_positive[3,0] = Theta_positive[0,3]

    # H-Ca
    Theta_positive[0,4] = 0.0612 + 3.275*1E-4*(T-298.15)
    Theta_positive[4,0] = Theta_positive[0,4]

    # Na-K
    Theta_positive[1,2] = -5.02312111E-2 + 14.0213141/T
    Theta_positive[2,1] = Theta_positive[1,2]

    # Na-Mg
    Theta_positive[1,3] = 0.07
    Theta_positive[3,1] = 0.07

    # Na-Ca
    Theta_positive[1,4] = 0.05
    Theta_positive[4,1] = 0.05

    # K-Mg
    Theta_positive[2,3] = 0.0
    Theta_positive[3,2] = 0.0

    # K-Ca
    Theta_positive[2,4] = 0.1156
    Theta_positive[4,2] = 0.1156

    # Sr-Na
    Theta_positive[5,1] = 0.07
    Theta_positive[1,5] = 0.07

    # Sr-K
    Theta_positive[5,2] = 0.01
    Theta_positive[2,5] = 0.01

    # Mg-Ca
    Theta_positive[3,4] = 0.007
    Theta_positive[4,3] = 0.007
    #print 5.31274136 - 6.3424248e-3*T - 9.83113847e2/T, "ca-mg" #Spencer et al 1990


    # Theta of negative ions  OH-=0; Cl-=1; B(OH)4-=2; HCO3-=3; HSO4-=4; CO3-=5; SO4-=6;
    Theta_negative = numpy.zeros((7,7)) # Array to hold Theta values between ion two ions (for numbering see list above)

    # Cl-SO4
    Theta_negative[1,6] = 0.07
    Theta_negative[6,1] = 0.07

    # Cl-CO3
    Theta_negative[1,5] = -0.092 #corrected after Pitzer and Peiper 1982
    Theta_negative[5,1] = -0.092 #corrected after Pitzer and Peiper 1982

    # Cl-HCO3
    Theta_negative[1,3] = 0.0359
    Theta_negative[3,1] = 0.0359

    # Cl-BOH4
    Theta_negative[1,2] = -0.0323 - 0.42333*1E-4*(T-298.15) - 21.926*1E-6*(T-298.15)*(T-298.15)
    Theta_negative[2,1] = Theta_negative[1,2]

    # CO3-HCO3
    Theta_negative[3,5] = 0.0
    Theta_negative[5,3] = 0.0

    # SO4-HSO4
    Theta_negative[4,6] = 0.0
    Theta_negative[6,4] = 0.0

    # OH-Cl
    Theta_negative[0,1] = -0.05 + 3.125*1E-4*(T-298.15) - 8.362*1E-6*(T-298.15)*(T-298.15)
    Theta_negative[1,0] = Theta_negative[0,1]

    # SO4-CO3
    Theta_negative[5,6] = 0.02
    Theta_negative[6,5] = 0.02

    # SO4-HCO3
    Theta_negative[3,6] = 0.01
    Theta_negative[6,3] = 0.01

    # SO4-BOH4
    Theta_negative[2,6] = -0.012
    Theta_negative[6,2] = -0.012

    # HSO4-Cl
    Theta_negative[1,4] = -0.006
    Theta_negative[4,1] = -0.006

    # OH-SO4
    Theta_negative[0,6] = -0.013
    Theta_negative[6,0] = -0.013
    
    # CO3-OH #http://www.aim.env.uea.ac.uk/aim/accent4/parameters.html
    Theta_negative[3,0] = 0.1
    Theta_negative[0,3] = 0.1


    # Phi
    # positive ions H+=0; Na+=1; K+=2; Mg2+=3; Ca2+=4; Sr2+=5
    # negative ions  OH-=0; Cl-=1; B(OH)4-=2; HCO3-=3; HSO4-=4; CO3-=5; SO4-=6;

    # Phi_PPN holds the values for cation-cation-anion
    Phi_PPN = numpy.zeros((6,6,7)) # Array to hold Theta values between ion two ions (for numbering see list above)

    # Na-K-Cl
    Phi_PPN[1,2,1] = 1.34211308E-2 - 5.10212917/T
    Phi_PPN[2,1,1] = Phi_PPN[1,2,1]

    # Na-K-SO4
    Phi_PPN[1,2,6] = 3.48115174E-2 - 8.21656777/T
    Phi_PPN[2,1,6] = Phi_PPN[1,2,6]

    # Na-Mg-Cl
    Phi_PPN[1,3,1] = 0.0199 - 9.51/T
    Phi_PPN[3,1,1] = Phi_PPN[1,3,1]

    # Na-Ca-Cl
    Phi_PPN[1,4,1] = -7.6398 -1.2990e-2*T + 1.1060e-5*T*T + 1.8475*lnT #Spencer et al 1990 # -0.003
    Phi_PPN[4,1,1] = Phi_PPN[1,4,1]
    #print -7.6398 -1.2990e-2*T + 1.1060e-5*T*T + 1.8475*lnT

    # Na-Ca-SO4
    Phi_PPN[1,4,6] = -0.012
    Phi_PPN[4,1,6] = Phi_PPN[1,4,6]

    # K-Mg-Cl
    Phi_PPN[2,3,1] = 0.02586 - 14.27/T
    Phi_PPN[3,2,1] = Phi_PPN[2,3,1]

    # K-Ca-Cl
    Phi_PPN[2,4,1] = 0.047627877 - 27.0770507/T
    Phi_PPN[4,2,1] = Phi_PPN[2,4,1]

    # K-Ca-SO4
    Phi_PPN[2,4,6] = 0.0
    Phi_PPN[4,2,6] = 0.0

    # H-Sr-Cl
    Phi_PPN[0,5,1] = 0.0054 - 2.1*1E-4*(T-298.15)
    Phi_PPN[5,0,1] = Phi_PPN[0,5,1]

    # H-Mg-Cl
    Phi_PPN[0,3,1] = 0.001 - 7.325*1E-4*(T-298.15)
    Phi_PPN[3,0,1] = Phi_PPN[0,3,1]

    # H-Ca-Cl
    Phi_PPN[0,4,1] = 0.0008 - 7.25*1E-4*(T-298.15)
    Phi_PPN[4,0,1] = Phi_PPN[0,4,1]

    # Sr-Na-Cl
    Phi_PPN[5,1,1] = -0.015
    Phi_PPN[1,5,1] = -0.015

    # Sr-K-Cl
    Phi_PPN[5,2,1] = -0.015
    Phi_PPN[2,5,1] = -0.015

    # Na-Mg-SO4
    Phi_PPN[1,3,6] = -0.015
    Phi_PPN[3,1,6] = -0.015

    # K-Mg-SO4
    Phi_PPN[2,3,6] = -0.048
    Phi_PPN[3,2,6] = -0.048

    # Mg-Ca-Cl
    Phi_PPN[3,4,1] = 4.15790220e1 + 1.30377312e-2*T -9.81658526e2/T -7.4061986*lnT #Spencer et al 1990 #-0.012
    Phi_PPN[4,3,1] = Phi_PPN[3,4,1]
    #print 4.15790220e1 + 1.30377312e-2*T -9.81658526e2/T -7.4061986*lnT

    # Mg-Ca-SO4
    Phi_PPN[3,4,6] = 0.024
    Phi_PPN[4,3,6] = 0.024

    # H-Na-Cl
    Phi_PPN[0,1,1] = 0.0002
    Phi_PPN[1,0,1] = 0.0002

    # H-Na-SO4
    Phi_PPN[0,1,6] = 0.0
    Phi_PPN[1,0,6] = 0.0

    # H-K-Cl
    Phi_PPN[0,2,1] = -0.011
    Phi_PPN[2,0,1] = -0.011

    # H-K-SO4
    Phi_PPN[0,2,1] = 0.197
    Phi_PPN[2,0,1] = 0.197


    # Phi_PPN holds the values for anion-anion-cation
    Phi_NNP = numpy.zeros((7,7,6)) # Array to hold Theta values between ion two ions (for numbering see list above)

    # Cl-SO4-Na
    Phi_NNP[1,6,1] = -0.009
    Phi_NNP[6,1,1] = -0.009

    # Cl-SO4-K
    Phi_NNP[1,6,2] = -0.21248147 + 37.5619614/T + 2.8469833*1E-3*T
    Phi_NNP[6,1,2] = Phi_NNP[1,6,2]

    # Cl-SO4-Ca
    Phi_NNP[1,6,4] = -0.018
    Phi_NNP[6,1,4] = -0.018

    # Cl-CO3-Ca
    Phi_NNP[1,5,4] = 0.016
    Phi_NNP[5,1,4] = 0.016

    # Cl-HCO3-Na
    Phi_NNP[1,3,1] = -0.0143
    Phi_NNP[3,1,1] = -0.0143

    # Cl-BOH4-Na
    Phi_NNP[1,2,1] = -0.0132
    Phi_NNP[2,1,1] = -0.0132

    # Cl-BOH4-Mg
    Phi_NNP[1,2,3] = -0.235
    Phi_NNP[2,1,3] = -0.235

    # Cl-BOH4-Ca
    Phi_NNP[1,2,4] = -0.8
    Phi_NNP[2,1,4] = -0.8

    # HSO4-SO4-Na
    Phi_NNP[4,6,1] = 0.0
    Phi_NNP[6,4,1] = 0.0

    # CO3-HCO3-Na
    Phi_NNP[3,5,1] = 0.0
    Phi_NNP[5,3,1] = 0.0

    # CO3-HCO3-K
    Phi_NNP[3,5,2] = 0.0
    Phi_NNP[5,3,2] = 0.0

    # Cl-SO4-Mg
    Phi_NNP[1,6,3] = -0.004
    Phi_NNP[6,1,3] = -0.004

    # Cl-HCO3-Mg
    Phi_NNP[1,3,3] = -0.0196
    Phi_NNP[3,1,3] = -0.0196

    # SO4-CO3-Na
    Phi_NNP[6,5,1] = -0.005
    Phi_NNP[5,6,1] = -0.005

    # SO4-CO3-K
    Phi_NNP[6,5,2] = -0.009
    Phi_NNP[5,6,2] = -0.009

    # SO4-HCO3-Na
    Phi_NNP[6,3,1] = -0.005
    Phi_NNP[3,6,1] = -0.005

    # SO4-HCO3-Mg
    Phi_NNP[6,3,3] = -0.161
    Phi_NNP[3,6,3] = -0.161

    # HSO4-Cl-Na
    Phi_NNP[4,1,1] = -0.006
    Phi_NNP[1,4,1] = -0.006

    # HSO4-SO4-K
    Phi_NNP[4,6,2] = -0.0677
    Phi_NNP[6,4,2] = -0.0677

    # OH-Cl-Na
    Phi_NNP[0,1,1] = -0.006
    Phi_NNP[1,0,1] = -0.006

    # OH-Cl-K
    Phi_NNP[0,1,2] = -0.006
    Phi_NNP[1,0,2] = -0.006

    # OH-Cl-Ca
    Phi_NNP[0,1,4] = -0.025
    Phi_NNP[1,0,4] = -0.025

    # OH-SO4-Na
    Phi_NNP[0,6,1] = -0.009
    Phi_NNP[6,0,1] = -0.009

    # OH-SO4-K
    Phi_NNP[0,6,2] = -0.05
    Phi_NNP[6,0,2] = -0.05

    return [beta_0, beta_1, beta_2, C_phi, Theta_negative, Theta_positive, Phi_NNP, Phi_PPN, C1_HSO4]









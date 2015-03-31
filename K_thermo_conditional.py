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

# definition of the function that takes (Temp) as input and returns the K at that temp
def CalculateKcond(Tc, Sal):
    sqrtSal = math.sqrt(Sal)
    T = Tc+273.15
    lnT= math.log(T)
    I = 19.924*Sal/(1000-1.005*Sal) # Ionic strength after Dickson 1990a; see Dickson et al 2007

    KspCcond = math.pow(10,-171.9065 -0.077993*T +2839.319/T +71.595*math.log10(T) +(-0.77712 +0.0028426*T +178.34/T)*sqrtSal -0.07711*Sal +0.0041249*Sal*sqrtSal)
    
    K1cond = math.pow(10,-3633.86/T +61.2172 -9.67770*lnT +0.011555*Sal - 0.0001152*Sal*Sal) #Dickson
    #K1cond = math.exp(290.9097 - 14554.21/T - 45.0575*lnT + (-228.39774 + 9714.36839/T + 34.485796*lnT)*sqrtSal + (54.20871 - 2310.48919/T - 8.19515*lnT)*Sal + (-3.969101 + 170.22169/T + 0.603627*lnT)*Sal*sqrtSal - 0.00258768*Sal*Sal) #Millero95

    K2cond = math.pow(10,-471.78/T -25.9290 +3.16967*lnT +0.01781*Sal - 0.0001122*Sal*Sal)
    
    Kwcond = math.exp(148.9652 - 13847.26/T -23.6521*lnT +(118.67/T -5.977 +1.0495*lnT)*sqrtSal -0.01615*Sal)
    
    Kbcond = math.exp((-8966.90-2890.53*sqrtSal-77.942*Sal+1.728*Sal*sqrtSal-0.0996*Sal*Sal)/T + (148.0248+137.1942*sqrtSal+1.62142*Sal) + (-24.4344-25.085*sqrtSal-0.2474*Sal)*lnT + 0.053105*sqrtSal*T) #Dickson90b
   
    KspAcond = math.pow(10,-171.945 -0.077993*T +2903.293/T +71.595*math.log10(T) +(-0.068393 +0.0017276*T +88.135/T)*sqrtSal -0.10018*Sal +0.0059415*Sal*sqrtSal)

    K0cond = math.exp(-60.2409 + 93.4517*100/T + 23.3585*math.log(T/100) + Sal*(0.023517 - 0.023656*T/100 +0.0047036*(T/100)*(T/100))) #Weiss74
    
    param_HSO4_cond = numpy.array([141.328, -4276.1, -23.093,  324.57,-13856, -47.986,  -771.54, 35474, 114.723, -2698, 1776]) #Dickson 1990
    KHSO4cond = math.exp(param_HSO4_cond[0] + param_HSO4_cond[1]/T + param_HSO4_cond[2]*math.log(T) + math.sqrt(I)*(param_HSO4_cond[3] + param_HSO4_cond[4]/T + param_HSO4_cond[5]*math.log(T)) + I*(param_HSO4_cond[6] + param_HSO4_cond[7]/T + param_HSO4_cond[8]*math.log(T)) + param_HSO4_cond[9]/T*I*math.sqrt(I) + param_HSO4_cond[10]/T*I*I + math.log(1-0.001005*Sal))

    
    return KspCcond, K1cond, K2cond, Kwcond, Kbcond, KspAcond, K0cond, KHSO4cond




# definition of the function that takes (Temp amd Sal) as input and returns the K at that temp
def CalculateKthermo(Tc, Sal):
    sqrtSal = math.sqrt(Sal)
    T = Tc+273.15
    lnT= math.log(T)
    
    #KspC = math.pow(10,-171.9065 -0.077993*T +2839.319/T +71.595*math.log10(T))
    KspC = math.exp(303.1308 - 13348.09/T -48.7537*lnT)
    
    #K1 = math.exp(290.9097 - 14554.21/T -45.0575*lnT) # Millero
    #K1 = math.pow(10,61.2172 - 3633.86/T -9.6777*lnT) #Dickson
    K1 = math.pow(10, -356.3094 - 0.06091964*T +21834.37/T +126.8339*math.log10(T) -1684915/T/T) #Plummer ad Busenberg82

    K2 = math.exp(207.6548 - 11843.79/T -33.6485*lnT)
    #K2 = math.pow(10, -107.8871 - 0.03252849*T +5151.79/T +38.92561*math.log10(T) -563713.9/T/T) #Plummer ad Busenberg82
    
    Kw = math.exp(148.9802 - 13847.26/T -23.6521*lnT)
    
    Kb = math.exp(148.0248 - 8966.90/T -24.4344*lnT)
    
    #KspA = math.pow(10,-171.945 -0.077993*T +2903.293/T +71.595*math.log10(T))
    
    KspA = math.exp(303.5363 - 13348.09/T -48.7537*lnT)
    
    
    K0 = math.pow(10, 108.3865 + 0.01985076*T - 6919.53/T - 40.45154*math.log10(T) + 669365/T/T) #Plummer and Busenberg 1982
    #K0 = math.exp(-60.2409 + 93.4517*100/T + 23.3585*math.log(T/100)) #truncated Dickson
    #K0 = math.exp(-60.2409 + 9345.17/T +18.7533*math.log(T)) #Millero&Pierrot ... totally off
    
    
    
    #param_HSO4 = numpy.array([562.69486, -13273.75, -102.5154, 0.2477538, -1.117033e-4]) #Clegg et al. 1994
    #K_HSO4 = math.pow(10,param_HSO4[0] + param_HSO4[1]/T + param_HSO4[2]*math.log(T) + param_HSO4[3]*T + param_HSO4[4]*T*T)
    
    param_HSO4 = numpy.array([141.411, -4340.704, -23.4825, 0.016637]) #Campbell et al. 1993
    #param_HSO4 = numpy.array([141.328, -4276.1, -23.093, 0]) #Dickson 1990
    #param_HSO4 = numpy.array([141.411, -4340.704, -23.4825, 0.016637])
    KHSO4 = math.pow(10,param_HSO4[0] + param_HSO4[1]/T + param_HSO4[2]*math.log(T) + param_HSO4[3]*T)
    
    return (KspC,K1, K2, Kw, Kb, KspA, K0, KHSO4)
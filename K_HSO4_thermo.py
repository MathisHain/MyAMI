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
def supplyK(T,I,S):
    I = pow(I,1)
    #param_HSO4 = numpy.array([562.69486, -13273.75, -102.5154, 0.2477538, -1.117033e-4]) #Clegg et al. 1994
    #K_HSO4 = math.pow(10,param_HSO4[0] + param_HSO4[1]/T + param_HSO4[2]*math.log(T) + param_HSO4[3]*T + param_HSO4[4]*T*T)

    param_HSO4 = numpy.array([141.411, -4340.704, -23.4825, 0.016637]) #Campbell et al. 1993
    #param_HSO4 = numpy.array([141.328, -4276.1, -23.093, 0]) #Dickson 1990
    #param_HSO4 = numpy.array([141.411, -4340.704, -23.4825, 0.016637])
    K_HSO4 = math.pow(10,param_HSO4[0] + param_HSO4[1]/T + param_HSO4[2]*math.log(T) + param_HSO4[3]*T)

    param_HSO4_cond = numpy.array([141.328, -4276.1, -23.093,  324.57,-13856, -47.986,  -771.54, 35474, 114.723, -2698, 1776]) #Dickson 1990
    K_HSO4_cond = math.exp(param_HSO4_cond[0] + param_HSO4_cond[1]/T + param_HSO4_cond[2]*math.log(T) + math.sqrt(I)*(param_HSO4_cond[3] + param_HSO4_cond[4]/T + param_HSO4_cond[5]*math.log(T)) + I*(param_HSO4_cond[6] + param_HSO4_cond[7]/T + param_HSO4_cond[8]*math.log(T)) + param_HSO4_cond[9]/T*I*math.sqrt(I) + param_HSO4_cond[10]/T*I*I)
    return [K_HSO4_cond, K_HSO4]

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

# definition of the function that takes (Temp) as input and returns the K at that temp
def supplyK(T,sqrtI):
    return math.exp(1590.2/T -12.641 +1.525*sqrtI)


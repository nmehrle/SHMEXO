#! /usr/bin/env python2.7
from pylab import *

data_jamrt = genfromtxt('jamrt_atm_NH3_2.7_H2O_2.0.txt')
temp_jamrt = data_jamrt[:,0]
pres_jamrt = data_jamrt[:,1]*1.E-1
cp_jamrt = data_jamrt[:,2]*1.E-7  # erg -> J/

data_armada = genfromtxt('armada_atm_NH3_2.7_H2O_2.7.txt', skip_header = 2)
temp_armada = data_armada[::-2,2]
cp_armada = data_armada[::-2,-3]

cp_interp = interp(temp_armada, temp_jamrt[::-1], cp_jamrt[::-1])

print '#%11s%12s%12s%12s' % ('TEM', 'CP(ARMADA)', 'CP(JAMRT)', 'DIFF')
for i in range(len(temp_armada)):
  print '%12.3f%12.3f%12.3f%12.3f' % (temp_armada[i], cp_armada[i], cp_interp[i],
    cp_armada[i] - cp_interp[i])

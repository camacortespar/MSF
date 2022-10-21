import numpy as np
import spiceypy as spice

import planetary_data as pd
import tools as _t
import spice_tools as st
from OrbitPropagator import OrbitPropagator as OP
from OrbitPropagator import null_perts

if __name__ == '__main__':

        #initial conditions
        coes=_t.tle2coes('test.txt',degres=True)

        print('a=',coes[0]) 
        print('e=',coes[1]) 
        print('i=',coes[2]) 
        print('ta=',coes[3]) 
        print('aop=',coes[4]) 
        print('raan=',coes[5]) 
        
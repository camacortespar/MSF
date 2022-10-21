import numpy as np

G_meters=6.67408e-11
G=G_meters*10**-9

G1=10.0**8 #kg*km**3/s**2/m**s

#AU to km
AU=1.4959787e8

#days to secs
day2sec=24*3600.0

#Path value
A='/home/wind/Escritorio/Simulacion/Proyecto-Simulaci-n/Orbital_Mechanics_with_python/spice_data/'

sun={
    'name':'sun',
    'mass':1.989e30,
    'mu':1.32712e11,
    'radius':695700.0,
    'G1':G1,
    'spice_file':A+'de432s.bsp',
    'deorbit_altitude':2*695700.0

}
atm=np.array([[63.096,2.059e-4],[251.189,5.909e-11],[1000.0,3.561e-15]])
earth={
    'name':'earth',
    'mass':5.972e24,
    'mu':5.972e24*G,
    'radius':6378.0,
    'J2':1.082635854e-3,
    'zs':atm[:,0], #km
    'rhos': atm[:,1]*10**8, # kg/km³
    'atm_rot_vector':np.array([0.0,0.0,72.9211e-6]), #rad/s
    'deorbit_altitude':100,
    'spice_file':A+'de432s.bsp'
}
moon={
    'name':'moon',
    'mass':7.34797309e22,
    'mu':7.34797309e22*G,
    'radius':1737.1,
    'Orbit_T':29*day2sec+12*3600.0+44*60.0+2.8,
    'dist2earth':384400.0,
    'spice_file':A+'de432s.bsp'
}
moon['orbit_w']=2*np.pi/moon['Orbit_T']

jupiter={
    'name':'jupiter barycenter',
    'mass':1.89813e27,
    'mu':1.89813e27*G,
    'spice_file':A+'de432s.bsp',
}
saturn={
    'name':'saturn barycenter',
    'mass':5.683e26,
    'mu':5.683e26*G,
    'spice_file':A+'de432s.bsp',
}
venus={
    'name':'venus barycenter',
    'mass':4.867e24,
    'mu':4.867e24*G,
    'spice_file':A+'de432s.bsp',
}
mars={
    'name':'mars barycenter',
    'mass':6.39e23,
    'mu':6.39e23*G,
    'spice_file':A+'de432s.bsp',
}
mercury={
    'name':'mercury barycenter',
    'mass':3.285e23,
    'mu':3.285e23*G,
    'spice_file':A+'de432s.bsp',
}
neptune={
    'name':'neptune barycenter',
    'mass':1.024e26,
    'mu':1.024e26*G,
    'spice_file':A+'de432s.bsp',
}
uranus={
    'name':'uranus barycenter',
    'mass':8.681e25,
    'mu':8.681e25*G,
    'spice_file':A+'de432s.bsp',
}

bodies={
'earth':earth,
'sun':sun,
'moon':moon,
'jupiter':jupiter,
'saturn':saturn,
'venus':venus,
'mars':mars,
'mercury':mercury,
'neptune':neptune,
'uranus':uranus,
}

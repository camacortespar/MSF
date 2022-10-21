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

SUN={
    'name':'SUN',
    'mass':1.989e30,
    'mu':1.32712e11,
    'radius':695700.0,
    'G1':G1,
    'spice_file':A+'de432s.bsp',
    'deorbit_altitude':2*695700.0

}
atm=np.array([[63.096,2.059e-4],[251.189,5.909e-11],[1000.0,3.561e-15]])
EARTH={
    'name':'EARTH',
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
MOON={
    'name':'MOON',
    'mass':7.34797309e22,
    'mu':7.34797309e22*G,
    'radius':1737.1,
    'Orbit_T':29*day2sec+12*3600.0+44*60.0+2.8,
    'dist2earth':384400.0,
    'spice_file':A+'de432s.bsp'
}
MOON['orbit_w']=2*np.pi/MOON['Orbit_T']



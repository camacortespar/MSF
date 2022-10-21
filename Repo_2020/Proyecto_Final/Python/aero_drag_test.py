
import planetary_data as pd
import tools as t
from OrbitPropagator import OrbitPropagator as OP
from OrbitPropagator import null_perts 

#time parameters
tspan=3600*24*3
dt=100.0

#central body
cb=pd.earth

if __name__ == '__main__':

	perts=null_perts()
	perts['aero']=True
	perts['Cd']=2.2
	perts['A']=(1e-3)**2/4.0 #km²

	#initial mass
	mass0=10.0 #kg


	#perigee and apogee
	rp=215+cb['radius']
	ra=300+cb['radius']

	#orbital elements angles
	raan=340.0
	i=65.2
	aop=58.0
	ta=332.0

	#other orbitaal elements
	a=(rp+ra)/2.0 #km
	e=(ra-rp)/(ra+rp)

	#initial state value
	state0=[a,e,i,ta,aop,raan]


	op=OP(state0,tspan,dt,degres=True,coes=True,mass0=mass0,perts=perts)

	op.plot_alts(show_plot=True,hours=True)
	op.plot_3d(show_plot=True)
	op.calculate_coes()
	op.plot_coes(show_plot=True,hours=True)
	op.calculate_apoapse_periapse()
	op.plot_apoapse_periapse(show_plot=True,hours=True)
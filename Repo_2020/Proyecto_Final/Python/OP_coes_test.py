
import planetary_data as pd
import tools as t
from OrbitPropagator import OrbitPropagator as OP

tspan=3600*24*1.0
dt=100.0

cb=pd.earth

if __name__ == '__main__':

	#a,e,i,ta,aop,raan=coes

	#ISS
	c0=[cb['radius']+414.0,0.0006189,51.6393,0.0,234.1955,1056372]

	#GEO
	c1=[cb['radius']+35800.0,0.0,0.0,0.0,0.0,0.0]

	#random
	c2=[cb['radius']+3000.0,0.3,20,0.0,15.0,40.0]


	#create orbit propagator
	op0=OP(c0,tspan,dt, coes=True)
	op1=OP(c1,tspan,dt, coes=True)
	op2=OP(c2,tspan,dt, coes=True)

	op0.propagate_orbit()
	op1.propagate_orbit()
	op2.propagate_orbit()

	t.plot_n_orbits([op0.rs,op1.rs,op2.rs],labels=['ISS','GEO', 'random'], show_plot=True)
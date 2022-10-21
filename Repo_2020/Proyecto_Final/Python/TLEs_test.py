
import planetary_data as pd
import tools as t
from OrbitPropagator import OrbitPropagator as OP

#time parameters
tspan=3600*24*1.0
dt=100.0

#central body
cb=pd.earth

if __name__ == '__main__':
	
	#ISS
	op0=OP(t.tle2coes('iss.txt'),tspan,dt,coes=True,degres=False)

	#COSMOS2251
	op1=OP(t.tle2coes('COSMOS2251.txt'),tspan,dt,coes=True,degres=False)

	#HJ-2A
	op2=OP(t.tle2coes('HJ-2A.txt'),tspan,dt,coes=True,degres=False)
	

	t.plot_n_orbits([op0.rs,op1.rs,op2.rs],labels=['ISS','COSMOS2251','HJ-2A'],show_plot=True,title='Orbita de satélites por TLEs')
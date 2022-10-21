
import planetary_data as pd

BODIES=['earth','moon','sun']

if __name__ == '__main__':
	
	for name in BODIES:

		body=pd.bodies[name]

		#calculate scape velocity at 1.1 radius
		v_esc=(2*body['mu']/(1.1*body['radius']))**0.5

		print('%s:\t\t%.5f km/s' %(name,v_esc))




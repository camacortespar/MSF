import math
import numpy as np
import datetime
import spiceypy as spice
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
plt.style.use('dark_background')

import tools as t
import planetary_data as pd
import spice_tools as st

pi=math.pi

sec2day=1/(3600*24)

def lamberts_problem_universal_variables(r0,r1,deltat,tm=1,mu=pd.earth['mu'],tol=1e-6,max_steps=200,psi=0,psi_u=4*pi**2,psi_l=-4*pi):
	
	#calculate square root of mu parameter
	sqrt_mu=math.sqrt(mu)

	#calculate norm of position vectors
	r0_norm=t.norm(r0)
	r1_norm=t.norm(r1)

	#calculate gamma parameter
	gamma=np.dot(r0,r1)/r0_norm/r1_norm

	#calculate beta parameter
	beta=tm*math.sqrt(1-gamma**2)

	#calculate A parameter
	A=tm*math.sqrt(r0_norm*r1_norm*(1+gamma))

	#if A=0, solution cant be calculated
	if A==0:
		return np.array([0,0,0]),np.array([0,0,0])

	#initial values of c2 and c3
	c2=0.5
	c3=1/6.0

	#counter and solved variables
	step=0
	solved=False

	#while tolerance not met and not at max steps
	for n in range(max_steps):

		#calculate B paramenter
		B=r0_norm+r1_norm+A*(psi*c3-1)/math.sqrt(c2)

		#if A and B parameters out of range 
		if A>0 and B<0:

			#increase lower psi value
			psi_l+=pi

			#recalculate B parameter
			B*=-1

		#calculate universal variable cubed
		chi3=math.sqrt(B/c2)**3

		#calculate deltat_ variable
		deltat_=(chi3*c3+A*math.sqrt(B))/sqrt_mu

		#if difference between deltat variables is with in range
		if abs(deltat-deltat_)<tol:

			#set solved variable true
			solved=True

			#break out of the for loop
			break

		#check diferent between deltat and deltat_
		if deltat_<=deltat:

			#adjust lower psi value
			psi_l=psi

		else:

			#adjust upper psi value
			psi_u=psi

		#update psi, c2 and c3 values 
		psi=(psi_u+psi_l)/2.0
		c2=C2(psi)
		c3=C3(psi)
		#print(psi)

	#check if maximun steps reached
	if not solved:

		#algorithm did not converge on a psi value
		print('Lamberts problem universal variables dint converge')
		return np.array([0,0,0])+np.array([0,0,0])

	#calculate coeficients
	f=1-B/r0_norm
	g=A*math.sqrt(B/mu)
	gdot=1-B/r1_norm

	#calculate velocity vectors
	v0=(r1-f*r0)/g
	v1=(gdot*r1-r0)/g

	#print(t.norm(v0),t.norm(v1))
	return v0,v1

def C2(psi):
	return (1-math.cos(math.sqrt(psi)))/psi

def C3(psi):
	return (math.sqrt(psi)-math.sin(math.sqrt(psi)))/(psi*math.sqrt(psi))

def interplanetary_porkchop_plot(config):
	
	_config={

	'planet0':'earth',
	'planet1':'mars barycenter',
	'departure0':'2020-07-01',
	'departure1':'2020-09-01',
	'arrival0':'2020-11-01',
	'arrival1':'2022-01-24',
	'mu':pd.sun['mu'],
	'step':1/sec2day,
	'frame':'ECLIPJ2000',
	'observer':'SOLAR SYSTEM BARYCENTER',
	'cuttof_v':20.0,
	'c3_levels':None,
	'vinf_levels':None,
	'tof_levels':None,
	'dv_levels':None,
	'dv_cmap':'RdPu_r',
	'figsize':(10,20),
	'lw':1.5,
	'title':'Porkchop plot',
	'title_dv':'Porkchop plot',
	'fontsize':15,
	'show':False,
	'filename':None,
	'filename_dv':None,
	'dpi':300,
	'load':False,
	}

	for key in config.keys():
		_config[key] = config[key]

	cutoff_c3=_config['cuttof_v']**2

	#array of departur and arrival times
	et_departures = np.arange(
		spice.utc2et(_config['departure0']),
		spice.utc2et(_config['departure1'])+_config['step'],
		_config['step']	)
	et_arrivals = np.arange(spice.utc2et(_config['arrival0']),
		spice.utc2et(_config['arrival1'])+_config['step'],
		_config['step'])

	#number of days in each array and total combinations
	ds=len(et_departures)
	as_=len(et_arrivals)
	total=ds*as_

	print('Departure days: %i.' %ds)
	print('Arrival days: %i.' %as_)
	print('Total combinations: %i.' %total)

	#create an empty array for c3, v infinity, and tof
	C3_shorts=np.zeros((as_,ds))
	C3_longs=np.zeros((as_,ds))
	V_inf_shorts=np.zeros((as_,ds))
	V_inf_longs=np.zeros((as_,ds))
	tofs=np.zeros((as_,ds))

	#create arrays for indexing the meshgrid
	x=np.arange(ds)
	y=np.arange(as_)

	#for each combination
	for na in y:
		for nd in x:
			
			#state of planet 0 at departure
			state_depart =st.load_ephemeris(_config['planet0'],[et_departures[nd]],_config['frame'],_config['observer'])[0]

			#state of planet1 at arrival
			state_arrival=st.load_ephemeris(_config['planet1'],[et_arrivals[na]],_config['frame'],_config['observer'])[0]

			#calculate time of flight
			tof=et_arrivals[na]-et_departures[nd]

			try:
				v_sc_depart_short, v_sc_arrive_short=lamberts_problem_universal_variables(state_depart[:3],state_arrival[:3],tof,tm=1,mu=_config['mu'])
			except:
				v_sc_depart_short=np.array([1000,1000,1000])
				v_scarrive_short=np.array([1000,1000,1000])

			v_sc_arrive_long=0

			try:
				v_sc_depart_long, v_sc_arrive_long=lamberts_problem_universal_variables(state_depart[:3],state_arrival[:3],tof,tm=-1,mu=_config['mu'])
			except:
				v_sc_depart_long=np.array([1000,1000,1000])
				v_scarrive_long=np.array([1000,1000,1000])

			#calculate C3 values departing
			C3_short=t.norm(v_sc_depart_short - state_depart[3:])**2
			C3_long=t.norm(v_sc_depart_long - state_depart[3:])**2

			#check for out of range values
			if C3_short > cutoff_c3: C3_short = cutoff_c3
			if C3_long > cutoff_c3: C3_long = cutoff_c3

			#calculate v infinity values arriving
			V_inf_short=t.norm(v_sc_arrive_short - state_arrival[3:])
			V_inf_long=t.norm(v_sc_arrive_long - state_arrival[3:])

			#check for out of range values
			if V_inf_short>_config['cuttof_v']: V_inf_short = _config['cuttof_v']
			if V_inf_long>_config['cuttof_v']: V_inf_long = _config['cuttof_v']

			#add values to corresponding arrays
			C3_shorts[na,nd]=C3_short
			C3_longs[na,nd]=C3_long
			V_inf_shorts[na,nd]=V_inf_short
			V_inf_longs[na,nd]=V_inf_long
			tofs[na,nd]=tof

		#print('%i/%i.'%(na,as_))

	tofs*=sec2day

	#total delta v
	dv_shorts=V_inf_shorts
	+np.sqrt(C3_shorts)
	dv_longs=V_inf_longs+np.sqrt(C3_longs)

	#create levels arrays
	if _config['c3_levels'] is None:
		_config['c3_levels'] = np.arange(10,50,2)
	if _config['vinf_levels'] is None:
		_config['vinf_levels'] = np.arange(0,15,1)
	if _config['tof_levels'] is None:
		_config['tof_levels'] = np.arange(100,500,20)
	if _config['dv_levels'] is None:
		_config['dv_levels'] = np.arange(3,20,0.5)

	lw=_config['lw']

	fig,ax=plt.subplots(figsize=_config['figsize'])

	c0=ax.contour(C3_shorts,levels=_config['c3_levels'],colors='m',linewidths=lw)
	c1=ax.contour(C3_longs,levels=_config['c3_levels'],colors='m',linewidths=lw)
	c2=ax.contour(V_inf_shorts,levels=_config['vinf_levels'],colors='deepskyblue',linewidths=lw)
	c3=ax.contour(V_inf_longs,levels=_config['vinf_levels'],colors='deepskyblue',linewidths=lw)
	c4=ax.contour(tofs,levels=_config['tof_levels'],colors='white',linewidths=lw*0.6)

	plt.clabel(c0,fmt= '%i')
	plt.clabel(c1,fmt= '%i')
	plt.clabel(c2,fmt= '%i')
	plt.clabel(c3,fmt= '%i')
	plt.clabel(c4,fmt= '%i')
	plt.plot([0],[0],'m')
	plt.plot([0],[0],'c')
	plt.plot([0],[0],'w')
	plt.legend([r'C3 ($\dfrac{km^2}{s^2}$)',r'$V_{\infinity}\;(\dfrac{km}{s})$',r'Time of Flight (days)'],
		bbox_to_anchor=(1.005,1.01),fontsize=10)

	ax.set_title(_config['title'],fontsize=_config['fontsize'])
	ax.set_ylabel('Arrival (Days Past %s)' %_config['arrival0'],fontsize=_config['fontsize'])
	ax.set_xlabel('Departure (Days Past %s)' %_config['departure0'],fontsize=_config['fontsize'])

	if _config['show']:
		plt.show()
	if _config['filename'] is not None:
		
		plt.savefig(_config['filename'],dpi=_config['dpi'])
		print('Saved',_config['filename'])

	plt.close()

	'''
	Delta v plot
	'''

	fig,ax=plt.subplots(figsize=_config['figsize'])

	c0=ax.contour(dv_shorts,levels=_config['dv_levels'],cmap=_config['dv_cmap'],linewidths=lw)
	c1=ax.contour(dv_longs,levels=_config['dv_levels'],cmap=_config['dv_cmap'],linewidths=lw)
	c2=ax.contour(tofs,levels=_config['tof_levels'],colors='c',linewidths=lw*0.6)

	plt.clabel(c0,fmt= '%.1f')
	plt.clabel(c1,fmt= '%.1f')
	plt.clabel(c2,fmt= '%.1f')

	ax.set_title(_config['title_dv'],fontsize=_config['fontsize'])
	ax.set_ylabel('Arrival (Days Past %s)' %_config['arrival0'],fontsize=_config['fontsize'])
	ax.set_xlabel('Departure (Days Past %s)' %_config['departure0'],fontsize=_config['fontsize'])

	if _config['show']:
		plt.show()
	if _config['filename_dv'] is not None:
		
		plt.savefig(_config['filename_dv'],dpi=_config['dpi'])
		print('Saved',_config['filename_dv'])

	plt.close()

def print_delta_v(v0,vf,state_depart,state_arrival):

	C3=t.norm(v0- state_depart)**2
	Vinf=t.norm(vf- state_arrival)

	delta_v=Vinf+math.sqrt(C3)

	print(delta_v)





















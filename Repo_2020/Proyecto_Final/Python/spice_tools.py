
import spiceypy as spice
import numpy as np

#retrives IDs, name and time converages of all objects in spk file
def get_objects(filename,display=False):
	objects=spice.spkobj(filename)
	ids,names,tcs_sec,tcs_cal=[],[],[],[]
	n=0
	if display:
		print('\nObjects in %s:'% filename)

	for o in objects:
		#id
		ids.append(o)

		#time coverage in seconds since J2000
		tc_sec=spice.wnfetd(spice.spkcov(filename,ids[n]),n)

		#convert time coverage to human readeble
		tc_cal=[spice.timout(f,"YYYY MON DD HR:MN:SC.### (TDB) ::TDB") for f in tc_sec]

		#append time coverages to output list
		tcs_sec.append(tc_sec)
		tcs_cal.append(tc_cal)

		#get name of body
		try:
			#add asociate name to list
			names.append(id2body(o))

		except:
			#called if body name does not exist
			names.append('Unknown Name')

		if display:
			print('id %i\t\tname: %s\t\ttc: %s --> %s' % (ids[-1],names[-1],tc_cal[0],tc_cal[1],))
	return ids,names,tcs_sec,tcs_cal

#return name of body given spice ID
def id2body(id_):
	return spice.bodc2n(id_)

#create time array for given time coverage
def tc2array(tcs,steps):
	arr=np.zeros((steps,1))
	arr[:,0]=np.linspace(tcs[0],tcs[1],steps)
	return arr

#get ephemerides data
def get_ephemerides_data(target,times,frame,observer):
	return np.array(spice.spkezr(target,times,frame,'NONE',observer)[0])

def load_ephemeris(target, times, frame, observer):
	
	if type( target )==str:
		return np.array(spice.spkezr(target,times,frame,'NONE',observer)[0])
	else:
		return np.array(spice.spkez(target,times,frame,'NONE',observer)[0])




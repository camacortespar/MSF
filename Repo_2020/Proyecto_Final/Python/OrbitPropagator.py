
import numpy as np
import math as m
import matplotlib.pyplot as plt
from scipy.integrate import ode
import spiceypy as spice
from time import time as time
from multiprocessing import Pool
from mpl_toolkits.mplot3d import Axes3D
plt.style.use('dark_background')

import planetary_data as pd
import tools as _t
import spice_tools as st
import spice_data as sd

d2r=np.pi/180
r2d=180/np.pi

def null_perts():
    return {
    'J2':False,
    'aero':False,
    'Cd':0,
    'A':0,
    'mu':0,
    'moon gravity':False,
    'solar_gravity':False,
    'thrust':0,
    'thrust_direction':0,
    'isp':0,
    'n_bodies':[],
    'srp':False,
    'CR':0,
    'A_srp':0,
    'C20':0,
    'BSTAR':False,
    'B*':0
    }


class OrbitPropagator:
    
    def __init__(self,state0,tspan,dt,date0='2020-02-03',Frame='J2000',coes=False,degres=True,cb=pd.earth, perts=null_perts(),mass0=0,sc={},propagator='lsoda',start=True):
        
        if coes:
            self.r0,self.v0=_t.coes2rv(state0,degres=degres,mu=cb['mu'])
        else:
            self.r0 = np. array (state0 [: 3])
            self.v0 = np. array (state0 [3:])

        #store parameters pased in
        self.dt=dt
        self.cb=cb
        self.mass=mass0
        self.propgator=propagator
        self.frame=Frame
        self.date0=date0
        self.perts=perts
        self.tspan=tspan


        #miscelaneus parameters
        self.n_steps=int(np.ceil(self.tspan/self.dt))+1
        self.ts=np.zeros((self.n_steps+1,1))
        self.y=np.zeros((self.n_steps+1,7))
        self.alts=np.zeros((self.n_steps+1,1))
        self.step=0
        self.spice_file_loaded=[]

        #initial conditions 
        self.y[0,:]=self.r0.tolist()+self.v0.tolist()+[self.mass]
        self.alts[0]=_t.norm(self.r0)-self.cb['radius']
        
        #initiate solver
        self.solver=ode(self.diffy_q)
        self.solver.set_integrator(propagator)
        self.solver.set_initial_value(self.y[0,:],0)

        #define stop conditions dictionary
        self.stop_conditions_dict=sc

        #define dictionary to check internal method
        self.stop_conditions_map={'max_alt':self.check_max_alt,'min_alt':self.check_min_alt}

        #create stop conditions function list with deorbit check
        self.stop_conditions_functions=[self.check_deorbit]

        #fill in the rest of stop conditions
        for key in self.stop_conditions_dict.keys():
            self.stop_conditions_functions.append(self.stop_conditions_map[key])

        #define perturbations dictionary
        self.perts=perts

        #check if loading spice data
        if self.perts['n_bodies'] or self.perts['srp']:

            #load leap seconds kernel
            spice.furnsh(sd.leap_seconds_kernel)

            #add to list of loaded files
            self.spice_file_loaded.append(sd.leap_seconds_kernel)

            #convert start date to seconds after J2000
            self.start_time=spice.utc2et(self.date0)

            #create time span in seconds after J2000
            self.spice_tspan=np.linspace(self.start_time,self.start_time+self.tspan,self.n_steps)

            #if srpget states of the sun
            if self.perts['srp']:

                #load spice file for given body
                spice.furnsh(self.cb['spice_file'])

                #add spice file to list
                self.spice_file_loaded.append(self.cb['spice_file'])

                #calculate central body states throughout entire propagation WRT sun
                self.cb['states']=st.get_ephemerides_data(self.cb['name'],self.spice_tspan,self.frame,'SUN')

        #load kernels for each body
        for body in self.perts['n_bodies']:

            #if spice hasn't all ready been loaded
            if body['spice_file'] not in self.spice_file_loaded:

                #load spice file
                spice.furnsh(body['spice_file'])

                #add to spice file list
                self.spice_file_loaded.append(body['spice_file'])

            #calculate bodi states WRT central body
            body['states']=st.get_ephemerides_data(body['name'],self.spice_tspan,self.frame,self.cb['name'])

        if start:
            #propagate orbit
            self.propagate_orbit()

    def check_deorbit(self):
        if self.alts[self.step]<self.cb['deorbit_altitude']:
            print('Spacecraft has deorbited after %.1f second' % self.ts[self.step])
            return False
        return True

    def check_max_alt(self):
        if self.alts[self.step]<self.stop_conditions_dict['max_alt']:
            print('Spacecraft reached max altitude after %.1f second' %self.ts[self.step])
            return False
        return True

    def check_min_alt(self):
        if self.alts[self.step]<self.stop_conditions_dict['min_alt']:
            print('Spacecraft reached min altitude after %.1f second' %self.ts[self.step])
            return False
        return True

    def check_escape_velocity(self):
        if _t.esc_v(_t.norm(self.y[step,:3]),mu=self.cb['mu'])<_t.norm(self.y[self.step,3:6]):
            self.print_stop_conditions('scape_velocity')
            return False
        return True

    def print_stop_conditions(self,parameter):

        print('Spacecraft has reached %s after %2.f days (%2.f hours or %2.f seconds)' %(parameter,self.ts[self.step]*days))

    def check_stop_conditions(self):
        #for each stop condition
        for sc in self.stop_conditions_functions:

            #if return False
            if not sc():

                #stop condition reached, return false
                return False

        #if stop condition is not reached return true
        return  True
    
    def propagate_orbit(self):

        print('Propagating orbit...')
   
        #propagate the orbit
        while self.solver.successful() and self.step<self.n_steps and self.check_stop_conditions():

            #integrate step
            self.solver.integrate(self.solver.t+self.dt)
            self.step+=1

            #extract values from solver
            self.ts[self.step]=self.solver.t
            self.y[self.step]=self.solver.y

            #calculate altitude at these time step
            self.alts[self.step]=_t.norm(self.solver.y[:3])-self.cb['radius']

        #exctract array after the propagation
        self.ts=self.ts[:self.step]    
        self.rs=self.y[:self.step,:3]
        self.vs=self.y[:self.step,3:6]
        self.masses=self.y[:self.step,-1]
        self.alts=self.alts[:self.step]
        
        
    def diffy_q(self,t,y):
    
        #unpack state
        rx,ry,rz,vx,vy,vz,mass=y
        r=np.array([rx,ry,rz])
        v=np.array([vx,vy,vz])
        norm_r=_t.norm(r)
        dmdt=0
        
        #two body acceleration
        a=-r*self.cb['mu']/norm_r**3

        #J2 perturbation
        if self.perts['J2']:
            z2=r[2]**2
            r2=norm_r**2
            tx=r[0]/norm_r*(5*z2/r2-1)
            ty=r[1]/norm_r*(5*z2/r2-1)
            tz=r[2]/norm_r*(5*z2/r2-3)

            a_j2=1.5*self.cb['J2']*self.cb['mu']*self.cb['radius']**2/norm_r**4*np.array([tx,ty,tz])

            a+=a_j2

        #aero drag perturbation
        if self.perts['aero']:

            if self.perts['BSTAR']:
                #calculate altitude and air density
                z=norm_r-self.cb['radius']
                rho=_t.calc_atmospheric_density(z)

                #calculate motion of s/c with respect to rotating atmosphere
                v_rel=v-np.cross(self.cb['atm_rot_vector'],r)

                drag=-v_rel*np.linalg.norm(v_rel)*self.perts['B*']/self.cb['radius']

                a+=drag

            else:

                #calculate altitude and air density
                z=norm_r-self.cb['radius']
                rho=_t.calc_atmospheric_density(z)

                #calculate motion of s/c with respect to rotating atmosphere
                v_rel=v-np.cross(self.cb['atm_rot_vector'],r)

                drag=-v_rel*0.5*rho*np.linalg.norm(v_rel)*self.perts['Cd']*self.perts['A']/self.mass

                a+=drag


        #thust perturbation
        if self.perts['thrust']:

            #thrust vector
            a+=self.perts['thrust_direction']*_t.normed(v)*self.perts['thrust']/mass/1000.0  #km/s²

            #derivate of totl mass
            dmdt=-self.perts['thrust']/self.perts['isp']/9.81

        #n body pertsurbations
        for body in self.perts['n_bodies']:

            #vector pointing from satellite to body
            r_cb2nb=body['states'][self.step,:3]

            #vector pointing from satellite to body
            r_sat2body=r_cb2nb-r

            #nth body acceleration vector
            a+=body['mu']*(r_sat2body/_t.norm(r_sat2body)**3-r_cb2nb/_t.norm(r_cb2nb)**3)


        #solar radiation preasure
        if self.perts['srp']:

            #vector pointing from sun to satellite
            r_sun2sc=self.cb['states'][self.step,:3]+r

            #srp vector from given ephemerides data 
            srp=((1+self.perts['CR'])*pd.sun['G1']*self.perts['A_srp'])/(mass*_t.norm(r_sun2sc)**3*r_sun2sc)

            a+=srp
        
        return [vx,vy,vz,a[0],a[1],a[2],dmdt]

    def calculate_coes(self,degres=True,print_results=False,parallel=False):

        print('Calculating COEs....')

        if parallel:
            start=time()
            states=_t.Parallel_states(self.y)
            pool=Pool()
            self.coes=np.array(pool.starmap(_t.rv2coes,states))
            self.coes_rel=self.coes-self.coes[0,:]
            print(time()-start)

        else:
            start=time()

            #preallocate arrays for coes
            self.coes=np.zeros((self.step,6))
            self.coes_rel=np.zeros((self.step,6))

        #fill array
        for n in range(self.step):
            self.coes[n,:]=_t.rv2coes(self.y[n,:6],mu=self.cb['mu'],degres=degres,print_results=print_results)

        self.coes_rel=self.coes-self.coes[0,:]

    def plot_coes(self,hours=False,days=False,show_plot=False,save_plot=False,title='COEs',figsize=(20,12),dpi=500,rel=True):
        print('Ploting COEs ...')

        #create figure and axses
        fig,axs=plt.subplots(nrows=2,ncols=3,figsize=figsize)

        #figure titles
        fig.suptitle(title,fontsize=20)

        #x axis
        if hours: 
            ts=self.ts/3600
            xlabel='Time (hours)'
        elif days:
            ts=self.ts/3600/24
            xlabel='Time (days)'
        else:
            ts=self.ts
            xlabel='Time (seconds)'

        #check to plote relative coes
        if rel:
            coes=self.coes_rel
        else:
            coes=self.coes

        #plot true anomaly
        axs[0,0].plot(ts,self.coes[:,3])
        axs[0,0].set_title('True anomaly vs.Time')
        axs[0,0].grid(True)
        axs[0,0].set_ylabel('Angle (degrees)')

        #plot semi major axis
        axs[1,0].plot(ts,self.coes[:,0])
        axs[1,0].set_title('Semi Major Axis vs.Time')
        axs[1,0].grid(True)
        axs[1,0].set_ylabel('Semi Major Axis (km)')
        axs[1,0].set_xlabel(xlabel)

        #plot eccentricity
        axs[0,1].plot(ts,self.coes[:,1])
        axs[0,1].set_title('Eccentricity vs.Time')
        axs[0,1].grid(True)

        #plot argument of periage
        axs[0,2].plot(ts,self.coes[:,4])
        axs[0,2].set_title('Argument of Periapse vs.Time')
        axs[0,2].grid(True)

        #plot inclination
        axs[1,1].plot(ts,self.coes[:,2])
        axs[1,1].set_title('Inclination vs.Time')
        axs[1,1].grid(True)
        axs[1,1].set_ylabel('Angle (degrees)')
        axs[1,1].set_xlabel(xlabel)

        #plot RAAN
        axs[1,2].plot(ts,self.coes[:,5])
        axs[1,2].set_title('RAAN vs.Time')
        axs[1,2].grid(True)
        axs[1,2].set_xlabel(xlabel)

        plt.subplots_adjust(wspace=0.3)
        plt.subplots_adjust(hspace=0.4)

        if show_plot:
            plt.show()
        if save_plot:
            plt.savefig(title+'.png',dpi=dpi)

    def plot_3d(self,show_plot=False,save_plot=False,title='Test Title'):
        
        #3D plot
        fig=plt.figure(figsize=(16,8))
        ax=fig.add_subplot(111,projection='3d')
        
        #plot trayectory and starting point
        ax.plot(self.rs[:,0],self.rs[:,1],self.rs[:,2],'r',label='Trayectory',zorder=10)
        ax.plot([self.rs[0,0]],[self.rs[0,1]],[self.rs[0,2]],'wo',label='Initial position',zorder=10)
               
        #plot earth
        _u,_v=np.mgrid[0:2*np.pi:20j,0:np.pi:10j]
        _x=self.cb['radius']*np.cos(_u)*np.sin(_v)
        _y=self.cb['radius']*np.sin(_u)*np.sin(_v)
        _z=self.cb['radius']*np.cos(_v)
        ax.plot_surface(_x,_y,_z,cmap='Blues',zorder=0)
        
        l=self.cb['radius']*2.0
        x,y,z=[[0,0,0],[0,0,0],[0,0,0]]
        u,v,w=[[l,0,0],[0,l,0],[0,0,l]]
        ax.quiver(x,y,z,u,v,w,color='k')
        
        #check for custom axes limits
        max_val=np.max(np.abs(self.rs))
        
        #set lables and title
        ax.set_xlim([-max_val,max_val])
        ax.set_ylim([-max_val,max_val])
        ax.set_zlim([-max_val,max_val])
        ax.set_xlabel('X (km)')
        ax.set_ylabel('Y (km)')
        ax.set_zlabel('Z (km)')
        ax.set_title(title)
        
        plt.legend()

        if show_plot:
            plt.show()
        if save_plot:
            plt.savefig(title+'.png',dpi=300)

    def calculate_apoapse_periapse(self):
        self.apoapse=self.coes[:,0]*(1+self.coes[:,1])
        self.periapse=self.coes[:,0]*(1-self.coes[:,1])

    def plot_apoapse_periapse(self,hours=False,days=False,show_plot=False,save_plot=False,title='Apoapse and Periapse',dpi=500):
        #create figure
        plt.figure(figsize=(20,10))

        if hours:
            ts=self.ts/3600.0
            x_unit='Hours'
        elif days:
            ts=self.ts/3600/24
            x_unit='Days'
        else:
            ts=self.ts
            x_unit='Seconds'

        #plot each
        plt.plot(ts,self.apoapse,'b',label='Apoapse')
        plt.plot(ts,self.periapse,'r',label='Periapse')

        #labels
        plt.xlabel('Time(%s)'% x_unit)
        plt.ylabel('Altitude (km)')

        #other parameters
        plt.grid(True)
        plt.title(title)
        plt.legend()

        if show_plot:
            plt.show()
        if save_plot:
            plt.savefig(title+'.png',dpi=dpi)

    def plot_alts(self,show_plot=False,save_plot=False,hours=False,days=False,title='Radial Distance vs. Time',figsize=(16,8),dpi=500):
        #create figure
        if hours:
            ts=self.ts/3600.0
            x_unit='Hours'
        elif days:
            ts=self.ts/3600/24
            x_unit='Days'
        else:
            ts=self.ts
            x_unit='Seconds'

        plt.figure(figsize=figsize)
        plt.plot(ts,self.alts,'w')
        plt.grid(True)
        plt.xlabel('Time(%s)'% x_unit)
        plt.ylabel('Altitude (km)')
        plt.title(title)
        if show_plot:
            plt.show()
        if save_plot:
            plt.savefig(title+'.png',dpi=dpi)

    def printcoes(self):
        #x,y,z,vx,vy,vz,m=self.y[-1,:]

        #print(x,y,z,vx,vy,vz)
        rf=self.coes[-1,:]
        print('a=',rf[0]) 
        print('e=',rf[1]) 
        print('i=',rf[2]) 
        print('ta=',rf[3]) 
        print('aop=',rf[4]) 
        print('raan=',rf[5]) 
            

        
            
            

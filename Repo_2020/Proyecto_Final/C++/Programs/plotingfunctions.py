import numpy as np
import matplotlib.pyplot as plt
import math as m
from mpl_toolkits.mplot3d import Axes3D
plt.style.use('dark_background')

import planetary_data as pd

def plot_n_orbits(rs,labels,cb=pd.EARTH, show_plot=False,save_plot=False,title='Many Orbits'):
    print('ploting orbits...')
    #3D plot
    fig=plt.figure(figsize=(300,300))
    ax=fig.add_subplot(111,projection='3d')
    
    #plot trayectory and starting point
    
    n=0
    for r in rs:
        ax.plot(r[:,0],r[:,1],r[:,2],label=labels[n],zorder=5)
        #ax.plot([r[0,0]],[r[0,1]],[r[0,2]],label='Initial position')
        n+=1
    
    #plot earth
    _u,_v=np.mgrid[0:2*np.pi:20j,0:np.pi:10j]
    _x=cb['radius']*np.cos(_u)*np.sin(_v)
    _y=cb['radius']*np.sin(_u)*np.sin(_v)
    _z=cb['radius']*np.cos(_v)
    ax.plot_surface(_x,_y,_z,cmap='Blues',zorder=0)
    
    l=cb['radius']*2.0
    x,y,z=[[0,0,0],[0,0,0],[0,0,0]]
    u,v,w=[[l,0,0],[0,l,0],[0,0,l]]
    ax.quiver(x,y,z,u,v,w,color='k')
    
    #check for custom axes limits
    max_val=np.max(np.abs(rs))
        
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
        plt.savefig(title+'.png',dpi=200)

def plot_coes(coes,hours=False,days=False,show_plot=False,save_plot=False,title='COEs'):
        print('Ploting COEs ...')
        ts=coes[:,6]
        #create figure and axses
        fig,axs=plt.subplots(nrows=2,ncols=3,figsize=(20,12))

        #figure titles
        fig.suptitle(title,fontsize=20)

        #x axis
        if hours: 
            ts=ts/3600
            xlabel='Time (hours)'
        elif days:
            ts=ts/3600/24
            xlabel='Time (days)'
        else:
            ts=ts
            xlabel='Time (seconds)'

        #check to plote relative coes


        #plot true anomaly
        axs[0,0].plot(ts,coes[:,3])
        axs[0,0].set_title('True anomaly vs.Time')
        axs[0,0].grid(True)
        axs[0,0].set_ylabel('Angle (degrees)')

        #plot semi major axis
        axs[1,0].plot(ts,coes[:,0])
        axs[1,0].set_title('Semi Major Axis vs.Time')
        axs[1,0].grid(True)
        axs[1,0].set_ylabel('Semi Major Axis (km)')
        axs[1,0].set_xlabel(xlabel)

        #plot eccentricity
        axs[0,1].plot(ts,coes[:,1])
        axs[0,1].set_title('Eccentricity vs.Time')
        axs[0,1].grid(True)

        #plot argument of periage
        axs[0,2].plot(ts,coes[:,4])
        axs[0,2].set_title('Argument of Periapse vs.Time')
        axs[0,2].grid(True)

        #plot inclination
        axs[1,1].plot(ts,coes[:,2])
        axs[1,1].set_title('Inclination vs.Time')
        axs[1,1].grid(True)
        axs[1,1].set_ylabel('Angle (degrees)')
        axs[1,1].set_xlabel(xlabel)

        #plot RAAN
        axs[1,2].plot(ts,coes[:,5])
        axs[1,2].set_title('RAAN vs.Time')
        axs[1,2].grid(True)
        axs[1,2].set_xlabel(xlabel)

        plt.subplots_adjust(wspace=0.3)
        plt.subplots_adjust(hspace=0.4)

        if show_plot:
            plt.show()
        if save_plot:
            plt.savefig(title+'.png',dpi=500)

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
            ts=ts/3600.0
            x_unit='Hours'
        elif days:
            ts=ts/3600/24
            x_unit='Days'
        else:
            ts=ts
            x_unit='Seconds'

        plt.figure(figsize=figsize)
        plt.plot(ts,alts,'w')
        plt.grid(True)
        plt.xlabel('Time(%s)'% x_unit)
        plt.ylabel('Altitude (km)')
        plt.title(title)
        if show_plot:
            plt.show()
        if save_plot:
            plt.savefig(title+'.png',dpi=dpi)
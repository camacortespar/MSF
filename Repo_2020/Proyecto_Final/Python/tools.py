import numpy as np
import datetime
import spiceypy as spice
import matplotlib.pyplot as plt
import math as m
from mpl_toolkits.mplot3d import Axes3D
plt.style.use('dark_background')

import planetary_data as pd
from OrbitPropagator import OrbitPropagator as OP

d2r=np.pi/180
r2d=180/np.pi
pi=m.pi

def norm(v):
    return np.linalg.norm(v)

def normed(v):
    return np.array(v)/norm(v)

def plot_n_orbits(rs,labels,cb=pd.earth, show_plot=False,save_plot=False,title='Many Orbits'):

    print('ploting orbits...')
    
    #3D plot
    fig=plt.figure(figsize=(20,20))
    ax=fig.add_subplot(111,projection='3d')
    
    #plot trayectory and starting point
    
    n=0
    max_val=0

    for r in rs:
        ax.plot(r[:,0],r[:,1],r[:,2],label=labels[n],zorder=5)
        #ax.plot([r[0,0]],[r[0,1]],[r[0,2]],label='Initial position')
        max_val = max( [ r.max(), max_val ] )
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

#convert classical orbital elements to r and v vectors
def coes2rv(coes,degres=False,mu=pd.earth['mu']):
    if degres:
        a,e,i,ta,aop,raan=coes
        i*=d2r
        ta*=d2r
        aop*=d2r
        raan*=d2r
    else:
        a,e,i,ta,aop,raan=coes

    E=ecc_anomaly([ta,e],'tae')

    r_norm=a*(1-e**2)/(1+e*np.cos(ta))

    #calculate r and v vectors in perifocal frame
    r_perif=r_norm*np.array([m.cos(ta),m.sin(ta),0])
    v_perif=m.sqrt(mu*a)/r_norm*np.array([-m.sin(E),m.cos(E)*m.sqrt(1-e**2),0])

    #rotation matrix for perifocal to ECI
    perif2eci=np.transpose(eci2perif(raan,aop,i))

    #calculate r and v in inertial frame
    r=np.dot(perif2eci,r_perif)
    v=np.dot(perif2eci,v_perif)

    return r,v

#convert r and v vectors to coes(classical orbital elements) 
def rv2coes(state,et=0,mu=pd.earth['mu'],degres=False,print_results=False):
    #calculate orbital elements for given state
    rp,e,i,raan,aop,ma,t0,mu,ta,a,T=spice.oscltx(state,et,mu)


    if degres:
        i*=r2d
        ta*=r2d
        aop*=r2d
        raan*=r2d
    if print_results:
        print(a,e,i,ta,aop,raan)

    return a,e,i,ta,aop,raan
    
#inertial to perifocal rotation matrix
def eci2perif(raan,aop,i):
    row0=[-m.sin(raan)*m.cos(i)*m.sin(aop)+m.cos(raan)*m.cos(aop),m.cos(raan)*m.cos(i)*m.sin(aop)+m.sin(raan)*m.cos(aop),m.sin(i)*m.sin(aop)]
    row1=[-m.sin(raan)*m.cos(i)*m.cos(aop)-m.cos(raan)*m.sin(aop),m.cos(raan)*m.cos(i)*m.cos(aop)-m.sin(raan)*m.sin(aop),m.sin(i)*m.cos(aop)]
    row2=[m.sin(raan)*m.sin(i),-m.cos(raan)*m.sin(i),m.cos(i)]
    return np.array([row0,row1,row2])
    
#calculate excentricity anomaly
def ecc_anomaly(arr,method,tol=1e-8):
    if method=='newton':
        #newton's method for iteratively finding E
        Me,e=arr
        if Me<np.pi/2.0: E0=Me+e/2.0
        else: E0=Me-e
        for n in range(200): #arbitrary max numberof steps
            ratio=(E0-e*np.sin(E0)-Me)/(1-e*np.cos(E0));
            if abs(ratio)<tol:
                if n==0: return E0
                else: return E1
            else:
                E1=E0-ratio
                E0=E1
            #did not converge
            return False
    elif method=='tae':
        ta,e=arr
        return 2*m.atan(m.sqrt((1-e)/(1+e))*m.tan(ta/2.0))
    else:
        print('Invalid method for eccentric anomaly')

#takes text file containing TLEs and return classical orbital elements and other parameters  
def tle2coes(tle_filename,mu=pd.earth['mu'], degres=False):
    #read tle file
    with open(tle_filename,'r') as f:
        lines=f.readlines()

    #separate into 3 lines
    line0=lines[0].strip() #name of satelite
    line1=lines[1].strip().split()
    line2=lines[2].strip().split()



    #epoch (year and day)
    epoch=line1[3]
    year,month,day,hour=calc_epoch(epoch)

    #collect coes
    if degres:
        #inclination
        i=float(line2[2])
        #right ascention of ascending node
        raan=float(line2[2])
        #eccentricity
        e_string=line2[4]
        e=float('0.' +e_string)
        #argument of perigee
        aop=float(line2[5])
        #mean anomaly
        Me=float(line2[6])
        #mean motion
        mean_motion=float(line2[7]) #revs/day
        #period
        T=1/mean_motion*24*3600 #seconds
        #semi mjor axis
        a=(T**2*mu/4.0/np.pi**2)**(1/3.0)

        #calculate eccentric anomaly
        E=ecc_anomaly([Me,e],'newton')

        #calculate the true anomaly
        ta=true_anomaly([E,e])

        print(a,e,i,ta,aop,raan)

        return a,e,i,ta,aop,raan    

    else:

        #inclination
        i=float(line2[2])*d2r #rad
        #right ascention of ascending node
        raan=float(line2[2])*d2r #rad
        #eccentricity
        e_string=line2[4]
        e=float('0.' +e_string)
        #argument of perigee
        aop=float(line2[5])*d2r #rad
        #mean anomaly
        Me=float(line2[6])*d2r #rad
        #mean motion
        mean_motion=float(line2[7]) #revs/day
        #period
        T=1/mean_motion*24*3600 #seconds
        #semi mjor axis
        a=(T**2*mu/4.0/np.pi**2)**(1/3.0)

        #calculate eccentric anomaly
        E=ecc_anomaly([Me,e],'newton')

        #calculate the true anomaly
        ta=true_anomaly([E,e])
        
        print(a,e,i,ta,aop,raan)

        return a,e,i,ta,aop,raan

def calc_epoch(epoch):
    #year
    year=int('20'+epoch[:2])

    epoch=epoch[2:].split('.')

    #day of the year
    day_of_year=int(epoch[0])-1

    #decimal hour of day
    hour=float('0.'+epoch[1])*24.0

    #get year-month-day
    date=datetime.date(year,1,1)+datetime.timedelta(day_of_year)

    #exact month and date
    month=float(date.month)
    day=float(date.day)

    return year,month,day,hour

def true_anomaly(arr):
    E,e=arr
    return 2*np.arctan(np.sqrt((1+e)/(1-e)))*np.tan(E/2.0)

def tle2rv(tle_filename):

    return coes2rv(tle2coes(tle_filename))

#calculate atmosferic density from given altitud
def calc_atmospheric_density(z):
    rhos,zs=find_rho_z(z)
    if rhos[0]==0: return 0.0

    Hi=-(zs[1]-zs[0])/m.log(rhos[1]/rhos[0])

    return rhos[0]*m.exp(-(z-zs[0])/Hi)

#find endpoints of altitudes and density surrounding imput altitude
def find_rho_z(z,zs=pd.earth['zs'],rhos=pd.earth['rhos']):
    if not 1.0<z<1000.0:
        return [[0.0,0.0],[0.0,0.0]]

    #find the two point surrounding the given input altitude
    for n in range(len(rhos)-1):
        if zs[n]<z<zs[n+1]:
            return [[rhos[n],rhos[n+1]],[zs[n],zs[n+1]]]

    #if out of range return 0
    return [[0.0,0.0],[0.0,0.0]]

def esc_v(r,mu=pd.earth['mu']):
    return m.sqrt(2*mu/r)

def hohmann_transfer(r0=0.0,r1=0.0,cb=pd.earth,coes0=[],coes1=[],dt=1,altitude=True,output_dir='',names=['Initial','Final','Transfer'],write_output=False,propagate=False):

    #check if coes passed in
    if coes0:

        r0=coes0[0]
        r1=coes1[0]

    #if passin in altitude, not semi mayor axixs
    elif altitude:

        r0+=cb['radius']
        r1+=cb['radius']

    #calculate semi mayor axis of transfer orbit
    a_transfer=(r0+r1)/2.0

    #calculate velocities of circular orbit
    v_circ_init=m.sqrt(cb['mu']/r0)
    v_circ_final=m.sqrt(cb['mu']/r1)

    #calculate eliptical orvit velocities (vis-viva equation)
    v0_transfer=m.sqrt(cb['mu']*(2/r0 -1/a_transfer))
    v1_transfer=m.sqrt(cb['mu']*(2/r1 -1/a_transfer))

    #calculate transfer time (half period)
    t_transfer=pi*m.sqrt(a_transfer**3 /cb['mu'])

    #calculate delta v values
    delta_vs=[v0_transfer- v_circ_init, v_circ_final- v1_transfer]

    if propagate:

        #if coes not pased in
        if not coes0:
            coes0=[r0,0.0,0.0,0.0,0.0,0.0]
            coes1=[r1,0.0,0.0,0.0,0.0,0.0]

        #calculate excentricity of transfer orbit
        e_transfer=1-r0/a_transfer

        #coes for transfer orbit
        coes_transfer=[a_transfer,e_transfer,coes0[2],0.0,coes0[4],coes0[5]]

        #calculate periods of initial and final orbit
        T0=pi*(r0**3 /cb['mu'])**0.5
        T1=pi*(r1**3 /cb['mu'])**0.5

        #create space craft instance and propagate orbit
        sc0=OP(coes0,T0,dt,coes=True)
        sc1=OP(coes1,T1,dt,coes=True)
        sc2=OP(coes_transfer,t_transfer,dt,coes=True)

        return sc0,sc1,sc2,delta_vs,t_transfer

    return delta_vs, t_transfer







        
    
        

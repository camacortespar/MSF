import spiceypy as spice
import os
import datetime
from dateutil.relativedelta import relativedelta
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

#Define our kernels path
kernelsPath = os.path.join(os.path.dirname(os.getcwd()), "kernels")

#Define the path to the leap second kernel
leapSecondsKernelPath = os.path.join(kernelsPath, "naif0012.tls")

#Load the kernel
spice.furnsh(leapSecondsKernelPath)

#Define the path to the solar system ephemeris kernel
solarSystemEphermisKernelPath = os.path.join(kernelsPath, "de405.bsp")

#Load the kernel
spice.furnsh(solarSystemEphermisKernelPath)


#Set our 
referenceFrame = "J2000"
target = "VENUS"
observer = "EARTH"

import datetime
from dateutil.relativedelta import relativedelta

#Lets make an array of 52 times over the next year
nowTime = datetime.datetime.now()
times = [nowTime]
for i in range(1, 365):
    times.append(nowTime + relativedelta(days=i))

#Convert to et
etTimes = [spice.str2et(str(x)) for x in times]

#Here we can use the same target/referenceFrame/observer variables set above
#Note : I'm using [0] to remove the light-time as we don't require it.
positions = spice.spkpos(target, etTimes, referenceFrame,'NONE', observer)[0]

#We'll define a few helper functions here so we can (spoiler alert) use them later on!
def makePlanetPositionPlot(title, size=8, axisLimit=1.5e+08, sunSize=400):

    #Make the figure
    fig = plt.figure(figsize=(size, size))

    #Make sub plot
    axis3d = fig.add_subplot(111, projection='3d')

    #Set axis limits
    axis3d.set_xlim([-axisLimit, axisLimit])
    axis3d.set_ylim([-axisLimit, axisLimit])
    axis3d.set_zlim([-axisLimit, axisLimit])

    #Set axis labels
    axis3d.set_xlabel('X (km)')
    axis3d.set_ylabel('Y (km)')
    axis3d.set_zlabel('Z (km)')

    #Create the sun
    axis3d.scatter([0.0], [0.0], [0.0], s=sunSize, c="orange")

    #Add a title
    plt.title(title, y=1.025)
    
    #Return the plt
    return axis3d, fig

def addPlanetToPlot(axis3D, planetPosition, size, color):
    return axis3D.scatter([planetPosition[0]], [planetPosition[1]], [planetPosition[2]], s=size, c=color)

	#Make the plot
	axis3D, fig = makePlanetPositionPlot('Mercury position relative to Earth over 365 days')

	#Add the data
	for position in positions:
    addPlanetToPlot(axis3D, position, 40, 'blue')

	#Show the plot
	
	plt.show()

# Clean up the kernels
spice.kclear()



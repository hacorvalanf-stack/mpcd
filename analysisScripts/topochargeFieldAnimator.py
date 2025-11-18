"""
	NIAMH-MPCD
	Topological charge rendering script

	Uses shendrukGroupStyle for formatting (from https://github.com/Shendruk-Lab/MPCDDefectLoader)
	Must install it or remove calls to shendrukGroupFormat

	Originally from Tyler N. Shendruk
"""

from pylab import *
from numpy import ma
from subprocess import call
from mpl_toolkits.mplot3d import axes3d
from mpl_toolkits.axes_grid1 import make_axes_locatable
import os
import json
import argparse

###########################################################
### Set up argparse
###########################################################
parser = argparse.ArgumentParser(description='Topological charge field rendering script.')
parser.add_argument("directorData", type=str, help="Path to the director data")
parser.add_argument("topoData", type=str, help="Path to the topological charge data")
parser.add_argument("inputname", type=str, help="Path to input .json file")
parser.add_argument('-s',"--start", type=int, default=0, help="Starting timestep for averaging")
parser.add_argument('-f',"--finish", type=int, default=9999999, help="Finishing timestep for averaging")
parser.add_argument('-a',"--avdim", type=str, default='z', help="Dimension to 'slice' over")
parser.add_argument("--qx", type=int, help="Only show every qx arrow in x", default=1)
parser.add_argument("--qy", type=int, help="Only show every qy arrow in y", default=1)
parser.add_argument("-l", "--length", type=float, help="Length of director lines", default=0.95)
parser.add_argument("-d", "--deflength", type=float, help="Length of nematic pointer lines", default=1.25)
parser.add_argument("-A", "--myAspect", type=str, help="'auto' or 'equal'", default="equal")
parser.add_argument("-k", "--keepFrames", type=int, help="0=don't keep (delete) frames; 1=keep frames", default=0)
parser.add_argument("-c", "--colourbar", type=int, help="1 adds a colourbar to the plot, 0 for none", default=1)
parser.add_argument("-g", "--groupStyle", type=int,
                    help="0=don't use group style; 1=use group style",
                    default=1)
args = parser.parse_args()

###########################################################
### Read arguments
###########################################################
print("Arguments:")
for arg, value in vars(args).items():
	print(f"\t{arg}: {value}")
directorData = args.directorData
topoData = args.topoData
inputName = args.inputname
start = args.start
finish = args.finish
qx = args.qx
qy = args.qy
avdim = args.avdim
myLength = args.length
defLength = args.deflength
myAspect = args.myAspect
keepFrames = args.keepFrames
cbFlag = args.colourbar
groupStyle = args.groupStyle

makeTransparent = False # Transparent backgrounds make crappy videos, but look good on webpages

###########################################################
### Format and style
###########################################################
if groupStyle:
	# Use our custom style and colours
	plt.style.use('shendrukGroupStyle')
	import shendrukGroupFormat as ed
# Colour map to use
if groupStyle:
	myMap=ed.plasma
else:
	myMap=plt.cm.plasma
# Adjust line width
myLW=1.0
# adjust line length
myLength *= sqrt(min(qx, qy))
defLength *= sqrt(min(qx, qy))/2.0
#Animation stuff
bitrate=5000
framerate=4		#Number of frames per second in the output video
codec='libx264'		#Other options include mpeg4
suffix='.mp4'

###########################################################
### Initialize
###########################################################
if avdim=='x':
	dim=0
	d1=1
	d2=2
elif avdim=='y':
	dim=1
	d1=0
	d2=2
elif avdim=='z':
	dim=2
	d1=0
	d2=1
else:
	print( "avdim must be 'x', 'y' or 'z' - not %s"%avdim )
	exit()

###########################################################
### Read json
###########################################################
if not os.path.isfile(inputName):
	print("%s not found."%inputName)
	exit()
with open(inputName, 'r') as f:
  input = json.load(f)
xyzSize=array([30,30,1])
if "domain" in input:
	xyzSize[0]=input['domain'][0]
	if(len(input['domain'])>1):
		xyzSize[1]=input['domain'][1]
	else:
		xyzSize[1]=1
	if(len(input['domain'])>2):
		xyzSize[2]=input['domain'][2]
	else:
		xyzSize[2]=1
# Data
XYZ = zeros(shape=(3,xyzSize[0],xyzSize[1],xyzSize[2]),dtype=float)
DIR = zeros(shape=(3,xyzSize[0],xyzSize[1],xyzSize[2]),dtype=float)
MEAN = zeros(shape=(3,xyzSize[d1],xyzSize[d2]),dtype=float)
MAG = zeros(shape=(xyzSize[d1],xyzSize[d2]),dtype=float)
XY = zeros(shape=(2,xyzSize[d1],xyzSize[d2]),dtype=float)
S = zeros(shape=(xyzSize[0],xyzSize[1],xyzSize[2]),dtype=float)
AVS = zeros(shape=(xyzSize[d1],xyzSize[d2]),dtype=float)
CHARGE = zeros(shape=(xyzSize[d1],xyzSize[d2]),dtype=float)
ANGLE = zeros(shape=(xyzSize[d1],xyzSize[d2]),dtype=float)

### Setup the animation
# Figure
width=8
height=6
if not cbFlag:
	width=height
if myAspect == 'auto':
	shrink_factor = 1.0
else:
	shrink_factor = float(xyzSize[d2])/float(xyzSize[d1])
	if xyzSize[d1] > xyzSize[d2]:
		height *= shrink_factor
	else:
		width /= shrink_factor
		shrink_factor = 1.0
# Create the figure
fig1, ax = plt.subplots(figsize=(width, height))
#Create the colorbar
CS3 = imshow(AVS.T,cmap=myMap,vmin=-0.5, vmax=0.5,aspect=myAspect,extent=[0,xyzSize[d1],0,xyzSize[d2]])
if cbFlag:
	cb=colorbar(CS3,shrink=1,aspect=20*shrink_factor,pad=0.04)
	cb.ax.set_ylabel(r'Topo. charge, $q$')
# Make labels
if avdim=='x':
	labX='y'
	labY='z'
elif avdim=='y':
	labX='x'
	labY='z'
elif avdim=='z':
	labX='x'
	labY='y'
#TODO: make this work for desynced avs.dat and topocharge.dat files

###########################################################
### Read the data for animation
###########################################################
print( 'Reading data ...' )
if not os.path.isfile(directorData):
	print("%s not found."%directorData)
	exit()
datainfile = open(directorData,"r")
if not os.path.isfile(topoData):
	print("%s not found."%topoData)
	exit()
topoinfile = open(topoData, "r")
# Toss header
for i in range(13):
	line = datainfile.readline()
	line2 = topoinfile.readline()
i=0
j=0
n=-1
while datainfile:
	i=i+1
	line = datainfile.readline()
	topoLine = topoinfile.readline()
	if( not line or not topoLine ):
		break
	else:
		if j>=start:
			#avS field
			t,Qx,Qy,Qz,Vx,Vy,Vz,s = line.split("\t",8)
			XYZ[0][int(Qx)][int(Qy)][int(Qz)] = float(Qx) + 0.5
			XYZ[1][int(Qx)][int(Qy)][int(Qz)] = float(Qy) + 0.5
			XYZ[2][int(Qx)][int(Qy)][int(Qz)] = float(Qz) + 0.5
			DIR[0][int(Qx)][int(Qy)][int(Qz)] = float(Vx)
			DIR[1][int(Qx)][int(Qy)][int(Qz)] = float(Vy)
			DIR[2][int(Qx)][int(Qy)][int(Qz)] = float(Vz)
			S[int(Qx)][int(Qy)][int(Qz)] = float(s)
			#topo field
			t,Qx,Qy,Qz,C,angle = topoLine.split("\t", 6)
			CHARGE[int(Qx)][int(Qy)] = float(C)
			ANGLE[int(Qx)][int(Qy)] = float(angle)

	if i==xyzSize[0]*xyzSize[1]*xyzSize[2]:
		j=j+1
		if j>finish:
			break
		if j<start or j>finish:
			pass
		else:
			#Sum
			if avdim=='x':
				for x in range(xyzSize[0]):
					for y in range(xyzSize[1]):
						for z in range(xyzSize[2]):
							MEAN[0][y][z]=MEAN[0][y][z]+DIR[0][x][y][z]
							MEAN[1][y][z]=MEAN[1][y][z]+DIR[1][x][y][z]
							MEAN[2][y][z]=MEAN[2][y][z]+DIR[2][x][y][z]
							AVS[y][z]=AVS[y][z]+S[x][y][z]
				for y in range(xyzSize[1]):
					for z in range(xyzSize[2]):
						for i in range(3):
							MEAN[i][y][z]/=xyzSize[0]
						AVS[y][z]/=xyzSize[0]
			elif avdim=='y':
				for x in range(xyzSize[0]):
					for y in range(xyzSize[1]):
						for z in range(xyzSize[2]):
							MEAN[0][x][z]=MEAN[0][x][z]+DIR[0][x][y][z]
							MEAN[1][x][z]=MEAN[1][x][z]+DIR[1][x][y][z]
							MEAN[2][x][z]=MEAN[2][x][z]+DIR[2][x][y][z]
						AVS[x][z]=AVS[x][z]+S[x][y][z]
				for x in range(xyzSize[0]):
					for z in range(xyzSize[2]):
						for i in range(3):
							MEAN[i][x][z]/=xyzSize[1]
						AVS[x][z]/=xyzSize[1]
			elif avdim=='z':
				for x in range(xyzSize[0]):
					for y in range(xyzSize[1]):
						for z in range(xyzSize[2]):
							MEAN[0][x][y]=MEAN[0][x][y]+DIR[0][x][y][z]
							MEAN[1][x][y]=MEAN[1][x][y]+DIR[1][x][y][z]
							MEAN[2][x][y]=MEAN[2][x][y]+DIR[2][x][y][z]
							AVS[x][y]=AVS[x][y]+S[x][y][z]
				for x in range(xyzSize[0]):
					for y in range(xyzSize[1]):
						for i in range(3):
							MEAN[i][x][y]/=xyzSize[2]
						AVS[x][y]/=xyzSize[2]
			for x in range(xyzSize[d1]):
				for y in range(xyzSize[d2]):
					MAG[x][y]=sqrt( MEAN[0][x][y]**2+MEAN[1][x][y]**2+MEAN[2][x][y]**2 )
			#Save the instantaneous or current velocity field frame
			# Make Mesh
			if avdim=='x':
				for y in range(xyzSize[1]):
					for z in range(xyzSize[2]):
						XY[0][y][z]=XYZ[d1][0][y][z]
						XY[1][y][z]=XYZ[d2][0][y][z]
			elif avdim=='y':
				for x in range(xyzSize[0]):
					for z in range(xyzSize[2]):
						XY[0][x][z]=XYZ[d1][x][0][z]
						XY[1][x][z]=XYZ[d2][x][0][z]
			elif avdim=='z':
				for x in range(xyzSize[0]):
					for y in range(xyzSize[1]):
						XY[0][x][y]=XYZ[d1][x][y][0]
						XY[1][x][y]=XYZ[d2][x][y][0]

			# Save frame
			n=n+1
			fig1 = plt.figure(1)
			plt.clf()
			quiver( XY[0][::qx, ::qy], XY[1][::qx, ::qy], 
					myLength*MEAN[0][::qx, ::qy], myLength*MEAN[1][::qx, ::qy], 
					CHARGE[::qx, ::qy]+0.5, cmap=myMap, clim=(0, 1), scale=50/myLength,width=0.005*myLW, headlength=0, headwidth=0,headaxislength=0, pivot='middle')
			for x in range(xyzSize[d1]):
				for y in range(xyzSize[d2]):
					if( x%qx==0 and y%qy==0 ):
						#plot the field data
						#plot +1/2 defects as lines 
						if (CHARGE[x][y] > 1.0e-6):
							plt.arrow(XY[0][x][y], XY[1][x][y], 0.5*defLength*math.cos(ANGLE[x][y]), 0.5*defLength*math.sin(ANGLE[x][y]), color='k')
							plt.arrow(XY[0][x][y], XY[1][x][y], -0.5*defLength*math.cos(ANGLE[x][y]), -0.5*defLength*math.sin(ANGLE[x][y]), color='k')
						#plot -1/2 defects with three lines (per symmetry)
						if (CHARGE[x][y] < -1.0e-6):
							plt.arrow(XY[0][x][y], XY[1][x][y], 0.5*defLength*math.cos(ANGLE[x][y]), 0.5*defLength*math.sin(ANGLE[x][y]), color='k')
							plt.arrow(XY[0][x][y], XY[1][x][y], 0.5*defLength*math.cos(ANGLE[x][y] + math.pi*2/3), 0.5*defLength*math.sin(ANGLE[x][y] + math.pi*2/3), color='k')
							plt.arrow(XY[0][x][y], XY[1][x][y], 0.5*defLength*math.cos(ANGLE[x][y] - math.pi*2/3), 0.5*defLength*math.sin(ANGLE[x][y] - math.pi*2/3), color='k')
			if cbFlag:
				cb=colorbar(CS3,shrink=1,aspect=20*shrink_factor,pad=0.04)
				cb.ax.set_ylabel(r'Topo. charge, $q$')
			xlabel(r'$%s$'%labX)
			ylabel(r'$%s$'%labY)
			plt.axis(xmax=xyzSize[d1], xmin=0, ymax=xyzSize[d2], ymin=0)
			name='frame%04d.png'%(n)
			plt.axis('off') 
			savefig( name,bbox_inches='tight',pad_inches=0, transparent=makeTransparent )
		#Zero matrix
		DIR= zeros( (3,xyzSize[0],xyzSize[1],xyzSize[2]),dtype=float )
		MEAN = zeros(shape=(3,xyzSize[d1],xyzSize[d2]),dtype=float)
		AVS = zeros(shape=(xyzSize[d1],xyzSize[d2]),dtype=float)
		i=0
datainfile.close()

###########################################################
### Visualize data
###########################################################
#Animate
print( "Animating ..." )
name='topoCharge_%s%s'%(avdim,suffix)
myCommand="rm %s"%name
call(myCommand,shell=True)
myCommand = "ffmpeg -f image2 -r %d"%(framerate)+" -i frame%04d.png"+" -vcodec %s -b %dk -r %d %s"%(codec,bitrate,framerate,name)
call(myCommand,shell=True)
if not keepFrames:
	myCommand="rm frame*.png"
	call(myCommand,shell=True)

"""
	NIAMH-MPCD
	Probability density rendering script
	To average histograms of vector values since histogram range fluctuates

	Uses shendrukGroupStyle for formatting (from https://github.com/Shendruk-Lab/MPCDDefectLoader)
	Must install it or remove calls to shendrukGroupFormat

	Created by Tyler Shendruk
"""

from pylab import *
from subprocess import call
from scipy import integrate
import os
import argparse

# Use our custom style and colours
plt.style.use('shendrukGroupStyle')
import shendrukGroupFormat as ed

###########################################################
### Set up argparse
###########################################################
parser = argparse.ArgumentParser(description='Probability density rendering '
											 'script.')
parser.add_argument('dataName', type=str,
					help="Path to data (should be a histogram output for a "
						 "vector quantity)")
parser.add_argument('-X', '--xaxis', type=str,
					help="x axis label --- What is this the histogram of?", default="")
parser.add_argument('-s', "--start", type=int, default=0, help="Starting timestep for averaging")
parser.add_argument('-f', "--finish", type=int, default=9999999, help="Finishing timestep for averaging")
parser.add_argument("-m", "--makeMovie", type=int, default=1,
					help="Whether or not to animate temporal data (0 or 1)")
parser.add_argument("-a", "--plotAv", type=int,
					help="Whether or not to plot average (0 or 1)",
					default=1)
parser.add_argument("-x", "--includeX", type=int,
					help="Whether or not to include x (0 or 1)", default=1)
parser.add_argument("-y", "--includeY", type=int,
					help="Whether or not to include y (0 or 1)", default=1)
parser.add_argument("-z", "--includeZ", type=int,
					help="Whether or not to include z (0 or 1)", default=1)
parser.add_argument("-k", "--keepFrames", type=int,
                    help="0=don't keep (delete) frames; 1=keep frames",
                    default=0)
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
dataName = args.dataName
xaxis = args.xaxis
start = args.start
finish = args.finish
makeMovie = args.makeMovie
keepFrames = args.keepFrames
plotAv = args.plotAv
_x = args.includeX
_y = args.includeY
_z = args.includeZ
groupStyle = args.groupStyle

###########################################################
### Assumed arguments that a user could change
###########################################################
numBins=101
avOutName=dataName.split(".dat")[0]+"_av"

###########################################################
### Format and style
###########################################################
if groupStyle:
	# Use our custom style and colours
	plt.style.use('shendrukGroupStyle')
	import shendrukGroupFormat as ed
#Animation stuff
bitrate=5000
framerate=12		#Number of frames per second in the output video
	# Note that the ideal for this is _very_ variable so play around with it
codec='libx264'		#Other options include mpeg4
suffix='.mp4'

# Figure
fig1 = plt.figure(1)

###########################################################
### Animation
###########################################################
print( 'Reading file for min/max/movie ...' )
if not os.path.isfile(dataName):
	print("%s not found."%dataName)
	exit()
infile = open(dataName,"r")
#toss header
for i in range(13):
	line = infile.readline()
t=0
min=999999.
max=-999999.
maxC=0.0
X=zeros(shape=(numBins),dtype=float)	#Start of the bins
C=zeros(shape=(3,numBins),dtype=float)
while infile:
	line = infile.readline()
	if( not line ):
		break
	else:
		time=float(line)
		for i in range(numBins):
			line = infile.readline()
			toss,x,cx,cy,cz = line.split("\t")
			X[i]=float(x)
			C[0][i]=float(cx)
			C[1][i]=float(cy)
			C[2][i]=float(cz)
		line = infile.readline()
		if(t>=start):
			dx=X[1]-X[0]
			min=minimum(min,X[0])
			max=maximum(max,X[-1]+dx)
			if(_x):
				maxC=maximum(maxC,C[0].max())
			if(_y):
				maxC=maximum(maxC,C[1].max())
			if(_z):
				maxC=maximum(maxC,C[2].max())
			if(makeMovie):
				X+=0.5*dx
				fig1 = plt.figure(1)
				plt.clf()
				title(r'Time, %s$\tau$'%(time))
				if(_x):
					plot(X,C[0],'-o',label=r'$x$')
				if(_y):
					plot(X,C[1],'-s',label=r'$y$')
				if(_z):
					plot(X,C[2],'-^',label=r'$z$')
				xlabel(r'%s'%xaxis)
				ylabel(r'Counts')
				plt.axis(xmax=max, xmin=min, ymax=maxC, ymin=0)
				name='frame%04d.png'%(t)
				savefig(name, bbox_inches='tight')
		t+=1
infile.close()
if(makeMovie):
	#Animate
	print( "\tAnimating ..." )
	name='%s_animation%s'%(avOutName,suffix)
	myCommand="rm %s"%name
	call(myCommand,shell=True)
	myCommand = "ffmpeg -f image2 -r %d"%(framerate)+" -i frame%04d.png"+" -vcodec %s -b %dk -r %d %s"%(codec,bitrate,framerate,name)
	call(myCommand,shell=True)
	if not keepFrames:
		myCommand="rm frame*.png"
		call(myCommand,shell=True)

###########################################################
### Average
###########################################################
print( 'Reading data ...' )
infile = open(dataName,"r")
# Keep header
header=[]
for i in range(13):
	header.append(infile.readline())
t=0
MmM=max-min
avX=linspace(min,max,num=numBins,endpoint=True,dtype=float)	#Start of the bins
avP=zeros(shape=(3,numBins),dtype=float)
while infile:
	line = infile.readline()
	if( not line ):
		break
	else:
		time=float(line)
		for i in range(numBins):
			line = infile.readline()
			toss,x,cx,cy,cz = line.split("\t")
			X[i]=float(x)
			C[0][i]=float(cx)
			C[1][i]=float(cy)
			C[2][i]=float(cz)
		line = infile.readline()
		if(t>=start):
			# Integrate under the curve
			dx=X[1]-X[0]
			for d in range(3):
				norm=integrate.simpson(C[d],dx=dx)
				C[d]/=norm
			X+=0.5*dx
			for i in range(numBins):
				I=int(numBins*(X[i]-min)/MmM)
				for d in range(3):
					avP[d][I]+=C[d][i]
		t+=1
infile.close()
dx=avX[1]-avX[0]
for d in range(3):
	norm=integrate.simpson(avP[d],dx=dx)
	avP[d]/=norm

outfile=open(avOutName+".dat", "w")
cutTime=header[-1].split("\t")[1:]
header[-1]="%s\t%s\t%s\t%s"%(cutTime[0],cutTime[3],cutTime[5],cutTime[7])
for h in header:
	outfile.write(h)
for i in range(numBins):
	outfile.write("%e\t%e\t%e\t%e\n"%(avX[i],avP[0][i],avP[1][i],avP[2][i]))
outfile.close()

avX+=0.5*dx
if(plotAv):
	fig1 = plt.figure(1)
	plt.clf()
	if(_x):
		plot(X,avP[0],'-o',label=r'$x$')
	if(_y):
		plot(X,avP[1],'-s',label=r'$y$')
	if(_z):
		plot(X,avP[2],'-^',label=r'$z$')
	xlabel(r'%s'%xaxis)
	ylabel(r'PDF')
	plt.axis(xmax=max, xmin=min, ymin=0)
	plt.legend(loc='best',frameon=False)
	savefig('%s.png'%(avOutName), bbox_inches='tight', transparent=True)
	savefig('%s.pdf'%(avOutName), bbox_inches='tight', transparent=True)
	show()

exit()
"""
	NIAMH-MPCD
	Density field rendering script

	Uses shendrukGroupStyle for formatting (from https://github.com/Shendruk-Lab/MPCDDefectLoader)
	Must install it or remove calls to shendrukGroupFormat

	Created by Timofey Kozhukhov
"""

from pylab import *
from subprocess import call
import os
import json
import argparse

from defectHandler import getDefectData

###########################################################
### Set up argsparse
###########################################################
parser = argparse.ArgumentParser(description='Density field rendering script.')
parser.add_argument("dataname", type=str, help="Path to the data (should be densityfield.dat)")
parser.add_argument("inputname", type=str, help="Path to input .json file")
parser.add_argument('-s',"--start", type=int, default=0, help="Starting timestep for averaging")
parser.add_argument('-f',"--finish", type=int, default=9999999, help="Finishing timestep for averaging")
parser.add_argument('-a',"--avdim", type=str, default='z', help="Dimension to 'average'/'slice' over for averaging")
parser.add_argument('-i',"--sliceIndex", type=int, default=0, help="Index of the plane to be sliced for averaging")
parser.add_argument("-A", "--myAspect", type=str, help="'auto' or 'equal'", default='equal')
parser.add_argument("-k", "--keepFrames", type=int, help="0=don't keep (delete) frames; 1=keep frames", default=0)
parser.add_argument("-N", "--normalise", type=int, help="Normalise the colourbar (1) or not (0)", default=0)
parser.add_argument("-p", "--savePDF", type=int, help="1 saves transparent pdfs for papers, 0 for none", default=0)
parser.add_argument("-c", "--colourbar", type=int, help="1 adds a colourbar to the plot, 0 for none", default=1)
parser.add_argument("-d", "--defectData", type=str, help="Path to defect data (if any)", default="")
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
dataName = args.dataname
inputName = args.inputname
start = args.start
finish = args.finish
avdim = args.avdim
sliceIndex = args.sliceIndex
myAspect = args.myAspect
keepFrames = args.keepFrames
normalise = args.normalise
savePDF = args.savePDF
cbFlag = args.colourbar
defectData = args.defectData
groupStyle = args.groupStyle

makeTransparent = False # Transparent backgrounds make crappy videos, but look good on webpages

###########################################################
### Format and style
###########################################################
if groupStyle:
	# Use our custom style and colours
	plt.style.use('shendrukGroupStyle')
	import shendrukGroupFormat as ed
# colour map used when species=1
if groupStyle:
	myMap=ed.plasma 
else:
	myMap=plt.cm.plasma

pop0Col = np.array([0, 1, 0]) # colour to show when ONLY pop0 is present (RGB)
pop1Col = np.array([1, 0, 0]) # colour to show when ONLY pop1 is present (RGB)
def getRGBAColField(totalPop, maxTotalPop, popProp):
	"""
	Helper method to get colours.

	Args:
		totalPop: total number of particles in given cell as a MxN np array
		maxTotalPop: the maximum total number of particles in any cell globally (scalar)
		popProp: proportion of pop0 particles as a MxN np array
	Returns:
		RGBA colour field as a MxNx4 np array
	"""
	# compute alpha
	a = totalPop/maxTotalPop 
	# compute RGB
	r = popProp*pop0Col[0] + (1-popProp)*pop1Col[0]
	g = popProp*pop0Col[1] + (1-popProp)*pop1Col[1]
	b = popProp*pop0Col[2] + (1-popProp)*pop1Col[2]
	return np.stack((r, g, b, a), axis=2)

#Animation stuff
bitrate=5000
framerate=12		#Number of frames per second in the output video
# Note that the ideal for this is _very_ variable so play around with it
codec='libx264'		#Other options include mpeg4
suffix='.mp4'

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
species=len(input['species'])		#1 = do total pop only; 2 = two species; will fail if more than 1
if species < 1:
	print( "Species must be 1 or 2 - not %s"%species )
	exit()
if species > 2:
	print( "More than 2 species are not supported!")
	print( "Species must be 1 or 2 - not %s"%species )
	exit()

###########################################################
### Initialize
###########################################################
if avdim=='x':
	print("Note, for now we assume avdim=z! Edit the script to fix this as necessary")
	dim=0
	d1=1
	d2=2
elif avdim=='y':
	print("Note, for now we assume avdim=z! Edit the script to fix this as necessary")
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
if sliceIndex>=xyzSize[dim]:
	print( "sliceIndex must be within the system: sliceIndex=%d but %s system size=%d"%(sliceIndex,avdim,xyzSize[dim]) )
	exit()

# Data
TOTPOP = zeros(shape=(xyzSize[0],xyzSize[1],xyzSize[2]),dtype=float)
sliceTOTPOP = zeros(shape=(xyzSize[d1],xyzSize[d2]),dtype=float)

POP0, slicePop0, POP1, slicePOP1, slicePOP0prop, sliceRGBA = None, None, None, None, None, None
if species == 2: # lets be memory efficient
	POP0 = zeros(shape=(xyzSize[0],xyzSize[1],xyzSize[2]),dtype=float)
	slicePOP0 = zeros(shape=(xyzSize[d1],xyzSize[d2]),dtype=float)
	slicePOP0prop = zeros(shape=(xyzSize[d1],xyzSize[d2]),dtype=float) # proportion of population 0 (POP0/TOTPOP)
	sliceRGBA = zeros(shape=(xyzSize[d1],xyzSize[d2],4),dtype=float) # RGBA colour for each cell

###########################################################
### defect handling if needed
###########################################################
LOADDEFECTS = False
defects = []
if defectData != "":
	print("Loading defects for rendering")
	LOADDEFECTS = True

	defContainer = getDefectData(defectData, np.array([xyzSize[0], xyzSize[1], xyzSize[2]]))
	for defList in defContainer:
		defects.append(defList.defectList)
	print("Finished loading defects")

# Figure
width=8
height=6
if not cbFlag:
	width=height
fig1, ax = plt.subplots(figsize=(width, height))
if myAspect == 'auto':
    shrink_factor = 1.0
elif float(xyzSize[d2])<float(xyzSize[d1]):
    shrink_factor = float(xyzSize[d2])/float(xyzSize[d1])
else:
    shrink_factor = 1.0

###########################################################
### Read the data for animation
###########################################################

### Setup the animation
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

print( 'Reading data ...' )
if not os.path.isfile(dataName):
	print("%s not found."%dataName)
	exit()
infile = open(dataName,"r")
# Toss header
for i in range(13):
	line = infile.readline()
i=0 # line counter (within timestep)
currTStep=0
n=-1
tStepMaxPop = 0
while infile:
	i=i+1
	line = infile.readline()
	if( not line ):
		break
	else:
		if currTStep>=start:
			currLine = line.split("\t")
			t,Qx,Qy,Qz,pop,mass,popSRD,popMD,popSW = [currLine[i] for i in range(9)] # get core data
			TOTPOP[int(Qx)][int(Qy)][int(Qz)] = float(pop)

			if species == 2: # now handle multispecies info
				sp0 = currLine[8]
				POP0[int(Qx)][int(Qy)][int(Qz)] = float(sp0) # this is all we need for now

	if i==xyzSize[0]*xyzSize[1]*xyzSize[2]: # end of a timestep
		currTStep=currTStep+1
		if currTStep>finish:
			break
		if currTStep<start or currTStep>finish:
			aaa=0
		else:
			# compute the information for the slice we want to render
			##NOTE: for now, assuming avdim=z and z=0 for simplicity
			sliceTOTPOP = TOTPOP[:,:,sliceIndex]
			sliceMax = np.max(sliceTOTPOP) # maximum population in the slice
			if species == 2: # for multispecies
				slicePop0 = POP0[:,:,sliceIndex]
				slicePop0prop = slicePop0/sliceTOTPOP
			# now can render!
			n = n + 1
			fig1 = plt.figure(1)
			plt.clf()

			# handle normalisation if required
			renderTOTPOP = sliceTOTPOP
			if normalise:
				renderTOTPOP = renderTOTPOP/sliceMax

			#handle single species render first
			imshow(renderTOTPOP,cmap=myMap,vmin=0, aspect=myAspect,extent=[0,xyzSize[d1],0,xyzSize[d2]])

			#Create the colorbar
			CS3 = imshow(renderTOTPOP.T,vmin=0,cmap=myMap,aspect=myAspect,extent=[0,xyzSize[d1],0,xyzSize[d2]])
			if cbFlag:
				cb=colorbar(CS3,shrink=shrink_factor,aspect=20*shrink_factor,pad=0.04)
				if normalise:
					cb.ax.set_ylabel(r'Number, $N_C / N_C^\mathrm{max}$')
				else:
					cb.ax.set_ylabel(r'Number, $N_C$')

			# perform matplotlib bits and save
			xlabel(r'$%s$'%labX)
			ylabel(r'$%s$'%labY)
			plt.axis(xmax=xyzSize[d1], xmin=0, ymax=xyzSize[d2], ymin=0)

			name='frame_s=1_%04d.png'%(n)
			namepdf='frame_s=1_%04d.pdf'%(n)
			
			## uncomment below for snapshots!
			plt.axis('off') 
			plt.subplots_adjust(left=0, right=1, top=1, bottom=0)
			savefig(name, bbox_inches='tight', transparent=makeTransparent)
			# save trans pdf
			if savePDF: savefig(namepdf, transparent=True, bbox_inches='tight')
				
			# also handle multispecies render if necessary
			if species == 2: # for multispecies
				plt.clf()
				sliceRGBA = getRGBAColField(sliceTOTPOP.T, sliceMax, slicePop0prop.T)
				imshow(sliceRGBA, aspect=myAspect, extent=[0,xyzSize[d1],0,xyzSize[d2]])
				# perform matplotlib bits and save
				xlabel(r'$%s$'%labX)
				ylabel(r'$%s$'%labY)
				plt.axis(xmax=xyzSize[d1], xmin=0, ymax=xyzSize[d2], ymin=0)
				name='frame_s=2_%04d.png'%(n)
				# savefig( name, bbox_inches='tight' )
				savefig( name, bbox_inches='tight', transparent=True )
				if savePDF: savefig(namepdf, transparent=True, bbox_inches='tight')
			i=0 # reset tStep line counter
infile.close()

#Animate
print( "Animating ..." )
# single species
name='density_%s%d%s'%(avdim,sliceIndex,suffix)
myCommand="rm %s"%name
call(myCommand,shell=True)
myCommand = "ffmpeg -f image2 -r %d"%(framerate)+" -i frame_s=1_%04d.png"+" -vcodec %s -b %dk -r %d %s"%(codec,bitrate,framerate,name)
call(myCommand,shell=True)

if species == 2: # for multispecies
	name='density_s=%d_%s%d%s'%(species,avdim,sliceIndex,suffix)
	myCommand="rm %s"%name
	call(myCommand,shell=True)
	myCommand = "ffmpeg -f image2 -r %d"%(framerate)+" -i frame_s=2_%04d.png"+" -vcodec %s -b %dk -r %d %s"%(codec,bitrate,framerate,name)
	call(myCommand,shell=True)

if not keepFrames:
	myCommand="rm frame*.png"
	call(myCommand,shell=True)

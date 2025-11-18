
//================================================================================
//
// name:   md.h
// author: ftessier
// date:   2005-05-03 @ 11:04:36
//
//================================================================================

#ifndef MD_H
#define MD_H

//================================================================================
// Prototypes
//================================================================================

// inline int	MDstep			    (simptr sim);
// inline int MDstep 			    (simptr sim,int MDmode,struct particleMPC *pSRD,struct bc *pBC,struct spec *SP,int GPOP,int NBC,struct cell ***CL);
int MDstep 			                (simptr sim,int MDmode,struct particleMPC *pSRD,struct bc *pBC,struct spec *SP,int GPOP,int NBC,struct cell ***CL);
// void VelocityVerletStep			(simptr sim);
void VelocityVerletStep 		    (simptr sim,int MDmode,struct particleMPC *pSRD,struct bc *pBC,struct spec *SP,int GPOP,int NBC,struct cell ***CL);
void ComputeElectrostaticForcesSRD	(simptr sim,struct particleMPC *pSRD,struct spec *SP,int GPOP, int steps,struct cell ***CL);

void ComputeForces			        (simptr sim);
void ComputeForcesSRD			    (simptr sim,int MDmode,particleMPC *pSRD,spec *SP,int GPOP,struct cell ***CL);
void ComputeCapForces			    (simptr sim);
void ComputeDispersionForces		(simptr sim);
void ComputeDispersionForcesSRD		(simptr sim,particleMPC *pSRD,spec *SP,int GPOP);
void ComputeDispersionForcesSRDCell	(simptr sim,int MDmode,struct particleMPC *pSRD,struct spec *SP,int GPOP,struct cell ***CL);
void ComputeCapDispersionForces		(simptr sim);
void ComputeElectrostaticForces		(simptr sim);
void ComputeAnchorForces	        (simptr sim);
void ComputeFeneForces  		    (simptr sim);
void ComputeBendForces              (simptr sim);
void ComputeDihedralForces          (simptr sim);
void ComputeNemForces               (simptr sim,struct spec *SP,struct cell ***CL);
void ComputeSqueezeForces   	    (simptr sim);

real LennardJones    	    	    (particleMD *p1, particleMD *p2, real dx, real dy, real dz, real rCut2);
real LennardJonesSRD        		(particleMD *p1,particleMPC *p2,spec *SP,real dx,real dy,real dz,real dt,real rCut2);
real LennardJonesCap	        	(particleMD *p1, particleMD *p2, real dx, real dy, real dz, real rCut2, simptr sim);
real Coulomb 			            (particleMD *p1, particleMD *p2, real dx, real dy, real dz, real rCutCoul2, real bjerrum, real lambda_Ds);
real HarmonicSpring 	        	(particleMD *p1, real dx, real dy, real dz, real k);
real HarmonicSpringTube     		(particleMD *p1, real dy, real dz, real k);
real QuadraticSpringTube	        (particleMD *p1, real dy, real dz, real k);
real SexticSpringTube       		(particleMD *p1, real dy, real dz, real k);
real LennardJonesTube		        (particleMD *p1, real dy, real dz, real k, simptr sim);
real FENE        			        (particleMD *p1, particleMD *p2, real dx, real dy, real dz, real k, real r0);
real bendHarmonic                   (particleMD *p1, particleMD *p2, particleMD *p3, real dx12, real dy12, real dz12, real dx23, real dy23, real dz23, real k, real equi);
real dihedralHarmonic               (particleMD *p1, particleMD *p2, particleMD *p3, particleMD *p4, real dx12, real dy12, real dz12, real dx23, real dy23, real dz23, real dx34, real dy34, real dz34, real k, real equi);
real bendNematic                    (particleMD *p1, particleMD *p2, real dx12, real dy12, real dz12, real dx23, real dy23, real dz23, real k, real S, real dt, struct spec *SP, struct cell *CL);

void ApplyPBC        				(simptr sim, real *px, real *py, real *pz);
void ApplyPBC_dr        			(simptr sim, real *dx, real *dy, real *dz);
void PeriodicBoundaries	        	(simptr sim);
void CheckNebrList      			(simptr sim);
void RefreshNebrList 	        	(simptr sim);
void RefreshCellList         		(simptr sim);

simptr launchMD 	        		(int argc, char *argv[]);
// void integrateMD		        	(simptr sim, int steps);
void integrateMD			        (simptr sim,int MDmode, int steps,struct particleMPC *pSRD,struct bc *pBC,struct spec *SP,int GPOP,int NBC,struct cell ***CL);
void integrateMD_Pinned	            (simptr sim,int MDmode, int steps,struct particleMPC *pSRD,struct bc *pBC,struct spec *SP,int GPOP,int NBC,struct cell ***CL);

#endif

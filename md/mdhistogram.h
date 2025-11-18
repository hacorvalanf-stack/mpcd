
//================================================================================
//
// name:   mdhistogram.h
// author: ftessier
// date:   2005-05-03 @ 11:04:36
//
//================================================================================


#ifndef MDHISTOGRAM_H
#define MDHISTOGRAM_H


//================================================================================
// Prototypes
//================================================================================

void 		SetupHistograms 			(simptr sim);
void 		HistogramSetupFiles 		(simptr sim);
void 		HistogramSetupFuncs 		(simptr sim);

void 		HistogramCollect			(simptr sim, histptr h);
void		HistogramAdd				(histptr h, real x, real y, real z, real value);
int 		HistogramIndexValue			(histptr h, real x, real y, real z);

void 		HistogramReset				(histptr h);
void 		HistogramNormalize 	 		(histptr h);
real 		HistogramGetBinSize	 		(histptr h);

histptr 	HistogramNew 		 		(histptr *h);
histptr 	HistogramDuplicate 	 		(histptr h, char *newlabel);
histptr 	HistogramGet		 		(histptr h, char *label);

void		HistogramPrint 		 		(simptr sim, histptr h);
void		HistogramPrintData	 		(simptr sim, histptr h, FILE *stream);

void 		HistogramAllocateBins		(histptr h);
void 		HistogramSetBinSize     	(real *min, real *max, int *n, real *bin);
void 		HistogramFreeList 	 		(histptr h);
int 		HistogramSpaceDomain 		(histptr h,  real s1,   real s2,   real s3);

void		HistogramFunction	 		(void *simvoid, histptr h, int action);
int 		HistogramFunction_bond_orientation	(void *simvoid, histptr h, int action);
int 		HistogramFunction_Tkin  	(void *simvoid, histptr h, int action);
int 		HistogramFunction_Tconf  	(void *simvoid, histptr h, int action);
int 		HistogramFunction_Tdivf 	(void *simvoid, histptr h, int action);
int 		HistogramFunction_q			(void *simvoid, histptr h, int action);
int 		HistogramFunction_fv		(void *simvoid, histptr h, int action);
int 		HistogramFunction_vx		(void *simvoid, histptr h, int action);
int 		HistogramFunction_dr_yz		(void *simvoid, histptr h, int action);
int 		HistogramFunction_bond		(void *simvoid, histptr h, int action);
int 		HistogramFunction_fenex_N	(void *simvoid, histptr h, int action);
int 		HistogramFunction_qeff_N (void *simvoid, histptr h, int action);
int 		HistogramFunction_hwall_N	(void *simvoid, histptr h, int action);
int 		HistogramFunction_ellipsoid	(void *simvoid, histptr h, int action);
int 		HistogramFunction_persist	(void *simvoid, histptr h, int action);

int HistogramFunction_bond_orientation_rho_y (void *simvoid, histptr h, int action);
int HistogramFunction_density_rho_y (void *simvoid, histptr h, int action);


void transformCoord(simptr sim,bc WALL);
void untransformCoord(simptr sim,bc WALL);

#endif

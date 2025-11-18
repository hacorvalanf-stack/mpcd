
//================================================================================
//
// name:   mdvmd.h
// author: tshen
// date:   2011-09-22 @ 13:59:36
//
//================================================================================


#ifndef MDVMD_H
#define MDVMD_H


//================================================================================
// Prototypes
//================================================================================

void 		SetupVMD		 			(simptr sim);
void 		CopyStepCounters			(int *dest, int *source);
void 		VMDSetupFiles 			(simptr sim);
void 		VMDSetupFuncs 			(simptr sim);
sceneptr 	VMDNew 		 			(sceneptr *s);
sceneptr 	VMDDuplicate 	 			(sceneptr s, char *newlabel);
sceneptr 	VMDGet		 			(sceneptr s, char *label);
void		VMDPrint	 		 		(simptr sim, sceneptr s);
void 		VMDFreeList 	 			(sceneptr s);
int 		VMDSpaceDomain 			(sceneptr s, real s1, real s2, real s3);
void		VMDFunction	 			(void *simvoid, sceneptr s, int action);

#endif

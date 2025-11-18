
//================================================================================
//
// name:   mdscene.h
// author: ftessier
// date:   2005-05-03 @ 11:04:36
//
//================================================================================


#ifndef MDSCENE_H
#define MDSCENE_H


//================================================================================
// Prototypes
//================================================================================

void 		SetupScenes		 			(simptr sim);
void 		CopyStepCounters			(int *dest, int *source);
void 		SceneSetupFiles 			(simptr sim);
void 		SceneSetupFuncs 			(simptr sim);
sceneptr 	SceneNew 		 			(sceneptr *s);
sceneptr 	SceneDuplicate 	 			(sceneptr s, char *newlabel);
sceneptr 	SceneGet		 			(sceneptr s, char *label);
void		ScenePrint	 		 		(simptr sim, sceneptr s);
void 		SceneFreeList 	 			(sceneptr s);
int 		SceneSpaceDomain 			(sceneptr s, real s1, real s2, real s3);
void		SceneFunction	 			(void *simvoid, sceneptr s, int action);
void 		ScenePrintObjects 			(simptr sim, FILE *stream);

#endif

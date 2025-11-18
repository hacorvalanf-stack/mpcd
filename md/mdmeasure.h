
//================================================================================
//
// name:   mdmeasure.h
// author: ftessier
// date:   2005-05-03 @ 11:04:36
//
//================================================================================


#ifndef MDMEASURE_H
#define MDMEASURE_H


#define	STEP_PROGRESS_MARK 100


//================================================================================
// Prototypes
//================================================================================

void		Measure					(simptr sim);
void		IncrementStepCounters 	(simptr sim);
int			SimulationPhaseCheck	(simptr sim);
void		SimulationPhaseReset	(simptr sim);
void		EvalProperties 			(simptr sim);
void		AccumProperties 		(simptr sim, int action);
void 		AverageProperties 		(void *simvoid);
void		FinalCalculations 		(simptr sim);
void		ZetaPotential 			(simptr sim);


#endif


//================================================================================
//
// name:   mdthermostat.h
// author: ftessier
// date:   2005-05-03 @ 11:04:36
//
//================================================================================


#ifndef MDTHERMOSTAT_H
#define MDTHERMOSTAT_H


//================================================================================
// Prototypes
//================================================================================

void ThermostatRescale	(void *simvoid);
void ThermostatDPD 		(particleMD *p1, particleMD *p2, real dx, real dy, real dz, real rCut2, real eta, real sigma, real dtSqrti, int groupThermDPD);

#endif

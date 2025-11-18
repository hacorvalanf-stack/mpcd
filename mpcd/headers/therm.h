#ifndef THERM_H
#define THERM_H

/* ****************************************** */
/* ****************************************** */
/* ****************************************** */
/* ************ DECLARE FUNCTIONS *********** */
/* ****************************************** */
/* ****************************************** */
/* ****************************************** */
/*
    Subroutines that relate to the MPCD fluid's
    thermal energy and the thermostat
*/
double thermostat( double KBT,double KBTNOW,double t,double RELAX,int TSTECH );
void scaleT( double KBT,double KBTNOW,double t,double RELAX,double VEL[],double VELNOW[],int TSTECH,spec SP[],int LC,bc *WALL,particleMPC *p,cell ***CL );
void heyes_cell( cell CL,double KBT,double KBTNOW,double RELAX );
double TEMP( particleMPC *pp,spec SP[],bc WALL[],double VEL[] );
double calcE( particleMPC *p,spec *pSP,bc WALL[] );
double calcE_MPC( particleMPC *p,spec *pSP );
double calcE_BC( bc *WALL );
double calcE_LC( cell ***CL,int LC,spec *pSP );
void avVel( cell ***CL,double AVVEL[] );
void avOri( particleMPC *p,double AVVEL[] );
double avEnstrophy( cell ***CL );

#endif

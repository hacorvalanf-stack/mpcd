#ifndef LC_H
#define LC_H

/* ****************************************** */
/* ****************************************** */
/* ****************************************** */
/* ************ DECLARE FUNCTIONS *********** */
/* ****************************************** */
/* ****************************************** */
/* ****************************************** */
/*
   These functions are the routines that deal
   with the liquid crystal component of the NIAMH-MPCD algorithm
*/

void LCcollision( cell *CL,spec *SP,double KBT,int zeroMFPot,double dt,double SG,int LC );
void addToTensOrderParam( particleMPC *pMPC,double **S );
void addToTensOrderParamVel( particleMPC *pMPC,double **S );
void tensOrderParam( cell *CL,double **S,int LC );
void tensOrderParamNNN( cell ***CL,double **S,int LC,int a,int b,int c );
double avOrderParam( particleMPC *p,int LC,double avDIR[] );
double avS4( particleMPC *p,int LC,double DIR[] );
double binderCumulant( cell ***CL,int L,int LC );

void magTorque( particleMPC *pMPC,spec *SP,double dt,double MAG[] );
void magTorque_all( particleMPC *pp,spec *SP,double dt,double MAG[] );
void magTorque_CL( cell *CL,spec *SP,double dt,double MAG[] );
void jefferysTorque( cell *CL,spec *SP,double dt );
void localVelGrad( cell ***CL );
void velGrad3D( cell ***CL );
void velGrad2D( cell ***CL );
void velGrad1D( cell ***CL );
void velGradD3Q15( cell ***CL );
void velGradD2Q9( cell ***CL );
double topoChargeLocal( cell ***CL, int i, int j, int k);
double topoAngleLocal( cell ***CL, int i, int j, int k, double charge);
double topoSmallestAngle( double u[], double v[]);
void computeQ(cell CL, double output[_2D][_2D]);

void larsonRotRate(double dudt[],double w[],double u[],double E[_3D][_3D],double tumbleParam);
void larsonRotRateOLD_AND_SLOW(double dudt[],double w[],double u[],double E[_3D][_3D],double tumbleParam);
void brielsRotRate(double dudt[],double w[],double u[],double E[_3D][_3D]);
void saintillanRotRate(double dudt[],double w[],double u[],double E[_3D][_3D],double tumbleParam);

void oriBC( particleMPC *pp,spec *SP,bc *WALL,double n[] );
void torqueLCBC( bc *WALL,double n[], double U0[], double torqueMPC[],double rodlength, double posColl[] );

void genrand_maierSaupe( double DIR[],double rotAx[],double rotAngle,double U[],double KBT,double S,double effM );
void genrand_maierSaupeGAUSS_3D( double rotAx[],double rotAngle,double U[],double KBT,double S,double effM );
void genrand_maierSaupeGAUSS_2D( double rotAx[],double rotAngle,double U[],double KBT,double S,double effM );
void genrand_maierSaupeEXP_3D( double rotAx[],double rotAngle,double U[],double KBT,double S,double effM );
void genrand_maierSaupeMetropolis_3D( double DIR[],double rotAx[],double rotAngle,double U[],double KBT,double S,double effM );
void genrand_maierSaupeMetropolis_2D( double DIR[],double rotAx[],double rotAngle,double U[],double KBT,double S,double effM );

void andersenROT_LC( cell *CL,spec *SP,specSwimmer SS,double KBT,double dt,double *CLQ,int outP );
void dipoleAndersenROT_LC( cell *CL,spec *SP,specSwimmer SS,double KBT,double RELAX,double dt,int RTECH,double *CLQ,int outP );

#endif


//================================================================================
//
// name:   mdutil.h
// author: ftessier
// date:   2005-05-03 @ 11:04:36
//
//================================================================================


#ifndef MDUTIL_H
#define MDUTIL_H


// Mersenne twister period parameters
#define mtN 624
#define mtM 397
#define MATRIX_A   0x9908b0dfUL	// constant vector a
#define UPPER_MASK 0x80000000UL // most significant w-r bits
#define LOWER_MASK 0x7fffffffUL // least significant r bits


// Mersenne twister tempering parameters
#define TEMPERING_MASK_B 		0x9d2c5680UL
#define TEMPERING_MASK_C 		0xefc60000UL
#define TEMPERING_SHIFT_U(y)  	(y >> 11)
#define TEMPERING_SHIFT_S(y)  	(y << 7)
#define TEMPERING_SHIFT_T(y)  	(y << 15)
#define TEMPERING_SHIFT_L(y)  	(y >> 18)


//================================================================================
// Prototypes
//================================================================================
unsigned long   RandomSeed 				(unsigned long seed);
unsigned long	RandomInteger			();
real 			RandomReal				();
real  			*RandomVector3D			(real *v);
void 			CoordinateOrder 	 	(int domain, real **p1, real **p2, real **p3);
// inline void		SwapRealPointers	 	(real **a, real **b);
void		SwapRealPointers	 	(real **a, real **b);
void 			CoordinateTransform  	(point *from, point *to, real capR);
real 			DistanceSquared 		(particleMD *p1, particleMD *p2, int coord, real capR);
void 			CopyStepCounters		(int *dest, int *source);
void 			FillStepCounters		(int *dest, int count);
int 			GroupToIndex 			(int group);
int 			GroupFromIndex 			(int index);
char 			*GroupToStr 			(int group, char *str);
char			*TypeToStr 				(int type, char *str);


#endif

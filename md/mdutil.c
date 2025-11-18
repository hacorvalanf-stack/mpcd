
//================================================================================
//
// name:   mdutil.c
// author: ftessier
// date:   2005-05-03 @ 11:04:36
//
// Utility routines for molecular dynamics program, including a random number
// generator in the form of a Mersenne Twister.
//
//================================================================================


#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <unistd.h>
#include "mdtypes.h"
#include "mdutil.h"


//================================================================================
// Global variables
//================================================================================

// Mersenne twister static variables
static unsigned long mt[mtN];		// the array for the state vector
static int mti=mtN+1;			// mti=mtN+1 means mt[mtN] is not initialized


//================================================================================
unsigned long RandomSeed (unsigned long seed)
//================================================================================
{
	// Initializing the array with a *NONZERO* seed. Setting initial seeds to
	// mt[mtN] using the generator Line 25 of Table 1 in [KNUTH 1981, The Art of
	// Computer Programming Vol. 2 (2nd Ed.), pp102]

	// Get a seed from time*pid if seed=0
	if (!seed) seed = time(0)*getpid();

	// Initialize mersenne twister array
	mt[0]= seed & 0xffffffff;
	for (mti=1; mti<mtN; mti++) {
		mt[mti] = (1812433253UL * (mt[mti-1] ^ (mt[mti-1] >> 30)) + mti);
		// See Knuth TAOCP Vol2. 3rd Ed. P.106 for multiplier.  In the previous
		// versions, MSBs of the seed affect only MSBs of the array mt[].
		// PREVIOUS: mt[mti] = (69069 * mt[mti-1]) & 0xffffffff;
		// 2002/01/09 modified by Makoto Matsumoto
		mt[mti] &= 0xffffffffUL;
		// for 32-bit machines
	}

	return (seed);
}


//================================================================================
unsigned long RandomInteger()
//================================================================================
{
	// Mersenne twister to generate random 32-bit integers in the interval
	// [0,4294967295].  [ http://www.math.keio.ac.jp/~matumoto/emt.html ]

	unsigned long y;
	static unsigned long mag01[2]={0x0UL, MATRIX_A};

	if (mti >= mtN) { 							// generate N words at one time
		int kk;
		if (mti == mtN+1)   					// if randomSeed() has not been called,
			RandomSeed(time(0)*getpid());		// a default initial seed is used.

		for (kk=0; kk<mtN-mtM; kk++) {
			y = (mt[kk]&UPPER_MASK)|(mt[kk+1]&LOWER_MASK);
			mt[kk] = mt[kk+mtM] ^ (y >> 1) ^ mag01[y & 0x1UL];
		}
		for (; kk<mtN-1; kk++) {
			y = (mt[kk]&UPPER_MASK)|(mt[kk+1]&LOWER_MASK);
			mt[kk] = mt[kk+(mtM-mtN)] ^ (y >> 1) ^ mag01[y & 0x1UL];
		}
		y = (mt[mtN-1]&UPPER_MASK)|(mt[0]&LOWER_MASK);
		mt[mtN-1] = mt[mtM-1] ^ (y >> 1) ^ mag01[y & 0x1UL];
		mti = 0;
	}

	y = mt[mti++];
	y ^= TEMPERING_SHIFT_U(y);
	y ^= TEMPERING_SHIFT_S(y) & TEMPERING_MASK_B;
	y ^= TEMPERING_SHIFT_T(y) & TEMPERING_MASK_C;
	y ^= TEMPERING_SHIFT_L(y);

	return y;
}


//================================================================================
real RandomReal()
//================================================================================
{
	// Scales a random integer to get a real number in the interval [0,1].

	return RandomInteger()*(2.3283064370807974e-10);	// 1.0 / 4294967295.0
}


//================================================================================
real *RandomVector3D (real *v)
//================================================================================
{
	// Produces unit vectors in 3 dimensions with uniformly distributed random
	// orientations, using a standard rejection method [rap95].

	real x, y, s;

	do {
		x = 2 * RandomReal() - 1;
		y = 2 * RandomReal() - 1;
		s = x*x + y*y;
	} while (s > 1);

	v[2] = 1 - 2*s;
	s = 2 * sqrt(1-s);
	v[0] = s*x;
	v[1] = s*y;

	return v;
}


//================================================================================
void CoordinateOrder (int domain, real **p1, real **p2, real **p3)
//================================================================================
{
	// Changes the order of coordinate pointers to match the choices of domain
	// coordinates.

	switch (domain) {

	case y_:  SwapRealPointers (p1, p2); break;
	case z_:  SwapRealPointers (p1, p3); SwapRealPointers (p2, p3); break;

	case xz_: SwapRealPointers (p2, p3); break;
	case yz_: SwapRealPointers (p1, p2); SwapRealPointers (p2, p3); break;

	}
}


//================================================================================
// inline void SwapRealPointers (real **a, real **b)
void SwapRealPointers (real **a, real **b)
//================================================================================
{
	// Swap two real pointers pointed to by **a and **b.

	real *tmp;

	tmp = *a;
	*a  = *b;
	*b  = tmp;
}


//================================================================================
void CoordinateTransform (point *from, point *to, real capR)
//================================================================================
{
	// Transform the coordinates of a point, acccording to the coord field of
	// the from and the to points. We always first transform to cartesian
	// coordinates, then to the final coordinate system.

	real	f1, f2, f3;
	real	t1, t2, t3, tmp;
	int		axis;

	// local coordinate variables
	t1 = f1 = to->c1 = from->c1;
	t2 = f2 = to->c2 = from->c2;
	t3 = f3 = to->c3 = from->c3;

	// return if points are identical
	if (memcmp(from, to, sizeof(point)) == 0) return;

	// 1. first convert everything to cartesian coordinates

	// cylindrical -> cartesian
	if (from->coord & CYLINDRICAL) {
		t1 = (f1-capR) * cos(f2);
		t2 = (f1-capR) * sin(f2);
		t3 = f3;
		axis = (from->coord & ~CYLINDRICAL);
		if (axis == x_) {
			tmp = t1;
			t1  = t3;
			t2  = tmp;
			t3  = t2;
		}
		else if (axis == y_) {
			tmp = t1;
			t1  = t2;
			t2  = t3;
			t3  = tmp;
		}
	}

	// spherical -> cartesian
	else if (from->coord & SPHERICAL) {
		t1 = f1 * sin(f3) * cos(f2);
		t2 = f1 * sin(f3) * sin(f2);
		t3 = f1 * cos(f3);
	}

	// compute new cartesian coordinates
	f1 = t1 + (from->x0 - to->x0);
	f2 = t2 + (from->y0 - to->y0);
	f3 = t3 + (from->z0 - to->z0);

	// 2. now convert from cartesian to new system

	// cartesian -> cylindrical
	if (to->coord & CYLINDRICAL) {
		axis = (to->coord & ~CYLINDRICAL);
		if (axis == x_) {
			tmp = f1;
			f1  = f2;
			f2  = f3;
			f3  = tmp;
		}
		else if (axis == y_) {
			tmp = f2;
			f1  = f3;
			f2  = f1;
			f3  = tmp;
		}

		f2 = f2-capR;
		f1 = f1-capR;

		t1 = sqrt (f1*f1 + f2*f2);
		t2 = atan (f2/f1);
		t3 = f3;
		if (f1 < 0) t2  = pi + t2;
		if (t2 < 0) t2 += 2*pi;
	}

	// cartesian -> spherical
	else if (to->coord & SPHERICAL) {
		t1 = sqrt (f1*f1 + f2*f2 + f3*f3);
		t2 = atan (f2/f1);
		if (f1 < 0) t2  = pi + t2;
		if (t2 < 0) t2 += 2*pi;
		if (t1 > 0)
			t3 = acos (f3/t1);
		else
			t3 = 0;
	}

	// return transformed coordinates
	to->c1 = t1;
	to->c2 = t2;
	to->c3 = t3;
}


/// Calculate the square of the distance between two particles. This function can
/// also handle difference metrics, for example calculate the curvilinear distance
/// in cylindrical coordinates of same radial distance. By default, the function
/// computes the CARTESIAN distance between the particles.
///
/// @param		p1 pointer to the first particle;
/// @param		p2 pointer to the second particle;
/// @param		coord the type of distance to compute, CARTESIAN (default), CYLINDRICAL or SPHERICAL
/// @warning	This function DOES NOT compute the real curvilinear distance,
///				since that would require a numerical integration. For distances
///				on the cylinder, we assume that the two points lie at the SAME r.
/// @todo		Implement the spherical shell distance calculation.

//================================================================================
real DistanceSquared (particleMD *p1, particleMD *p2, int coord, real capR)
//================================================================================
{
	point	r0, r1, r2;
	real	dx, dy, dz;
	real	R, d2=0;

	// setup points
	memset (&r0, 0, sizeof(point));
	memset (&r1, 0, sizeof(point));
	memset (&r2, 0, sizeof(point));
	r0.coord = CARTESIAN;
	r1.coord = CARTESIAN;
	r2.coord = CARTESIAN;
	if (coord & CYLINDRICAL || coord & SPHERICAL) {
		r1.coord = coord;
		r2.coord = coord;
	}

	// transform first point
	r0.c1 = p1->rx;
	r0.c2 = p1->ry;
	r0.c3 = p1->rz;
	CoordinateTransform (&r0, &r1, capR);

	// transform second point
	r0.c1 = p2->rx;
	r0.c2 = p2->ry;
	r0.c3 = p2->rz;
	CoordinateTransform (&r0, &r2, capR);

	// calculate distance
	if (coord & CYLINDRICAL) {
		R  = r1.c1;
		dx = 0;
		dy = r2.c2 - r1.c2;
		dz = r2.c3 - r1.c3;
		if (dy>pi)  dy = 2*pi-dy;
		if (dy<-pi) dy = 2*pi+dy;
		dy = R*dy;
		d2 = dx*dx + dy*dy + dz*dz;
	}
	else if (coord & SPHERICAL) {
	}
	else {
		dx = r2.c1 - r1.c1;
		dy = r2.c2 - r1.c2;
		dz = r2.c3 - r1.c3;
		d2 = dx*dx + dy*dy + dz*dz;
	}

	// return distance squared
	return (d2);
}


//================================================================================
void CopyStepCounters (int *dest, int *source)
//================================================================================
{
	// Copy the step counter information from source to destination using the
	// memcpy library call (source and destination may not overlap). Useful to
	// copy step counter parameters to other locations (e.g. histograms and scenes).

	memcpy (dest, source, (PHASE_COUNT+1)*sizeof(int));
}


//================================================================================
void FillStepCounters (int *dest, int count)
//================================================================================
{
	int i;

	for (i=1; i<=PHASE_COUNT+1; i++) {
		dest[i] = count;
	}
}


//================================================================================
int GroupToIndex (int group)
//================================================================================
{
	int index=0;

	while (group > 0) {
		index++;
		group = group >> 1;
	}
	return index;
}


//================================================================================
int GroupFromIndex (int index)
//================================================================================
{
	if (index <= 0) return 0;
	else return 1<<(index-1);
}


//================================================================================
char *GroupToStr (int group, char *str)
//================================================================================
{
	// Convert INDIVIDUAL group bits into a literal string

	if 		(group == GROUP_NONE)		snprintf (str, STRMAX, "GROUP_NONE");
	else if (group == GROUP_FLUID)		snprintf (str, STRMAX, "GROUP_FLUID");
	else if (group == GROUP_WALL)		snprintf (str, STRMAX, "GROUP_WALL");
	else if (group == GROUP_WALL_IN)	snprintf (str, STRMAX, "GROUP_WALL_IN");
	else if (group == GROUP_WALL_OUT)	snprintf (str, STRMAX, "GROUP_WALL_OUT");
	else if (group == GROUP_ION)		snprintf (str, STRMAX, "GROUP_ION");
	else if (group == GROUP_ION_POS)	snprintf (str, STRMAX, "GROUP_ION_POS");
	else if (group == GROUP_ION_NEG)	snprintf (str, STRMAX, "GROUP_ION_NEG");
	else if (group == GROUP_MONOMER)	snprintf (str, STRMAX, "GROUP_MONOMER");
	else if (group == GROUP_GRAFT)		snprintf (str, STRMAX, "GROUP_GRAFT");
	else								snprintf (str, STRMAX, "GROUP_UNKNOWN");
	return str;
}


//================================================================================
char *TypeToStr (int type, char *str)
//================================================================================
{
	// Convert a type number into a literal type string

	if 		(type == TYPE_WALL)			snprintf (str, STRMAX, "TYPE_WALL");
	else if	(type == TYPE_FLUID)		snprintf (str, STRMAX, "TYPE_FLUID");
	else if (type == TYPE_MONOMER)		snprintf (str, STRMAX, "TYPE_MONOMER");
	else								snprintf (str, STRMAX, "TYPE_UNKNOWN");
	return str;
}

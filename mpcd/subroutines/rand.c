///
/// @file
///
/// @brief Methods for handling random number generation.
///
/// Methods for handling random number generation in NIAMH-MPCD. Most methods are interfaces that call one of the two
/// underlying random number generators:
/// - The Mersenne Twister, taken from http://www.math.sci.hiroshima-u.ac.jp/~m-mat/MT/MT2002/CODES/mt19937ar.c.
/// - The xoshiro128++ algorithm, taken from https://prng.di.unimi.it/xoshiro128plusplus.c.
///
/// The choice of random number generator is controlled by whether `RNG_MERSENNE` is defined or not. If not (the
/// default), then xoroshiro is used. Xoroshiro is faster, but has a much shorter period than the Mersenne Twister,
/// albeit still sufficient for the purposes of MPCD.
///
/// Base methods using the Mersenne Twister are prefixed with `MT_`, while those using xoroshiro are prefixed with `X_`.
/// Methods with no prefix are those freely available to the user, and will call the appropriate underlying method.
///

# include <math.h>
# include <sys/time.h>
# include <stdio.h>
# include <unistd.h>
# include <stdint.h>

# include "../headers/definitions.h"
# include "../headers/SRDclss.h"
# include "../headers/mtools.h"
# include "../headers/pout.h"

/* ****************************************** */
/* ****************************************** */
/* ****************************************** */
/* ************* MERSENNE TWISTER *********** */
/* ****************************************** */
/* ****************************************** */
/* ****************************************** */

/// @brief MT variable: The state vector of the RNG.
static unsigned long mt[NN];
/// @brief MT variable: MT counter variable. `mti==NN+1` means mt[NN] is not initialized.
static int mti=NN+1;

///
/// @brief Generate a random seed, if necessary, and initialize the Mersenne Twister random number generator.
///
/// Generates a random seed (if necessary), then uses this to initialise the state of the Mersenne Twister random number
/// generator. Initialisation is near identical to the method MT_init_genrand(), but for legacy reasons does not call
/// it explicitly.
///
/// @param seed Input seed. If zero, a random seed is generated using the time in microseconds add the process id.
/// @see MT_init_genrand()
/// @return The seed used to initialise the random number generator. Primarily to allow the user to check what seed was generated.
///
unsigned long MT_RandomSeedSRD (unsigned long seed)
{
  // Adapted from Fred Tessier code

  // Get a seed from time+pid if seed=0
  struct timeval tv;
  gettimeofday(&tv, NULL); // Get the time to use the microseconds as an "random" seed
  if (!seed) seed = tv.tv_usec+getpid();

  // Initialize mersenne twister array
  mt[0]= seed & 0xffffffff;
  for (mti=1; mti<NN; mti++) {
      mt[mti] = (1812433253UL * (mt[mti-1] ^ (mt[mti-1] >> 30)) + mti);
      mt[mti] &= 0xffffffffUL;
  }

  return (seed);
}
///
/// @brief Initialize the Mersenne Twister random number generator with a given seed.
///
/// Initialize the state of the Mersenne Twister random number generator with a given seed. Note that this explicitly
/// uses the given seed, unlike MT_RandomSeedSRD() which will generate a random seed if the input is zero.
///
/// @param s The seed to use to initialise the random number generator.
/// @see MT_RandomSeedSRD()
///
void MT_init_genrand(unsigned long s){
  // Mersenne twister
  // http://www.math.sci.hiroshima-u.ac.jp/~m-mat/MT/MT2002/CODES/mt19937ar.c.
  // Initializes mt[NN] with a seed.

  mt[0]= s & 0xffffffffUL;
  for (mti=1; mti<NN; mti++) {
    mt[mti] =
            (1812433253UL * (mt[mti-1] ^ (mt[mti-1] >> 30)) + mti);
    /* See Knuth TAOCP Vol2. 3rd Ed. P.106 for multiplier. */
    /* In the previous versions, MSBs of the seed affect   */
    /* only MSBs of the array mt[].                        */
    /* 2002/01/09 modified by Makoto Matsumoto             */
    mt[mti] &= 0xffffffffUL;
    /* for >32 bit machines */
  }
}
///
/// @brief Generates a random int32 number on the [0,0xffffffff]-interval using Mersenne Twister.
///
/// This performs a twist transformation to the state vector to generate a new random number. It first performs a sanity
/// check to see if the seed has not been properly set (indicated by the state variable `mti`), before performing the
/// twist. It finally uses the state vector to generate a random number.
///
/// @return The pseudo-randomly generated number.
///
unsigned long MT_genrand_int32(void){
  // Mersenne twister
  // http://www.math.sci.hiroshima-u.ac.jp/~m-mat/MT/MT2002/CODES/mt19937ar.c
  // Generates a random number on [0,0xffffffff]-interval

  unsigned long y;
  static unsigned long mag01[2]={0x0UL, MATRIX_A};
  struct timeval tv;

  /* mag01[x] = x * MATRIX_A  for x=0,1 */
  if (mti >= NN) { /* generate NN words at one time */
      int kk;

      if (mti == NN+1) {  /* if init_genrand() has not been called, */
          gettimeofday(&tv, NULL); // Get the time to use the microseconds as an "random" seed
          MT_init_genrand(tv.tv_usec);
          //init_genrand(5489UL); /* a default initial seed is used */
      }
      for (kk=0;kk<NN-MM;kk++) {
          y = (mt[kk]&UPPER_MASK)|(mt[kk+1]&LOWER_MASK);
          mt[kk] = mt[kk+MM] ^ (y >> 1) ^ mag01[y & 0x1UL];
      }
      for (;kk<NN-1;kk++) {
          y = (mt[kk]&UPPER_MASK)|(mt[kk+1]&LOWER_MASK);
          mt[kk] = mt[kk+(MM-NN)] ^ (y >> 1) ^ mag01[y & 0x1UL];
      }
      y = (mt[NN-1]&UPPER_MASK)|(mt[0]&LOWER_MASK);
      mt[NN-1] = mt[MM-1] ^ (y >> 1) ^ mag01[y & 0x1UL];

      mti = 0;
  }
  y = mt[mti++];
  /* Tempering */
  y ^= (y >> 11);
  y ^= (y << 7) & 0x9d2c5680UL;
  y ^= (y << 15) & 0xefc60000UL;
  y ^= (y >> 18);

  return y;
}

/* ****************************************** */
/* ****************************************** */
/* ****************************************** */
/* *************** xoshiro128++ ************* */
/* ****************************************** */
/* ****************************************** */
/* ****************************************** */

/*  Taken from: https://prng.di.unimi.it/xoshiro128plusplus.c

  Written in 2019 by David Blackman and Sebastiano Vigna (vigna@acm.org)

  To the extent possible under law, the author has dedicated all copyright
  and related and neighboring rights to this software to the public domain
  worldwide. This software is distributed without any warranty.

  See <http://creativecommons.org/publicdomain/zero/1.0/>. */

/// @brief Xoroshiro variable: the state vector of the RNG.
static unsigned long X_state[4];
/// @brief Xoroshiro variable: flag to show whether this has been seeded or not.
int X_seeded = 0;

///
/// @brief Rotates a 32-bit integer left by a given number of bits.
///
/// This is a helper function for the xoshiro128++ RNG. It rotates a 32-bit integer left by a given number of bits,
/// which is used several times in the generation method.
///
/// @param x The 32-bit integer to rotate.
/// @param k The number of bits to rotate by.
/// @return The rotated32 integer.
///
static inline unsigned long X_rotl(const long int x, int k) {
  return (x << k) | (x >> (32 - k));
}

///
/// @brief Initialise the xoshiro128++ RNG with a given seed using SplitMix64.
///
/// Initialises the xoroshiro RNG. In order to do this, an instance of the SplitMix64 algorithm is initialised with the
/// given seed. 4 usages of SplitMix64 are then applied to generate the initial state of xoroshiro from the seed.
///
/// This additional step is necessary because xoroshiro requires 4 pseudo-random 32-bit integers to initialise the
/// state, rather than a single state variable (as in Mersenne Twister).
///
/// @param s The seed to initialise with.
///
void X_init_genrand(unsigned long s) {
  int i; // counting variable

  /* SplitMix64 code taken from: https://github.com/svaarala/duktape/blob/master/misc/splitmix64.c
   * Written in 2015 by Sebastiano Vigna (vigna@acm.org)
      To the extent possible under law, the author has dedicated all copyright
      and related and neighboring rights to this software to the public domain
      worldwide. This software is distributed without any warranty.
      See <http://creativecommons.org/publicdomain/zero/1.0/>. */

  unsigned long sm_state = s; // splitmix64 state

  for (i = 0; i < 4; i++) {
      // generate the next value in splitmix
      unsigned long z = (sm_state += UINT64_C(0x9E3779B97F4A7C15));
      z = (z ^ (z >> 30)) * UINT64_C(0xBF58476D1CE4E5B9);
      z = (z ^ (z >> 27)) * UINT64_C(0x94D049BB133111EB);

      X_state[i] = z ^ (z >> 31); // output value from splitmix
  }

  X_seeded = 1; // mark as seeded
}

///
/// @brief Performs a pre-processing step to initialise the Xoroshiro RNG.
///
/// If the seed is set as 0 then a new seed is generated using the current time and process id. A safety check is then
/// performed to ensure that xoroshiro has not been properly seeded in X_init_genrand(), and if it has not then it
/// is initialised with the new seed.
///
/// @param seed A single 32-bit integer to use as the seed. If set as zero, a new seed is generated.
/// @see X_init_genrand()
/// @return The seed used to initialise the random number generator. Primarily to allow the user to check what seed was generated.
///
unsigned long X_RandomSeedSRD (unsigned long seed) {
  // Perform the same seeding as in the MT

  struct timeval tv;
  gettimeofday(&tv, NULL); // Get the time to use the microseconds as an "random" seed
  if (!seed) seed = tv.tv_usec+getpid();

  if (X_seeded == 0) {
      X_init_genrand(seed);
  }

  return seed;
}

///
/// @brief Generates a random 32-bit integer using the xoshiro128++ RNG.
///
/// Generates a random 32-bit integer using the xoshiro128++ RNG. This is a 32-bit RNG with a period of 2^128-1. First,
/// it ensures the RNG has been properly initialised, and if it has not then it is initialised with a new seed. It
/// proceeds to perform a rotation using the state, giving the next random number. The state is then updated per the
/// algorithm, doing a bitwise <b>XO</b>R, a <b>ro</b>tation, a <b>sh</b>ift, and a <b>ro</b>tation - Hence the name
/// Xoroshiro.
///
/// @return A random 32-bit integer.
///
unsigned long X_genrand_int32(void) {
  if (X_seeded == 0) { // ensure seed is properly set, if not then seed with a random value
      X_RandomSeedSRD(0);
  }

  const unsigned long result = X_rotl(X_state[0] + X_state[3], 7) + X_state[0];
  const unsigned long t = X_state[1] << 9;

  X_state[2] ^= X_state[0];
  X_state[3] ^= X_state[1];
  X_state[1] ^= X_state[2];
  X_state[0] ^= X_state[3];

  X_state[2] ^= t;

  X_state[3] = X_rotl(X_state[3], 11);

  return result;
}

/* ****************************************** */
/* ****************************************** */
/* ****************************************** */
/* ************ interface methods *********** */
/* ****************************************** */
/* ****************************************** */
/* ****************************************** */
///
/// @brief Takes the input seed and checks to see if a random seed needs to be generated. Proceeds to initialise the RNG.
///
/// Interface method. Takes the input seed and checks to see if a random seed needs to be generated. Proceeds to
/// initialise the RNG.
///
/// @param seed Seed used to initialise the generators. If set to zero, a new seed is generated.
/// @return The seed used to initialise the RNG.
///
unsigned long RandomSeedSRD (unsigned long seed)
{
  #ifdef RNG_MERSENNE
    return MT_RandomSeedSRD(seed);
  #else
    return X_RandomSeedSRD(seed);
  #endif
}

///
/// @brief Generates a random 32-bit integer using the RNG.
///
/// Interface method. Initialises the RNG without adjusting the seed. See RandomSeedSRD() for a similar method but
/// also generates a seed.
///
/// @param s The seed to use to initialise the RNG.
/// @see RandomSeedSRD()
///
void init_genrand(unsigned long s){
  #ifdef RNG_MERSENNE
    MT_init_genrand(s);
  #else
    X_init_genrand(s);
  #endif
}

///
/// @brief Generates a random 32-bit integer using the given RNG.
///
/// Interface method. Generates a random 32-bit integer using the given RNG. Used as the base method for all other
/// random number generation methods.
///
/// @return The generated random 32-bit integer.
///
unsigned long genrand_int32(void){
  //Base RNG method. Returns a random unsigned long.
  #ifdef RNG_MERSENNE
    return MT_genrand_int32();
  #else
    return X_genrand_int32();
  #endif
}

///
/// @brief Generates a random 31-bit integer by generating a random 32-bit integer and then shifting it.
///
/// Performs a bit-shift on a randomly generated 32-bit integer, to create a 31-bit integer.
///
/// @see genrand_int32()
/// @return A random 31-bit integer.
///
long genrand_int31(void){
  /*
  Mersenne twister
  http://www.math.sci.hiroshima-u.ac.jp/~m-mat/MT/MT2002/CODES/mt19937ar.c
  Generates a random number on [0,0x7fffffff]-interval
  */
  return (long)(genrand_int32()>>1);
}

///
/// @brief Generates a real number on the [0, 1) interval.
///
/// Generates a random 32-bit integer on the [0, 0xffffffff] interval, and then divides it by 2^32 to create a real
/// number on the [0, 1) interval.
///
/// Uses a constant variable `divisor` to control the precision real number - A higher number will give a higher
/// precision double.
///
/// @see genrand_int32()
/// @return The generated real number.
///
double genrand_real(void){
  const double divisor = 4294967296.0; // 2^32 by default
  return (genrand_int32() % (unsigned long) divisor) * (1.0/divisor); // modulo ensures this is always bounded
}

///
/// @brief Generates +1.0 or -1.0 with equal probability.
///
/// Generates a random real and maps this to either +1.0 or -1.0 with equal probability.
///
/// @see genrand_real()
/// @return The generated +1.0 or -1.0.
///
double genrand_pmOne(void){
	double rand = genrand_real();
	if(rand<=0.5) return -1.0;
	else return 1.0;
}

///
/// @brief Generate a random unit vector uniformly distributed about a cone around the x-axis.
///
/// Pivots off the `dimension` to uniformly generate a unit vector within a cone distribution:
/// - 3D: Sample the z component on [cos(theta), 1] and then sample the azimuthal angle on [0, 2pi).
/// - 2D: Simply sample the azimuthal angle homogeneously.
/// - 1D: Checks if cos(phi) would be parallel or anti-parallel, returning +1.0 or -1.0.
///
/// @param vec The vector to store the result in. Used as a return variable. Must be the same size as `dimension`.
/// @param theta The angle of the cone in radians.
/// @param dimension The dimension of the cone to generate in. Must match the dimension of `vec`.
/// @see genrand_real()
///
void genrand_coneNP( double vec[],double theta,int dimension ) {
	// Generate a random, uniformly distributed normalized vector for a direction within cone around "north pole"/x-axis

	double z,phi,ct;

	ct=cos(theta);
	if( dimension==_3D ) {
		// Sample z on [cos(theta),1] and phi on [0,2*pi]
		z=genrand_real()*(1.0-ct)+ct;
		phi=genrand_real()*2.0*pi;
		vec[0]=z;
		vec[1]=sqrt(1.0-z*z)*cos(phi);
		vec[2]=sqrt(1.0-z*z)*sin(phi);
	}
	else if( dimension==_2D ) {
		// Sample angle homogeneously
		phi=theta*(1.-2.*genrand_real());
		vec[0]=cos(phi);
		vec[1]=sin(phi);
	}
    else if( dimension==_1D ) {
		// A bit of a funny definition of a cone in 1D.
        // Checks if cos(phi) 2D would be parallel or antiparallel
		phi=theta*0.5*pi*(1.-2.*genrand_real());
        if(phi>1.0 || phi<-1.0) vec[0]=-1.0;
        else vec[0]=1.0;
	}
	else printf("Warning: genrand_coneNP() only programmed for DIM={3,2,1}, not DIM=%d\n",dimension);
}

///
/// @brief Generate a random unit vector uniformly distributed about a cone around a specified axis.
///
/// Generates a cone around the x-axis using genrand_coneNP(), then perform a rodrigues rotation to rotate the cone.
///
/// @param axis The axis for the cone to move around. Must be 3D.
/// @param vecOut The vector to store the result in. Used as a return variable. Must be 3D.
/// @param theta The angle of the cone in radians.
/// @param dimension The dimension we're working in.
/// @see genrand_coneNP()
///
void genrand_cone( double axis[],double vecOut[],double theta,int dimension ) {
	double rotAx[_3D],randVec[_3D]={0.0},xaxis[_3D]={0.0};
	double angle;
	int i;
	xaxis[0]=1.0;

	//Generate random vec about z-axis/north pole
	genrand_coneNP( randVec,theta,dimension );
	//Rotate north pole to align with the axis direction
	//Find the axis that the x-axis must be rotated about
	crossprod( axis,xaxis,rotAx );
	angle = absAngle( axis,xaxis,_3D );
	//Even if 2D, this must be 3D cuz then rotAX will be in 3rd dimension.
	//This rotation is why 3D vectors were required
	rodriguesRotation( randVec,rotAx,angle );
	for( i=0; i<dimension; i++ ) vecOut[i]=randVec[i];
}

/* ****************************************** */
/* ****************************************** */
/* ****************************************** */
/* ************** GAUSSIAN DIST ************* */
/* ****************************************** */
/* ****************************************** */
/* ****************************************** */

///
/// @brief Generate a random number from a Gaussian distribution of mean 0 and standard deviation 1.
///
/// Uses a Box-Muller transformation to turn a uniform random number on [0,1) into a Gaussian random number. Taken from
/// http://www.taygeta.com/random/gaussian.html.
///
/// @see genrand_real()
/// @return The randomly generated Gaussian.
///
float genrand_gauss( void ) {
/*
   Box-Muller transformation to turn a uniform random number
   between 0-1 into a gaussian distribution of mean 0 and
   StDev of 1. From
   http://www.taygeta.com/random/gaussian.html
*/
	float x1,x2,w,y1;
  // 	float y2;
	do{
		x1 = 2. * genrand_real() - 1.;
		x2 = 2. * genrand_real() - 1.;
		w = x1*x1+x2*x2;
	} while( w >= 1. );
	w = sqrt( (-2. * log( w )) / w );
	y1 = x1*w;
	return y1;
}

///
/// @brief Generate a Maxwell-Boltzmann distributed random number (by scaling a Gaussian by `KBT/M`).
///
/// Generates a Gaussian random number and scales it by `KBT/M` to get a Maxwell-Boltzmann distributed random number.
///
/// @param KBT The temperature in MPCD units of energy.
/// @param M The mass in MPCD units of mass.
/// @see genrand_gauss()
/// @return The Maxwell-Boltzmann generated random number.
///
double genrand_gaussMB(double KBT,double M) {
	double sigma=sqrt(KBT/M);
	return sigma*genrand_gauss();
}

///
/// @brief Generate a random number from a Gaussian distribution of mean `mu` and standard deviation `sigma`.
///
/// Rescales a standard Gaussian distributed number by `sigma` and adds `mu` to it, replicating generating from the
/// relevent Gaussian.
///
/// @param mu The mean of the Gaussian distribution.
/// @param sigma The standard deviation of the Gaussian distribution.
/// @see genrand_gauss()
/// @return The randomly generated Gaussian.
///
double genrand_gaussGen(double mu,double sigma) {
	return sigma*genrand_gauss()+mu;
}

///
/// @brief Generate a random number from an exponential distribution with mean `lambda` using Box-Muller.
///
/// Generates a random real number, and performs a Box-Muller transformation to get a Gaussian distributed random
/// number.
///
/// @param lambda The mean of the exponential distribution.
/// @return The randomly generated number from the exponential distribution.
///
double genrand_exp(double lambda) {
	return -lambda*log( genrand_real() );
}

///
/// @brief Generate a random integer from a Poisson distribution with mean `lambda`.
///
/// Pivots off small or large `lambda` to generate a poisson distributed random number:
/// - For small `lambda`, repeatedly multiply random numbers until they reach `exp(-lambda)`, returning the amount of
/// numbers generated.
/// - For large `lambda`, uses a similar `STEP` method so that `exp(-STEP)` doesn't underflow.
///
/// Algorithm taken from Junhao, based on Knuth. Explanation available at
/// https://en.wikipedia.org/wiki/Poisson_distribution#Computational_methods.
///
/// @param lambda The mean of the Poisson distribution.
/// @see genrand_real()
/// @return A random integer generated from the Poisson distribution.
///
int genrand_poisson(double lambda) {
	double L,eSTEP,STEP=50.0;	//STEP is chosen for double precision
	double r,myExp=M_E,p=1.0;
	int k=0;

	// For "small" values of lambda
	if( lambda<STEP ) {
		L=exp(-lambda);
		do{
			k+=1;
			r = genrand_real();
			p*=r;
		}while( p>L );
	}
	// For "large" values of lambda the exponent is not well known
	else if( lambda<10*STEP ) {
		L=lambda;
		eSTEP=exp(STEP);
		do{
			k+=1;
			r = genrand_real();
			p*=r;
			if( p<myExp && L>0.0 ) {
				if( L>STEP ){
					p*=eSTEP;
					L-=STEP;
				}
				else{
					p*=L;
					L=-1.0;
				}
			}
		}while( p>1.0 );
	}
	// For large values of lambda just use a normal distribution approximation
	else {
		L=genrand_gaussGen( lambda,sqrt(lambda) );
		k=1+(int)L;
	}
	return k-1;
}
///
/// @brief Generate a random number from a Rayleigh (`xe^{-x^2}`) distribution with standard deviation `std`.
///
/// Found from the inverse-transform sampling method from a uniformly generated random number. More info here:
/// https://en.wikipedia.org/wiki/Rayleigh_distribution#Generating_random_variates.
///
/// @param std The standard deviation of the distribution.
/// @see genrand_real()
/// @return The randomly generated number from the Rayleigh distribution.
///
float genrand_rayleigh( float std ) {
	float r,lnr;	// Input uniform random number
	float x;	// Distributed random number

	// Uniform Random number
	do{
		r = genrand_real();
	} while ( feq(r,0.0) );
	lnr = log( r );
	x = std * sqrt( -2.*lnr );

	return x;
}

///
/// @brief Generate a random normalised vector uniformly distributed on a sphere.
///
/// Pivots off the `dimension` to uniformly generate a unit vector within a sphere:
/// - 3D: Generate the z component randomly between -1 and 1, and then generate a random angle in [0, pi). Use the angle
/// to find the x and y components.
/// - 2D: Simply generate a random angle and use that to find the x and y components.
/// - 1D: Wraps around genrand_pmOne() which is equivalent.
///
/// @param vec The vector to return the random normalised vector in. Must have dimension of `dimension`.
/// @param dimension The dimension to be calculated on. Must match dimension of `vec`.
/// @see genrand_real()
/// @see genrand_pmOne()
///
void genrand_sphere( double vec[],int dimension ) {
// 	int i;
// 	for( i=0; i<dimension; i++ ) vec[i]=genrand_gauss();
// 	norm( vec,dimension );

	double u,theta,t1;

	if( dimension==_3D ) {
		u=genrand_pmOne()*genrand_real();
		theta=genrand_real()*pi;
		t1=sqrt(1.-u*u);
		vec[0] = t1*cos(theta);
		vec[1] = t1*sin(theta);
		vec[2] = u;
	}
	else if( dimension==_2D ) {
		theta=2*genrand_real()*pi;
		vec[0]=cos(theta);
		vec[1]=sin(theta);
	}
	else if( dimension==_1D ) vec[0]=genrand_pmOne();
}
/* ****************************************** */
/* ****************************************** */
/* ****************************************** */
/* ************* RANDOM VECTORS ************* */
/* ****************************************** */
/* ****************************************** */
/* ****************************************** */

///
/// @brief Generate a non-normalised vector where the components are randomly uniformly distributed in the range [0,1).
///
/// Generates a random vector based off `dimension`. Each component is generated using genrand_real().
///
/// @param v The vector to return the random vector in. Must have dimension of `dimension`.
/// @param doShift Integer flag on whether to do the randomness operation. Set to 0 to disable and return a vector of 0s.
/// @param dimension The dimension to be calculated on. Must match dimension of `v`.
/// @see genrand_real()
/// @return The randomly generated vector, identical to `v`.
///
double *ranshift( double *v,int doShift,int dimension ) {
	double x,y,z;
	if( doShift ) {
		if( dimension >= _3D ) z = genrand_real();
		else z = 0.0;
		if( dimension >= _2D ) y = genrand_real();
		else y=0.0;
		x = genrand_real();
		v[0] = x;
		v[1] = y;
		v[2] = z;
	}
	else {
		v[0] = 0.0;
		v[1] = 0.0;
		v[2] = 0.0;
	}
	return v;
}

///
/// @brief Generate a normalised vector in 3D the components are randomly uniformly distributed in the range [-1,1).
///
/// Generates pairs of random numbers in the range [0,1) until they have a magnitude of less than 1. Uses this to create
/// the third vector component, re-normalising the x and y components before returning.
///
/// @param v The randomly generated vector to be returned.
/// @see genrand_real()
/// @return The randomly generated vector. Identical to `v`.
///
double *ranvec3D( double *v ) {
	double x,y,s,r;
	do{
		r = genrand_real();
		x = 2.*r - 1.;
		r = genrand_real();
		y = 2.*r - 1.;
		s = x*x + y*y;
	}while( s > 1. );
	v[2] = 1. - 2.*s;
	s = 2.*sqrt( 1.-s );
	v[0] = s * x;
	v[1] = s * y;
	return v;
}

///
/// @brief Generate a normalised vector in 2D the components are randomly uniformly distributed in the range [-1,1).
///
/// Generates a random number in the range [0,1), then uses this to create the second vector component. Similar to
/// ranvec3D().
///
/// @param v The randomly generated vector to be returned.
/// @see genrand_real()
/// @see ranvec3D()
/// @return The randomly generated vector. Identical to `v`.
///
double *ranvec2D( double *v ) {
	double x,y;
	x = genrand_real();
	y = sqrt( 1.-x*x );
	v[0] = x;
	v[1] = y;
	return v;
}

///
/// @brief Generates a random normalised vector in 2D or 3D. Uses ranvec2D() or ranvec3D() respectively.
///
/// Generates a random normalised vector in 2D or 3D. Uses ranvec2D() or ranvec3D() respectively.
///
/// @param v The randomly generated vector to be returned. Must have dimension of `dimension`.
/// @param dimension The dimension to be calculated on. Must match dimension of `v`, and should be either 2 or 3.
/// @see ranvec2D()
/// @see ranvec3D()
/// @return The randomly generated vector. Identical to `v`.
///
double *ranvec( double *v,int dimension ) {
	if( dimension==_3D ) ranvec3D( v );
	if( dimension==_2D ) ranvec2D( v );
	return v;
}

///
/// @brief Randomly picks an integer corresponding to one of all MPCD particles.
///
/// Randomly picks an integer corresponding to one of all MPCD particles.
///
/// @param POP Total population of MPCD particles.
/// @see genrand_real()
/// @return The randomly generated integer corresponding to one of all MPCD particles.
///
int rand_particle( int POP ) {
	return (int)( genrand_real()*(double)POP );
}

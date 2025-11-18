
//================================================================================
//
// name:   mdpolymer.c
// author: ftessier
// date:   2005-05-03 @ 11:04:36
//
//================================================================================
///
/// @file md.c
/// Polymer properties measurements for the Molecular Dynamics program.
///
//================================================================================


#include <math.h>
#include <string.h>
#include "mderror.h"
#include "mdtypes.h"
#include "mdfiles.h"
#include "mdpolymer.h"


/// Measure the geometrical properties of all polymers in the system (position of
/// the center of mass, radius of gyration and end to end distance).
///
/// @param		sim a pointer to a simulation structure
/// @return		void
/// @warning	The grafted chain angle calculations in this function are hardwired
/// 			for the case of a cylindrical capillary
/// @todo		Change the angle calculations to be geometry dependent

//================================================================================
void PolymerGeometry (simptr sim)
//================================================================================
{
	int			i, j, k, nPolymer;
	itemPoly 	*polymer;
	particleMD	*p1, *p2;
	real		dx, dy, dz;
	real		m, mtot;
	real		sumx, sumy, sumz;
	real		cmx, cmy, cmz;
	real		rg, rgx, rgy, rgz, rh;
	real		h,  hx,  hy,  hz, hr;
	real		x0, y0, z0;
	real		yn, zn, rn;
	real		tx, ty, tz;
	real		a, b, c2;
	real 		theta, thetamean, thetasig, phi, phimean, phisig;
	int			nangle;

	// local sim variables
	polymer  = sim->polymer.items;
	nPolymer = sim->polymer.n;

	// write simulation time in data files
	fprintf (GetSimStream(sim->files,"POLYMER-rh"),  				"%.3f ", sim->tNow);
	fprintf (GetSimStream(sim->files,"POLYMER-cmx"),  				"%.3f ", sim->tNow);
	fprintf (GetSimStream(sim->files,"POLYMER-cmy"),  				"%.3f ", sim->tNow);
	fprintf (GetSimStream(sim->files,"POLYMER-cmz"),  				"%.3f ", sim->tNow);
	fprintf (GetSimStream(sim->files,"POLYMER-rg"),   				"%.3f ", sim->tNow);
	fprintf (GetSimStream(sim->files,"POLYMER-rgx"),  				"%.3f ", sim->tNow);
	fprintf (GetSimStream(sim->files,"POLYMER-rgy"),  				"%.3f ", sim->tNow);
	fprintf (GetSimStream(sim->files,"POLYMER-rgz"),  				"%.3f ", sim->tNow);
	fprintf (GetSimStream(sim->files,"POLYMER-h"),    				"%.3f ", sim->tNow);
	fprintf (GetSimStream(sim->files,"POLYMER-hx"),   				"%.3f ", sim->tNow);
	fprintf (GetSimStream(sim->files,"POLYMER-hy"),   				"%.3f ", sim->tNow);
	fprintf (GetSimStream(sim->files,"POLYMER-hz"),   				"%.3f ", sim->tNow);
	fprintf (GetSimStream(sim->files,"POLYMER-hr"),   				"%.3f ", sim->tNow);
	fprintf (GetSimStream(sim->files,"POLYMER-angle-theta"),		"%.3f ", sim->tNow);
	fprintf (GetSimStream(sim->files,"POLYMER-angle-theta-mean"),  	"%.3f ", sim->tNow);
	fprintf (GetSimStream(sim->files,"POLYMER-angle-phi"), 			"%.3f ", sim->tNow);
	fprintf (GetSimStream(sim->files,"POLYMER-angle-phi-mean"), 	"%.3f ", sim->tNow);
	fprintf (GetSimStream(sim->files,"POLYMER-Velocity"),       "\ntimestep indexed\n");
	fprintf (GetSimStream(sim->files,"POLYMER-Acceleration"),   "\ntimestep indexed\n");
	// initialize variables
	thetamean = 0;
	thetasig  = 0;
	phimean	  = 0;
	phisig 	  = 0;
	nangle    = 0;

	// loop over polymers
	for (i=0; i<nPolymer; i++) {

		// pointer to first monomer
		p1 = polymer[i].p1;
		if (!p1) continue;

		// skip grafting point
		if (polymer[i].grafted) p1 = p1->next;
		if (!p1) continue;

		// Print velocity and acceleration
		j = 0;
		while (p1) {
			fprintf (GetSimStream(sim->files,"POLYMER-Velocity"),  "%d %.6G %.6G %.6G", j, p1->vx, p1->vy, p1->vz);
			fprintf (GetSimStream(sim->files,"POLYMER-Acceleration"),  "%d %.6G %.6G %.6G", j, p1->ax, p1->ay, p1->az);
			fprintf (GetSimStream(sim->files,"POLYMER-Velocity"), "\n");
			fprintf (GetSimStream(sim->files,"POLYMER-Acceleration"), "\n");
			p1 = p1->next;
			j = j+1;
		}

		// pointer to first monomer
		p1 = polymer[i].p1;
		if (!p1) continue;

		// skip grafting point
		if (polymer[i].grafted) p1 = p1->next;
		if (!p1) continue;

		// reset center of mass variables
		mtot = 0;
		sumx = 0;
		sumy = 0;
		sumz = 0;

		// calculate center of mass
		while (p1) {
			m = p1->mass;
			mtot += m;
			sumx += p1->wx * m;
			sumy += p1->wy * m;
			sumz += p1->wz * m;
			p1 = p1->next;
		}
		cmx = sumx / mtot;
		cmy = sumy / mtot;
		cmz = sumz / mtot;

		// write cm values
		fprintf (GetSimStream(sim->files,"POLYMER-cmx"),  "%.8G ", cmx);
		fprintf (GetSimStream(sim->files,"POLYMER-cmy"),  "%.8G ", cmy);
		fprintf (GetSimStream(sim->files,"POLYMER-cmz"),  "%.8G ", cmz);

		// reset radius of gyration and end-to-end distance variables
		mtot = 0;
		sumx = 0;
		sumy = 0;
		sumz = 0;

		// get pointer to first monomer
		p1 = polymer[i].p1;
		if (!p1) continue;

		// skip grafting point
		if (polymer[i].grafted) p1 = p1->next;
		if (!p1) continue;

		// init end-to end-distance
		hx = p1->wx;
		hy = p1->wy;
		hz = p1->wz;

		// calculate radius of gyration
		while (p1) {
			m  = p1->mass;
			mtot += m;
			dx = p1->wx - cmx;
			dy = p1->wy - cmy;
			dz = p1->wz - cmz;
			sumx += dx*dx * m;
			sumy += dy*dy * m;
			sumz += dz*dz * m;
			// end-to end-distance (for last monomer)
			if (!(p1->next)) {
				hx = fabs(p1->wx-hx);
				hy = fabs(p1->wy-hy);
				hz = fabs(p1->wz-hz);
			}
			p1 = p1->next;
		}
		rgx = sqrt(sumx/mtot);
		rgy = sqrt(sumy/mtot);
		rgz = sqrt(sumz/mtot);
		rg  = sqrt((sumx+sumy+sumz)/mtot);
		h   = sqrt (hx*hx+hy*hy+hz*hz);

		// special radial value of h for the capillary case
		hr = sqrt (hy*hy + hz*hz);

		// write rg values
		fprintf (GetSimStream(sim->files,"POLYMER-rg"),   "%.4G ", rg);
		fprintf (GetSimStream(sim->files,"POLYMER-rgx"),  "%.4G ", rgx);
		fprintf (GetSimStream(sim->files,"POLYMER-rgy"),  "%.4G ", rgy);
		fprintf (GetSimStream(sim->files,"POLYMER-rgz"),  "%.4G ", rgz);

		// write h values
		fprintf (GetSimStream(sim->files,"POLYMER-h"),    "%.4G ", h);
		fprintf (GetSimStream(sim->files,"POLYMER-hx"),   "%.4G ", hx);
		fprintf (GetSimStream(sim->files,"POLYMER-hy"),   "%.4G ", hy);
		fprintf (GetSimStream(sim->files,"POLYMER-hz"),   "%.4G ", hz);
		fprintf (GetSimStream(sim->files,"POLYMER-hr"),   "%.4G ", hr);

		// calculate tilt and azimutal angle (grafted polymers)
		if (polymer[i].grafted) {

			// pointer to first monomer
			p1 = polymer[i].p1;
			if (!p1) continue;

			// skip grafting point
			p1 = p1->next;
			if (!p1) continue;

			// first monomer position and normal unit vector (capillary)
			x0 = p1->rx;
			y0 = p1->ry;
			z0 = p1->rz;
			rn = sqrt(y0*y0+z0*z0);
			yn = -y0/rn;
			zn = -z0/rn;

			// cm vector
			tx = cmx-x0;
			ty = cmy-y0;
			tz = cmz-z0;
			c2 = (tx*tx+ty*ty+tz*tz);

			// normal and surface lengths (b and a)
			b  = ty*yn + tz*zn;
			a  = sqrt(c2-b*b);

			// tilt and azimutal angle
			theta = atan (a/b) * 180/pi;
			if (fabs(tx) > 0)
				phi = acos (tx/a) * 180/pi;
			else
				phi = 0.0;

			thetamean += theta;
			thetasig  += theta*theta;
			phimean   += phi;
			phisig 	  += phi*phi;
			nangle++;

			// write tilt values
			fprintf (GetSimStream(sim->files,"POLYMER-angle-theta"),	"%.4G ", theta);
			fprintf (GetSimStream(sim->files,"POLYMER-angle-phi"),		"%.4G ", phi);
		}

		// get pointer to first monomer
		p1 = polymer[i].p1;
		if (!p1) continue;

		// skip grafting point
		if (polymer[i].grafted) p1 = p1->next;
		if (!p1) continue;

		// calculate hydrodynamic radius (EQUAL MASSES ONLY!)
		rh = 0;
		k  = 0;
		while (p1) {
			p2 = p1->next;
			while (p2) {
				dx  = p2->wx - p1->wx;
				dy  = p2->wy - p1->wy;
				dz  = p2->wz - p1->wz;
				rh += 2 / sqrt (dx*dx + dy*dy + dz*dz);
				k++;
				p2 = p2->next;
			}
			p1 = p1->next;
		}
		rh = 1/(rh/k);

		// write rh values
		fprintf (GetSimStream(sim->files,"POLYMER-rh"),   "%.4G ", rh);
	}

	// average tilt and azimutal angles
	if (nangle > 0) {
		thetamean /= nangle;
		thetasig   = sqrt(thetasig/nangle - thetamean*thetamean);
		phimean	  /= nangle;
		phisig     = sqrt(phisig/nangle - phimean*phimean);
		fprintf (GetSimStream(sim->files,"POLYMER-angle-theta-mean"),   "%.4G %.4G ", thetamean, thetasig);
		fprintf (GetSimStream(sim->files,"POLYMER-angle-phi-mean"), 	"%.4G %.4G ", phimean, phisig);
	}

	// finish lines in data files
	fprintf (GetSimStream(sim->files,"POLYMER-rh"),  				"\n");
	fprintf (GetSimStream(sim->files,"POLYMER-cmx"),  				"\n");
	fprintf (GetSimStream(sim->files,"POLYMER-cmy"),  				"\n");
	fprintf (GetSimStream(sim->files,"POLYMER-cmz"),  				"\n");
	fprintf (GetSimStream(sim->files,"POLYMER-rg"),   				"\n");
	fprintf (GetSimStream(sim->files,"POLYMER-rgx"),  				"\n");
	fprintf (GetSimStream(sim->files,"POLYMER-rgy"),  				"\n");
	fprintf (GetSimStream(sim->files,"POLYMER-rgz"),  				"\n");
	fprintf (GetSimStream(sim->files,"POLYMER-h"),    				"\n");
	fprintf (GetSimStream(sim->files,"POLYMER-hx"),   				"\n");
	fprintf (GetSimStream(sim->files,"POLYMER-hy"),   				"\n");
	fprintf (GetSimStream(sim->files,"POLYMER-hz"),   				"\n");
	fprintf (GetSimStream(sim->files,"POLYMER-angle-theta"),   		"\n");
	fprintf (GetSimStream(sim->files,"POLYMER-angle-theta-mean"), 	"\n");
	fprintf (GetSimStream(sim->files,"POLYMER-angle-phi"),   		"\n");
	fprintf (GetSimStream(sim->files,"POLYMER-angle-phi-mean"), 	"\n");
}


/// Compute force exerted by the fisrt FENE bond in all polymer molecule. This can be
/// useful to calculate the drag force on grafted polymer coil, for example.
///
/// @param		sim a pointer to a simulation structure
/// @return		void
/// @warning	This function assumes that the monomers are linked with FENE springs

//================================================================================
void PolymerForce (simptr sim)
//================================================================================
{
	int			i, nPolymer;
	itemPoly 	*polymer;
	particleMD	*p1, *p2;
	real		dx, dy, dz, r2;
	real		fx, fy, fz, fMag;
	real		k, r0, r02;

	// local sim variables
	polymer  = sim->polymer.items;
	nPolymer = sim->polymer.n;
	k		 = sim->kFene;
	r0		 = sim->r0Fene;
	r02  	 = r0*r0;

	// write simulation time in data files
	fprintf (GetSimStream(sim->files,"POLYMER-f"),  "%.3f ", sim->tNow);
	fprintf (GetSimStream(sim->files,"POLYMER-fx"), "%.3f ", sim->tNow);
	fprintf (GetSimStream(sim->files,"POLYMER-fy"), "%.3f ", sim->tNow);
	fprintf (GetSimStream(sim->files,"POLYMER-fz"), "%.3f ", sim->tNow);

	// loop over polymers
	for (i=0; i<nPolymer; i++) {

		// get pointer to first monomer (grafted polymers only)
		if (!polymer[i].grafted) continue;
		p1 = polymer[i].p1;
		if (!p1) continue;
		p2 = p1->next;
		if (!p2) continue;

		// calculate FENE force
		dx = p2->wx - p1->wx;
		dy = p2->wy - p1->wy;
		dz = p2->wz - p1->wz;
		r2 = dx*dx + dy*dy + dz*dz;
		fMag = k/(1-r2/r02);
		fx	= fMag * dx;
		fy	= fMag * dy;
		fz	= fMag * dz;

		// write f values
		fprintf (GetSimStream(sim->files,"POLYMER-f"),   "%.4G ", fMag);
		fprintf (GetSimStream(sim->files,"POLYMER-fx"),  "%.4G ", fx);
		fprintf (GetSimStream(sim->files,"POLYMER-fy"),  "%.4G ", fy);
		fprintf (GetSimStream(sim->files,"POLYMER-fz"),  "%.4G ", fz);
	}

	// finish lines in data files
	fprintf (GetSimStream(sim->files,"POLYMER-f"),  "\n");
	fprintf (GetSimStream(sim->files,"POLYMER-fx"), "\n");
	fprintf (GetSimStream(sim->files,"POLYMER-fy"), "\n");
	fprintf (GetSimStream(sim->files,"POLYMER-fz"), "\n");
}


//================================================================================
void PolymerPersistence (simptr sim)
//================================================================================
{
	int			i, nPolymer;
	real		r1, r2, s;
	real		v1x, v1y, v1z;
	real		v2x, v2y, v2z;
	itemPoly 	*polymer;
	particleMD	*p1, *p2;

	// local sim variables
	polymer  = sim->polymer.items;
	nPolymer = sim->polymer.n;

	// persistence length
	for (i=0; i<nPolymer; i++) {
		p1 = polymer[i].p1;
		p2 = p1->next;
		s = 0;
		while (p1 && p2) {
			// vector 1
			v1x = p2->wx - p1->wx;
			v1y = p2->wy - p1->wy;
			v1z = p2->wz - p1->wz;
			r1  = sqrt(v1x*v1x + v1y*v1y + v1z*v1z);

			// arc length
			s = r1;
			while (p2 && p2->next) {

				// vector 2
				v2x = p2->next->wx - p2->wx;
				v2y = p2->next->wy - p2->wy;
				v2z = p2->next->wz - p2->wz;
				r2  = sqrt(v2x*v2x + v2y*v2y + v2z*v2z);

				// next bond
				p2 = p2->next;
				s += r2;
			}
			p1 = p1->next;
			p2 = p1->next;
		}
		printf("Polymer %d persistence length %lf\n",i,s);
	}
}

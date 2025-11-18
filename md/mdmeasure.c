
//================================================================================
//
// name:   mdmeasure.c
// author: ftessier
// date:   2005-05-03 @ 11:04:36
//
//================================================================================
///
/// @file mdmeasure.c
/// Data gathering functions. The major function Measure is called at every MD step
/// and is respondible for dispatching data gathering functions. Other functions
/// serve to increment step counters and adjust to the different simulation phases.
/// Basic properties (e.g. energy) can be handled through the AccumProperties
/// function, while more advanced data collection is probably best handled by
/// data histograms. Finally, other actions can be registered in the simulation
/// actions array to be triggered at specified intervals.
///
//================================================================================


#include <math.h>
#include <string.h>
#include <unistd.h>
#include <time.h>
#include "mderror.h"
#include "mdtypes.h"
#include "mdfiles.h"
#include "mdhistogram.h"
#include "mdscene.h"
#include "mdpolymer.h"
#include "mdmeasure.h"
#include "mdsetup.h"
#include "mdutil.h"


/// Performs measurements periodically during the simulation. First, a mark is
/// written to the log file periodically to monitor the progress of the simulation.
/// Then, the frequency of the various actions is determined by step counters
/// pertaining to the current simulation phase sim->phase. Basic measurements
/// (such as energy) are collected at each step with the Eval and Accum functions.
/// Hitograms and Scenes contain their own step counters, while other actions
/// are triggered by step counters in the sim structure. Basically
/// a counter is an array in which the element 0 is the current count, and where
/// the subsequent element i correspond to the period of the action trigger during
/// simulation phase i. When the current count is equal to the value of the
/// period, the action is triggered and the current count is reset to 0.
///
/// @param		sim a pointer to a simulation structure
/// @see		the function SetupStepCounters for how the step counter list is built
/// @return		void

//================================================================================
void Measure (simptr sim)
// void Measure (simptr sim,bc WALL)
//================================================================================
{
	int		 i;
	itemStep *step;

	// progress report
	if (sim->step[count_] % STEP_PROGRESS_MARK == 0) {
		LOG ("mark\n");
		fflush(0);
	}

	// evaluate and accumulate basic quantities
	EvalProperties  (sim);
	AccumProperties (sim, ACTION_COLLECT);

	// histograms
	HistogramCollect (sim, sim->histograms);
// 	HistogramCollect (sim, sim->histograms,WALL);
	HistogramPrint   (sim, sim->histograms);

	// scenes
	ScenePrint (sim, sim->scenes);

	// perform actions based on step counters (skip global step counter)
	for (i=1; i<sim->stepCounter.n; i++) {
		step = &(sim->stepCounter.items[i]);
		if (step->counter[count_] == step->counter[sim->phase]) {
			if(sim->warmupMD==POS_WARMUP){
				step->action ((void *) sim);
			}
			step->counter[count_] = 0;
		}
	}
}


/// Increment the step count for all step counters in the step list. There is a
/// list of step counters in the sim structure, as well as step counters
/// local to each histogram (collect and print) and scene (print only). This
/// function simply increments the element 0 of each of the counters, to
/// represent the number of steps taken since the last counter reset.
///
/// @param		sim a pointer to a simulation structure
/// @return		void
/// @warning	This function MUST be called for each MD step.

//================================================================================
void IncrementStepCounters (simptr sim)
//================================================================================
{
	int		 i;
	histptr  h;
	sceneptr s;

	// increment step counters
	for (i=0; i < sim->stepCounter.n; i++) {
		sim->stepCounter.items[i].counter[count_]++;
	}

	// increment histogram step counters
	h = sim->histograms;
	while (h) {
		if (h->active) {
			h->stepCollect [count_]++;
			h->stepPrint   [count_]++;
		}
		h = h->next;
	}

	// increment scene step counters
	s = sim->scenes;
	while (s) {
		if (s->active) {
			s->stepPrint   [count_]++;
		}
		s = s->next;
	}

	// update simulation md time
	sim->tNow = sim->step[count_] * sim->dt;

}


/// Checks if a new simulation phase is about to begin. The simulation is divided
/// in different phases, e.g., equilibration, production, etc. New phases can be
/// added by simply adding fields in the phase-related quantities in the simulation
/// input file, and adjusting the PHASE_COUNT define in the Makefile. This function
/// simply checks the main step counter to see if it has reached the number of steps
/// in the phase. If so, then a function is called to prepare the simulation for
/// the next phase. Note that the simulation will simply skip over phases that have
/// a stepcount <= 0, so that phases can easily be turned on or off in the input
/// file. This function is also responsible for setting the DONE flag to signal
/// the end of the simulation.
///
/// @param		sim a pointer to a simulation structure
/// @return		a flag to specify if the simulation is over or not
/// @warning	This function relies on a correct value of PHASE_COUNT, defined
/// 			in the Makefile.

//================================================================================
int SimulationPhaseCheck (simptr sim)
//================================================================================
{
	if (sim->step[count_] >= sim->step[sim->phase]) {
		do {
			sim->phase++;
			if (sim->phase > PHASE_COUNT) return DONE;
		} while (sim->nStep[sim->phase] <= 0);
		SimulationPhaseReset(sim);
	}

	return !DONE;
}


/// Prepares a new simulation phase. First, the end step of the new phase is calculated
/// by adding the requested phase steps in the input file to the current step count.
/// Next, all the step counters are reset to 0 (except the overall step counter). Then
/// the histogram and scene step counters are likewise reset, as weel as the data
/// accumulators. Finally, the thermostat groups are adjusted and information about
/// the new phase is printed in the log file.
///
/// @param		sim a pointer to a simulation structure
/// @return		void

//================================================================================
void SimulationPhaseReset (simptr sim)
//================================================================================
{
	int			i, nAtom;
	int			groupThermRescale, groupThermDPD;
	int			phase;
	particleMD	*atom;
	histptr		h;
	sceneptr	s;

	// local sim variables
	nAtom = sim->atom.n;
	atom  = sim->atom.items;
	phase = sim->phase;
	h	  = sim->histograms;
	s	  = sim->scenes;

	// calculate step parameter for this phase
	sim->step[phase] = sim->step[count_] + sim->nStep[phase];

	// reset simulation step counters (except global counter 0)
	for (i=1; i < sim->stepCounter.n; i++) {
		sim->stepCounter.items[i].counter[count_] = 0;
	}

	// reset histogram step counters
	while (h) {
		h->stepCollect [count_] = 0;
		h->stepPrint   [count_] = 0;
		h = h->next;
	}

	// reset scene step counters
	while (s) {
		s->stepPrint   [count_] = 0;
		s = s->next;
	}

	// reset accumulators
	AccumProperties (sim, ACTION_RESET);

	// update thermostat rescale group
	groupThermRescale 		= sim->groupThermRescale[phase];
	sim->nAtomThermRescale 	= 0;
	sim->kinETherm  		= 0;
	if (groupThermRescale > 0) {
		for (i=0; i<nAtom; i++) {
			if (atom[i].group & groupThermRescale) sim->nAtomThermRescale++;
		}
	}

	// update thermostat DPD group
	groupThermDPD			= sim->groupThermDPD[phase];
	sim->nAtomThermDPD		= 0;
	if (groupThermDPD > 0) {
		for (i=0; i<nAtom; i++) {
			if (atom[i].group & groupThermDPD) sim->nAtomThermDPD++;
		}
	}

	// report
	LOG ("--------------------------------------------------------------\n");
	LOG ("Starting simulation phase %d\n",	phase);
	LOG ("  %-30s = %010x\n", 	"groupThermRescale",	groupThermRescale);
	LOG ("  %-30s = %d\n", 		"nAtomThermRescale",	sim->nAtomThermRescale);
	LOG ("  %-30s = %010x\n", 	"groupThermDPD",		groupThermDPD);
	LOG ("  %-30s = %d\n", 		"nAtomThermDPD",		sim->nAtomThermDPD);
	LOG ("--------------------------------------------------------------\n");
}


/// Perform some additional transformation on the values collected during the
/// MD integration, in order to evaluate physical properties, e.g., divide
/// the calculated energies by the number of atoms to obtain the energy per
/// atom, etc. Keep this efficient as it is called every MD step.
///
/// @param		sim a pointer to a simulation structure
/// @return		void

//================================================================================
void EvalProperties (simptr sim)
//================================================================================
{
	// energies per atom
	sim->kinE /= sim->atom.n;
	sim->potE /= sim->atom.n;
	sim->totE  = sim->kinE + sim->potE;
	if (sim->nAtomThermRescale > 0) sim->kinETherm /= sim->nAtomThermRescale;

	// breakdown of potential energy
	if (sim->atom.n   > 0) sim->ljE   /= sim->atom.n;
	if (sim->anchor.n > 0) sim->harmE /= sim->anchor.n;
	if (sim->charge.n > 0) sim->coulE /= sim->charge.n;
	if (sim->fene.n   > 0) sim->feneE /= sim->fene.n;
	if (sim->fene.n   > 0) sim->nemE /= sim->fene.n;	//These high-jacked the fene pairs
	if (sim->bend.n   > 0) sim->bendE /= sim->bend.n;
}


/// Accumulate basic physical properties in running sums. The action parameter
/// specifies whether we want to reset the running sums, add one sample or
/// perform the average of the accumulated samples.
///
/// @param		sim a pointer to a simulation structure
/// @param		action the manipulation to perform on the running sums
/// @return		void

//================================================================================
void AccumProperties (simptr sim, int action)
//================================================================================
{
	int 	stepAvg;
	real 	a, b, c;

	// local sim variables
	stepAvg = sim->stepAvg[sim->phase];

	// perform action
	switch (action) {

		case (ACTION_RESET):
			// reset accumulators
			sim->s_kinE = sim->ss_kinE = 0.0;
			sim->s_potE = sim->ss_potE = 0.0;
			sim->s_totE = sim->ss_totE = 0.0;
			break;

		case (ACTION_COLLECT):
			// add current values to running sums
			sim->s_kinE	  += sim->kinE;
			sim->ss_kinE  += sim->kinE*sim->kinE;
			sim->s_potE	  += sim->potE;
			sim->ss_potE  += sim->potE*sim->potE;
			sim->s_totE	  += sim->totE;
			sim->ss_totE  += sim->totE*sim->totE;
			break;

		case (ACTION_AVERAGE):
			// average values over previous stepAvg segment
			sim->s_kinE /= stepAvg;
			sim->s_potE /= stepAvg;
			sim->s_totE /= stepAvg;
			a = sim->ss_kinE/stepAvg - sim->s_kinE*sim->s_kinE;
			b = sim->ss_potE/stepAvg - sim->s_potE*sim->s_potE;
			c = sim->ss_totE/stepAvg - sim->s_totE*sim->s_totE;
			if (a > 0) sim->ss_kinE = sqrt(a); else sim->ss_kinE = 0;
			if (b > 0) sim->ss_potE = sqrt(b); else sim->ss_potE = 0;
			if (c > 0) sim->ss_totE = sqrt(c); else sim->ss_totE = 0;
			break;

		default:
			break;
	}
}


/// Prints average physical properties accumulated over a certain simulation time
/// interval. The AccumProperties is called with the AVERAGE action and the
/// quantities are printed to specific simulation files. The measured properties
/// are then reset and other functions are called to handle more specialized
/// measurements, e.g., polymer properties etc. Histograms and scenes have their
/// own private counters, so they are not handled by this generic function.
///
/// @param		simvoid	void pointer cast of a pointer to a simulation structure
/// @return		void
/// @warning	This function is not automatically called by the program, and
///             has to be registered as an action in a counter to be called at all.
///	@warning	Note also that the function must be called with a void* version of
/// 			the sim pointer; it is typecasted back right away.

//================================================================================
void AverageProperties (void *simvoid)
//================================================================================
{
	simptr sim = (simptr) simvoid;

	// average accumulated measurements
	AccumProperties (sim, ACTION_AVERAGE);

	// print average system energies
	fprintf ( GetSimStream (sim->files,"energy"),
		  		"%f %.8G %.8G %.8G %.8G %.8G %.8G %.8G %.8G %.8G\n",
		  		sim->tNow, sim->s_kinE, sim->s_potE, sim->s_totE, sim->ljE,
		  		sim->harmE, sim->coulE, sim->feneE, sim->bendE, sim->nemE);

	// reset measurements accumulators
	AccumProperties (sim, ACTION_RESET);

	// compute electrostatic properties
//	ZetaPotential	(sim);

	// compute polymer properties
	PolymerGeometry (sim);
	PolymerForce	(sim);
}


/// Performs final calculations at the end of the main simulation loop. Anything
/// that should be calculated after the MD integration is complete should go in
/// in here. A good example would be the calculation of correlation functions or
/// frame averages. For now, this function simply prints a message in the log, as
/// further calculations are carried out by other programs from the generated data.
///
/// @param		sim a pointer to a simulation structure
/// @return		void

//================================================================================
void FinalCalculations (simptr sim)
//================================================================================
{
	LOG ("Final calculations completed\n");
}


/// Calculates the value of the zeta potential directly by summing the contribution
/// of all charges. The zeta potential is calculated at r=a, the internal radius
/// of the capillary (caprInTrue), and averaged from a number of different points
///
/// @param		sim a pointer to a simulation structure
/// @return		void

//================================================================================
void ZetaPotential (simptr sim)
//================================================================================
{
	int			i, n, k;
// 	real		a;
// 	real		x0, y0, z0;
	real		x0;
	real		dx, dy, dz;
	real		Nx, Nq, q, q0;
	real		*box, *boxHalf;
	particleMD	*p1;
	item1STD	*charge;
	real		r, r2, rCutCoul2=0, bjerrumkT=0;
	real		phi, zeta, dzeta;

	// local sim variables
	charge  	= sim->charge.items;
	bjerrumkT   = sim->bjerrum * sim->kT[sim->phase];
// 	a			= sim->caprInTrue;
	box			= sim->box;
	boxHalf 	= sim->boxHalf;

	// calculate coulomb cutoff or return if there are no charges
	if (sim->charge.n > 0) {
		rCutCoul2  = sim->rCutCoul * sim->rCutCoul;
	} else return;

	// set number of calculation points
	Nx = 30;
	Nq = 30;

	// loop over points on the inner capillary surface
	k=0;
	zeta=0;
	dzeta=0;
	for (x0=-boxHalf[x_]; x0<=boxHalf[x_]*(1-1.0/Nx); x0 += 2.0*boxHalf[x_]/Nx) {
		q0 = RandomReal()*2.0*pi/Nq;
		for (q=q0; q<=q0+2.0*pi*(1-0.5/Nq); q+=2.0*pi/Nq) {
// 			y0 = a*sin(q);
// 			z0 = a*cos(q);

			// sum potential contribution of all charges
			phi = 0;
			n = sim->charge.n;
			for (i=0; i<n; i++) {
				// extract particle pointer
				p1 = charge[i].p1;
				// calculate dr
				dx = (p1->rx - x0);
				dy = (p1->ry - x0);
				dz = (p1->rz - x0);
				// periodic boundary condition (don't call ApplyPBD to save time)
				if (sim->pbcond & PBC_COND_x) {	if 		(dx >= boxHalf[x_]) dx -= box[x_];
												else if	(dx < -boxHalf[x_]) dx += box[x_]; }
				if (sim->pbcond & PBC_COND_y) {	if 		(dy >= boxHalf[y_]) dy -= box[y_];
												else if	(dy < -boxHalf[y_]) dy += box[y_]; }
				if (sim->pbcond & PBC_COND_z) { if 		(dz >= boxHalf[z_]) dx -= box[z_];
												else if	(dz < -boxHalf[z_]) dz += box[z_]; }
				// calculate the distance
				r2 = dx*dx + dy*dy + dz*dz;
				// compute the potential contribution of this charge
				if (r2 <= rCutCoul2) {
					r = sqrt(r2);
					phi += bjerrumkT * p1->q * (1/r);
				}
			}
			zeta  += phi;
			dzeta += phi*phi;
			k++;
		}
	}

	// calculate average zeta potential
	zeta /= k;
	dzeta = dzeta/k - zeta*zeta;
	fprintf (GetSimStream(sim->files,"zeta"),  "%.3f %.3f %.3f\n", sim->tNow, zeta, dzeta/sqrt(k));

}

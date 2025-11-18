
//================================================================================
//
// name:   mdscene.c
// author: ftessier
// date:   2005-05-03 @ 11:04:36
//
// Scene functions to print out visualization files
//
//================================================================================


#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "mdtypes.h"
#include "mderror.h"
#include "mdmemory.h"
#include "mdfiles.h"
#include "mdutil.h"
#include "mdsrd.h"
#include "mdscene.h"


//================================================================================
void SetupScenes (simptr sim)
//================================================================================
{
	// Initializes all scenes (used to output visualization data).

	sceneptr s;

	// report
	LOG ("Setting up scenes\n");

	// include additioinal user histogram definitions
	#include "mdscenedef.h"

	// setup scene stream pointers and function pointers
	SceneSetupFiles (sim);
	SceneSetupFuncs (sim);

	// report
	s = sim->scenes;
	while (s) {
		if (s->file) {
			if (s->groupInc) LOG ("  %-25s<- included %#010X\n", s->label, s->groupInc);
			if (s->groupExc) LOG ("  %-25s-> excluded %#010X\n", s->label, s->groupExc);
		}
		s = s->next;
	}
}


//================================================================================
void SceneSetupFuncs (simptr sim)
//================================================================================
{
	// Links scenes to their handling functions. For now, we point everything to
	// the a default scene function, but eventually, we could do something more
	// complex, like with the histograms.

	sceneptr s;

	// local sim variables
	s = sim->scenes;

	// loop over all scenes
	while (s) {

		// scenefunc pointer
		s->scenefunc = &SceneFunction;

		// next scene
		s = s->next;
	}
}


//================================================================================
void SceneSetupFiles (simptr sim)
//================================================================================
{
	// Links scenes to their respective output file by saving the appropriate
	// stream into the stream member of the histogram structure. Data files MUST
	// be open and bear the same label as the scene for this setup to
	// work, otherwise the program will report an error and abort.

	sceneptr s = sim->scenes;

	while (s) {
		if (s->file) s->stream = GetSimStream (sim->files, s->label);
		s = s->next;
	}
}


//================================================================================
void ScenePrint (simptr sim, sceneptr s)
//================================================================================
{
	// Loop through the supplied scene list and check if they are scheduled
	// for data printing. If so, call the histfunc with PRINT if func pointer is
	// not null.

	int	phase;

	// local simulation variables
	phase = sim->phase;

	// check every histogram
	while (s) {
		if (s->active && s->stream) {
			if (sim->step[count_] == 0 || s->stepPrint[count_] == s->stepPrint[phase]) {
				s->stepPrint[count_] = 0;
				if(sim->warmupMD==POS_WARMUP){
					if (s->scenefunc) s->scenefunc ((void *) sim, s, ACTION_PRINT);
				}
			}
		}
		s = s->next;
	}
}


//================================================================================
void ScenePrintObjects (simptr sim, FILE *stream)
//================================================================================
{
	// Write object properties to file. These object properties are used to
	// determine the visual characteristics of the particles in the scenic viewer.

	// write object properties
	fprintf (stream, "\n\n");

	// particles
	fprintf (stream, "<object label=\"fluid\" 			id=\"%d\" size=\"1\" color=\"0,0,1\"	   detail=\"10\" visible=\"1\" /> \n", OBJ_FLUID);
	fprintf (stream, "<object label=\"wall\" 				id=\"%d\" size=\"1\" color=\"0.5,0.6,0.6\" detail=\"10\" visible=\"1\" /> \n", OBJ_WALL);
	fprintf (stream, "<object label=\"ion+\" 				id=\"%d\" size=\"1\" color=\"0.8,0,0\"	 detail=\"10\" visible=\"1\" /> \n", OBJ_ION_POS);
	fprintf (stream, "<object label=\"ion-\" 				id=\"%d\" size=\"1\" color=\"0.2,0.2,0.2\" detail=\"10\" visible=\"1\" /> \n", OBJ_ION_NEG);
	fprintf (stream, "<object label=\"monomer\" 			id=\"%d\" size=\"1\" color=\"0.0,0.6,0.2\" detail=\"10\" visible=\"1\" /> \n", OBJ_MONOMER);
	fprintf (stream, "<object label=\"monomer-first\" 		id=\"%d\" size=\"1\" color=\"0.0,0.3,0.1\" detail=\"10\" visible=\"1\" /> \n", OBJ_MONOMER_1);
	fprintf (stream, "<object label=\"monomerGraft\" 		id=\"%d\" size=\"1\" color=\"0.9,0.8,0.0\" detail=\"10\" visible=\"1\" /> \n", OBJ_MONOMER_GRAFT);
	fprintf (stream, "<object label=\"monomerGraft-first\" 	id=\"%d\" size=\"1\" color=\"0.9,0.6,0.0\" detail=\"10\" visible=\"1\" /> \n", OBJ_MONOMER_GRAFT_1);

	// periodic bounding box
	fprintf (stream, "<object label=\"PBCbox\" 				id=\"%d\" size=\"%f,%f,%f\" color=\"0.7,0.7,0.7\" render=\"wireframe\" visible=\"1\" shape=\"box\" linewidth=\"2\" />\n",OBJ_PBCBOX, sim->box[x_], sim->box[y_], sim->box[z_]);

	fprintf (stream, "\n\n");
}


//================================================================================
int SceneSpaceDomain (sceneptr s, real s1, real s2, real s3)
//================================================================================
{
	// Check if the coordinates of the point lie within the scene's space
	// domain. Note that space domain has its own dimension and coordinate
	// system. Returns 0 if the point is outside the space domain, 1 otherwise.

	// points outside space domain (reject)
	if ( (s->sd >= 1) && (s1 < s->s1min || s1 > s->s1max) ) return 0;
	if ( (s->sd >= 2) && (s2 < s->s2min || s2 > s->s2max) ) return 0;
	if ( (s->sd >= 3) && (s3 < s->s3min || s3 > s->s3max) ) return 0;

	// points inside space domain (accept)
	return 1;
}


//================================================================================
sceneptr SceneNew (sceneptr *s)
//================================================================================
{
	// Adds a scene to the beginning of the linked list of histograms
	// pointed to by *s.

	sceneptr new;

	// allocate memory for the new simscene structure
	new = (sceneptr) mycalloc (1, sizeof(simhist));

	// insert the new scene at beginning of the list
	if (s) {
		new->next = *s;
		*s = new;
	}

	// return a pointer to the newly created scene
	return new;
}


//================================================================================
sceneptr SceneDuplicate (sceneptr s, char *newlabel)
//================================================================================
{
	// Duplicates the scene pointed to by s and inserts it AFTER s, giving
	// it the new label.

	sceneptr new;

	if (!s) return 0;

	// allocate memory for the new simscene structure
	new = (sceneptr) mycalloc (1, sizeof(simhist));

	// copy scene information and attribute new label
	memcpy   (new, s, sizeof(simscene));
	snprintf (new->label, STRMAX, "%s", newlabel);

	// insert new scene at beginning of the list
	new->next = s->next;
	s->next   = new;

	// return a pointer to the newly created scene
	return new;
}


//================================================================================
sceneptr SceneGet (sceneptr s, char *label)
//================================================================================
{
	// Finds the scene coresponding to the label in the linked list of scenes
	// pointed to by s. Returns a pointer to the scene, or generates a fatal
	// error if the requested scene is not found.

	int found=0;

	while (!found && s) {
		if (strcmp(label, s->label) == 0) found = 1;
		else s = s->next;
	}

	if (!found) {
		error (ESCENELIST);
	}

	return s;
}

//================================================================================
void SceneFunction (void *simvoid, sceneptr s, int action)
//================================================================================
{
	// General scene operation function. The requested action is specified
	// by the action flag. For the moment, only the PRINT action is supported.

	simptr		sim;
	int 		i, nAtom;
	particleMD	*atom, *p1;
	point		r, rs;
	real		*sc1=0, *sc2=0, *sc3=0;

	// local simulation variables
	sim	  	= (simptr) simvoid;
	atom	= sim->atom.items;
	nAtom 	= sim->atom.n;

	// execute requested action
	switch (action) {

		case ACTION_PRINT:

			// write objects before the first frame
			if (s->frame == 0) ScenePrintObjects (sim, s->stream);

			// increase frame count
			s->frame++;

			// print frame header
			fprintf (s->stream, "<frame>\n");

			// initialize local point variables
			memset (&r,  0, sizeof(point));
			memset (&rs, 0, sizeof(point));
			r.coord  = CARTESIAN;
			rs.coord = s->scoord;

			// select space domain if applicable
			if (s->sdomain >= 0) {
				sc1 = &(rs.c1);
				sc2 = &(rs.c2);
				sc3 = &(rs.c3);
				CoordinateOrder (s->sdomain, &sc1, &sc2, &sc3);
			}

			// loop over all atoms
			for (i=0; i<nAtom; i++) {

				// particle pointer
				p1 = atom+i;

				// apply include-exclude rules
				if (!(p1->group & s->groupInc)) continue;
				if (  p1->group & s->groupExc) 	continue;


				// get particle coordinates and adjust for periodic boundaries
				r.c1 = p1->rx;
				r.c2 = p1->ry;
				r.c3 = p1->rz;
				ApplyPBC (sim, &r.c1, &r.c2, &r.c3);


				// filter points outside space domain if applicable
				if (s->sdomain >= 0) {
					CoordinateTransform (&r, &rs,0);
					if (!SceneSpaceDomain (s, *sc1, *sc2, *sc3)) continue;
				}

				// print the atom
				fprintf (s->stream, "%d %.6G %.6G %.6G %d", i, r.c1, r.c2, r.c3, p1->object);
				// print link if there is one
				if (p1->next) {
					p1 = p1->next;
					// get particle coordinates and adjust for periodic boundaries
					r.c1 = p1->rx;
					r.c2 = p1->ry;
					r.c3 = p1->rz;
					ApplyPBC (sim, &r.c1, &r.c2, &r.c3);
					fprintf (s->stream, " 1 %.6G %.6G %.6G", r.c1, r.c2, r.c3);
				}
				fprintf (s->stream, "\n");
			}

			// write pbc box
			fprintf (s->stream, "%d 0 0 0 %d\n", i, OBJ_PBCBOX);

			// print frame trailer
			fprintf (s->stream, "</frame>\n\n");

			break;

		default:
			break;
	}
}


//================================================================================
void SceneFreeList (sceneptr s)
//================================================================================
{
	// Recursively frees all the memory allocated for the linked list of scenes.

	if (!s) return;
	SceneFreeList (s->next);

	// free histogram structure
	free (s);
}





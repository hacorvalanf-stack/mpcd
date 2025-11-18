
//================================================================================
//
// name:   mderror.c
// author: ftessier
// date:   2005-05-03 @ 11:04:36
//
//================================================================================
///
/// @file mderror.c
/// Error handling code. The single function in this file prints out error
/// messages based on the supplied error code. Admittedly, this program has
/// very basic error checking and reporting.
///
//================================================================================


#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <errno.h>
#include "mderror.h"

#ifdef MPI
#include <mpi.h>
#endif


/// Error function to report errors, either standard library ones or custom
/// defined errors. The actual exit call occurs at the very end, so if the
/// error is not fatal, we can return to the calling function from within the
/// case statement (before the break). The function also prints the source
/// file name and line where the error occured, to assist in debugging.
///
/// @param 		code error code (defined in the header file)
/// @param		filename name of the code file in which the error occured
/// @param		line line on which the error occured
/// @return		void
/// @warning	This program does not have very strong error checking (i.e., it's possible to
///				fool it quite easily), only basic errors are reported.

//================================================================================
void ReportError (int code, char *filename, int line)
//================================================================================
{
	FILE *stream = stderr;

	// print error location
	fprintf (stream, "Error occured in %s at line %d\n", filename, line);

	// print error message
    switch (code) {

		case (EABORT):
		break;

		case (EPARSE):
			fprintf (stream, "Error parsing command-line arguments.\n");
			break;

		case (ESETUP):
			fprintf (stream, "Error in setting up the simulation.\n");
			break;

		case (EMACRO):
			fprintf (stream, "Error translating macros while reading parameters\n");
			break;

		case (EFILEEXCL):
			fprintf (stream, "Ouput files already exist.\n");
			break;

		case (EFILE):
			perror ("Cannot acces file");
			break;

		case (EINPUT):
			perror ("Cannot open input file");
			break;

		case (EFILELIST):
			fprintf (stream, "Cannot find requested file in file list.\n");
			break;

		case (ECHKWRITE):
			fprintf (stream, "Error writing checkpoint file!\n");
			break;

		case (ECHKREAD):
			fprintf (stream, "Error reading checkpoint file!\n");
			break;

		case (EALLOC):
			perror ("Memory allocation error");
			break;

		case (ELATTICE):
			fprintf (stream, "Cannot handle this crystal lattice type in the current context.\n");
			break;

		case (EGEOMETRY):
			fprintf (stream, "Cannot handle this geometry in the current context.\n");
			break;

		case (ETYPE):
			fprintf (stream, "Cannot handle this particle type in the current context.\n");
			break;

		case (ECAPILLARY):
			fprintf (stream, "Cannot build capillary as specified. Check input file.\n");
			break;

		case (EHIST):
			fprintf (stream, "Histogram operation error.\n");
			break;

		case (EHISTDIM):
			fprintf (stream, "Unsupporter histogram dimension.\n");
			break;

		case (EHISTLIST):
			fprintf (stream, "Cannot find requested histogram in histogram list.\n");
			break;

		case (ESCENELIST):
			fprintf (stream, "Cannot find requested scene in histogram list.\n");
			break;

		case (ECELLSORT):
			fprintf (stream, "Cell sorting error; atom is out of bounds!\n");
			break;

		case (ELOOPMAX):
			fprintf (stream, "Exceeded maximum loop count!\n");
			break;

		case (EGROWTH):
			fprintf (stream, "Recursive polymer growth has failed!\n");
			break;

		case (EGRAFT):
			fprintf (stream, "Failed to setup grafted polymers (distance contraint is probably to large).\n");
			break;

		case (ECHARGE):
			fprintf (stream, "Failed to setup charges (distance contraint is probably to large).\n");
			break;

		case (ELAYOUT):
			fprintf (stream, "Cannot handle this layout in current context.\n");
			break;

		case (EDENSITY):
			fprintf (stream, "Cannot handle this kind of density in current context.\n");
			break;

		case (ENAN):
			fprintf (stream, "Got a NaN value!\n");
			break;

		default:
			fprintf (stream, "%s.\n", strerror(errno));

    }

    // exit program
    fprintf (stream, "Aborting.\n\n");
    exit (code);
}


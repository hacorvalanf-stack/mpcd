
//================================================================================
//
// name:   mdmemory.c
// author: ftessier
// date:   2005-05-03 @ 11:04:36
//
// memory allocation and deallocation routines for molecular dynamics program
//
//================================================================================


#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "mdtypes.h"
#include "mderror.h"
#include "mdhistogram.h"
// #include "mdscene.h"
#include "mdvmd.h"
#include "mdfiles.h"
#include "mdmemory.h"


/// Allocates contiguous cleared memory, with out of memory error trapping.
/// Simply redefines calloc to handle allocation errors, such as out of memory error.
/// Note that calloc is a "cleared" alloc, meaning that the allocated space is
/// filled with 0's.
///
/// @param		n number of blocks to allocate
/// @param		size size of each block (in bytes)
/// @return		void
/// @warning	Make sure that n is greater than 0, otherwise on some platforms
///				this will crash with an EINVAL error (e.g. on Westgrid machines).
/// @see		calloc man page

//================================================================================
void *mycalloc (size_t n, size_t size)
//================================================================================
{
	void *ptr;

	if (n<=0 || size<=0) return (NULL);
	ptr = (void *) calloc (n, size);
	if (!ptr) error (EALLOC);
	return (ptr);
}


/// Reallocates memory to segment of a different size. This is used to change array
/// sizes dynamically, e.g., when we increase the size of a list. Simply redefines
/// realloc to handle allocation errors.
///
/// @param		ptr pointer to the allocated space to modify
/// @param		size new size of the memory block
/// @return 	void
/// @see		realloc man page

//================================================================================
void *myrealloc (void *ptr, size_t size)
//================================================================================
{
	if (size <= 0) {
		if (ptr) free (ptr);
		return (NULL);
	}
	ptr = (void *) realloc (ptr, size);
	if (!ptr) error (EALLOC);
	return (ptr);
}


/// Allocates additional memory to store particle information. Everytime the
/// function is called it allocates space for BLOCK_LARGE additional particles.
/// Therefore the array may be partly empty, but this is better than calling
/// realloc many times.
///
///	@param		list a pointer to a list structure for atoms
/// @return		void
/// @warning	Don't pass the structure, but a pointer to it, because this
///				function modifies members of the list structure.

//================================================================================
void GrowListAtom (listAtom *list)
//================================================================================
{
	int n = list->n + BLOCK_LARGE;

	list->items = (particleMD *) myrealloc (list->items, (size_t) (n*sizeof(particleMD)));
	list->max   = n;
}


//================================================================================
void GrowListPoly (listPoly *list)
//================================================================================
{
	// Allocates additional memory for the list of pointers to polymers.

	int n = list->n + BLOCK_LARGE;

	list->items = (itemPoly *) myrealloc (list->items, (size_t) (n*sizeof(itemPoly)));
	list->max   = n;
}


//================================================================================
void GrowList1STD (list1STD *list)
//================================================================================
{
	// Allocates additional memory for the list of pointers to all anchored
	// atoms.

	int	n = list->n + BLOCK_LARGE;

	list->items = (item1STD *) myrealloc (list->items, (size_t) (n*sizeof(item1STD)));
	list->max   = n;
}


//================================================================================
void GrowList2STD (list2STD *list)
//================================================================================
{
	// Allocates additional memory to hold the neighbor pointer pairs. Each
	// element is a structure of p1, p2.

	int	n = list->n + BLOCK_LARGE;

	list->items = (item2STD *) myrealloc (list->items, (size_t) (n*sizeof(item2STD)));
	list->max   = n;
}


//================================================================================
void GrowList3STD (list3STD *list)
//================================================================================
{
	// Allocates additional memory to hold the neighbor pointer triplets. Each
	// element is a structure of p1, p2, p3.

	int	n = list->n + BLOCK_LARGE;

	list->items = (item3STD *) myrealloc (list->items, (size_t) (n*sizeof(item3STD)));
	list->max   = n;
}


//================================================================================
void GrowList4STD (list4STD *list)
//================================================================================
{
	// Allocates additional memory to hold the neighbor pointer quadruplets. Each
	// element is a structure of p1, p2, p3, p4.

	int	n = list->n + BLOCK_LARGE;

	list->items = (item4STD *) myrealloc (list->items, (size_t) (n*sizeof(item4STD)));
	list->max   = n;
}


//================================================================================
void GrowList2PBC (list2PBC *list)
//================================================================================
{
	// Allocates additional memory for pointers to periodic boundary neighbor
	// pairs (see explanation above). Each element is a structure of p1, p2,
	// plus a pbc pointer to the appropriate periodic boundary offset in
	// calculating the distance between the pair. Whether a pair is STD or PBC
	// is decided upon building the neighbor list.

	int	n = list->n + BLOCK_LARGE;

	list->items = (item2PBC *) myrealloc (list->items, (size_t) (n*sizeof(item2PBC)));
	list->max   = n;
}


//================================================================================
void GrowListStep (listStep *list)
//================================================================================
{
	// Allocates additional memory for pointers to simulation step counters. Note
	// that the counter values for each phase are NOT in the counter items. Rather
	// the list just points to the actual counters.

	int	n = list->n + BLOCK_SMALL;

	list->items = (itemStep *) myrealloc (list->items, (size_t) (n*sizeof(itemStep)));
	list->max   = n;
}


//================================================================================
int AddItemPoly (listPoly *list, particleMD *p1, int grafted)
//================================================================================
{
	int	k = list->n;

	if (k >= list->max) GrowListPoly (list);
	list->items[k].p1 		= p1;
	list->items[k].grafted  = grafted;
	list->n++;
	return k;
}


//================================================================================
int AddItem1STD (list1STD *list, particleMD *p1)
//================================================================================
{
	int	k = list->n;

	if (k >= list->max) GrowList1STD (list);
	list->items[k].p1 = p1;
	list->n++;
	return k;
}


//================================================================================
int AddItem2STD (list2STD *list, particleMD *p1, particleMD *p2)
//================================================================================
{
	int	k = list->n;

	if (k >= list->max) GrowList2STD (list);
	list->items[k].p1 = p1;
	list->items[k].p2 = p2;
	list->n++;
	return k;
}


//================================================================================
int AddItem3STD (list3STD *list, particleMD *p1, particleMD *p2, particleMD *p3)
//================================================================================
{
	int	k = list->n;

	if (k >= list->max) GrowList3STD (list);
	list->items[k].p1 = p1;
	list->items[k].p2 = p2;
	list->items[k].p3 = p3;
	list->n++;
	return k;
}


//================================================================================
int AddItem4STD (list4STD *list, particleMD *p1, particleMD *p2, particleMD *p3, particleMD *p4)
//================================================================================
{
	int	k = list->n;

	if (k >= list->max) GrowList4STD (list);
	list->items[k].p1 = p1;
	list->items[k].p2 = p2;
	list->items[k].p3 = p3;
	list->items[k].p4 = p4;
	list->n++;
	return k;
}


//================================================================================
int AddItem2PBC (list2PBC *list, particleMD *p1, particleMD *p2, real *pbc)
//================================================================================
{
	int	k = list->n;

	if (k >= list->max) GrowList2PBC (list);
	list->items[k].p1  = p1;
	list->items[k].p2  = p2;
	list->items[k].pbc = pbc;
	list->n++;
	return k;
}


//================================================================================
int AddItemStep (listStep *list, int *counter, void (*action)())
//================================================================================
{
	int	k = list->n;

	if (k >= list->max) GrowListStep (list);
	list->items[k].counter = counter;
	list->items[k].action  = action;
	list->n++;
	return k;
}


//================================================================================
void RemoveItem1STD (list1STD *list, int i)
//================================================================================
{
	if (list->n > 0) {
		if (i < list->n-1)
			memmove (list->items+i, list->items+i+1, (list->n-1-i)*sizeof(item1STD));
		list->n--;
	}
}


//================================================================================
void ResetListAtom (listAtom *list)
//================================================================================
{
	list->items = NULL;
	list->n     = 0;
	list->max   = 0;
}


//================================================================================
void ResetList1STD (list1STD *list)
//================================================================================
{
	list->items = NULL;
	list->n     = 0;
	list->max   = 0;
}


//================================================================================
void ResetList2STD (list2STD *list)
//================================================================================
{
	list->items = NULL;
	list->n     = 0;
	list->max   = 0;
}


//================================================================================
void ResetList3STD (list3STD *list)
//================================================================================
{
	list->items = NULL;
	list->n     = 0;
	list->max   = 0;
}


//================================================================================
void ResetList4STD (list4STD *list)
//================================================================================
{
	list->items = NULL;
	list->n     = 0;
	list->max   = 0;
}


//================================================================================
void ResetList2PBC (list2PBC *list)
//================================================================================
{
	list->items = NULL;
	list->n     = 0;
	list->max   = 0;
}


//================================================================================
void ResetListPoly (listPoly *list)
//================================================================================
{
	list->items = NULL;
	list->n     = 0;
	list->max   = 0;
}


//================================================================================
void AllocateCellList (simptr sim)
//================================================================================
{
	// Allocate space for the cells and cellList which allow us to update the
	// neighbor list much more efficiently (see Allen and Tildesley).

	int i, j, cx, cy, cz;
	int nAtom, ***cell, *cellList;

	// local simulation variables
	nAtom = sim->atom.n;
	cx	= sim->nCellAxis[x_];
	cy	= sim->nCellAxis[y_];
	cz	= sim->nCellAxis[z_];

	// allocate an array of pointers (X axis)
	cell = (int ***) mycalloc ((size_t) cx, sizeof(int **));

	// allocate a contiguous block for (XY plane)
	cell[0] = (int **) mycalloc ((size_t) (cx*cy), sizeof(int *));
	for (i=1; i<cx; i++) {
		cell[i] = cell[i-1]+cy;
	}

	// allocate a contiguous block for actual elements (XYZ volume)
	cell[0][0] = (int *) mycalloc ((size_t) (cx*cy*cz), sizeof(int));
	for (i=0; i<cx; i++) {
		for (j=0; j<cy; j++) {
			if (i==0 && j==0) continue;
			cell[i][j] = cell[i][j-1]+cz;
		}
	}

	// allocate "linked list" of cell atoms
	cellList = (int *) mycalloc (nAtom, sizeof(int));

	// save pointers in sim structure
	sim->cell	  = cell;
	sim->cellList = cellList;
}


//================================================================================
void AllocatePBCList (simptr sim)
//================================================================================
{
	// Allocate memory for the permanent list of PBC offsets along the three
	// axes of the simulation box. Indices run from 0 to 26, i.e. there are 3^3
	// possiblities of combining x, y, and z. Index 13 is the PBC_INDEX_NONE
	// case, which corresponds to an offset of (0,0,0).

	int 	i, j, k, n;
	real	**pbc;

	// allocate an array of real pointers
	pbc = (real **) mycalloc ((size_t) 27, sizeof(real *));

	// allocate a contiguous block for the whole array
	pbc[0] = (real *) mycalloc ((size_t) (3*27), sizeof(real));

	// adjust pointers
	for (i=1; i<27; i++) {
		pbc[i] = pbc[i-1]+3;
	}

	// set initial values
	n=0;
	for (i=-1; i<=1; i++) {
		for (j=-1; j<=1; j++) {
			for (k=-1; k<=1; k++) {
				pbc[n][x_] = i*sim->box[x_];
				pbc[n][y_] = j*sim->box[y_];
				pbc[n][z_] = k*sim->box[z_];
				n++;
			}
		}
	}

	// save pointers in sim structure
	sim->pbc = pbc;
}


//================================================================================
void FreeFileList (simfile *fptr)
//================================================================================
{
	// Iteratively frees the memory allocated for the linked list of all simfile
	// structures that contain information about all simulation output files.

	if (!fptr) return;
	FreeFileList (fptr->next);
	free (fptr);
}


//================================================================================
void FreeMemory (simptr sim)
//================================================================================
{
	// Frees all the memory allocated by the simulation program. This is not
	// strictly necessary, since the system will free memory when the program
	// exits, but it is good practice!

	// free lists
	free (sim->atom.items);
	free (sim->nebrSTD.items);
	free (sim->nebrPBC.items);
	free (sim->anchor.items);
	free (sim->charge.items);
	free (sim->fene.items);
	free (sim->bend.items);
	free (sim->polymer.items);
	free (sim->stepCounter.items);
	free (sim->wall.items);
	free (sim->fluid.items);
	free (sim->surface.items);

	// free arrays
	free (sim->pbc[0]);
	free (sim->pbc);
	free (sim->cell[0][0]);
	free (sim->cell[0]);
	free (sim->cell);
	free (sim->cellList);

	// free parameter list
	free (sim->param);

	// free histogram list
	HistogramFreeList (sim->histograms);

	// free scene list
// 	SceneFreeList (sim->scenes);
	VMDFreeList (sim->scenes);

	// free simfile list
	FreeFileList (sim->files);

	// free sim pointer
	free (sim);
}

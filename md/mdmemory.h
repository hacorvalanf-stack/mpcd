
//================================================================================
//
// name:   mdmemory.h
// author: ftessier
// date:   2005-05-03 @ 11:04:36
//
// Header file for mdmem.c
//
//================================================================================


#ifndef MDMEMORY_H
#define MDMEMORY_H


#define BLOCK_LARGE	1024
#define	BLOCK_SMALL	16

//================================================================================
// Prototypes
//================================================================================
void		*mycalloc				(size_t n, size_t size);
void		*myrealloc				(void *ptr, size_t size);
void		GrowListAtom			(listAtom *list);
void		GrowListPoly			(listPoly *list);
void		GrowList1STD	 		(list1STD *list);
void		GrowList2STD	 		(list2STD *list);
void		GrowList3STD	 		(list3STD *list);
void		GrowList4STD	 		(list4STD *list);
void		GrowList2PBC 			(list2PBC *list);
void 		GrowListStep			(listStep *list);
int 		AddItemPoly				(listPoly *list, particleMD *p1, int grafted);
int 		AddItem1STD 			(list1STD *list, particleMD *p1);
int 		AddItem2STD 			(list2STD *list, particleMD *p1, particleMD *p2);
int 		AddItem3STD 			(list3STD *list, particleMD *p1, particleMD *p2, particleMD *p3);
int         AddItem4STD             (list4STD *list, particleMD *p1, particleMD *p2, particleMD *p3, particleMD *p4);
int 		AddItem2PBC 			(list2PBC *list, particleMD *p1, particleMD *p2, real *pbc);
int 		AddItemStep 			(listStep *list, int *counter, void (*action)());
void 		RemoveItem1STD 			(list1STD *list, int i);
void		ResetListAtom 			(listAtom *list);
void		ResetListPoly 			(listPoly *list);
void		ResetList1STD 			(list1STD *list);
void		ResetList2STD 			(list2STD *list);
void		ResetList3STD 			(list3STD *list);
void		ResetList4STD 			(list4STD *list);
void		ResetList2PBC 			(list2PBC *list);
void		AllocateCellList		(simptr sim);
void		AllocatePBCList 		(simptr sim);
void 		AllocateHistogramBins 	(histptr h);
void 		FreeHistogramList 		(histptr h);
void 		FreeFileList 			(simfile *fptr);
void		FreeMemory				(simptr sim);


#endif


//================================================================================
//
// name:   mderror.h
// author: ftessier
// date:   2005-05-03 @ 11:04:36
//
//================================================================================


#ifndef MDERROR_H
#define MDERROR_H


// Error codes

#define EABORT		1000
#define EPARSE		1001
#define EMACRO		1002
#define	EFILEEXCL	1004
#define	EFILE		1005
#define	EINPUT		1006
#define	EFILELIST	1007
#define	ECHKWRITE	1008
#define	ECHKREAD	1009
#define	ESETUP		1010
#define	ELAYOUT		1011
#define	EDENSITY	1012

#define	EALLOC		2000

#define ENAN		3000

#define ELATTICE	6000
#define EGEOMETRY	6001
#define ETYPE		6002
#define ELOOPMAX	6003
#define EGROWTH		6004
#define	EGRAFT		6005
#define	ECHARGE		6006

#define ECAPILLARY	7000

#define EHIST		8000
#define EHISTDIM	8001
#define EHISTLIST	8002
#define ESCENELIST	8003

#define ECELLSORT	9000

#define ERRNO		-1


//================================================================================
// Prototypes
//================================================================================
void ReportError (int code, char *filename, int line);


#endif

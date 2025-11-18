
//================================================================================
//
// name:   mdfiles.h
// author: ftessier
// date:   2005-05-03 @ 11:04:36
//
// Header file for mdfiles.c
//
//================================================================================


#ifndef MDFILES_H
#define MDFILES_H


//================================================================================
// Macro definitions
//================================================================================

#define LINEBLOCK 256
#define	LINESIZE  256
#define KEYSIZE   256

#define BINARY_DESC_SIM	 			1
#define BINARY_DESC_ATOM			2
#define BINARY_DESC_ATOM_START		3
#define BINARY_DESC_POLYMER			4
#define BINARY_DESC_FILE			5
#define BINARY_DESC_HISTOGRAM		6
#define BINARY_DESC_HISTBIN			7
#define BINARY_DESC_HISTCOUNT		8
#define BINARY_DESC_SCENE			9


//================================================================================
// Prototypes
//================================================================================

void	SetupDataFiles		(simptr sim);
void	CloseDataFiles		(simptr sim);
void	ReopenDataFiles 	(simptr sim);
void	CheckpointWrite		(void *simvoid);
void	CheckpointRead		(simptr sim);
void	WriteBinary 		(int desc, int count, int size, void *ptr, FILE *f);
void	WriteHeader			(simfile *output, simptr sim);
void	Report				(simptr sim, FILE *stream, char *msg);

void	WriteObjects 		(simptr sim, FILE *label);
void	WriteScene 			(void *simvoid);

void	ReadParameters		(char *inputFile, paramptr param, int nParam, char *label, char *stop);
void	WriteParameters 	(FILE *output, paramptr param, int nParam, char *label);

char	**LinesToArray 		(FILE *input, const char *stop);
void	FilterString		(char **stringptr, const char *filter);
void	TranslateMacros		(char **stringptr);
void	ReplaceMacro		(char **stringptr, char *macro, char *value);

void	ParseOptions		(int argc, char *argv[], simptr sim, simoptions *options);
void	SetSimLabel 		(char *label, int pid);
void	SetWorkingDir		(char *dir, char *simdir);

simfile *FileNew 			(simfile *fptr);
simfile *GetSimFile 		(simfile *fptr, char *label);
FILE	*GetSimStream 	    (simfile *fptr, char *label);
char	*GetSimFileName 	(simfile *fptr, char *label);


#endif

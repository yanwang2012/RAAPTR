/*! \file 

\brief Utilities used by MaxPhase codes.

\author Soumya D. Mohanty
*/
#include "maxphase.h"
#include <gsl/gsl_vector.h>
#include <stdio.h>
#include <sys/types.h>
#include <dirent.h>
#include <fnmatch.h>
#include <string.h>

/*! 
Allocate a fitness function parameter structure for a given dimensionality of the fitness function.
*/
struct fitFuncParams * ffparam_alloc(size_t nDim /*!> Number of dimensions */){
	struct fitFuncParams *ffp = (struct fitFuncParams *)malloc(sizeof(struct fitFuncParams));
	ffp->nDim = nDim;
	ffp->rmin = gsl_vector_calloc(nDim);
	ffp->rangeVec = gsl_vector_calloc(nDim);
	ffp->realCoord = gsl_vector_calloc(nDim);
	return ffp;
}

/*!
Deallocate a fitness function parameter structure.
*/
void ffparam_free(struct fitFuncParams *ffp){
	if (ffp == NULL){
		printf("Invalid FitFuncParams supplied\n");
		return;
	}
	gsl_vector_free(ffp->rmin);
	gsl_vector_free(ffp->rangeVec);
	gsl_vector_free(ffp->realCoord);
	free(ffp);
}

/*!					   			   
Takes standardized coordinates and
returns real coordinates using the supplied range and minimum limits. 
The range and limits can be different for the different coordinates. 

Notes:
 -  Derived from sr2vector.c
 -  Shifted to using gsl_vector
*/
void s2rvector(const gsl_vector *xVec, /*!< Standardized coordinates of the point.*/
 			   const gsl_vector *rmin,       /*!< Minimum value of each coordinate.*/
			   const gsl_vector *rangeVec,   /*!< Range of each coordinate. */
               gsl_vector *realCoord  /*!< Returns real coordinate values.*/)
{
			   
	size_t lpc;
	double x, rv, rmi, rc;
	size_t ncols = xVec->size;
	
    for(lpc=0; lpc < ncols; lpc++){
		x = gsl_vector_get(xVec,lpc);
		rv = gsl_vector_get(rangeVec,lpc);
		if(rv<0){
			printf("Invalid range\n");
			abort();
		}
		rmi = gsl_vector_get(rmin,lpc);
		rc = x*rv+rmi;
    	gsl_vector_set(realCoord,lpc,rc);
    }
}

/*! 
Returns 0 or 1 corresponding to invalid/valid
coordinates. The invalid flag is set if any of the coordinates
 fall outside the closed interval [0,1].

Notes:
 - Derived from chkstdsrchrng.c
*/
size_t chkstdsrchrng(const gsl_vector *xVec/*!< Standardized coordinates of the point.*/)
{
	double x;
    unsigned int validPt,lpc;
	size_t ncols = xVec->size;;
	
	/*printf("from chkstdsrchrng\n");*/
    validPt=1;
	for (lpc=0; lpc <ncols; lpc++){
		x = gsl_vector_get(xVec,lpc);
		if (x<0 || x>1){
			validPt = 0;
			break;
		}						
	}
	
	return validPt;
	
}

/*! Limit each compononent of a vector to a specified range */
void limitVecComponent(gsl_vector *xVec, double min, double max){
	size_t lpc;
	size_t nDim = xVec->size;
	double x;
	for (lpc = 0; lpc < nDim; lpc++){
		x = gsl_vector_get(xVec,lpc);
		if (x < min)
			gsl_vector_set(xVec,lpc,min);
		else if (x > max)
			gsl_vector_set(xVec,lpc,max);
	}
}

/*! List all files in a directory with a given extension. 
   - First parameter: Extension (without a '.'). Example: "mat" and not ".mat".
   - Second parameter: Directory name
   - Third parameter:  number of files found.
   - Fourth parameter: length of longest file name string.
   - Output: List of filenames. Dynamically allocated memory: must be freed by 
           calling function.
 */
char ** listfileswext (const char *ext, const char *dirName, size_t *nFiles, size_t *maxFileNameLen)
{
	/*printf("------- %s -------\n",ext);
	printf("Extension string length %zu\n",strlen(ext));*/
  char *pattern = (char *)malloc((3+strlen(ext))*sizeof(char));
  DIR *dp;
  struct dirent *ep;
  char **fileList;
  size_t lpc;
 
  pattern[0] = '*';
  pattern[1] = '.';
  pattern[2] = '\0';
  strcat(pattern,ext);
  
  /*printf("Pattern %s\n",pattern);*/
  /* Step 1: Count the number of files with the required 
     extension.
  */
  size_t countValidFiles = 0;
  dp = opendir(dirName);
  if (dp != NULL){
      while ((ep = readdir(dp))){
	        /*printf("Checking file %s\n",ep->d_name);*/
			if(!fnmatch(pattern, ep->d_name, (FNM_FILE_NAME|FNM_PERIOD))){
				/* Match. Increment counter */
				/*printf("found match with pattern %s\n",pattern);*/
				countValidFiles++;
			}
	  }
      (void) closedir (dp);
	  /*
	    Apparently, there is no need to free ep as it is declared to be 'static' in the readdir function
	  */
  }
  else{
    printf ("Couldn't open the directory %s\n",dirName);
	free(pattern);
	return NULL;
  }
  *nFiles = countValidFiles;
  /* Step 2: Create storage for list of filenames */
  fileList = (char **)malloc(countValidFiles*sizeof(char *));
  /* Can't find a better way than to repeat the whole loop again */
  countValidFiles = 0;
  dp = opendir(dirName);
  if (dp != NULL){
      while ((ep = readdir(dp))){
			if(!fnmatch(pattern, ep->d_name, (FNM_FILE_NAME|FNM_PERIOD))){
				fileList[countValidFiles] = (char *)malloc((strlen(ep->d_name)+1)*sizeof(char));
				strcpy(fileList[countValidFiles],ep->d_name);
				/* Match. Increment counter */
				countValidFiles++;
			}
	  }
      (void) closedir (dp);
  }
  else{
    printf ("Couldn't open the directory %s\n",dirName);
	return NULL;
  }
  
	/*Find longest filename */
	size_t fileNameLen;
	*maxFileNameLen = 0;
	size_t lpc1;
	for (lpc1 = 0; lpc1 < *nFiles; lpc1++){
		fileNameLen = strlen(fileList[lpc1]);
		if ( fileNameLen > *maxFileNameLen)
			*maxFileNameLen = fileNameLen;
	}
  
  /* Wrap up */
  free(pattern);
  return fileList;
}
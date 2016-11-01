/* Test suite for chkstdsrchrng.c */
#include <stdio.h>
void chkstdsrchrng(double *[], unsigned int *, int, int);
int main(){
   /* Test matrix */
   int nrows = 5, ncols = 3;
   int lpr, lpc;
   double testMat[5][3]={{.1,.2,.5},
                        {1.0,.1,.2},
					    {.1,.2,0.0},
					    {2.0,.1,.2},
					    {.1,.2,2.0}};
   double *xVec[5];
   unsigned int validPts[5];
   /* Display test matrix and result*/   
   for (lpr=0; lpr<nrows; lpr++){
	   xVec[lpr] = &testMat[lpr][0];
   		for (lpc=0; lpc<ncols; lpc++){
			printf("%f ",*(xVec[lpr]+lpc));
		}
		printf("\n");
   }   
   chkstdsrchrng(xVec,validPts,nrows,ncols);
   printf("-------------\n");
   for (lpr=0; lpr<nrows; lpr++){
   		for (lpc=0; lpc<ncols; lpc++){
			printf("%f ",*(xVec[lpr]+lpc));
		}
		printf(": %u",validPts[lpr]);
		printf("\n");
   }     
}
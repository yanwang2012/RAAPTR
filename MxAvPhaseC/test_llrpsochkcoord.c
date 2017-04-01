#include "LLR_PSO.h"
#include "maxphase.h"
#include <gsl/gsl_vector.h>
#include <gsl/gsl_math.h>
/* Test
size_t  llrpsochkcoord(const size_t wrapAngles, const gsl_vector *xVec, 
                       const gsl_vector *rmin, const gsl_vector *rangeVec,
                       gsl_vector *realCoord)

gcc -c test_llrpsochkcoord.c -I/usr/local/include -I/Applications/MATLAB_R2016a.app/extern/include

gcc test_llrpsochkcoord.o LLR_PSO.o maxphaseutils.o -L/usr/local/lib 
    -L/Applications/MATLAB_R2016a.app/bin/maci64 -lmat -lmx -Xlinker -rpath -Xlinker 
    /Applications/MATLAB_R2016a.app/bin/maci64 -lgsl -lgslcblas -lm -o test_llrpsochkcoord.out
*/

int main(){

    /* Test matrix */
    int nrows = 5, ncols = 7;
    int lpr, lpc;
	/* 
	   coord #0: alpha; default range is [0, 2*pi]
	         #1: delta; default range is [pi/2,-pi/2]
		Coord #2: Frequency
	         #3: phi0; default range is [0,pi]
		Coord #4: Amplitude
	         #5: iota; default range is [0,pi] 
	         #6: thetaN; default range is [0,pi] 
	*/
	//Minimum of coordinate range for Angular parameters (in degrees)
	double rminAngleDeg[5] = {50.0000 , -20.0000,  0,  0, 0};
	//Maximum of coordinate range for Angular parameters (in degrees)
	double rmaxAngleDeg[5] = {117.0000, 20.0000,  180.0000,  180.0000, 180.0000};
	
	printf("Minimum and Maximum values of angular coordinates:\n");
	printf("alpha: %f, %f\n",rminAngleDeg[0], rmaxAngleDeg[0]);
	printf("delta: %f, %f\n",rminAngleDeg[1], rmaxAngleDeg[1]);
	printf("phi0: %f, %f\n",rminAngleDeg[2], rmaxAngleDeg[2]);
	printf("iota: %f, %f\n",rminAngleDeg[3], rmaxAngleDeg[3]);
	printf("thetaN: %f, %f\n",rminAngleDeg[4], rmaxAngleDeg[4]);
	
	//Minimum for non-angular coordinates
	double rminNonAngles[2] = {0.0100, -5.0000};
	//Maximum for non-angular coordinates
	double rmaxNonAngles[2] = {0.0200, -1.0000};
	
	//Test coordinates for only Angular coordinates (in degrees)
	double realCoordAngles[5][5] = {{117.000000, 0.000000, 91.000000, 32.000000, 139.000000},
									{831.000000, 0.000000, 91.000000, 32.000000, 139.000000},
									{51.000000, 0.000000, 91.000000, 32.000000, 139.000000},
									{143.000000, -12.000000, 190.000000, 32.000000, 139.000000},
									{831.000000, -375.0, 91.000000, 32.000000, 139.000000}};
	//Test coordinates for non-angular coordinates
	double realCoordNonAngles[2] = {12e-3, -3.5};
    /*double stdCoord[5][7]= { { 0.397222,1.166667,0.200000,0.505556,0.375000,0.177778,0.772222},
							{ 2.397222,0.322222,0.200000,0.505556,0.375000,0.177778,0.772222},
							{ 2.397222,1.166667,0.200000,0.505556,0.375000,0.177778,0.772222},
							{ 0.397222,0.322222,0.200000,1.055556,0.375000,0.177778,0.772222},
							{ 0.397222,0.322222,0.200000,0.505556,0.375000,0.177778,0.772222} };
	*/

  	/*double rmin[7] = { 0.000000000000000,-1.570796326794897,0.010000000000000,0.000000000000000,-5.000000000000000,0.000000000000000,0.000000000000000 }
;
  	double rmax[7] = { 6.283185307179586,1.570796326794897,0.020000000000000,3.141592653589793,-1.000000000000000,3.141592653589793,3.141592653589793 }
;*/
	//Degrees to radians conversion factor
	double deg2rad = M_PI/180;
	double rad2deg = 180/M_PI;
	//Merge minima of coordinate ranges; convert angular coordinates to radians
	double rmin[7];
	rmin[0] = rminAngleDeg[0]*deg2rad;
	rmin[1] = rminAngleDeg[1]*deg2rad;
	rmin[2] = rminNonAngles[0];
	rmin[3] = rminAngleDeg[2]*deg2rad;
	rmin[4] = rminNonAngles[1];
	rmin[5] = rminAngleDeg[3]*deg2rad;
	rmin[6] = rminAngleDeg[4]*deg2rad;
	
	//Merge maxima of coordinate ranges; convert angular coordinates to radians
	double rmax[7];
	rmax[0] = rmaxAngleDeg[0]*deg2rad;
	rmax[1] = rmaxAngleDeg[1]*deg2rad;
	rmax[2] = rmaxNonAngles[0];
	rmax[3] = rmaxAngleDeg[2]*deg2rad;
	rmax[4] = rmaxNonAngles[1];
	rmax[5] = rmaxAngleDeg[3]*deg2rad;
	rmax[6] = rmaxAngleDeg[4]*deg2rad;
	
	//Merge coordinates; convert angular coordinates to radians
	double realCoord[5][7];
	for(lpr = 0; lpr < nrows; lpr++){
		realCoord[lpr][0] = realCoordAngles[lpr][0]*deg2rad;
		realCoord[lpr][1] = realCoordAngles[lpr][1]*deg2rad;
		realCoord[lpr][2] = realCoordNonAngles[0];
		realCoord[lpr][3] = realCoordAngles[lpr][2]*deg2rad;
		realCoord[lpr][4] = realCoordNonAngles[1];
		realCoord[lpr][5] = realCoordAngles[lpr][3]*deg2rad;
		realCoord[lpr][6] = realCoordAngles[lpr][4]*deg2rad;
	}
	
	//Standardize coordinates
	double stdCoord[5][7];
	for(lpr = 0; lpr < nrows; lpr++){
		for(lpc = 0; lpc < ncols; lpc++){
			stdCoord[lpr][lpc] = (realCoord[lpr][lpc] - rmin[lpc])/(rmax[lpc]-rmin[lpc]);
		}
	}
	
	//Standardize;
	/*
	double rmin[7] = { 0.872664625997165,-0.349065850398866,0.010000000000000,0.000000000000000,-5.000000000000000,0.000000000000000,0.000000000000000 };
	double rmax[7] = { 2.042035224833366,0.349065850398866,0.020000000000000,3.141592653589793,-1.000000000000000,3.141592653589793,3.141592653589793 }
	;
	double stdCoord[5][7]= { { 1.000000,0.500000,0.200000,0.505556,0.375000,0.177778,0.772222},
							{ 11.656716,0.500000,0.200000,0.505556,0.375000,0.177778,0.772222},
							{ 0.014925,0.500000,0.200000,0.505556,0.375000,0.177778,0.772222},
							{ 1.388060,0.200000,0.200000,1.055556,0.375000,0.177778,0.772222},
							{ 1.388060,0.200000,0.200000,0.505556,0.375000,0.177778,0.772222} };
	*/

	
	gsl_vector *xVec = gsl_vector_alloc(ncols);
	gsl_vector *rcOut = gsl_vector_alloc(ncols);
	gsl_vector *rminVec = gsl_vector_alloc(ncols);
	gsl_vector *rngVec = gsl_vector_alloc(ncols);
	for(lpc = 0; lpc < ncols; lpc++){
		gsl_vector_set(rminVec, lpc, rmin[lpc]);
		gsl_vector_set(rngVec, lpc, rmax[lpc] - rmin[lpc]);
	}
	
	size_t validPt;
	for (lpr=0; lpr < nrows; lpr++){
		
		printf("Row # %zu\n",lpr+1);
		
		printf("Unit vector components from input alpha and delta:\n");
		double alpha = realCoord[lpr][0];
		double delta = realCoord[lpr][1];
		
		printf("%f, %f, %f\n", sin(M_PI_2 - delta)*cos(alpha),
		                       sin(M_PI_2 - delta)*sin(alpha),
							   cos(M_PI_2 - delta));
		
		for (lpc = 0; lpc < ncols; lpc++){
			gsl_vector_set(xVec, lpc, stdCoord[lpr][lpc]);			
		}
		
		printf("Input std coords: \n");
		for (lpc = 0; lpc < ncols; lpc++){
			printf("%f, ", gsl_vector_get(xVec, lpc));			
		}
		printf("\n");
		
		printf("Input real coords: \n");
		for (lpc = 0; lpc < ncols; lpc++){
			printf("%f, ", realCoord[lpr][lpc]);			
		}
		printf("\n");	
		
		printf("Input real angles (degrees):\n");
		printf("%f, %f, %f, %f, %f\n", realCoord[lpr][0]*rad2deg,
		        realCoord[lpr][1]*rad2deg, realCoord[lpr][3]*rad2deg, 
				realCoord[lpr][5]*rad2deg, realCoord[lpr][6]*rad2deg);	
				
		validPt = llrpsochkcoord(1, rminVec, rngVec, xVec, rcOut);
		
		printf("Output real angles (degrees):\n");
		printf("%f, %f, %f, %f, %f\n", gsl_vector_get(rcOut,0)*rad2deg,
		        gsl_vector_get(rcOut,1)*rad2deg, gsl_vector_get(rcOut,3)*rad2deg, 
				gsl_vector_get(rcOut,5)*rad2deg, gsl_vector_get(rcOut,6)*rad2deg);	
		
		printf("Output real coords: \n");
		for (lpc = 0; lpc < ncols; lpc++){
			printf("%f ",gsl_vector_get(rcOut, lpc));			
		}
		printf("\n");
		
		printf("Output std coords:\n");
		for (lpc = 0; lpc < ncols; lpc++){
			printf("%f, ", gsl_vector_get(xVec, lpc));			
		}
		printf("\n");		
		
		
		alpha = gsl_vector_get(rcOut,0);
		delta = gsl_vector_get(rcOut,1);
		printf("Unit vector components from output alpha and delta:\n");
		printf("%f, %f, %f\n", sin(M_PI_2 - delta)*cos(alpha),
		                       sin(M_PI_2 - delta)*sin(alpha),
							   cos(M_PI_2 - delta));
		printf("ValidPt: %zu\n --------\n", validPt);
	}
	
	//Wrap up
	gsl_vector_free(xVec);
	gsl_vector_free(rcOut);
	gsl_vector_free(rminVec);
	gsl_vector_free(rngVec);
	
}
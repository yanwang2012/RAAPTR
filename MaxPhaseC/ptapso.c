#include "ptapso.h"
#include "maxphase.h"
#include <stdio.h>
#include <gsl/gsl_multimin.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_math.h>
/*! \file
\brief Particle Swarm Optimization (PSO) and support functions.

The function \ref ptapso is the main one. The rest are support functions 
that do jobs such as initialization, memory allocation/deallocation, 
and handling output.

\author Soumya D. Mohanty
*/
/*! 
This function accepts a pointer to a fitness function and searches 
for its global optimum using Particle Swarm Optimization. The fitness
function must accept coordinates in the range [0,1]. See \ref fitfunc_example for
an example of the interface required for a fitness function.

Notes on the PSO implementation used:
   - Follows the prescription of Bratton, Kennedy, 2007.
   - Local best (lbest) PSO with three nearest neighbors in a ring topology.
   - Linearly deacreasing inertia weight.
   - Velocity clamping
*/
void ptapso(size_t nDim, /*!< Number of search dimensions */
            double (*fitfunc)(gsl_vector *, void *), /*!< Pointer to Fitness function */
			void *ffParams, /*!< Fitness function parameter structure */
            struct psoParamStruct *psoParams, /*!< PSO parameter structure */
			struct returnData *psoResults /*!< Output structure */){
				
	/* 
	  Random numbers can be generated on the fly or read from a file. If a
	  valid file is specified, generation of random numbers is overriden.
	*/
	gsl_rng *rngGen = psoParams->rngGen;
	
	/* Initialize local minimizer of gbest */
	gsl_multimin_function func2minimz;
	func2minimz.n = nDim; /*dimensionality of function to minimize */
	func2minimz.f = dummyfitfunc;//fitfunc; /* Name of function to minimize */
	struct dummyFitFuncParam dffp;
	dffp.trufuncPr = fitfunc;
	dffp.trufuncParam = ffParams;
	func2minimz.params = &dffp;//ffParams; /* Parameters needed by this function */
	/* Local Minimization method: Nelder Mead */
	gsl_multimin_fminimizer *minimzrState = gsl_multimin_fminimizer_alloc(gsl_multimin_fminimizer_nmsimplex2, nDim);
	/* Initial step vector of local minimization method */
	gsl_vector *locMinStp = gsl_vector_alloc(nDim);
	gsl_vector_set_all(locMinStp,psoParams->locMinStpSz);
	size_t lpLocMin;/* Local minimization iteration counter */
	int status;
	
	/* PSO loop counters */
	size_t lpParticles, lpPsoIter;
	/* Number of particles */
	const size_t popsize = psoParams->popsize;
	/* Number of iterations */
	const size_t maxSteps = psoParams->maxSteps;
	/* Information about a particles is stored in a struct array.
	*/
    struct particleInfo pop[popsize];
	/* initialize particles */
	for (lpParticles = 0; lpParticles < popsize; lpParticles++){
	      initPsoParticles(&pop[lpParticles], nDim, rngGen);
    }	
	/* Variables needed to find and track gbest */
	double gbestFitVal = GSL_POSINF;
	gsl_vector *gbestCoord = gsl_vector_alloc(nDim);
	gsl_vector *partSnrCurrCol = gsl_vector_alloc(popsize);
	size_t bestfitParticle;
	double currBestFitVal;
	/* Variables needed to keep track of number of fitness function evaluations */
	unsigned char computeOK;
	size_t funcCount;
	/* Variables needed in PSO dynamical equation update */
	size_t lpNbrs; /* Loop counter over nearest neighbors */
	size_t nNbrs = 3;
	size_t ringNbrs[nNbrs];
	double nbrFitVal; /*Fitness of a neighbor */
	size_t lbestPart; /* local best particle */
	double lbestFit; /* Fitness of local best particle */
	/* Variables needed for PSO dynamical equations */
	size_t lpCoord;
	gsl_vector *accVecPbest = gsl_vector_alloc(nDim);
	gsl_vector *accVecLbest = gsl_vector_alloc(nDim);
	gsl_vector *chi1Vec = gsl_vector_alloc(nDim);
	gsl_vector *chi2Vec = gsl_vector_alloc(nDim);
	
	/* 
	   Start PSO iterations from the second iteration since the first is used
	   above for initialization.
	*/
	for (lpPsoIter = 1; lpPsoIter <= maxSteps-1; lpPsoIter++){
		if (psoParams->debugDumpFile != NULL){
			fprintf(psoParams->debugDumpFile,"Loop %zu \n",lpPsoIter);
			particleInfoDump(psoParams->debugDumpFile,pop,popsize);
		}		
        /* Calculate fitness values */
		for (lpParticles = 0; lpParticles < popsize; lpParticles++){
			/* Evaluate fitness */
			pop[lpParticles].partSnrCurr = fitfunc(pop[lpParticles].partCoord,ffParams);
			/* Separately store all fitness values -- needed to find best particle */
			gsl_vector_set(partSnrCurrCol,lpParticles,pop[lpParticles].partSnrCurr);
			/* Check if fitness function was actually evaluated or not */
	        computeOK = ((struct fitFuncParams *)ffParams)->fitEvalFlag;
	        funcCount = 0;
	        if (computeOK){
			    /* Increment fitness function evaluation count */
	            funcCount = 1;
			}
	        pop[lpParticles].partFitEvals += funcCount;
			/* Update pbest fitness and coordinates if needed */
	        if (pop[lpParticles].partSnrPbest > pop[lpParticles].partSnrCurr){
	            pop[lpParticles].partSnrPbest = pop[lpParticles].partSnrCurr;
	            gsl_vector_memcpy(pop[lpParticles].partPbest,pop[lpParticles].partCoord);
	        }
	    }
		
		/* Find the best particle in the current iteration */
		bestfitParticle = gsl_vector_min_index(partSnrCurrCol);
	    currBestFitVal = pop[bestfitParticle].partSnrCurr; 
	    if (gbestFitVal > currBestFitVal){
		/* 
		   Do local minimization iterations since gbest has changed.
		*/
		   	gsl_multimin_fminimizer_set(minimzrState,&func2minimz,
			                            pop[bestfitParticle].partCoord,
										locMinStp);
			funcCount = 0;
			
			for (lpLocMin = 0; lpLocMin < psoParams->locMinIter; lpLocMin++){
				status = gsl_multimin_fminimizer_iterate(minimzrState);
				  //A non-zero value of status indicates some type of failure  
				if (status)
					break;
				 /*Note that the function evaluation count is only an approximate
				   one for the nmsimplex2 algorithm as GSL routines 
				   do not return this information.*/
				funcCount += nDim+1;
				pop[bestfitParticle].partSnrCurr = gsl_multimin_fminimizer_minimum(minimzrState);
                gsl_vector_memcpy(pop[bestfitParticle].partCoord, minimzrState->x);
			}
			
	        pop[bestfitParticle].partFitEvals += funcCount;
			/* Update particle pbest */
			pop[bestfitParticle].partSnrPbest = pop[bestfitParticle].partSnrCurr;
			gsl_vector_memcpy(pop[bestfitParticle].partPbest,pop[bestfitParticle].partCoord);
			/* Update gbest */
			gbestFitVal = pop[bestfitParticle].partSnrCurr;
			gsl_vector_memcpy(gbestCoord,pop[bestfitParticle].partCoord);
			/* Update list of fitness Values */
			gsl_vector_set(partSnrCurrCol,bestfitParticle,pop[bestfitParticle].partSnrCurr);
		}
		
		/* Get lbest */
	    for (lpParticles = 0; lpParticles < popsize; lpParticles++){
	           if (lpParticles == 0){
	                   ringNbrs[0]=popsize-1; ringNbrs[1]=lpParticles; ringNbrs[2]=lpParticles+1;
			   		}
			   else if (lpParticles == popsize -1){
			       ringNbrs[0]=lpParticles-1; ringNbrs[1]=lpParticles; ringNbrs[2]=0;
			   		}
			   else{
				   ringNbrs[0]=lpParticles-1; ringNbrs[1]=lpParticles; ringNbrs[2]=lpParticles+1;
			   }					  
			   /* Get best particle in neighborhood */
			   lbestPart = ringNbrs[0];
			   lbestFit = gsl_vector_get(partSnrCurrCol,lbestPart);
			   for (lpNbrs = 1; lpNbrs < nNbrs; lpNbrs++){
				   nbrFitVal = gsl_vector_get(partSnrCurrCol,ringNbrs[lpNbrs]);
				   if (nbrFitVal < lbestFit){
					   lbestPart = ringNbrs[lpNbrs];
					   lbestFit = nbrFitVal;
				   }
			   }
	           if (lbestFit < pop[lpParticles].partSnrLbest){
	               pop[lpParticles].partSnrLbest = lbestFit;
	               gsl_vector_memcpy(pop[lpParticles].partLocalBest,
				                     pop[lbestPart].partCoord);
	           }
		}
        

	    for (lpParticles = 0; lpParticles < popsize; lpParticles++){
			/* Update inertia Weight */
		    pop[lpParticles].partInertia = psoParams->dcLaw_a-(psoParams->dcLaw_b/psoParams->dcLaw_c)*lpPsoIter;
			if (pop[lpParticles].partInertia < psoParams->dcLaw_d)
				pop[lpParticles].partInertia = psoParams->dcLaw_d;
		    /* Acceleration vector to pbest */
			gsl_vector_memcpy(accVecPbest,pop[lpParticles].partPbest);
			gsl_vector_sub(accVecPbest,pop[lpParticles].partCoord);
			/* Acceleration vector to lbest */
			gsl_vector_memcpy(accVecLbest,pop[lpParticles].partLocalBest);
			gsl_vector_sub(accVecLbest,pop[lpParticles].partCoord);
			/* Random weights for acceleration components */
			for (lpCoord = 0; lpCoord < nDim; lpCoord++){
				gsl_vector_set(chi1Vec,lpCoord,gsl_rng_uniform(rngGen));
			}
			for (lpCoord = 0; lpCoord < nDim; lpCoord++){
				gsl_vector_set(chi2Vec,lpCoord,gsl_rng_uniform(rngGen));
			}
			/* Multiply random weights */
			gsl_vector_mul(accVecPbest,chi1Vec);
			gsl_vector_mul(accVecLbest,chi2Vec);
			/* Scale with acceleration constants */
			gsl_vector_scale(accVecPbest,psoParams->c1);
			gsl_vector_scale(accVecLbest,psoParams->c2);
            /* Velocity update 			
		        pop(k,partVelCols)=partInertia*pop(k,partVelCols)+...
		                           c1*(pop(k,partPbestCols)-pop(k,partCoordCols))*chi1+...
		                           c2*(pop(k,partLocalBestCols)-pop(k,partCoordCols))*chi2;
			*/
			gsl_vector_scale(pop[lpParticles].partVel,pop[lpParticles].partInertia);
			gsl_vector_add(pop[lpParticles].partVel,accVecPbest);
			gsl_vector_add(pop[lpParticles].partVel,accVecLbest);
			/*
			  Apply max. velocity threshold
		        maxvBustCompPos = find(pop(k,partVelCols) > max_velocity);
		        maxvBustCompNeg = find(pop(k,partVelCols) < -max_velocity);
		        if ~isempty(maxvBustCompPos)
		            pop(k,partVelCols(maxvBustCompPos))= max_velocity;
		        end
		        if ~isempty(maxvBustCompNeg)
		            pop(k,partVelCols(maxvBustCompNeg))= -max_velocity(1);
		        end
			*/
			limitVecComponent(pop[lpParticles].partVel, - psoParams->max_velocity, psoParams->max_velocity);
			/* 
				Position Update
			*/
	        gsl_vector_add(pop[lpParticles].partCoord,pop[lpParticles].partVel);        
	    }
		
		if (psoParams->debugDumpFile != NULL){
			fprintf(psoParams->debugDumpFile,"After dynamical update\n");   
			particleInfoDump(psoParams->debugDumpFile,pop,popsize);
			fprintf(psoParams->debugDumpFile,"--------\n");			      
	    }
	}
	
	/* Prepare output */
	psoResults->totalIterations = lpPsoIter-1;
	/* 	actualEvaluations = sum(pop(:,partFitEvalsCols)); */
	psoResults->totalFuncEvals = 0;
	for (lpParticles = 0; lpParticles < popsize; lpParticles ++){
		psoResults->totalFuncEvals += pop[lpParticles].partFitEvals;
	}
	gsl_vector_memcpy(psoResults->bestLocation, gbestCoord);
	psoResults->bestFitVal = gbestFitVal;
	
	/* Free function minimizer state */
	gsl_multimin_fminimizer_free(minimzrState);
	/* Deallocate vectors */
	gsl_vector_free(locMinStp);
	gsl_vector_free(gbestCoord);
	gsl_vector_free(partSnrCurrCol);
	gsl_vector_free(accVecPbest);
    gsl_vector_free(accVecLbest); 
	gsl_vector_free(chi1Vec);
	gsl_vector_free(chi2Vec); 
	/* Deallocate members of pop */
	for(lpParticles = 0; lpParticles < popsize; lpParticles++){
		particleinfo_free(&pop[lpParticles]);
	}
}

/*Dummy function to wrap supplied function so that the 
  call to the GSL local optimizer is compatible */
double dummyfitfunc(const gsl_vector *xVec, void *dffParams){
	struct dummyFitFuncParam *dfp = (struct dummyFitFuncParam *)dffParams;
	double (*trueFitFunc)(gsl_vector *, void *) = dfp->trufuncPr;
	gsl_vector *xVec2 = (gsl_vector *)xVec;
	double funcVal = trueFitFunc(xVec2,dfp->trufuncParam);
}

/*! Initializer of particle position, velocity, and other properties. */
void initPsoParticles(struct particleInfo *p, size_t nDim, gsl_rng *rngGen){
	
	double rngNum;
	size_t lpCoord;
	
	particleinfo_alloc(p,nDim);
	
	for (lpCoord = 0; lpCoord < nDim; lpCoord++){
		rngNum = gsl_rng_uniform(rngGen);
		gsl_vector_set(p->partCoord,lpCoord,rngNum);
	}
	
	for (lpCoord = 0; lpCoord < nDim; lpCoord++){
		rngNum = gsl_rng_uniform(rngGen);
		rngNum = - gsl_vector_get(p->partCoord,lpCoord)
			     + rngNum;
		gsl_vector_set(p->partVel,lpCoord,rngNum);
	}
	
	if(gsl_vector_memcpy(p->partPbest, p->partCoord))
		printf("Error in copying vectors\n");	
	
	p->partSnrPbest = GSL_POSINF;
	p->partSnrCurr = 0;
	p->partSnrLbest = GSL_POSINF;
	p->partInertia = 0;
	p->partFitEvals = 0;
}


/*! Allocate storage for members of particleInfo structure */
void particleinfo_alloc(struct particleInfo *p, size_t nDim){
	p->partCoord = gsl_vector_alloc(nDim); /* Current coordinates */
	p->partVel = gsl_vector_alloc(nDim);  /* Current velocity */
	p->partPbest = gsl_vector_alloc(nDim); /* Coordinates of pbest */
	p->partLocalBest = gsl_vector_alloc(nDim); /* Coordinates of neighborhood best */
}

/*! Free the storage assigned to members of particleInfo structure */
void particleinfo_free(struct particleInfo *p){
	gsl_vector_free(p->partCoord);
	gsl_vector_free(p->partVel);
	gsl_vector_free(p->partPbest);
	gsl_vector_free(p->partLocalBest);	
}

/*! Allocate storage for returnData struct members */
struct returnData * returnData_alloc(size_t nDim){
	struct returnData *psoResults = (struct returnData *)malloc(sizeof(struct returnData));
	gsl_vector *bestLocation = gsl_vector_alloc(nDim);
	psoResults->bestLocation = bestLocation;
	return psoResults;
};

/*! Free storage assigned to returnData struct members */
void returnData_free(struct returnData *psoResults){
	gsl_vector_free(psoResults->bestLocation);
	free(psoResults);
};


/*! Dump information stored in particleInfo struct */
void particleinfo_fwrite(FILE *outF, struct particleInfo *p){
	
	size_t nDim = p->partCoord->size;
	size_t lpc;
	
	fprintf(outF,"Particle locations in standardized coordinates\n");
	for (lpc = 0; lpc < nDim; lpc++){
		fprintf(outF,"%f ",gsl_vector_get(p->partCoord,lpc));
	}
	fprintf(outF,"\n");
	fprintf(outF,"Particle velocities in standardized coordinates\n");
	for (lpc = 0; lpc < nDim; lpc++){
		fprintf(outF,"%f ",gsl_vector_get(p->partVel,lpc));
	}
	fprintf(outF,"\n -------- \n");
}

/*! Dump particleInfo struct array information as a matrix
with all information pertaining to one particle in a row.
*/
void particleInfoDump(FILE *outF, struct particleInfo *p, size_t popsize){	
	size_t nDim = p[0].partCoord->size;
	size_t lpParticles, lpCoord;
	
	for (lpParticles = 0; lpParticles < popsize; lpParticles++){
		for(lpCoord = 0; lpCoord < nDim; lpCoord++){
			fprintf(outF,"%lf ",gsl_vector_get(p[lpParticles].partCoord,lpCoord));
		}
		for(lpCoord = 0; lpCoord < nDim; lpCoord++){
			fprintf(outF,"%lf ",gsl_vector_get(p[lpParticles].partVel,lpCoord));
		}
		for(lpCoord = 0; lpCoord < nDim; lpCoord++){
			fprintf(outF,"%lf ",gsl_vector_get(p[lpParticles].partPbest,lpCoord));
		}
		fprintf(outF,"%lf ",p[lpParticles].partSnrPbest);
		fprintf(outF,"%lf ",p[lpParticles].partSnrCurr);	
		fprintf(outF,"%lf ",p[lpParticles].partSnrLbest);
		fprintf(outF,"%lf ",p[lpParticles].partInertia);		
		for(lpCoord = 0; lpCoord < nDim; lpCoord++){
			fprintf(outF,"%lf ",gsl_vector_get(p[lpParticles].partLocalBest,lpCoord));
		}
		fprintf(outF,"X ");	
		fprintf(outF,"%zu ",p[lpParticles].partFitEvals);		
		fprintf(outF,"\n");				
	}
}


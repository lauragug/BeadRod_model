/*
 * File: BeadSpring.h
 *
 * Contains the structure definition and functions declaration
 *
 */
#ifndef _beadSpring_h
#define _beadSpring_h

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <ctype.h>
#include <string.h>
#include "timer.h"
#include "math.h"
#include "rand.h"
#include "mpi.h"
//#include "omp.h"
//#include <atlas_enum.h>
//#include "clapack.h"

// GSL Headers
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
//#define HAVE_INLINE

#define DIM 3

double sqrt2;
double var;
double b2;
double smallDt;
double bigDt;

#define MAG(vec) sqrt(vec[0]*vec[0]+vec[1]*vec[1]+vec[2]*vec[2])

/* Simulation Paramter structure */
typedef struct _parameter{

	int NChains;				// Number of chains
	int N;					// Number of beads per chain
	int stepEquil;				// Equilibration time steps
	int stepFinal;				// Final run time steps
	double dt;				// Time step size
	double kShear[3][3];			// Velocity gradient tensor : Shear Flow
	double kExt[3][3];			// Velocity gradient tensor : Extensional Flow
	double k[3][3];
	double Pe;				// Peclet Number
	double bk;				// Kuhn Step size
	double rodLength;			// Rod length
	double beadRad;				// Bead Radius
	double nu;				//
	double tol;				// Tolerance for contraints
	int restart;
	int sampleFreq;
	int sampleMin;
	int sample;

	int startDyn;

}parameterStruct;


// InitializeNode: Intializes the Bead coordiantes using random walk 
// unCoiled = 0	Random Walk
// unCoiled = 1 x-stretch
// unCoiled = 2 y-stretch
// unCoiled = 3 z-stretch
void Initialize(parameterStruct *param, double beadCoord[][DIM], double rodQ[][DIM], double *rodT, int unCoiled, int procId);

// Initialize from restart file
void InitializeFromFile(parameterStruct *param, double beadCoord[][DIM], double rodQ[][DIM], double *rodT, int procId);

// Generates Gaussian Random Varibles
inline void GaussianRandom(double *y1, double *y2);

// Updates Weiner Process
void UpdateWeinerProcess(parameterStruct *param, double beadDW[][DIM]);
void UpdateWeinerProcessGSL(parameterStruct *param, double beadDW[][DIM], const gsl_rng *r);

// Predictor Step
void PredictorStep(parameterStruct *param, double beadCoord[][DIM], double beadFBend[][DIM], double beadFWall[][DIM], double beadDW[][DIM], double *D_RPY, double *Alpha, double rStar[][DIM]);

// Calculates tensions in the rod
void CalculateTension(parameterStruct *param, double *D_RPY, double *Alpha,double rodQ[][DIM], double rStar[][DIM], double *rodT);

// Corrector Step 
void CorrectorStep(parameterStruct *param, double rodQ[][DIM], double beadCoord[][DIM], double beadFC[][DIM], double *D_RPY, double *Alpha, double rStar[][DIM], double *rodT);

// Checks for convergence
int CheckConvergence(parameterStruct *param, double *QOld, double *QNew);

// Computes End End Length 
void ComputeEndLength(parameterStruct *param, int n, double beadCoord[][DIM], double *time, double *endendLength, double *endendProject);

// Computes Coeff for BireFringence
void ComputeBireFringenceCoeff(parameterStruct *param, int n, double rodQ[][DIM], double *var1, double *var2);

// Writes coordinates
void WriteBeadCoordinates(parameterStruct *param, int n, double beadCoord[][DIM], FILE *coord);

// Writes stresses
//void WriteStresses(parameterStruct *param, int n, rodStruct *rod, FILE *stress);

// Writes output to files
void WriteToFile(parameterStruct *param, double *endendLength, double *endendProject, double *varBF1, double *varBF2);

// Writes the restart file
void WriteRestartFile(parameterStruct *param, int n, double beadCoord[][DIM], double rodQ[][DIM], double *rodT, int procId, int chain);

// Writes output 
void WriteOutputFile(parameterStruct *param, double *time, double *endendLength, double *endendProject, FILE *restartOutPut, FILE *outPut);
#endif

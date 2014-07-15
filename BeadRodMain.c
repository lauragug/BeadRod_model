#include "BeadRod.h"

int main(int argc, char *argv[]) {

	int i,n,chain;
	parameterStruct param = {		// Contains all simulation paramters
		.NChains = 1,			// Number of chains
		.N=201,				// Number of beads per chain
		.stepEquil = 000,		// Number of time steps for equilibration
		.stepFinal = 50e+7,		// Number of time steps for production run
		.dt=0.1e-6,			// Time step size
		.kShear = {{0.,1.,0.},
			{0.,0.,0.},
			{0.,0.,0.} }, 	// Velocity gradient tensor : Shear Flow
		.kExt = {{1.,0.,0.},
			{0.,-0.5,0.},
			{0.,0.,-0.5}}, 	// Velocity gradient tensor : Extensional Flow
		.k = {{1.,0.,0.},
			{0.,-1.,0.},
			{0.,0.,0}}, 		// Velocity gradient tensor : Extensional Flow
		.Pe = 0.10,			// Peclet Number
		.bk = .132,			// Kuhn length
		.beadRad = .101,		// Bead Radius
		.rodLength = 0.,		// Rod length
		.nu = 0.0034,
		.tol = 1e-6,			// Tolerance
		.restart = 1, 
		.sampleFreq = 1000,
		.sampleMin = 0*5e7,
		.sample = 1e3,

		.startDyn = 0,
	};
	int N = param.N;
	double Pe = param.Pe;
	int sample = param.sample;//param.stepFinal/param.sampleFreq;

	double *D_RPY;				// Tensor D
	double *Alpha;				// Coefficient for weiner process	
	double endendLength[sample];		// End to end length
	double endendProject[sample];	        // Projected end to end length
	double varBF1[sample];			// For calculating Birefringence var1 = <u1*u1-u2*u2>
	double varBF2[sample];			// For calculating Birefringence var2 = <u1*u2>
	double rStar[N][3];			// Value of bead position after predictor Step
	double rodT[N-1];			// Value of rod tension
	double beadCoord[N][3];			// Bead coordinates
	double beadDW[N][3];			// Bead Weiner 
	double beadFC[N][3];			// Bead constraint force
	double beadFBend[N][3];			// Bead bending force 
	double beadFWall[N][3];			// Bead wall force 
	double rodQ[N-1][3];			// Rod orientation

	char fileName[100] = "Coords_procId";
	timeType start,end;
	FILE *histor,*coord,*stress,*outPut,*restartOutPut;

	double global_endendLength[sample],global_endendProject[sample];
	double global_varBF1[sample],global_varBF2[sample];
	double time[sample];
	int startEquil;

	int tid,nthreads;
	int numProcs,procId;

	// Global variables: to avoid computaion everytime
	sqrt2 = sqrt(2.);
	var = sqrt(param.dt);
	b2 = 1. + pow(param.tol,2.);  
	smallDt = param.dt;//*.01;
	bigDt = param.dt;

	MPI_Init(&argc,&argv);

	/* Star the timer */
	MPI_Barrier(MPI_COMM_WORLD);
	start = MPI_Wtime();

	MPI_Comm_rank(MPI_COMM_WORLD, &procId);
	MPI_Comm_size(MPI_COMM_WORLD, &numProcs);

	// Seeds for random number generator
	// Using rand function generator from GSL
	const gsl_rng_type *T;
	gsl_rng *r;

	gsl_rng_env_setup();
	T = gsl_rng_default;
	r = gsl_rng_alloc(T);
	gsl_rng_set (r,procId);

	// Open the histor.dat file
	if(procId == 0){
		histor = fopen("histor.dat","w+");
//		sprintf(fileName,"%s%d",fileName,procId);
//		coord = fopen(fileName,"w+");
		// To make sure OutPut.dat is not already present
		if(!param.restart)
			rename("OutPut.dat","OutPut_Old.dat");
	}

	// Check which chain to start with 
	// Simply check whether the RestartFile is there
	char fileRes[100] = "Restart_chain";
	char fileChain[100];
	FILE *resFile;
	int startChain = procId;

	if(param.restart > 0){
		for(i=startChain;i<param.NChains;i += numProcs){
			sprintf(fileChain,"%s%d",fileRes,i);
			resFile = fopen(fileChain,"r");
			if(resFile == NULL){	
				startChain = i-numProcs;
				break;
			}
			fclose(resFile);
		}
	}

	for(chain=startChain;chain<param.NChains;chain += numProcs){

		for(i=0;i<sample;i++){
			varBF1[i] 	= 0.;
			varBF2[i] 	= 0.;
			endendLength[i] = 0.;
			endendProject[i] = 0.;
		}

		// Initialize the chain
		Initialize(&param,beadCoord,rodQ,rodT,0,chain);

		// Check for the OutPut_Restart.dat and set other restart parameters
		char line[100];
		if(procId == 0){
			outPut = fopen("OutPut.dat","r");
			if(outPut != NULL){
				restartOutPut = fopen("OutPuttemp.dat","w+");
				for(i=0;i<param.startDyn/param.sampleFreq;i++){
					fgets(line,sizeof(line),outPut);
					fputs(line,restartOutPut);
				}
			}
			else
				restartOutPut = fopen("OutPuttemp.dat","a+");
		}

		MPI_Barrier(MPI_COMM_WORLD);

		// Initialize
		param.Pe = Pe;	
		// Production run
		for(n=param.startDyn;n<param.stepFinal;n++){
			//ComputeBireFringenceCoeff(&param,n,rodQ,varBF1,varBF2);
			ComputeEndLength(&param,n,beadCoord,time,endendLength,endendProject);
			if((n+1)%(sample*param.sampleFreq)==0)
				WriteOutputFile(&param,time,endendLength,endendProject,restartOutPut,outPut);

			CalculateWallForce(&param,beadCoord,beadFWall,beadFC);
			CalculateBendingForce(&param,beadCoord,beadFBend);
			UpdateWeinerProcessGSL(&param,beadDW,r);
			PredictorStep(&param,beadCoord,beadFBend,beadFWall,beadDW,D_RPY,Alpha,rStar);
			CalculateTension(&param,D_RPY,Alpha,rodQ,rStar,rodT);
			CorrectorStep(&param,rodQ,beadCoord,beadFC,D_RPY,Alpha,rStar,rodT);
			//WriteStresses(&param,n,rods,stress);
			if((n+1)%(sample*param.sampleFreq)==0)
				WriteRestartFile(&param,n,beadCoord,rodQ,rodT,procId,chain);

			if((procId == 0) && (n%1000000 == 0)){
				fprintf(histor,"Relaxaing Run Done: Step: %d  Chain: %d\n",n,chain);
				fflush(histor);
			}
		}

		if(procId == 0){
			if(outPut)
				fclose(outPut);
			if(restartOutPut)
				fclose(restartOutPut);
			rename("OutPuttemp.dat","OutPut.dat");
			fprintf(histor,"Done: Chain: %d  ProcId: %d\n",chain,procId);
			fflush(histor);
		}
	}

	end = MPI_Wtime();

	if(procId == 0){
		fprintf(histor,"Time (in secs): %e\n",end-start);
//		WriteToFile(&param,global_endendLength,global_endendProject,global_varBF1,global_varBF2);
		fclose(histor);
//		fclose(coord);
	}	

	gsl_rng_free (r);
	MPI_Finalize();
	return 0;
}

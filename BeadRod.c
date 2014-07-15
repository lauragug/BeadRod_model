/*
 * File: quad.c
 * --------------
 *
 */
#include "BeadRod.h"


void Initialize(parameterStruct *param, double beadCoord[][DIM], double rodQ[][DIM], double *rodT, int unCoiled, int chain){

	int i,j;
	float stepL = 1.0;
	double mag;
	int N = param->N;
	double com[3]={0,0,0};
	double length[3] = {0.,0.,0.};

	// Checking if restart file exists
	char fileName[100] = "Restart_chain";
	FILE *restartFile;
	sprintf(fileName,"%s%d",fileName,chain);
	restartFile = fopen(fileName,"r");
	if(param->restart && restartFile){
		InitializeFromFile(param,beadCoord,rodQ,rodT,chain);
		return;
	}

	param->startDyn = 0;
	// First Rod is Straight
	rodQ[0][0] = 0.;
	rodQ[0][1] = 1.;
	rodQ[0][2] = 0.; 
	for(j=0;j<3;j++){
		beadCoord[0][j] = 0.0;
		beadCoord[1][j] = beadCoord[0][j] + rodQ[0][j];
	}

	if (unCoiled == 0){
		for(i=1;i<N-1;i++){

initialize:
			rodQ[i][0] = 0.;//frand(-stepL,stepL);
			rodQ[i][1] = 1.;//frand(-stepL,stepL);
			rodQ[i][2] = 0.;//frand(-stepL,stepL); 

			mag = pow(pow(rodQ[i][0],2) + pow(rodQ[i][1],2) + pow(rodQ[i][2],2),0.5);	

			rodQ[i][0] /= mag;
			rodQ[i][1] /= mag;
			rodQ[i][2] /= mag;
			rodT[i] = 0.;

			beadCoord[i+1][0] = beadCoord[i][0] + rodQ[i][0];
			beadCoord[i+1][1] = beadCoord[i][1] + rodQ[i][1];
			beadCoord[i+1][2] = beadCoord[i][2] + rodQ[i][2];

			if(beadCoord[i+1][1] < 0.)
				goto initialize;
		}
	}
	else
	{
		length[unCoiled-1] = 1.;
		for(i=0;i<N-1;i++){
			rodT[i] = 0.;
			for(j=0;j<3;j++){
				rodQ[i][j] = length[j];
				beadCoord[i+1][j] = beadCoord[i][j] + rodQ[i][j];
				com[j] += beadCoord[i+1][j];	
			}
		}
	}
}


void InitializeFromFile(parameterStruct *param, double beadCoord[][DIM], double rodQ[][DIM], double *rodT, int chain){

	double time;
	int beadCntr,segCntr;
	char fileName[100] = "Restart_chain";
	char line[100];
	FILE *restartFile;

	sprintf(fileName,"%s%d",fileName,chain);
	restartFile = fopen(fileName,"r");

	fgets(line,sizeof(line),restartFile);
	sscanf(line,"%ld",&param->startDyn);

	fgets(line,sizeof(line),restartFile);
	sscanf(line,"%le",&time);

	beadCntr = 0;
	while(beadCntr < param->N){
		fgets(line,sizeof(line),restartFile);
		sscanf(line,"%le%le%le",&beadCoord[beadCntr][0],&beadCoord[beadCntr][1],&beadCoord[beadCntr][2]);
		beadCntr++;
	}

	segCntr = 0;
	while(segCntr < (param->N-1)){
		fgets(line,sizeof(line),restartFile);
		sscanf(line,"%le%le%le",&rodQ[segCntr][0],&rodQ[segCntr][1],&rodQ[segCntr][2]);
		segCntr++;
	}

	segCntr = 0;
	while(segCntr < (param->N-1)){
		fgets(line,sizeof(line),restartFile);
		sscanf(line,"%le",&rodT[segCntr]);
		segCntr++;
	}
}

void UpdateWeinerProcessGSL(parameterStruct *param, double beadDW[][DIM], const gsl_rng *r){

	int i=0;
	int N = param->N;
	var = sqrt(param->dt);

	for(i=2;i<N;i++){
		beadDW[i][0] = gsl_ran_gaussian_ziggurat(r,var);
		beadDW[i][1] = gsl_ran_gaussian_ziggurat(r,var);
		beadDW[i][2] = gsl_ran_gaussian_ziggurat(r,var);
	}

	beadDW[0][0] = 0.; 
	beadDW[0][1] = 0.;
	beadDW[0][2] = 0.;
	beadDW[1][0] = 0.; 
	beadDW[1][1] = 0.;
	beadDW[1][2] = 0.;
}

void UpdateWeinerProcess(parameterStruct *param, double beadDW[][DIM]){

	int i=0;
	double b[3];
	int N = param->N;
	var = sqrt(param->dt);

	if (N%2 != 0){
		GaussianRandom(&beadDW[0][0],&b[0]);
		GaussianRandom(&beadDW[0][1],&b[1]);
		GaussianRandom(&beadDW[0][2],&b[2]);

		for(i=1;i<N;i+=2){
			GaussianRandom(&beadDW[i][0],&beadDW[i+1][0]);
			GaussianRandom(&beadDW[i][1],&beadDW[i+1][1]);
			GaussianRandom(&beadDW[i][2],&beadDW[i+1][2]);
		}
	}
	else{
		for(i=0;i<N;i+=2){
			GaussianRandom(&beadDW[i][0],&beadDW[i+1][0]);
			GaussianRandom(&beadDW[i][1],&beadDW[i+1][1]);
			GaussianRandom(&beadDW[i][2],&beadDW[i+1][2]);
		}
	}

	beadDW[0][0] = 0.; 
	beadDW[0][1] = 0.;
	beadDW[0][2] = 0.;
	for(i=1;i<N;i++){
		beadDW[i][0] *= var; 
		beadDW[i][1] *= var;
		beadDW[i][2] *= var;
	}
}

inline void GaussianRandom(double *y1, double *y2){

	double x1, x2, w;

	do {
		x1 = 2.0 * frand(0.,1.) - 1.0;
		x2 = 2.0 * frand(0.,1.) - 1.0;
		w = x1 * x1 + x2 * x2;
	} while ( w >= 1.0 );

	w = sqrt( (-2.0 * log( w ) ) / w );
	*y1 = x1 * w;
	*y2 = x2 * w;
}


static long dgtsv(long N, long NRHS, double *DL, double *D, double *DU, double *B, long LDB)
{
	extern void dgtsv_(const long *Np, const long *NRHSp, double *DL,
			double *D, double *DU, double *B, const long *LDBp,
			long *INFOp);
	long info;
	dgtsv_(&N, &NRHS, DL, D, DU, B, &LDB, &info);
	return info;
}


static long dgttrf(long N, double *DL, double *D, double *DU, double *DU2, int *IPIV)
{
	extern void dgttrf_(const long *Np, double *DL, double *D, double *DU, double *DU2, int *IPIV, long *INFOp);
	long info;
	dgttrf_(&N, DL, D, DU, DU2, IPIV, &info);
	return info;
}


static long dgttrs(char C, long N, long NRHS, double *DL, double *D, double *DU, 
		double *DU2, int *IPIV, double *B, long LDB)
{
	extern void dgttrs_(const char *Cp, const long *Np, const long *NRHSp, 
			double *DL, double *D, double *DU, double *DU2, 
			int *IPIV, double *B, const long *LDBp, long *INFOp);
	long info;
	dgttrs_(&C, &N, &NRHS, DL, D, DU, DU2, IPIV, B, &LDB, &info);
	return info;
}


void CalculateBendingForce(parameterStruct *param, double beadCoord[][DIM], double beadFBend[][DIM]){

	int i,j;
	int N = param->N;
	double lpbyl = 50000;

	if(N==2){
		for(i=0;i<=N;i++)
			for(j=0;j<3;j++)
				beadFBend[i][j] = 0.;

		return;
	}	

	if(N == 3){
		for(j=0;j<3;j++){
			i = 0;
			beadFBend[i][j] = lpbyl*(-beadCoord[i][j] + 2.*beadCoord[i+1][j] - beadCoord[i+2][j]);

			i = 2;
			beadFBend[i][j] = lpbyl*(-beadCoord[i-1][j] + 2.*beadCoord[i][j] - beadCoord[i+1][j]);

			i = 2;
			beadFBend[i][j] = lpbyl*(-beadCoord[i][j] + 2.*beadCoord[i-1][j] - beadCoord[i-2][j]);
		}
		return;
	}

	for(i=2;i<=N-3;i++){
		for(j=0;j<3;j++)
			beadFBend[i][j] = lpbyl*(-beadCoord[i-2][j] + 4.*beadCoord[i-1][j] 
					- 6.*beadCoord[i][j] + 4.*beadCoord[i+1][j] - beadCoord[i+2][j]);
	}

	for(j=0;j<3;j++){
		i = 1;
		beadFBend[i][j] = lpbyl*(2*beadCoord[i-1][j] - 5.*beadCoord[i][j] 
				+ 4.*beadCoord[i+1][j] - beadCoord[i+2][j]);

		i = N-2; 
		beadFBend[i][j] = lpbyl*(2*beadCoord[i+1][j] - 5.*beadCoord[i][j] 
				+ 4.*beadCoord[i-1][j] - beadCoord[i-2][j]);

		i = 0;
		beadFBend[i][j] = lpbyl*(-beadCoord[i][j] + 2.*beadCoord[i+1][j] - beadCoord[i+2][j]);

		i = N-1;
		beadFBend[i][j] = lpbyl*(-beadCoord[i][j] + 2.*beadCoord[i-1][j] - beadCoord[i-2][j]);
	}
	// Bead 0 is tethered
	beadFBend[0][0] = 0.;
	beadFBend[0][1] = 0.;
	beadFBend[0][2] = 0.;
	beadFBend[1][0] = 0.;
	beadFBend[1][1] = 0.;
	beadFBend[1][2] = 0.;
}


void CalculateWallForce(parameterStruct *param, double beadCoord[][DIM], double *beadFWall){

	int i,j;
	int N = param->N;
	double hjB,hjT; 
	double rc = 0.3;
	double Cw = 35;		
	double wallBOT = 0.;

	// Update bead force from spring force	
	for(i=2;i<N;i++){
		beadFWall[3*i+0] = 0.;
		beadFWall[3*i+2] = 0.;
		
		hjB = +(beadCoord[i][1])-(wallBOT);
		if(hjB < 0){
			fprintf(stderr,"hj < 0\n");
			abort();
		}

		if(hjB < 0.25)
			beadFWall[3*i+1] = Cw*pow(rc/0.25,6);
		else if(hjB > 0.25 && hjB < 0.45)
			beadFWall[3*i+1] = Cw*pow(rc/hjB,6);
		else
			beadFWall[3*i+1] = 0;
	}

	i = 0;
	beadFWall[3*i+0] = 0;
	beadFWall[3*i+1] = 0;
	beadFWall[3*i+2] = 0;

	i = 1;
	beadFWall[3*i+0] = 0;
	beadFWall[3*i+1] = 0;
	beadFWall[3*i+2] = 0;
}


void PredictorStep(parameterStruct *param, double beadCoord[][DIM], double beadFBend[][DIM], double beadFWall[][DIM], double beadDW[][DIM], double *D_RPY, double *Alpha, double rStar[][DIM]){

	int i,j;
	int N = param->N;
	double Pe = param->Pe;
	double dt = param->dt;

	for(i=0;i<N;i++){
		for(j=0;j<3;j++){

			rStar[i][j] = beadCoord[i][j] + Pe*(param->k[j][0]*beadCoord[i][0] +
				param->k[j][1]*beadCoord[i][1]+param->k[j][2]*beadCoord[i][2])*dt +
				(beadFBend[i][j] + beadFWall[i][j])*dt + sqrt2*beadDW[i][j];	
		}
	}
}


void CalculateTension(parameterStruct *param, double *D_RPY, double *Alpha, double rodQ[][DIM], double rStar[][DIM], double *T){

	int i,j;
	int N = param->N;
	double dt = param->dt;
	double B[N-1];
	double D[N-1], DL[N-2], DU[N-2], DU2[N-3];
	int IPIV[N-1];
	double bstar[N-1][3];
	double nonLinear;
	double residual;
	double TOL = 1e-6;	
	int info;
	int tid;

	// Diagonal of A for AX = B
	for(i=0;i<N-1;i++){
		D[i] = 0.;
		B[i] = T[i];
		for(j=0;j<3;j++){
			bstar[i][j] = rStar[i+1][j]-rStar[i][j];
			D[i] -= bstar[i][j]*rodQ[i][j];
		}
		D[i] *= 4*dt;
	}

	// Lower Upper Diagonal of A for AX = B 
	for(i=0;i<N-2;i++){
		DL[i] = 0.;
		DU[i] = 0.;
		for(j=0;j<3;j++){
			DL[i] += bstar[i+1][j]*rodQ[i][j];
			DU[i] += bstar[i][j]*rodQ[i+1][j];
		}
		DL[i] *= 2*dt;
		DU[i] *= 2*dt;
	}

	// Compute the LU decomposition of A
	info = dgttrf(N-1, DL, D, DU, DU2, IPIV);
	if (info != 0) fprintf(stderr, "failure with error %d\n", info);

	do{
		for(i=0;i<N-1;i++){
			T[i] = B[i];
		}

		i = 0;
		B[i] = b2;
		for(j=0;j<3;j++){
			nonLinear = - 2.*T[i]*rodQ[i][j] + T[i+1]*rodQ[i+1][j];
			B[i] -= (bstar[i][j]*bstar[i][j] + nonLinear*nonLinear*dt*dt);
		}


		i = N-2;
		B[i] = b2;
		for(j=0;j<3;j++){
			nonLinear = T[i-1]*rodQ[i-1][j] - 2.*T[i]*rodQ[i][j];
			B[i] -= (bstar[i][j]*bstar[i][j] + nonLinear*nonLinear*dt*dt);
		}

		for(i=1;i<N-2;i++){
			B[i] = b2;
			for(j=0;j<3;j++){
				nonLinear = T[i-1]*rodQ[i-1][j] - 2.*T[i]*rodQ[i][j] + T[i+1]*rodQ[i+1][j];
				B[i] -= (bstar[i][j]*bstar[i][j] + nonLinear*nonLinear*dt*dt);
			}
		}

		// Solve the tridiagonal system
		info = dgttrs('N', N-1, 1, DL, D, DU, DU2, IPIV, B, N-1);
		if (info != 0) fprintf(stderr, "failure with error %d\n", info);

		residual = 0.;
		for(i=0;i<N-1;i++){
			residual += pow((B[i]-T[i]),2);
		}
		residual = sqrt(residual);
		//printf("%e\n",residual);

	}while (residual > TOL);

	for(i=0;i<N-1;i++)
		T[i] = B[i];
}


void CorrectorStep(parameterStruct *param, double rodQ[][DIM], double beadCoord[][DIM], double beadFC[][DIM], double *D_RPY, double *Alpha, double rStar[][DIM], double *rodT){

	int i,j;
	int N = param->N;
	double dt = param->dt;
	double com[3]={0,0,0};
	int tid;

	// First Bead is tethered
	// F_tethered = -F_All
//	i = 0;
//	for(j=0;j<3;j++){
//		beadFC[i][j] = rodT[i]*rodQ[i][j];
//		beadCoord[i][j] = rStar[i][j] + beadFC[i][j]*dt;
//		beadFTeth[j] = -beadCoord[i][j]/dt;
//		beadCoord[i][j] = 0.;
//	}

	for(i=2;i<N-1;i++)
		for(j=0;j<3;j++){
			beadFC[i][j] = rodT[i]*rodQ[i][j] - rodT[i-1]*rodQ[i-1][j];
			beadCoord[i][j] = rStar[i][j] + beadFC[i][j]*dt;
			rodQ[i-1][j] = beadCoord[i][j] - beadCoord[i-1][j];
		}

	i = N-1;	
	for(j=0;j<3;j++){
		beadFC[i][j] = -rodT[i-1]*rodQ[i-1][j];
		beadCoord[i][j] = rStar[i][j] + beadFC[i][j]*dt;
		rodQ[i-1][j] = beadCoord[i][j] - beadCoord[i-1][j];
	}
}


/*int CheckConvergence(parameterStruct *param, double *QOld, double *QNew){

	int i,converged = 1;
	double residual=0.0;
	double TOL = 1e-8;
	int N = param->N;

	// Residual computation for interior nodes
	#pragma omp for private(i) reduction(+:residual)
	for(i=0;i<N-1;i++){
		residual += pow((QOld[i]-QNew[i]),2);
	}
	residual = pow(residual,0.5);
	printf("%e\n",residual);
	// Convergence criterion
	if(residual < TOL)
		converged = 0; 

	return converged;
}
*/

void ComputeEndLength(parameterStruct *param, int n, double beadCoord[][DIM], double *time, double *endendLength, double *endendProject){

	int i,j;
	int N = param->N;
	int nchains = param->NChains;
	double Rend2 = 0.;
	double Rmax, Rmin, temp; 
	int index;

	if((n%param->sampleFreq) != 0)
		return;

	// Implimented only for 0 degree flow direction
	Rmax = Rmin = beadCoord[0][0];
	
	for(i=1;i<N;i++){
	
		temp = beadCoord[i][0];	
		if(temp > Rmax)
			Rmax = temp;
		else if(temp < Rmin)
			Rmin = temp;
	}

	for(j=0;j<3;j++){
		Rend2 += pow(beadCoord[N-1][j]-beadCoord[0][j],2.);
	}

	index = n/param->sampleFreq;
	index = index%param->sample;
	endendProject[index] += pow(Rmax-Rmin,2)/nchains;
	endendLength[index] += Rend2/nchains;
	time[index] = n*param->dt;
}


void ComputeBireFringenceCoeff(parameterStruct *param, int n, double rodQ[][DIM], double *var1, double *var2){

	int i;
	int N = param->N;
	int nchains = param->NChains;
	double variable1=0., variable2=0.;
	int index;

	if((n%param->sampleFreq) != 0)
		return;

	for(i=0;i<N-1;i++){
		variable1 += (pow(rodQ[i][0],2.)-pow(rodQ[i][1],2.));	
		variable2 += rodQ[i][0]*rodQ[i][1];
	}
	variable1 /= nchains;
	variable2 /= nchains;

	index = n/param->sampleFreq;
	{
		var1[index] += variable1;
		var2[index] += variable2;
	}
}


void WriteBeadCoordinates(parameterStruct *param, int n, double beadCoord[][DIM], FILE *coord){

	int i;
	double dt = param->dt;

	if(n < param->sampleMin)
		return;

	if((n%param->sampleFreq) != 0)
		return;

	fprintf(coord,"%e",n*dt);

	for(i=0;i<param->N;i+=5){
		fprintf(coord,"\t%e\t%e\t%e",beadCoord[i][0],beadCoord[i][1],beadCoord[i][2]);
	}
	fprintf(coord,"\n");
}

/*
void WriteStresses(parameterStruct *param, int n, rodStruct *rod, FILE *stress){

	int i;
	double dt = param->dt;
	double Txy = 0.;

	if((n%1000) != 0)
		return;

	for(i=0;i<param->N-1;i++){
		Txy += rod[i].Q[0]*rod[i].Q[1]*rod[i].tension;
	}
	
	fprintf(stress,"%e\t%e\n",n*dt,Txy);
}
*/

void WriteToFile(parameterStruct *param, double *endendLength, double *endendProject, double *varBF1, double *varBF2){

	int i;
	double sample_dt = param->dt*param->sampleFreq;
	double bireFringence;
	FILE *file1;
	int sampleSize = param->stepFinal/param->sampleFreq;

	file1 = fopen("OutPut.dat","a+");

	for(i=0;i<sampleSize;i++){
				
			bireFringence = pow(pow(varBF1[i],2.) + 4.*pow(varBF2[i],2.),0.5);
			fprintf(file1,"%e\t%e\t%e\t%e\n",i*sample_dt,endendLength[i],endendProject[i],bireFringence);	
	}

	fclose(file1);
}

void WriteRestartFile(parameterStruct *param, int n, double beadCoord[][DIM], double rodQ[][DIM], double *rodT, int procId, int chain){

	int i;
	FILE *restartFile;
	char fileName[100] = "Restart_chain";	
		
	sprintf(fileName,"%s%d",fileName,chain);
	restartFile = fopen(fileName,"w+");	

	fprintf(restartFile,"%d\n",n+1);
	fprintf(restartFile,"%le\n",n*param->dt);

	for(i=0;i<param->N;i++)
		fprintf(restartFile,"%.10f\t%.10f\t%.10f\n",beadCoord[i][0],beadCoord[i][1],beadCoord[i][2]);

	for(i=0;i<(param->N-1);i++)
		fprintf(restartFile,"%.10f\t%.10f\t%.10f\n",rodQ[i][0],rodQ[i][1],rodQ[i][2]);

	for(i=0;i<(param->N-1);i++)
		fprintf(restartFile,"%.10f\n",rodT[i]);

	fclose(restartFile);

	MPI_Barrier(MPI_COMM_WORLD);
}

void WriteOutputFile(parameterStruct *param, double *time, double *endendLength, double *endendProject, FILE *restartOutPut, FILE *outPut){

	int i;
	double sample_dt = param->dt*param->sampleFreq;
	FILE *outputFile;
	int sample = param->sample;
	double global_endendLength[sample],global_endendProject[sample];
	double tmp_time,tmp_endendLength,tmp_endendProject;
	int procId;
	char line[100];

	MPI_Comm_rank(MPI_COMM_WORLD, &procId);

	MPI_Reduce(endendLength,global_endendLength,sample,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);	
	MPI_Reduce(endendProject,global_endendProject,sample,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);	

	for(i=0;i<sample;i++)
		endendLength[i] = endendProject[i] = 0.;

	if(procId == 0){
		if(outPut != NULL){
			fgets(line,sizeof(line),outPut);
			// Read the old OutPut.dat file
			if(!feof(outPut)){
				for(i=0;i<sample;i++){
					if(i>0)
						fgets(line,sizeof(line),outPut);
					sscanf(line,"%le%le%le",&tmp_time,&tmp_endendLength,&tmp_endendProject);

					global_endendLength[i] += tmp_endendLength;
					global_endendProject[i] += tmp_endendProject;

					fprintf(restartOutPut,"%.10f\t%.10f\t%.10f\n",time[i],global_endendLength[i],global_endendProject[i]);
				}
			}
			else{
				for(i=0;i<sample;i++){
					fprintf(restartOutPut,"%.10f\t%.10f\t%.10f\n",time[i],global_endendLength[i],global_endendProject[i]);
				}
			}
		}
		else{
			for(i=0;i<sample;i++){
				fprintf(restartOutPut,"%.10f\t%.10f\t%.10f\n",time[i],global_endendLength[i],global_endendProject[i]);
			}
		}
	fflush(restartOutPut);
	}
}

#include <ilcplex/ilocplex.h>
#include <cmath>
#include <unordered_set>
#include <string>
#include <sstream>
#include <vector>
#include <cstddef>  // std::size_t
#include <numeric>
#include <set>
#include <limits>
#include <algorithm>
#include <omp.h>
#include <chrono>
#include <time.h>
#include <sys/time.h>
#include "global.h"
#include "functions.h"
#include "mt64.h"


using namespace std;

double get_wall_time()
{
	struct timeval time;
	if (gettimeofday(&time,NULL)){
		return 0;
	}
	return (double)time.tv_sec + (double)time.tv_usec*.000001;
}

double get_cpu_time()
{
	return (double)clock() / CLOCKS_PER_SEC;
}

void readData (FormatData * fData_p)
{
	//! create data environment
	fData_p->dataEnv = IloEnv();
	IloEnv * env = & fData_p->dataEnv;

	cout << "Start reading data from files... " << endl;

	//! read number of stages
	readArray<IloInt> (fData_p->numStage, "../data/numStage.dat");
	IloInt numStage = fData_p->numStage;

	readArray<IloInt> (fData_p->breakstage, "../data/breakstage.dat");

	fData_p->initState = IloNumArray(*env);
	readArray<IloNumArray> (fData_p->initState, "../data/initState.dat");
	fData_p->initStateBin = IloNumArray(*env);
	readArray<IloNumArray> (fData_p->initStateBin, "../data/initStateBin.dat");

	readArray<IloInt> (fData_p->numFWsample, "../data/numFWsample.dat");
	
	fData_p->numScen = IloIntArray(*env);
	readArray<IloIntArray> (fData_p->numScen, "../data/numScen.dat");
	
	fData_p->totalScen = fData_p->numScen[0];
	for (int t = 1; t < numStage; ++t)
		fData_p->totalScen *= fData_p->numScen[t];

	readArray<IloInt> (fData_p->intX, "../data/intX.dat");

	fData_p->xBound = IloNumArray2(*env);
	readArray<IloNumArray2> (fData_p->xBound, "../data/xBound.dat");

	fData_p->thetaBound = IloNumArray2(*env);
	readArray<IloNumArray2> (fData_p->thetaBound, "../data/thetaBound.dat");

	fData_p->T = IloNumArray2(*env);
	readArray<IloNumArray2> (fData_p->T, "../data/T.dat");

	fData_p->TT = IloNumArray2(*env);
	readArray<IloNumArray2> (fData_p->TT, "../data/TT.dat");

	cout << "Basic SDDIP data read." << endl;
	
	//! read objective coefficients for x (original) variables 
	fData_p->x = IloNumArray2(*env);
	readArray<IloNumArray2> (fData_p->x, "../data/x.dat");

	fData_p->xBin = IloNumArray2(*env);
	readArray<IloNumArray2> (fData_p->xBin, "../data/xBin.dat");

	fData_p->y1 = IloNumArray2(*env);
	readArray<IloNumArray2> (fData_p->y1, "../data/y1.dat");

	fData_p->y2 = IloNumArray2(*env);
	readArray<IloNumArray2> (fData_p->y2, "../data/y2.dat");

	//! read matrices in constraint A_tx_t + W1_ty1_t + W2_ty2_t + B_tz_t >= b_t
	fData_p->A = IloNumArray2(*env);
	readArray<IloNumArray2> (fData_p->A, "../data/A.dat");

	fData_p->ABin = IloNumArray2(*env);
	readArray<IloNumArray2> (fData_p->ABin, "../data/ABin.dat");

	fData_p->B = IloNumArray2(*env);
	readArray<IloNumArray2> (fData_p->B, "../data/B.dat");

	fData_p->BBin = IloNumArray2(*env);
	readArray<IloNumArray2> (fData_p->BBin, "../data/BBin.dat");

	fData_p->W1 = IloNumArray2(*env);
	readArray<IloNumArray2> (fData_p->W1, "../data/W1.dat");

	fData_p->W2 = IloNumArray2(*env);
	readArray<IloNumArray2> (fData_p->W2, "../data/W2.dat");
	
	fData_p->b = IloNumArray2(*env);
	readArray<IloNumArray2> (fData_p->b, "../data/rhs.dat");

	cout << "Initialization data read." << endl;

	fData_p->uncertaintySource = IloIntArray(*env);
	readArray<IloIntArray> (fData_p->uncertaintySource, "../data/uncertaintySource.dat");

	if ( fData_p->uncertaintySource[0] )
	{
		fData_p->scen.x = IloNumArray3(*env);
		fData_p->scen.xBin = IloNumArray3(*env);
		readArray<IloNumArray3> (fData_p->scen.x, "../data/scenX.dat");
		readArray<IloNumArray3> (fData_p->scen.xBin, "../data/scenXBin.dat");
	}
	if ( fData_p->uncertaintySource[1] )
	{
		fData_p->scen.y1 = IloNumArray3(*env);
		readArray<IloNumArray3> (fData_p->scen.y1, "../data/scenY1.dat");
	}
	if ( fData_p->uncertaintySource[2] )
	{
		fData_p->scen.y2 = IloNumArray3(*env);
		readArray<IloNumArray3> (fData_p->scen.y2, "../data/scenY2.dat");
	}
	if ( fData_p->uncertaintySource[3] )
	{
		fData_p->scen.A = IloNumArray4(*env);
		fData_p->scen.ABin = IloNumArray4(*env);
		readArray<IloNumArray4> (fData_p->scen.A, "../data/scenA.dat");
		readArray<IloNumArray4> (fData_p->scen.ABin, "../data/scenABin.dat");
	}
	if ( fData_p->uncertaintySource[4] )
	{
		fData_p->scen.B = IloNumArray4(*env);
		fData_p->scen.BBin = IloNumArray4(*env);
		readArray<IloNumArray4> (fData_p->scen.B, "../data/scenB.dat");
		readArray<IloNumArray4> (fData_p->scen.BBin, "../data/scenBBin.dat");
	}
	if ( fData_p->uncertaintySource[5] )
	{
		fData_p->scen.W1 = IloNumArray4(*env);
		readArray<IloNumArray4> (fData_p->scen.W1, "../data/scenW1.dat");
	}
	if ( fData_p->uncertaintySource[6] )
	{
		fData_p->scen.W2 = IloNumArray4(*env);
		readArray<IloNumArray4> (fData_p->scen.W2, "../data/scenW2.dat");
	}
	if ( fData_p->uncertaintySource[7] )
	{
		fData_p->scen.b = IloNumArray3(*env);
		readArray<IloNumArray3> (fData_p->scen.b, "../data/scenRHS.dat");
	}

	return;
}

void buildModel (Model * models, FormatData * fData_p)
{
	cout << "Start to build one model for each stage..." << endl;

	int t, i, j, k;
	
	IloInt dimX, dimZ;
	IloInt dimY1 = fData_p->y1[0].getSize();  // current stage integral variables
	IloInt dimY2 = fData_p->y2[0].getSize();  // current stage continuous variables
	IloInt nconstr1 = fData_p->b[0].getSize();  // number of constraints


	/*!
		Verify data dimensions
	*/
	cout << "number of constraints 1: " << nconstr1 << endl;
	cout << "x(z) dim:" << fData_p->x[0].getSize() << endl;
	cout << "x(z) dim (discretized):" << fData_p->xBin[0].getSize() << endl;
	cout << "y1 dim:" << dimY1 << endl;
	cout << "y2 dim:" << dimY2 << endl;
	cout << "A dim: " << fData_p->A.getSize() << 
			"x" << fData_p->A[0].getSize() << endl;
	cout << "A dim (discretized): " << fData_p->ABin.getSize() << 
			"x" << fData_p->ABin[0].getSize() << endl;
	cout << "B dim: " << fData_p->B.getSize() <<
			"x" << fData_p->B[0].getSize() << endl;
	cout << "B dim (discretized): " << fData_p->BBin.getSize() << 
			"x" << fData_p->BBin[0].getSize() << endl;
	cout << "W1 dim: " << fData_p->W1.getSize() << 
			"x" << fData_p->W1[0].getSize() << endl;
	cout << "W2 dim: " << fData_p->W2.getSize() << 
			"x" << fData_p->W2[0].getSize() << endl;
	
	/*!
	  Create subproblem for each stage
	*/
	for (t = 0; t < fData_p->numStage; ++t)
	{
		//! create cplex environment and initialize model
		IloEnv newEnv = IloEnv();
		models[t].mod = IloModel(newEnv);

		/*!
			construct state variables (x), local integer variables (y1),
			and local copy of state variables (z).

			Breakstage dictates the state space change. Range of breakstage: [0, T].
			Once breakstage is chosen:
				a. t < breakstage: stage problem t uses binary state variables
				b. t >= breakstage: stage problem t uses mixed integer state variables
			breakstage = 0, all stage problems use original state variables
			breakstage = T, all stage problems use binary state variables
		*/

		if ( t < fData_p->breakstage )
		{
			dimX = fData_p->xBin[0].getSize();
			dimZ = dimX;
			models[t].x = IloNumVarArray(newEnv, dimX, 0.0, 1.0, ILOINT);
			models[t].z = IloNumVarArray(newEnv, dimZ, 0.0, 1.0, ILOFLOAT);
		}
		else
		{
			dimX = fData_p->x[0].getSize();
			dimZ = dimX;

			//! split bounds for integer and continuous state variables
			IloNumArray ix_lb(newEnv), ix_ub(newEnv), cx_lb(newEnv), cx_ub(newEnv);
			for ( i = 0; i < dimX; i++ )
			{
				if ( i < fData_p->intX )
				{
					ix_lb.add(fData_p->xBound[0][i]);
					ix_ub.add(fData_p->xBound[1][i]);
				}
				else
				{
					cx_lb.add(fData_p->xBound[0][i]);
					cx_ub.add(fData_p->xBound[1][i]);
				}
			}
		
			models[t].x = IloNumVarArray(newEnv);
			models[t].x.add(IloNumVarArray(newEnv, ix_lb, ix_ub, ILOINT));
			models[t].x.add(IloNumVarArray(newEnv, cx_lb, cx_ub, ILOFLOAT));
			models[t].z = IloNumVarArray(newEnv, fData_p->xBound[0], fData_p->xBound[1], ILOFLOAT);
		}

		//! construct local variables
		models[t].y1 = IloNumVarArray(newEnv, dimY1, 0, IloInfinity, ILOINT);
		models[t].y2 = IloNumVarArray(newEnv, dimY2, 0, IloInfinity, ILOFLOAT);
		// TODO: include variable bounds when creating problems.
		// models[t].y1 = IloNumVarArray(newEnv, dimY1, fData_p->y1Bound[0], fData_p->y1Bound[1], ILOINT);
		// models[t].y2 = IloNumVarArray(newEnv, dimY2, fData_p->y2Bound[0], fData_p->y2Bound[1], ILOFLOAT);

		//! construct theta variables
		models[t].theta = IloNumVar(newEnv, fData_p->thetaBound[0][t], fData_p->thetaBound[1][t]);

		//! set variable names
		char varName[100];
		for ( i = 0; i < dimX; ++i )
		{
			sprintf(varName, "x_%d", i+1);
			models[t].x[i].setName(varName);
		}
		
		for ( i = 0; i < dimY1; ++i )
		{
			sprintf(varName, "y1_%d", i+1);
			models[t].y1[i].setName(varName);
		}
		
		for ( i = 0; i < dimY2; ++i )
		{
			sprintf(varName, "y2_%d", i+1);
			models[t].y2[i].setName(varName);
		}

		for ( i = 0; i < dimZ; ++i )
		{
			sprintf(varName, "z_%d", i+1);
			models[t].z[i].setName(varName);
		}

		sprintf(varName, "theta");
		models[t].theta.setName(varName);

		//! add variables to model
		models[t].mod.add(models[t].x);
		models[t].mod.add(models[t].y1);
		models[t].mod.add(models[t].y2);		
		models[t].mod.add(models[t].z);
		models[t].mod.add(models[t].theta);

		//! construct continuous relaxation
		models[t].xLP = IloConversion(newEnv, models[t].x, ILOFLOAT);
		models[t].yLP = IloConversion(newEnv, models[t].y1, ILOFLOAT);

		// cout << "variables construction completed." << endl;

		//! create objective function and add to model
		models[t].obj = IloObjective(newEnv);
		IloExpr objExpr(newEnv);
		if ( t < fData_p->breakstage )
			objExpr = IloScalProd(models[t].x, fData_p->xBin[t]); 
		else 
			objExpr = IloScalProd(models[t].x, fData_p->x[t]);
		objExpr += IloScalProd(models[t].y1, fData_p->y1[t]) + 
				IloScalProd(models[t].y2, fData_p->y2[t]);
		if ( t < fData_p->numStage )
			objExpr += models[t].theta;
		models[t].obj.setExpr(objExpr);
		models[t].obj.setSense(IloObjective::Minimize);
		char objName[100];
		sprintf(objName, "objective_%d", t);
		models[t].obj.setName(objName);
		models[t].mod.add(models[t].obj);
		objExpr.end();

		// cout << "objective construction completed." << endl;

		//! Add constraints A_tx_t + B_tz_t + W1_ty1_t + W2_ty2_t >= b_t
		models[t].constr1 = IloRangeArray(newEnv);
		if ( t < fData_p->breakstage )
		{
			for (i = 0; i < nconstr1; ++i)
			{
				IloExpr expr = IloScalProd(models[t].x, fData_p->ABin[i]) + 
								IloScalProd(models[t].z, fData_p->BBin[i]) +
								IloScalProd(models[t].y1, fData_p->W1[i]) +
								IloScalProd(models[t].y2, fData_p->W2[i]);
				models[t].constr1.add(expr >= fData_p->b[t][i]);
				expr.end();
			}
		}
		else
		{
			for (i = 0; i < nconstr1; ++i)
			{
				IloExpr expr = IloScalProd(models[t].x, fData_p->A[i]) + 
								IloScalProd(models[t].z, fData_p->B[i]) +
								IloScalProd(models[t].y1, fData_p->W1[i]) +
								IloScalProd(models[t].y2, fData_p->W2[i]);
				models[t].constr1.add(expr >= fData_p->b[t][i]);
				expr.end();
			}
		}
		models[t].mod.add(models[t].constr1);

		// cout << "A_tx_t + B_tz_t + W1_ty1_t + W2_ty2_t >= b_t constructed." << endl;

		//! Add constraints z_t = x_{t-1}, rhs initialized as 0
		models[t].constr2 = IloRangeArray(newEnv);
		for (i = 0; i < dimZ; ++i)
		{
			IloExpr expr = models[t].z[i];
			//! initialize the state variable for models[0]
			if ( t == 0 )
			{
				if ( t < fData_p->breakstage )
					models[t].constr2.add(expr == fData_p->initStateBin[i]);
				else
					models[t].constr2.add(expr == fData_p->initState[i]);
			}
			else  //! all others left as zero for future update
				models[t].constr2.add(expr == 0.0);
			expr.end();
		}
		models[t].mod.add(models[t].constr2);

		// cout << "z_t = x_{t-1} constructed." << endl;

		//! Initialize cut constraints
		models[t].cuts = IloRangeArray(newEnv);

		//! create cplex algorithms
		models[t].cplex = IloCplex(models[t].mod);

		//! cplex parameters setting
		models[t].cplex.setParam(IloCplex::Param::MIP::Tolerances::MIPGap, MIPTOL);
		models[t].cplex.setParam(IloCplex::Param::RandomSeed, rand() % CPX_BIGINT);
		models[t].cplex.setParam(IloCplex::Param::TimeLimit, TIME_LIMIT_SM);
		models[t].cplex.setParam(IloCplex::Param::Simplex::Tolerances::Feasibility, 1e-8);
		models[t].cplex.setParam(IloCplex::Param::Simplex::Tolerances::Optimality, 1e-8);
		models[t].cplex.setParam(IloCplex::Param::Threads, MIPTHREAD);
		models[t].cplex.setOut(models[t].mod.getEnv().getNullStream());
		// models[t].cplex.setWarning(models[t].mod.getEnv().getNullStream());

		// write model to lp file
		// char fileName[100];
		// sprintf(fileName, "model_%d.lp", t);
		// models[t].cplex.exportModel(fileName);
	}
	return;
} // End of function buildModel


void getSamplePaths (SamplePath & sample, FormatData * fData_p)
{
	// cout << "Sampling forward paths from scenarios..." << endl;
	
	sample.clear();
	IloEnv * env = &(fData_p->dataEnv);

	for (int p = 0; p < fData_p->numFWsample; ++p)
	{			
		//! initialize each sample path with first stage data
		if ( fData_p->uncertaintySource[0] )
		{
			sample.x.add(IloNumArray2(*env));
			if ( fData_p->breakstage <= 0 )
				sample.x[p].add(fData_p->x[0]);
			else
				sample.x[p].add(fData_p->xBin[0]);
		}
		if ( fData_p->uncertaintySource[1] )
		{
			sample.y1.add(IloNumArray2(*env));
			sample.y1[p].add(fData_p->y1[0]);
		}
		if ( fData_p->uncertaintySource[2] )
		{
			sample.y2.add(IloNumArray2(*env));
			sample.y2[p].add(fData_p->y2[0]);
		}
		if ( fData_p->uncertaintySource[3] )
		{
			sample.A.add(IloNumArray3(*env));
			if ( fData_p->breakstage <= 0 )
				sample.A[p].add(fData_p->A);
			else
				sample.A[p].add(fData_p->ABin);
		}
		if ( fData_p->uncertaintySource[4] )
		{
			sample.B.add(IloNumArray3(*env));
			if ( fData_p->breakstage <= 0 )
				sample.B[p].add(fData_p->B);
			else
				sample.B[p].add(fData_p->BBin);
		}
		if ( fData_p->uncertaintySource[5] )
		{
			sample.W1.add(IloNumArray3(*env));
			sample.W1[p].add(fData_p->W1);
		}
		if ( fData_p->uncertaintySource[6] )
		{
			sample.W2.add(IloNumArray3(*env));
			sample.W2[p].add(fData_p->W2);
		}
		if ( fData_p->uncertaintySource[7] )
		{
			sample.b.add(IloNumArray2(*env));
			sample.b[p].add(fData_p->b[0]);
		}

		//! sample forward path for stage 1 to T-1
		for (int t = 1; t < fData_p->numStage; ++t)
		{
			int chosen;

			double pdf = 1.0/fData_p->numScen[t];  // probability of each outcome
			double U = genrand64_real1();  // generate a random number between 0 and 1
			chosen = int(U/pdf);  // the index of the selected outcome

			if ( fData_p->uncertaintySource[0] )
			{
				if ( t < fData_p->breakstage )
					sample.x[p].add(fData_p->scen.xBin[t][chosen]);
				else
					sample.x[p].add(fData_p->scen.x[t][chosen]);
			}
			if ( fData_p->uncertaintySource[1] )
				sample.y1[p].add(fData_p->scen.y1[t][chosen]);
			if ( fData_p->uncertaintySource[2] )
				sample.y2[p].add(fData_p->scen.y2[t][chosen]);
			if ( fData_p->uncertaintySource[3] )
			{
				if ( t < fData_p->breakstage )
					sample.A[p].add(fData_p->scen.ABin[t][chosen]);
				else
					sample.A[p].add(fData_p->scen.A[t][chosen]);
			}
			if ( fData_p->uncertaintySource[4] )
			{
				if ( t < fData_p->breakstage )
					sample.B[p].add(fData_p->scen.BBin[t][chosen]);
				else
					sample.B[p].add(fData_p->scen.B[t][chosen]);
			}
			if ( fData_p->uncertaintySource[5] )
				sample.W1[p].add(fData_p->scen.W1[t][chosen]);
			if ( fData_p->uncertaintySource[6] )
				sample.W1[p].add(fData_p->scen.W2[t][chosen]);
			if ( fData_p->uncertaintySource[7] )
				sample.b[p].add(fData_p->scen.b[t][chosen]);
		} // End of stage for-loop
	} // End of sample paths for loop
	return;
} // End of getSamplePaths

void forward (Model * models, FormatData * fData_p, int start, int end,
	const SamplePath sample, IloNumArray3 & fwdSoln, int numSelect, IloNumArray2 & ub)
/*!
	Forward step of SDDiP algorithm.
	User can specify start and end stage in the forward step,
	to accommodate the "bouncing" technique.
*/
{
	// cout << "Start the forward process..." << endl;

	int p, t, i;
	IloInt nconstr1 = models[0].constr1.getSize();
	IloInt sampleSize = fData_p->numFWsample;

	//! create an array to record the objective function value for each sample path
	IloNumArray sampleObj(fData_p->dataEnv, sampleSize);
	//! create a 3d array to record the candidate solutions for each sample path
	IloNumArray3 fwdSolnAll(fData_p->dataEnv, fData_p->numStage);

	for ( t = start; t < end; ++t)
	{
		fwdSolnAll[t] = IloNumArray2(fData_p->dataEnv);
		fwdSoln[t] = IloNumArray2(fData_p->dataEnv);
	}
	
	//! Find candidate solutions for each sample path
	for ( p = 0; p < sampleSize; ++p )
	{
		// cout << "================================" << endl;
		// if (sampleSize > 100)
		// {
		// 	cout << "Compute solution for sample path p =  " << p+1 << endl;
		// 	if ( p+1 != sampleSize )
		// 		cout << "\e[A";
		// }

		for ( t = start; t < end; ++t )
		{
			// cout << "Forward: current stage: " << t << endl;
			// cin.get();

			if ( t > 0 )
			{
				//! update stage problem parameter
				if ( fData_p->uncertaintySource[0] )
					models[t].obj.setLinearCoefs(models[t].x, sample.x[p][t]);
				if ( fData_p->uncertaintySource[1] )
					models[t].obj.setLinearCoefs(models[t].y1, sample.y1[p][t]);
				if ( fData_p->uncertaintySource[2] )
					models[t].obj.setLinearCoefs(models[t].y2, sample.y2[p][t]);

				for ( i = 0; i < nconstr1; ++i )
				{
					if ( fData_p->uncertaintySource[3] )
						models[t].constr1[i].setLinearCoefs(models[t].x, sample.A[p][t][i]);
					if ( fData_p->uncertaintySource[4] )
						models[t].constr1[i].setLinearCoefs(models[t].z, sample.B[p][t][i]);
					if ( fData_p->uncertaintySource[5] )
						models[t].constr1[i].setLinearCoefs(models[t].y1, sample.W1[p][t][i]);
					if ( fData_p->uncertaintySource[6] )
						models[t].constr1[i].setLinearCoefs(models[t].y2, sample.W2[p][t][i]);
					if ( fData_p->uncertaintySource[7] )
						models[t].constr1[i].setLB(sample.b[p][t][i]);
				}
				// cout << "Update previous stage solution." << endl;

				//! pass the latest candidate solution to next stage 
				if ( t > start )
				{
					if ( t == fData_p->breakstage )  // state space change, need to convert candidate solutions
					{
						IloNumArray convertedSol(fData_p->dataEnv, models[t].x.getSize());
						for ( i = 0; i < models[t].x.getSize(); ++i )
							convertedSol[i] = IloScalProd(fData_p->T[i], fwdSolnAll[t-1][p]);
						models[t].constr2.setBounds(convertedSol, convertedSol);
					}
					else
						models[t].constr2.setBounds(fwdSolnAll[t-1][p], fwdSolnAll[t-1][p]);
				}
			}
			// End of update problem models[t]

			// models[t].cplex.exportModel("FW.lp");
			// cout << "good here." << endl;

			// solve the current MIP model
			if ( ! models[t].cplex.solve() )
			{
				// cout << fwdSolnAll[t-1][p] << endl;
				// cout << models[t-1].mod << endl;
				char fileName[100];
				sprintf(fileName, "fw_model_%d.lp", t);
				models[t].cplex.exportModel(fileName);
				cout << "Sample " << p << ", Stage " << t << endl;
				cout << "Solution status is " << models[t].cplex.getStatus() << endl;
				throw ("Forward problem is infeasible...");
			}
			else if ( models[t].cplex.getStatus() != IloAlgorithm::Optimal )
			{
				char fileName[100];
				sprintf(fileName, "fw_model_%d.lp", t);
				models[t].cplex.exportModel(fileName);
				cout << "Sample " << p << ", Stage " << t << endl;
				cout << "Solution status is " << models[t].cplex.getStatus() << endl;
				throw ("Forward problem is not optimal...");
			}
			else
			{
				IloNumArray vals(fData_p->dataEnv);	
				models[t].cplex.getValues(vals, models[t].x);
				//! deal with numerical issues
				for ( i = 0; i < vals.getSize(); ++i )
				{
					if ( (vals[i] < 0) || (abs(vals[i]) < 9e-7))
						vals[i] = 0;
					if ( models[t].cplex.isMIP() )
						vals[i] = round(vals[i]*100)/100;
				}
				//! record candidate solution
				fwdSolnAll[t].add(vals);
				//! update objective function value for this sample path
				IloNum stageCost = models[t].cplex.getObjValue();
				if ( t < fData_p->numStage )
					stageCost -= models[t].cplex.getValue(models[t].theta);
				sampleObj[p] += stageCost;
			} // End of solving stage problem
		} // End of loop over stages
	} //End of loop over sample paths

	/*!
		Select fwdSoln based on user's preference:
		choose the forSoln that has the biggest objective function value
	*/
	IloIntArray idx(fData_p->dataEnv);
	getLargestK(sampleObj, numSelect, idx);
	for (i = 0; i < idx.getSize(); ++i)
		for (t = start; t < end; ++t)
			fwdSoln[t].add( fwdSolnAll[t][idx[i]] );

	//! calculate CI for upper bounds only if a full FW step is done
	if ( (start == 0) && (end == fData_p->numStage) )
	{
		double center = IloSum(sampleObj)/sampleSize,
			sumSq = 0.0,
			halfLength = 0.0,
			sd = 0.0;
		for ( p = 0; p < sampleSize; ++p)
			sumSq += pow( (sampleObj[p] - center), 2);

		if ( sampleSize > 1 )
		{
			sd = sqrt(sumSq / (sampleSize - 1));
			halfLength = ZALPHA * sd / sqrt(sampleSize);
		}
		ub[0].add(center - halfLength);
		ub[1].add(center);
		ub[2].add(center + halfLength);
		ub[3].add(IloMin(sampleObj));
		ub[4].add(IloMax(sampleObj));
		ub[5].add(sd);

		ofstream output ("sampleObj.csv", ios::out | ios::app);
		if ( output.is_open() && (sampleObj.getSize() > 100) )
		{
			for (i = 0; i < sampleObj.getSize(); ++i)
				output << sampleObj[i] << ",";
			output << endl;
		}
	}
	sampleObj.end();
	fwdSolnAll.end();
}  //! End of forward pass to generate candidate solutions

bool sortByValue(const Pair &lhs, const Pair &rhs) { return lhs.value > rhs.value; }

void getLargestK(IloNumArray arr, int k, IloIntArray & idx)
{
	unsigned size = arr.getSize();
	vector<Pair> target(size);
	for (vector<Pair>::size_type i = 0; i != size; ++i)
    {
        target[i].value = arr[i];
        target[i].index = i;
    }

    sort(target.begin(), target.end(), sortByValue);

    for (int i = 0; i < k; ++i)
    	idx.add(target[i].index);

    return;
}

void backward (Model * models, FormatData * fData_p, int start, int end, 
	const IloNumArray3 fwdSoln, IloNumArray & lb, const CutSwitch cut)
{
	// cout << "Start the backward process...." << endl;

	int t, k, j;
	IloNumArray2 x, xBin, y1, y2, b;
	IloNumArray3 A, ABin, B, BBin, W1, W2;
	IloInt dimY1 = fData_p->y1[0].getSize();
	IloInt dimY2 = fData_p->y2[0].getSize();
	IloInt breakstage = fData_p->breakstage;
	
	for ( t = end - 1; t > start; --t )
	{
		// cout << "Backward: current stage t =  " << t << " " << endl;
		int nconstr1 = models[t].constr1.getSize();
		int nconstr2 = models[t].constr2.getSize();
		int nscen = fData_p->numScen[t];
		int nsoln = fwdSoln[t-1].getSize();
		
		IloEnv env = models[t].mod.getEnv();

		for ( j = 0; j < nsoln; ++j )
		{
			// cout << "evaluate solution: " << j << endl;
			// cout << "*************************************" << endl;
			IloNumArray MIPobj(env, nscen);
			IloNumArray SBobj(env, nscen);
			IloNumArray LGobj(env, nscen);
			IloNumArray LPobj(env, nscen);			
			IloNumArray LPdualAvg(env, nconstr2);
			IloNumArray LGdualAvg(env, nconstr2);
			IloNumArray convertedSol(env);

			// update constraints z = x, check if need to convert variable space
			if ( t == fData_p->breakstage ) // model change, need to convert candidate solutions
			{
				for (int i = 0; i < models[t].x.getSize(); ++i )
					convertedSol.add(IloScalProd(fData_p->T[i], fwdSoln[t-1][j]));
				models[t].constr2.setBounds(convertedSol, convertedSol);
			}
			else
				models[t].constr2.setBounds(fwdSoln[t-1][j], fwdSoln[t-1][j]);

			IloModel modcp = IloGetClone(env, models[t].mod);

			// Start parallelization
			int nProcessors = omp_get_max_threads();

			// set the number of threads
			// omp_set_num_threads(nProcessors)
			omp_set_num_threads(1);

			#pragma omp parallel
			{
				double MIPobj_private[nscen];
				double SBobj_private[nscen];
				double LGobj_private[nscen];
				double LPobj_private[nscen];
				double LGdualAvg_private[nconstr2];
				double LPdualAvg_private[nconstr2];
				int i, m;

				for (i = 0; i < nscen; ++i)
				{
					MIPobj_private[i] = 0;
					SBobj_private[i] = 0;
					LGobj_private[i] = 0;
					LPobj_private[i] = 0;
				}
				for (m = 0; m < nconstr2; ++m)
				{
					LGdualAvg_private[m] = 0;
					LPdualAvg_private[m] = 0;
				}
				
				double LPobjScen, SBobjScen, LGobjScen, MIPobjScen;

				#pragma omp for private(k,i,m,LPobjScen,SBobjScen,LGobjScen,MIPobjScen)\
				schedule(static)
				for ( k = 0; k < nscen; ++k )
				{
					int tid = omp_get_thread_num();
					// cout << "thread " << tid << " is assigned with scenario " << k << endl; 	
					// cout << "Solving scenario " << k << " in stage " << t << endl;

					// ========================================
					// copy model[t] for scenario k
					// ========================================					
					Model scenModel;
					scenModel.env = IloEnv();		
					scenModel.mod = IloModel(scenModel.env);
					scenModel.mod = IloGetClone(scenModel.env, modcp);
					scenModel.cplex = IloCplex(scenModel.env);
					scenModel.cplex.setOut(scenModel.env.getNullStream());
					scenModel.cplex.extract(scenModel.mod);
					scenModel.cplex.setParam(IloCplex::Param::Threads, MIPTHREAD);
					
					scenModel.x = IloNumVarArray(scenModel.env);
					scenModel.y1 = IloNumVarArray(scenModel.env);
					scenModel.y2 = IloNumVarArray(scenModel.env);
					scenModel.z = IloNumVarArray(scenModel.env);
					scenModel.constr1 = IloRangeArray(scenModel.env);
					scenModel.constr2 = IloRangeArray(scenModel.env);
					scenModel.cuts = IloRangeArray(scenModel.env);

					IloNumArray dual(scenModel.env, nconstr2);
					IloRange range = IloRange();
					// cout << "thread " << tid << " initialized model " << k << endl; 

					for (IloModel::Iterator it(scenModel.mod); it.ok(); ++it)
					{
						if ((*it).isVariable())
						{
							string name((*it).getName());
							if ( name.find("x") != string::npos)
	   							scenModel.x.add((*it).asVariable());
	   					
		   					if ( name.find("y1") != string::npos)
		   						scenModel.y1.add((*it).asVariable());

		   					if ( name.find("y2") != string::npos)
		   						scenModel.y2.add((*it).asVariable());

		   					if ( name.find("z") != string::npos)
		   						scenModel.z.add((*it).asVariable());
		   				}
	   					else if ((*it).isConstraint())
	   					{
	   						if (scenModel.constr1.getSize() < nconstr1)
	   						{	
	   							IloRangeI * impl = dynamic_cast<IloRangeI *>((*it).asConstraint().getImpl());
	   							range = IloRange(impl);
	   							scenModel.constr1.add(range);
	   						}
	   						else if (scenModel.constr2.getSize() < nconstr2)
	   						{
	   							IloRangeI * impl = dynamic_cast<IloRangeI *>((*it).asConstraint().getImpl());
	   							range = IloRange(impl);
	   							scenModel.constr2.add(range);
	   						}
	   						else
	   						{
	   							IloRangeI * impl = dynamic_cast<IloRangeI *>((*it).asConstraint().getImpl());
	   							range = IloRange(impl);
	   							scenModel.cuts.add(range);
	   						}
	   					}
	   					else if ((*it).isObjective())
	   						scenModel.obj = (*it).asObjective();
	   				}
	   				// end of model copy
	   				// cout << "thread " << tid << " finished copy model." << endl; 
	   
	   				// =============================================================
					// update objective coefficients and constraints if necessary
					// =============================================================
					#pragma omp critical
					{
						if ( fData_p->uncertaintySource[0] )
						{
							if ( t >= fData_p->breakstage )
								scenModel.obj.setLinearCoefs(scenModel.x, fData_p->scen.x[t][k]);
							else
								scenModel.obj.setLinearCoefs(scenModel.x, fData_p->scen.xBin[t][k]);
						}
						if ( fData_p->uncertaintySource[1] )
							scenModel.obj.setLinearCoefs(scenModel.y1, fData_p->scen.y1[t][k]);
						if ( fData_p->uncertaintySource[2] )
							scenModel.obj.setLinearCoefs(scenModel.y2, fData_p->scen.y2[t][k]);

						for ( i = 0; i < nconstr1; ++i )
						{
							if ( fData_p->uncertaintySource[3] )
							{
								if ( t >= fData_p->breakstage )
									scenModel.constr1[i].setLinearCoefs(scenModel.x, fData_p->scen.A[t][k][i]);
								else
									scenModel.constr1[i].setLinearCoefs(scenModel.x, fData_p->scen.ABin[t][k][i]); 
							}
							if ( fData_p->uncertaintySource[4] )
							{
								if ( t >= fData_p->breakstage )
									scenModel.constr1[i].setLinearCoefs(scenModel.z, fData_p->scen.B[t][k][i]);
								else
									scenModel.constr1[i].setLinearCoefs(scenModel.z, fData_p->scen.BBin[t][k][i]); 
							}
							if ( fData_p->uncertaintySource[5] )
								scenModel.constr1[i].setLinearCoefs(scenModel.y1, fData_p->scen.W1[t][k][i]);
							if ( fData_p->uncertaintySource[6] )
								scenModel.constr1[i].setLinearCoefs(scenModel.y2, fData_p->scen.W2[t][k][i]);
							if ( fData_p->uncertaintySource[7] )
								scenModel.constr1[i].setLB(fData_p->scen.b[t][k][i]);
						}
					}
					// end of updates
					// cout << "thread " << tid << " finish update model" << endl;

					// ============================================================================
					// The following section computes necessary cut coefficients from each scenario
					// ============================================================================
					// cout << "Computing cut coefficients... " << endl;
					
					// cout << LPobj_private << endl;

					if ( t >= breakstage ) // only add Benders' cuts
					{
						getBcutCoef(scenModel, LPobjScen, dual);
						LPobj_private[k] += LPobjScen;
						
						for (i = 0; i < nconstr2; ++i)
							LPdualAvg_private[i] += dual[i];

						if ( (t == fData_p->breakstage) && (cut.SB) )
						{
							getSBcutCoef(scenModel, dual, SBobjScen);
							SBobj_private[k] += SBobjScen;
						}
					}
					else // models before breakstage, discretized formulation
					{
						if ( cut.LG + cut.I )
						{
							// cout << "compute Int coeff" << endl;
							getIcutCoef(scenModel, MIPobjScen);
							MIPobj_private[k] += MIPobjScen;
						}

						if ( cut.B + cut.SB + cut.LG )
						{
							getBcutCoef(scenModel, LPobjScen, dual);
							LPobj_private[k] += LPobjScen;
							for (i = 0; i < nconstr2; ++i)
								LPdualAvg_private[i] += dual[i];

							if ( cut.SB )
							{
								// cout << "compute SB coeff" << endl;
								getSBcutCoef(scenModel, dual, SBobjScen);
								SBobj_private[k] += SBobjScen;
							}

							if ( cut.LG )
							{
								IloNumArray pi_opt(scenModel.env);
								// cout << "solve LG" << endl;
								// cout << "LP dual: " << dual << endl;
								if ( cut.LEVEL )
									pi_opt = getLGcutCoefLevel(scenModel, dual, fwdSoln[t-1][j], LGobjScen);
								else
									pi_opt = getLGcutCoefSubgrad(scenModel, dual, fwdSoln[t-1][j], LGobjScen);
								// cout << pi_opt << endl;
								
								LGobj_private[k] += LGobjScen;
								for (i = 0; i < nconstr2; ++i)
									LGdualAvg_private[i] += pi_opt[i];
								// #pragma critical (LG)
								// {
								// 	for (i = 0; i < nconstr2; ++i)
								// 		LGdualAvg[i] += dual[i] / nscen;
								// }
								// cout << "optimal LG dual:" << pi_opt << endl;
								// cout << "LG solved." << endl;
								pi_opt.end();
							}
						}
					}
					
					scenModel.env.end();

					// cout << "thread " << tid << " finished computation scenario " << k << endl;

				}// End of loop over all scenarios in stage t

				// cout << LPobj_private << endl;

				#pragma omp critical (collectInfo)
				{
					if (cut.B + cut.SB)
					{
						for ( i = 0; i < nconstr2; ++i )
							LPdualAvg[i] += LPdualAvg_private[i]/nscen;
						if (cut.B)
						{	
							for ( m = 0; m < nscen; ++m )
								LPobj[m] += LPobj_private[m];
						}
						if (cut.SB)
						{
							for ( m = 0; m < nscen; ++m )
								SBobj[m] += SBobj_private[m];
						}
					}
					if (cut.LG)
					{
						for ( m = 0; m < nscen; ++m )
							LGobj[m] += LGobj_private[m];
						for ( i = 0; i < nconstr2; ++i )
							LGdualAvg[i] += LGdualAvg_private[i]/nscen;
					}
					if (cut.I)
					{
						for ( m = 0; m < nscen; ++m )
							MIPobj[m] += MIPobj_private[m];
					}
				}
				// cout << "final LPdualAvg: " << LPdualAvg << endl;

				// for ( int m = 0; m < nscen; ++m )
				// {
				// 	for (int n = 0; n < nconstr2; ++n)
				// 		cout << LPdualAvg_private[m][n] << " ";
				// 	cout << endl;
				// }				
			}
			// cout << "final LPobj: " << LPobj << endl;
			// cout << "+++++++++++++++++++++++++++++++++++++++++" << endl;
			// cout << "final LPdualAvg: " << LPdualAvg << endl;
			// for (int i = 0; i < nconstr2; ++i )
			// {
			// 	if (cut.B + cut.SB)
			// 		LPdualAvg[i] = LPdualAvg[i] / nscen;
			// 	if (cut.LG)
			// 		LGdualAvg[i] = LGdualAvg[i] / nscen;
			// }
			// cout << LPdualAvg << endl;

			

			// ====================================================
			// The following section add cuts to models[t-1]
			// ====================================================
			// cout << "Adding cuts to previous stage model... " << endl;

			if ( t > breakstage )
				addBcut(models[t-1], LPobj, LPdualAvg, fwdSoln[t-1][j]);
			else if ( t == breakstage )
			{
				IloNumArray LPdualAvgBin(fData_p->dataEnv);
				for (int i = 0; i < models[t-1].x.getSize(); ++i)
					LPdualAvgBin.add(IloScalProd(LPdualAvg, fData_p->TT[i]));
				if ( cut.SB )
					addSBcut(models[t-1], SBobj, LPdualAvgBin, fwdSoln[t-1][j]);
				else
					addBcut(models[t-1], LPobj, LPdualAvgBin, fwdSoln[t-1][j]);
			}
			else
			{
				if ( cut.B )
					addBcut(models[t-1], LPobj, LPdualAvg, fwdSoln[t-1][j]);

				if ( cut.SB )
					addSBcut(models[t-1], SBobj, LPdualAvg, fwdSoln[t-1][j]);
				 
				if ( cut.LG )
					addLGcut(models[t-1], LGobj, LGdualAvg, fwdSoln[t-1][j]);
	
				if ( cut.I )
					addIcut(models[t-1], MIPobj, fwdSoln[t-1][j], fData_p->thetaBound[0][t-1]);
			}

			// free momery
			MIPobj.end();
			SBobj.end();
			LGobj.end();
			LPobj.end();
			LPdualAvg.end();
			LGdualAvg.end();
			convertedSol.end();
		} // End of loop over unique candidate solutions
	} // End of loop over stages

	// cout << "================================" << endl;
	// cout << "Solving the first stage problem." << endl;

	// solve models[0] as a MIP
	if ( (start == 0) && (end == fData_p->numStage) )
	{
		if ( ! models[0].cplex.solve() ) // infeasible
		{
			cout << "First stage problem status: " << models[0].cplex.getStatus() << endl;
			throw ("First-stage problem infeasible...");
		}
		else if ( models[0].cplex.getStatus() != IloAlgorithm::Optimal )
		{
			cout << "First-stage problem status: " << models[0].cplex.getStatus() << endl;
			throw ("First-stage problem not optimal...");
		}
		else
		{
			lb.add(models[0].cplex.getObjValue());
		}
	}
	
	return;
}

void getIcutCoef(Model & model, double & MIPobjScen)
{
	IloAlgorithm::Status solStatus;
	if ( model.cplex.solve() ) // feasible
	{	
		solStatus = model.cplex.getStatus();
		if ( solStatus == IloAlgorithm::Optimal )
			MIPobjScen = model.cplex.getObjValue();
		else // not optimal
		{
			cout << "Current MIP is: " << solStatus << endl;
			throw ("MIP not optimal...");
		}
	}
	else // infeasible
	{
		char fileName[100];
		sprintf(fileName, "inf_MIP.lp");
		model.cplex.exportModel(fileName);
		cout << "Current MIP is " << model.cplex.getStatus() << endl;
		throw ("MIP infeasible...");
	}
	return;
}

void getBcutCoef(Model & model, double & LPobjScen, IloNumArray & dual)
{
	IloEnv env = model.mod.getEnv();
	IloModel LP(env);
	LP.add(model.mod);
	IloCplex cplexLP(LP);
	if (cplexLP.isMIP())
	{
		LP.add(IloConversion(env, model.x, ILOFLOAT));
		LP.add(IloConversion(env, model.y1, ILOFLOAT));
	}
	cplexLP.setOut(env.getNullStream());
	// cplexLP.setParam(IloCplex::Param::Threads, 2);

	IloAlgorithm::Status solStatus;

	if (cplexLP.solve())
	{
		solStatus = cplexLP.getStatus();
		if ( solStatus == IloAlgorithm::Optimal )
		{
			LPobjScen = cplexLP.getObjValue();
			// cout << "Objective Value: " << LPobjScen << endl;
			cplexLP.getDuals(dual, model.constr2);
		}
		else // not optimal
		{
			cout << "Current LP is " << solStatus << endl;
			throw ("LP not optimal...");
		}
	}
	else // infeasible
	{
		char fileName[100];
		sprintf(fileName, "inf_LP.lp");
		cplexLP.exportModel(fileName);
		cout << "Current LP is " << cplexLP.getStatus() << endl;
		throw ("LP infeasible...");
	}

	LP.end();
	cplexLP.end();
	return;
}

void getSBcutCoef(Model & model, const IloNumArray LPdual, double & SBobjScen)
{
	IloEnv env = model.mod.getEnv();
	IloModel LG(env);
	LG.add(model.constr1);
	LG.add(model.cuts);
	// LG.add(model.mod);
	// LG.remove(model.constr2);
	// LG.remove(model.obj);
	// create new objective function with additional term -\pi_LP' z
	IloObjective objLG = IloMinimize(env);
	IloExpr expr = model.obj.getExpr() - IloScalProd(LPdual, model.z);
	objLG.setExpr(expr);
	LG.add(objLG);
	IloCplex cplexLG(LG);
	cplexLG.setParam(IloCplex::Param::MIP::Tolerances::MIPGap, MIPTOL);
	cplexLG.setParam(IloCplex::Param::Threads, MIPTHREAD);
	cplexLG.setOut(env.getNullStream());

	IloAlgorithm::Status solStatus;
	if ( cplexLG.solve() )
	{
		solStatus = cplexLG.getStatus();
		if ( solStatus == IloAlgorithm::Optimal )
			SBobjScen = cplexLG.getObjValue();
		else // not optimal
		{
			cout << "Current LG is " << solStatus << endl;
			throw ("LG not optimal...");
		}
	}
	else // infeasible
	{
		char fileName[100];
		sprintf(fileName, "inf_LG.lp");
		cplexLG.exportModel(fileName);
		cout << "Current LG is " << cplexLG.getStatus() << endl;
		throw ("LG infeasible...");
	}

	LG.end();
	cplexLG.end();
	return;
}

IloNumArray getLGcutCoefLevel(Model & model, const IloNumArray LPdual,
	const IloNumArray fwdSoln, double & LGobjScen)
{
	IloEnv env = model.mod.getEnv();
	IloModel LG(env);
	LG.add(model.z);
	LG.add(model.constr1);
	LG.add(model.cuts);
	IloObjective objLG = IloMinimize(env);
	IloExpr expr = model.obj.getExpr() - IloScalProd(LPdual, model.z);
	objLG.setExpr(expr);
	LG.add(objLG);
	IloCplex cplexLG(LG);
	cplexLG.setParam(IloCplex::Param::MIP::Tolerances::MIPGap, MIPTOL);
	cplexLG.setParam(IloCplex::Param::Threads, MIPTHREAD);
	cplexLG.setOut(model.mod.getEnv().getNullStream());

	IloAlgorithm::Status solStatus;
	// IloNumArray pi_opt(model.mod.getEnv());

	IloNumArray2 pi(env);
	IloNumArray2 grad(env);
	IloNumArray z(env);
	IloNumArray vpi(env);
	IloNumArray lb(env);
	IloNumArray pi_new(env, LPdual.getSize());
	for ( int i = 0; i < pi_new.getSize(); ++i )
		pi_new[i] = LPdual[i];
	// z.add(IloNumArray(model.mod.getEnv()));
	
	double bestUB = model.cplex.getObjValue(),
		   bestLB = -numeric_limits<double>::infinity(),
		   delta = numeric_limits<double>::infinity(),
		   gnorm = numeric_limits<double>::infinity(),
		   lambda = MIPTOL, level;
	int iter = 0, ibest;
	while ( delta / bestUB > MIPTOL )// && (gnorm > EPSILON) )
	{
		if ( ! cplexLG.solve() )
		{
			char fileName[100];
			sprintf(fileName, "inf_LG.lp");
			cplexLG.exportModel(fileName);
			cout << "Current LG is " << cplexLG.getStatus() << endl;
			throw ("LG infeasible...");
		}
		else if ( cplexLG.getStatus() != IloAlgorithm::Optimal )
		{
			cout << "Current LG is " << cplexLG.getStatus() << endl;
			throw ("LG not optimal...");
		}
		else
		{
			pi.add(pi_new);
			// z.add(IloNumArray(model.mod.getEnv()));
			cplexLG.getValues(z, model.z);
			vpi.add(cplexLG.getObjValue());
			// cout << pi_new.getSize() << " " << fwdSoln[t-1][j].getSize() << endl;
			lb.add(vpi[iter] + IloScalProd(pi_new, fwdSoln));
			grad.add(IloNumArray(model.mod.getEnv()));
			for ( int i = 0; i < z.getSize(); ++i )
				grad[iter].add(fwdSoln[i] - z[i]);
			gnorm = sqrt(IloScalProd(grad[iter], grad[iter]));

			if ( lb[iter] > bestLB )
				ibest = iter;
			bestLB = lb[ibest];
			delta = bestUB - bestLB;
			level = bestUB - lambda * delta;
			// cout << " current Delta: " << delta << endl;
			
			if ( delta / bestUB > MIPTOL ) //&& (gnorm > EPSILON) )
			{
				pi_new = projection(pi[ibest], lb, grad, pi, level);
				expr = model.obj.getExpr() - IloScalProd(pi_new, model.z);
				objLG.setExpr(expr);
				// cout << pi_new << endl;		
				iter += 1;
				// cout << iter << endl;
			}
		}
	}
	LGobjScen = lb[iter];
	pi.end();
	grad.end();
	z.end();
	vpi.end();
	lb.end();
	LG.end();
	cplexLG.end();

	return pi_new;
}

IloNumArray projection (const IloNumArray xstar, IloNumArray fval,
	IloNumArray2 grad, IloNumArray2 x0, double level)
{
	IloEnv env;
	IloInt dim = xstar.getSize();

	IloModel m(env);

	IloNumVarArray x(env, dim, -IloInfinity, IloInfinity, ILOFLOAT);
	m.add(x);

	IloObjective obj(env);
	IloExpr objExpr = IloScalProd(x, x) - 2 * IloScalProd(x, xstar);
	obj.setExpr(objExpr);
	obj.setSense(IloObjective::Minimize);
	m.add(obj);

	IloRangeArray constr(env);
	int i, nrow = x0.getSize();
	for (i = 0; i < nrow; ++i)
	{
		IloExpr expr(env);
		expr = fval[i] + IloScalProd(grad[i], x) - IloScalProd(grad[i], x0[i]);
		// cout << fval[i] << " " << IloScalProd(grad[i], x0[i]) << " " << level << endl;
		constr.add(expr >= level);
	}
	m.add(constr);

	IloCplex cplex(m);
	cplex.setOut(env.getNullStream());

	IloNumArray opt(xstar.getEnv());
	if ( cplex.solve() )
	{	
		cplex.getValues(opt, x);
		// cplex.exportModel("projection.lp");
		// cin.get();
	}
	else
	{
		cout << " No solution exists for the projection problem. Please check!! " << endl;
		cout << " current level: " << level << endl;
		cplex.exportModel("projection.lp");
	}
	env.end();

	return opt;
}

IloNumArray getLGcutCoefSubgrad(Model & model, const IloNumArray LPdual,
	const IloNumArray fwdSoln, double & LGobjScen)
{
	IloEnv env = model.mod.getEnv();
	IloModel LG(env);
	LG.add(model.z);
	LG.add(model.constr1);
	LG.add(model.cuts);
	IloObjective objLG = IloMinimize(env);
	IloExpr expr = model.obj.getExpr() - IloScalProd(LPdual, model.z);
	objLG.setExpr(expr);
	LG.add(objLG);
	IloCplex cplexLG(LG);
	cplexLG.setParam(IloCplex::Param::MIP::Tolerances::MIPGap, MIPTOL);
	cplexLG.setParam(IloCplex::Param::Threads, MIPTHREAD);
	cplexLG.setOut(env.getNullStream());

	IloAlgorithm::Status solStatus;
	int iter = 0, maxIter = 2000;
	double bestUB = model.cplex.getObjValue();

	IloNumArray pi(env, LPdual.getSize());
	for ( int i = 0; i < pi.getSize(); ++i )
		pi[i] = LPdual[i];
	IloNumArray grad(env);
	double stepSize, vpi;
	double gap = numeric_limits<double>::infinity();
	double gnorm = numeric_limits<double>::infinity();
	
	while ( (iter < maxIter) && (gap > MIPTOL) && (gnorm > EPSILON) )
	{
		if ( ! cplexLG.solve() )
		{
			char fileName[100];
			sprintf(fileName, "inf_LG.lp");
			cplexLG.exportModel(fileName);
			cout << "Current LG is " << cplexLG.getStatus() << endl;
			throw ("LG infeasible...");
		}
		else if ( cplexLG.getStatus() != IloAlgorithm::Optimal )
		{
			cout << "Current LG is " << cplexLG.getStatus() << endl;
			throw ("LG not optimal...");
		}
		else
		{
			stepSize = 20 / sqrt(iter+1);
			vpi = cplexLG.getObjValue();
			double delta = bestUB - vpi - IloScalProd(pi, fwdSoln);
			// cout << "current delta: " << delta << endl;
			gap = delta / bestUB;
			cplexLG.getValues(grad, model.z);
			for (int i = 0; i < grad.getSize(); ++i )
				grad[i] -= fwdSoln[i];
			gnorm = sqrt(IloScalProd(grad, grad));
			if ( (gap > MIPTOL) && (gnorm > EPSILON) )
			{
				for (int i = 0; i < grad.getSize(); ++i )
				{
					pi[i] -= stepSize / (gnorm + EPSILON) * grad[i];
					objLG.setLinearCoef(model.z[i], -pi[i]);
				}
				iter += 1;
			}
		}
	}
	LGobjScen = vpi + IloScalProd(pi, fwdSoln);

	LG.end();
	cplexLG.end();
	grad.end();

	return pi;
}

void addIcut(Model & model, const IloNumArray MIPobj,
	const IloNumArray fwdSoln, const IloNum lb)
{
	// cout << "adding integer cuts..." << endl;
	IloExpr expr = model.theta;
	IloNum mipObjAvg = IloSum(MIPobj)/MIPobj.getSize();
	for (int i = 0; i < fwdSoln.getSize(); ++i )
	{
		if ( abs( fwdSoln[i] ) < EPSILON )
			expr += (mipObjAvg - lb) * model.x[i];
		else
			expr -= (mipObjAvg - lb) * model.x[i];
	}
	IloNum rhs = mipObjAvg - (mipObjAvg - lb) * IloSum(fwdSoln);
	model.cuts.add(expr >= rhs);
	model.mod.add(expr >= rhs);
	expr.end();
}

void addBcut(Model & model, const IloNumArray LPobj,
	const IloNumArray LPdualAvg, const IloNumArray fwdSoln)
{
	// cout << "adding Benders' cuts..." << endl;
	IloExpr expr = model.theta - IloScalProd(LPdualAvg, model.x);
	IloNum rhs = IloSum(LPobj) / LPobj.getSize() - IloScalProd(LPdualAvg, fwdSoln);
	model.cuts.add(expr >= rhs);
	model.mod.add(expr >= rhs);
	expr.end();
}

void addSBcut(Model & model, const IloNumArray SBobj,
	const IloNumArray LPdualAvg, const IloNumArray fwdSoln)
{
	// cout << "adding Strengthened Benders' cuts..." << endl;
	IloExpr expr = model.theta - IloScalProd(LPdualAvg, model.x);
	IloNum rhs = IloSum(SBobj) / SBobj.getSize();
	model.cuts.add(expr >= rhs);
	model.mod.add(expr >= rhs);
	expr.end();
}

void addLGcut(Model & model, const IloNumArray LGobj,
	const IloNumArray LGdualAvg, const IloNumArray fwdSoln)
{
	// cout << "adding Lagrangian cuts..." << endl;
	IloExpr expr = model.theta - IloScalProd(LGdualAvg, model.x);
	IloNum rhs = IloSum(LGobj) / LGobj.getSize() - IloScalProd(LGdualAvg, fwdSoln);
	model.cuts.add(expr >= rhs);
	model.mod.add(expr >= rhs);
	expr.end();
}

// double avg ( vector<float> & v )
// {
//     double sumV = 0.0;
//     for ( unsigned i=0; i < v.size(); i++)
//         sumV += v[i];
            
//     return sumV / v.size();
// }
// // End of mean funtion

// double std_dev ( vector<float> & v )
// {
//         double sumSquare = 0.0;
//         double var =0.0;
        
//         double mean = avg(v);
        
//         for ( unsigned j = 0; j < v.size(); ++j )
//             sumSquare += pow((v[j] - mean),2);
       
//         return var = sqrt(sumSquare / v.size());
// }
// End of standard deviation funtion

void usage (char *progname)
{
	cerr << "Usage:  " << progname << " arg1 arg2 arg3 arg4 arg5 [arg6]" << endl;
	cerr << "At least 6 parameters must be specified." << endl;
	cerr << "arg1: 0 -- turn off Benders' cuts;" << endl;
	cerr << "      1 -- turn on Benders' cuts." << endl;
	cerr << "arg2: 0 -- turn off Strengthened Benders' cuts;" << endl;
	cerr << "      1 -- turn on Strengthened Benders' cuts." << endl;
	cerr << "arg3: 0 -- turn off Lagragian cuts;" << endl;
	cerr << "      1 -- turn on Lagragian cuts." << endl;
	cerr << "arg4: 0 -- use subgradient method;" << endl;
	cerr << "      1 -- use Level method." << endl;
	cerr << "arg5: 0 -- turn off Integer L-shaped cuts;" << endl;
	cerr << "      1 -- turn on Integer L-shaped cuts." << endl;
	cerr << "arg6: [optional] used as the seed of the random number generator." << endl;
	cerr << "      If not provided, system will generate one automatically." << endl;
} // END usage

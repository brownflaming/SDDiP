#include <ilcplex/ilocplex.h>
#include <cmath>
#include <ctime>
#include <string>
#include <sstream>
#include <fstream>
#include <vector>
#include <chrono>
#include <algorithm> // std::max
#include "global.h"
#include "functions.h"
#include "mt64.h"


ILOSTLBEGIN

int main (int argc, char *argv[])
{
	try 
	{
		if ( argc != 6 && argc != 7 )
		{
			usage (argv[0]);
			throw (-1);
		}

		CutSwitch cut;

		unsigned long long seed;
		if ( argc == 7 )
			seed = atoi(argv[6]);
		else
			seed = chrono::system_clock::now().time_since_epoch().count();
		init_genrand64(seed);

		// create a new formated data structure and read data from file
		FormatData fData;
		FormatData * fData_p = & fData;
		readData (fData_p);

		cout << "All data has been read into fData." << endl;
		cout << "==================================================" << endl;
		cout << "Number of stage: " << fData.numStage << endl;
		cout << "Breakstage: " << fData.breakstage << endl;
		cout << "Number of scenarios at each stage: " << fData.numScen << endl;
		cout << "Total number of scenarios: " << fData.totalScen << endl;
		cout << "Number of FW samples: " << fData.numFWsample << endl;
		cout << "==================================================" << endl;

		//! construct models from fData
		Model * models = new Model[fData.numStage];
		buildModel(models, fData_p);
		cout << "Model construction completed..." << endl;

		//! start siddp method
		double runtime, cputime, simulation_time, fwdtime, bwdtime;
		double wall0, wall1, cpu0, cpu1;
		wall0 = get_wall_time();
		cpu0 = get_cpu_time();

		cout << "Starting SDDP procedure ... " << endl;

		/*! Specify if "bouncing" option is used
			If yes, segment indicates how to decompose planning horizon,
		*/
		bool bouncing = 0;
		int segment = 4;
		int repeat[4] = {2,2,2,2};

		//! initialization
		IloNumArray lb(fData.dataEnv); // double array to record lowerbound
		IloNumArray2 ub(fData.dataEnv, 6); // double array to record ub
		ub[0] = IloNumArray(fData.dataEnv); // CI left
		ub[1] = IloNumArray(fData.dataEnv); // mean
		ub[2] = IloNumArray(fData.dataEnv); // CI right
		ub[3] = IloNumArray(fData.dataEnv); // minimum
		ub[4] = IloNumArray(fData.dataEnv); // maximum
		ub[5] = IloNumArray(fData.dataEnv); // std_dev

		SamplePath sample;
		sample.x = IloNumArray3(fData.dataEnv);
		sample.y1 = IloNumArray3(fData.dataEnv);
		sample.y2 = IloNumArray3(fData.dataEnv);
		sample.A = IloNumArray4(fData.dataEnv);
		sample.W1 = IloNumArray4(fData.dataEnv);
		sample.W2 = IloNumArray4(fData.dataEnv);
		sample.B = IloNumArray4(fData.dataEnv);
		sample.b = IloNumArray3(fData.dataEnv);

		IloNumArray3 fwdSoln;
		
		//! start with SDDP by ignoring all integrality constraints
		models[t].mod.add(models[t].xLP);
		models[t].mod.add(models[t].yLP);
		//! Set cut generation configuration (Benders only)
		cut.B = 1; cut.SB = 0; cut.LG = 0; cut.LEVEL = 0;; cut.I = 0;

		double sddp_gap, simulation_gap;
		int iteration = 0; // iteration counter

		// forward step configuration
		fData.numFWsample = 3; // number of forwad sample paths
		int numSelect = 1; // number of fwdSoln chosen from all candidates
		int initSampleSize = fData.numFWsample;
		int numBranch = fData.numScen[1];

		// start the loop until some stopping creterior is satisfied
		do
		{
			cout << "================================" << endl;
			cout << "Iteration: " << iteration + 1 << endl;

			// if reaches MINITER, reinstate integrality cosntraints
			if ( iteration + 1 == MINITER )
			{	
				for ( int t = 0; t <= fData.breakstage; ++t )
				{
					models[t].xLP.end();
					models[t].yLP.end();
				}

				cut.B = atoi(argv[1]);
				cut.SB = atoi(argv[2]);
				cut.LG = atoi(argv[3]);
				cut.LEVEL = atoi(argv[4]);;
				cut.I = atoi(argv[5]);
				cout << "Change back to user setting." << endl;
			}

			getSamplePaths(sample, fData_p);
			fwdSoln = IloNumArray3(fData.dataEnv, fData.numStage);
			forward(models, fData_p, 0, fData.numStage, sample, fwdSoln, numSelect, ub);

			// cout << "Forward pass completed." << endl;

			if ( bouncing )
			{
				int first, last;
				for ( int m = 0; m < segment; m++ )
				{
					first = max(int(fData_p->numStage/segment * (segment-m-1)-1), 0);
					last = fData_p->numStage/segment * (segment-m);
					for ( int k = 0; k < repeat[segment-m-1]; ++k )
					{
						backward(models, fData_p, first, last, fwdSoln, lb, cut);
						getSamplePaths(sample, fData_p);
						if ( k < repeat[segment-m-1] - 1 )
							forward(models, fData_p, first, last, sample, fwdSoln, numSelect, ub);
					}
				}
			}
			
			backward(models, fData_p, 0, fData.numStage, fwdSoln, lb, cut);

			// cout << "Backward pass completed." << endl;

			cout << "Lower bound: " << lb[iteration] << "  " <<
				"CI for upper bound mean: [" << ub[0][iteration]  << ", " << 
				ub[2][iteration] << "]" << endl;

			iteration += 1;

		} while ( (iteration < MAXITER) && (runtime < TIME_LIMIT) );

		cout << "================================" << endl;
		cout << "Increasing sample size." << endl;
		fData.numFWsample = LARGE_SAMPLE;
		getSamplePaths(sample, fData_p);
		fwdSoln = IloNumArray3(fData.dataEnv, fData.numStage);
		forward(models, fData_p, 0, fData.numStage, sample, fwdSoln, 0, ub);
		lb.add(lb[iteration-1]);

		cout << "Lower bound: " << lb[iteration] << "   " <<
				"CI for upper bound: [" << ub[0][iteration]  <<
				", " <<  ub[2][iteration] << "]" << endl;
		sddp_gap = (ub[2][iteration] - lb[iteration])/ub[2][iteration];
		printf("SDDP gap: %.2f%%. \n", sddp_gap * 100);

		wall1 = get_wall_time();
		cpu1 = get_cpu_time();

		runtime = wall1 - wall0;
		cputime = cpu1 - cpu0;
		printf("Total running time %.2f seconds.\n", runtime);
		printf("Total cpu time %.2f seconds.\n", cputime);
		

		// *******************************
		// Similation
		// *******************************
		// cout << "================================" << endl;
		// cout << "Start simulation.... " << endl;
		
		// wall0 = get_wall_time();
		// // restore integrality constraints
		// for ( int t = fData.breakstage + 1; t < fData.numStage; ++t )
		// {
		// 	models[t].xLP.end();
		// 	models[t].yLP.end();
		// 	// for ( int i = 0; i < fData.intX; ++i )
		// 	// 	models[t].mod.add(IloConversion(models[t].mod.getEnv(), models[t].x[i], ILOINT));
		// 	// models[t].mod.add(IloConversion(models[t].mod.getEnv(), models[t].y1, ILOINT));
		// }

		// // forward simulation to test the performance of current approximation
		// fData.numFWsample = LARGE_SAMPLE;
		// fData.scen.b.clear();
		// fData.numScen.clear();
		// readArray<IloNumArray3> (fData.scen.b, "../data/scenRHS_test.dat");
		// fwdSoln = IloNumArray3(fData.dataEnv, fData.numStage);
		// getSamplePaths(sample, fData_p);
		// forward(models, fData_p, 0, fData.numStage, sample, fwdSoln, 0, ub);
		// simulation_gap = (ub[2][iteration+1] - lb[iteration])/ub[2][iteration+1];
		// wall1 = get_wall_time();
		// simulation_time = wall1 - wall0;
		// printf("Simulation gap: %.2f%%. \n", simulation_gap * 100);
		// printf("Total simulation time %.2f seconds.\n", simulation_time);

		// // *******************************
		// // Output
		// // *******************************
		// ofstream outputLB ("lb.csv", ios::out | ios::app);
		// if ( outputLB.is_open() )
		// {
		// 	for (int i = 0; i < lb.getSize(); ++i)
		// 		outputLB << lb[i] << ",";
		// 	outputLB << endl;
		// }
		// outputLB.close();


		// ofstream outputUB_l ("ub_l.csv", ios::out | ios::app);
		// if ( outputUB_l.is_open() )
		// {
		// 	for (int i = 0; i < ub[0].getSize() - 1; ++i)
		// 		outputUB_l << ub[0][i] << ",";
		// 	outputUB_l << endl;
		// }
		// outputUB_l.close();

		// ofstream outputUB_c ("ub_c.csv", ios::out | ios::app);
		// if ( outputUB_c.is_open() )
		// {
		// 	for (int i = 0; i < ub[1].getSize() - 1; ++i)
		// 		outputUB_c << ub[1][i] << ",";
		// 	outputUB_c << endl;
		// }
		// outputUB_c.close();

		// ofstream outputUB_r ("ub_r.csv", ios::out | ios::app);
		// if ( outputUB_r.is_open() )
		// {
		// 	for (int i = 0; i < ub[2].getSize() - 1; ++i)
		// 		outputUB_r << ub[2][i] << ",";
		// 	outputUB_r << endl;
		// }
		// outputUB_r.close();
		
		// ofstream summary ("summary.csv", ios::out | ios::app);
		// if ( summary.is_open() )
		// {
		// 	summary << cut.B << "," 
		// 		<< cut.SB << ","
		// 		<< cut.LG * cut.LEVEL << ","
		// 		<< cut.LG * (!cut.LEVEL) << ","
		// 		<< cut.I << ",,"
		// 		<< numBranch << ","
		// 		<< fData.breakstage << ","
		// 		<< initSampleSize << ","
		// 		<< numSelect << ","
		// 		<< bouncing << ","
		// 		<< iteration << ","
		// 		<< lb[iteration] << ","
		// 		<< ub[0][iteration] << ","
		// 		<< ub[2][iteration] << ","
		// 		<< sddp_gap * 100 << ","
		// 		<< runtime << ","
		// 		<< cputime << ","
		// 		<< LARGE_SAMPLE << ","
		// 		<< simulation_time << ","
		// 		<< ub[3][iteration+1] << ","
		// 		<< ub[4][iteration+1] << ","
		// 		<< ub[5][iteration+1] << ","
		// 		<< ub[1][iteration+1] << ","
		// 		<< ub[0][iteration+1] << ","
		// 		<< ub[2][iteration+1] << ","
		// 		<< simulation_gap * 100 << endl;		
		// }
		// summary.close();
		
		// free memory
		delete [] models;
		lb.end();
		ub.end();
		fwdSoln.end();
		sample.end();
	}
	catch (const IloException& e)
	{
		cerr << "Exception caught: " << e << endl;
  	}
   	catch (char const* status)
   	{
   		cerr << "Exception caught: " << status;
   	}
   	catch (...)
   	{
      	cerr << "Unknown exception caught!" << endl;
   	}

   	return 0;
}

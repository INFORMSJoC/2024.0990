package FluidModels;
import java.io.FileNotFoundException;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Iterator;

import ilog.concert.*;
import ilog.cplex.*;


public class AlgorithmHandler {



	IloCplex cplex;
	IloCplex cplexB;
	double[][] xsol;			//Solution fluid model
	double[][] xRsol;			//Solution target model
	double[][] xDsol;			//Solution target model

	double probRin;			//In training probabilities
	double probDin;

	double probRout;			//Out of sample probabilities
	double probDout;
	double flow;				// Total flow

	public AlgorithmHandler(DataHandler data) {
		probRout = 0;	
		probDout = 0;
		probRin = 0;	
		probDin = 0;
		flow = 0;		
		xsol = new double[data.numCounty][data.numCounty];
		xRsol = new double[data.numCounty][data.numCounty];
		xDsol = new double[data.numCounty][data.numCounty];
	}

	void modelFluidSAA(DataHandler data, double lambda, double limit, int peopleLimit) throws IloException, FileNotFoundException {
		int M = 100000000;
		int numSim = 100; //data.countyList.get(0).rVotesSim.length;

		cplex = new IloCplex();
		cplex.setParam(IloCplex.Param.TimeLimit, 300);
		//cplex.setParam(IloCplex.IntParam.MIPEmphasis, 4);
		//cplex.setOut(null);

		//Compute admissible pairs
		int countVars = 0;
		ArrayList<int[]> pairs = new ArrayList<int[]>();
		for (int i = 0; i < data.numCounty; i++) {
			for (int j = 0; j < data.numCounty; j++) {
				if(!data.countyList.get(i).state.equals(data.countyList.get(j).state) && data.distMatrix[i][j] <= limit
						&& data.countyList.get(i).pop >= 50000 && data.countyList.get(j).pop >= 50000 && (data.countyList.get(i).avgRfrac >= 0.5 || data.countyList.get(i).avgDfrac >= 0.5) ) { //Don't move inside the same state!
					int[] aux = new int[2];
					aux[0] = i;
					aux[1] = j;
					pairs.add(aux);
					countVars++;
				}
			}
		}
		System.out.println("TOTAL NUMBER OF FLOW VARS: "+countVars);

		// Define decision variables
		IloNumVar[] x = new IloNumVar[countVars];
		IloNumVar[][] zR = new IloNumVar[data.numCounty][numSim];				// Votes per county
		IloNumVar[][] zD = new IloNumVar[data.numCounty][numSim];
		IloNumVar[][] TR = new IloNumVar[data.numStates][numSim];				// Votes per state
		IloNumVar[][] TD = new IloNumVar[data.numStates][numSim];
		IloNumVar elecR[] = new IloNumVar[numSim];
		IloNumVar elecD[] = new IloNumVar[numSim];
		IloNumVar[][] vD = new IloNumVar[data.numStates][numSim];				//Binary aux
		IloNumVar[][] vR = new IloNumVar[data.numStates][numSim];
		// Special districts 0=Maine1, 1=Maine2, ... 4=Nebraska3
		IloNumVar[][] TRE = new IloNumVar[5][numSim];				// Votes per special district
		IloNumVar[][] TDE = new IloNumVar[5][numSim];
		IloNumVar[][] vDE = new IloNumVar[5][numSim];				//Binary aux
		IloNumVar[][] vRE = new IloNumVar[5][numSim];
		// Prefered party wins per simulation
		IloNumVar partyWins[] = new IloNumVar[numSim];

		for (int k = 0; k < countVars; k++) {
			int i = pairs.get(k)[0];
			x[k] = cplex.numVar(0, data.countyList.get(i).pop, IloNumVarType.Float);
		}

		for (int q = 0; q < numSim; q++) {

			for (int i = 0; i < data.numCounty; i++) {
				zR[i][q] = cplex.numVar(0, Double.MAX_VALUE, IloNumVarType.Float);
				zD[i][q] = cplex.numVar(0, Double.MAX_VALUE, IloNumVarType.Float);	
			}
			for (int s = 0; s < data.numStates; s++) {
				TR[s][q] = cplex.numVar(0, Double.MAX_VALUE, IloNumVarType.Float);
				TD[s][q] = cplex.numVar(0, Double.MAX_VALUE, IloNumVarType.Float);
				vD[s][q] = cplex.numVar(0, 1, IloNumVarType.Bool);
				vR[s][q] = cplex.numVar(0, 1, IloNumVarType.Bool);
			}
			for (int s = 0; s < 5; s++) {
				TRE[s][q] = cplex.numVar(0, Double.MAX_VALUE, IloNumVarType.Float);
				TDE[s][q] = cplex.numVar(0, Double.MAX_VALUE, IloNumVarType.Float);
				vDE[s][q] = cplex.numVar(0, 1, IloNumVarType.Bool);
				vRE[s][q] = cplex.numVar(0, 1, IloNumVarType.Bool);
			}
			elecR[q] = cplex.numVar(0, 9999, IloNumVarType.Int); 
			elecD[q] = cplex.numVar(0, 9999, IloNumVarType.Int); 
			partyWins[q] = cplex.numVar(0, 1, IloNumVarType.Bool); 
		}
		// Total people moved out of a county
		for (int i = 0; i < data.numCounty; i++) {
			IloLinearNumExpr expr = cplex.linearNumExpr();
			for (int k = 0; k < countVars; k++) {
				int origin = pairs.get(k)[0];
				if(origin == i) {
					expr.addTerm(1, x[k]);
				}
			}
			cplex.addLe(expr, lambda*data.countyList.get(i).pop);
		}
		// Total people moved into a county
		for (int i = 0; i < data.numCounty; i++) {
			IloLinearNumExpr expr = cplex.linearNumExpr();
			for (int k = 0; k < countVars; k++) {
				int dest = pairs.get(k)[1];
				if(dest == i) {
					expr.addTerm(1, x[k]);
				}
			}
			cplex.addLe(expr, lambda*data.countyList.get(i).pop);
		}
		// Total votes for a county
		for (int q = 0; q < numSim; q++) {
			IloLinearNumExpr[] exprR = new IloLinearNumExpr[data.numCounty];
			IloLinearNumExpr[] exprD = new IloLinearNumExpr[data.numCounty];
			// Initialize constraints
			for (int i = 0; i < data.numCounty; i++) {
				exprR[i] = cplex.linearNumExpr();
				exprD[i] = cplex.linearNumExpr();
			}
			// Consider movement
			for (int k = 0; k < countVars; k++) {
				int origin = pairs.get(k)[0];
				int destination = pairs.get(k)[1];
				double fracRorigin = data.countyList.get(origin).rVotesSim[q]/(data.countyList.get(origin).pop+0.0);
				double fracDorigin = data.countyList.get(origin).dVotesSim[q]/(data.countyList.get(origin).pop+0.0);
				exprR[origin].addTerm(fracRorigin, x[k]);
				exprR[destination].addTerm(-fracRorigin, x[k]);
				exprD[origin].addTerm(fracDorigin, x[k]);
				exprD[destination].addTerm(-fracDorigin, x[k]);
			}

			// Add coefficients and rhs for all counties
			for (int i = 0; i < data.numCounty; i++) {
				County c = data.countyList.get(i);
				exprR[i].addTerm(1, zR[i][q]);
				exprD[i].addTerm(1, zD[i][q]);
				cplex.addEq(exprR[i], c.rVotesSim[q]);
				cplex.addEq(exprD[i], c.dVotesSim[q]);

			}
		}

		// Total votes for a state
		for (int q = 0; q < numSim; q++) {
			for (int s = 0; s < data.numStates; s++) {
				IloLinearNumExpr exprR = cplex.linearNumExpr();
				exprR.addTerm(1, TR[s][q]);
				for (int i = 0; i < data.numCounty; i++) {
					if(data.countyList.get(i).state.equals(data.states[s])) {
						exprR.addTerm(-1, zR[i][q]);
						//System.out.println("MATCH COUNTY "+i+" "+data.states[s]);
					}
				}
				cplex.addEq(exprR, 0);

				IloLinearNumExpr exprD = cplex.linearNumExpr();
				exprD.addTerm(1, TD[s][q]);
				for (int i = 0; i < data.numCounty; i++) {
					if(data.countyList.get(i).state.equals(data.states[s])) {
						exprD.addTerm(-1, zD[i][q]);
						//System.out.println("MATCH COUNTY "+i+" "+data.states[s]);
					}
				}
				cplex.addEq(exprD, 0);

			}
			// Total votes for special districts
			for (int s = 0; s < 5; s++) {
				IloLinearNumExpr exprR = cplex.linearNumExpr();
				exprR.addTerm(1, TRE[s][q]);
				if(s==0) {
					for (int i = 0; i < data.maine1.size(); i++) {
						exprR.addTerm(-1, zR[data.maine1.get(i)][q]);
						//System.out.println("MAINE 1 COUNTY "+data.countyList.get(data.maine1.get(i)).name);
					}
					//Only half of the votes for shared county 23011
					exprR.addTerm(0.5, zR[data.countyIndex.get(23011)][q]);
				}
				else if(s==1) {
					for (int i = 0; i < data.maine2.size(); i++) {
						exprR.addTerm(-1, zR[data.maine2.get(i)][q]);
						//System.out.println("MAINE 2 COUNTY "+data.countyList.get(data.maine2.get(i)).name);
					}
					//Only half of the votes for shared county 23011
					exprR.addTerm(0.5, zR[data.countyIndex.get(23011)][q]);
				}
				else if(s==2) {
					for (int i = 0; i < data.nebraska1.size(); i++) {
						exprR.addTerm(-1, zR[data.nebraska1.get(i)][q]);
					}
					//Only half of the votes for shared county 31153
					exprR.addTerm(0.5, zR[data.countyIndex.get(31153)][q]);
				}
				else if(s==3) {
					for (int i = 0; i < data.nebraska2.size(); i++) {
						exprR.addTerm(-1, zR[data.nebraska2.get(i)][q]);
					}
					//Only half of the votes for shared county 31153
					exprR.addTerm(0.5, zR[data.countyIndex.get(31153)][q]);
				}
				else if(s==4) {
					for (int i = 0; i < data.nebraska3.size(); i++) {
						exprR.addTerm(-1, zR[data.nebraska3.get(i)][q]);
					}
				}

				cplex.addEq(exprR, 0);


				IloLinearNumExpr exprD = cplex.linearNumExpr();
				exprD.addTerm(1, TDE[s][q]);
				if(s==0) {	
					for (int i = 0; i < data.maine1.size(); i++) {
						exprD.addTerm(-1, zD[data.maine1.get(i)][q]);
						//System.out.println("MAINE 1 COUNTY "+data.countyList.get(data.maine1.get(i)).name);
					}
					//Only half of the votes for shared county 23011
					exprD.addTerm(0.5, zD[data.countyIndex.get(23011)][q]);
				}
				else if(s==1) {	
					for (int i = 0; i < data.maine2.size(); i++) {
						exprD.addTerm(-1, zD[data.maine2.get(i)][q]);
						//System.out.println("MAINE 2 COUNTY "+data.countyList.get(data.maine2.get(i)).name);
					}
					//Only half of the votes for shared county 23011
					exprD.addTerm(0.5, zD[data.countyIndex.get(23011)][q]);
				}
				else if(s==2) {	
					for (int i = 0; i < data.nebraska1.size(); i++) {
						exprD.addTerm(-1, zD[data.nebraska1.get(i)][q]);
					}
					//Only half of the votes for shared county 31153
					exprD.addTerm(0.5, zD[data.countyIndex.get(31153)][q]);
				}
				else if(s==3) {	
					for (int i = 0; i < data.nebraska2.size(); i++) {
						exprD.addTerm(-1, zD[data.nebraska2.get(i)][q]);
					}
					//Only half of the votes for shared county 31153
					exprD.addTerm(0.5, zD[data.countyIndex.get(31153)][q]);
				}
				else if(s==4) {	
					for (int i = 0; i < data.nebraska3.size(); i++) {
						exprD.addTerm(-1, zD[data.nebraska3.get(i)][q]);
					}
				}

				cplex.addEq(exprD, 0);

			}
		}
		// Elec votes count
		for (int q = 0; q < numSim; q++) {

			// Every state must be won
			for (int s = 0; s < data.numStates; s++) {
				IloLinearNumExpr expr = cplex.linearNumExpr();
				expr.addTerm(1, vR[s][q]);
				expr.addTerm(1, vD[s][q]);
				cplex.addEq(expr, 1);
			}
			// Every SD must be won
			for (int s = 0; s < 5; s++) {
				IloLinearNumExpr expr = cplex.linearNumExpr();
				expr.addTerm(1, vRE[s][q]);
				expr.addTerm(1, vDE[s][q]);
				cplex.addEq(expr, 1);
			}

			// Binary activations R
			for (int s = 0; s < data.numStates; s++) {
				IloLinearNumExpr expr = cplex.linearNumExpr();
				expr.addTerm(1, TR[s][q]);
				expr.addTerm(-1, TD[s][q]);
				expr.addTerm(-1, vR[s][q]);
				expr.addTerm(-M, vR[s][q]);
				cplex.addGe(expr, -M);
			}
			for (int s = 0; s < 5; s++) {
				IloLinearNumExpr expr = cplex.linearNumExpr();
				expr.addTerm(1, TRE[s][q]);
				expr.addTerm(-1, TDE[s][q]);
				expr.addTerm(-1, vRE[s][q]);
				expr.addTerm(-M, vRE[s][q]);
				cplex.addGe(expr, -M);
			}

			// Binary activations D
			for (int s = 0; s < data.numStates; s++) {
				IloLinearNumExpr expr = cplex.linearNumExpr();
				expr.addTerm(1, TD[s][q]);
				expr.addTerm(-1, TR[s][q]);
				expr.addTerm(-1, vD[s][q]);
				expr.addTerm(-M, vD[s][q]);
				cplex.addGe(expr, -M);
			}
			for (int s = 0; s < 5; s++) {
				IloLinearNumExpr expr = cplex.linearNumExpr();
				expr.addTerm(1, TDE[s][q]);
				expr.addTerm(-1, TRE[s][q]);
				expr.addTerm(-1, vDE[s][q]);
				expr.addTerm(-M, vDE[s][q]);
				cplex.addGe(expr, -M);
			}

			//Total Electoral votes for R
			IloLinearNumExpr exprA = cplex.linearNumExpr();
			exprA.addTerm(1, elecR[q]);
			for (int s = 0; s < data.numStates; s++) {
				exprA.addTerm(-data.eVotes[s], vR[s][q]);
			}
			for (int s = 0; s < 5; s++) {
				exprA.addTerm(-1, vRE[s][q]);
			}
			cplex.addEq(exprA, 0);

			//Total Electoral votes for D
			IloLinearNumExpr exprB = cplex.linearNumExpr();
			exprB.addTerm(1, elecD[q]);
			for (int s = 0; s < data.numStates; s++) {
				exprB.addTerm(-data.eVotes[s], vD[s][q]);
			}
			for (int s = 0; s < 5; s++) {
				exprB.addTerm(-1, vDE[s][q]);
			}
			cplex.addEq(exprB, 0);

		}
		// total number of people moved
		IloLinearNumExpr exprP = cplex.linearNumExpr();
		for (int k = 0; k < countVars; k++) {
			exprP.addTerm(1, x[k]);
		}
		cplex.addLe(exprP, peopleLimit);

		// Prefered party wins count
		for (int q = 0; q < numSim; q++) {
			IloLinearNumExpr expr = cplex.linearNumExpr();
			expr.addTerm(1, elecR[q]);
			expr.addTerm(-270, partyWins[q]);
			cplex.addGe(expr, 0);
		}

		// Objective prob	
		IloLinearNumExpr obj = cplex.linearNumExpr();
		for (int q = 0; q < numSim; q++) {
			obj.addTerm(1, partyWins[q]);
		}

		cplex.addMaximize(obj);
		cplex.solve();

		if(cplex.getStatus() == IloCplex.Status.Infeasible || cplex.getStatus() == IloCplex.Status.InfeasibleOrUnbounded){
			System.out.println("		MODEL FLUID SAA is infeasible!");
		}
		else if(cplex.getStatus() == IloCplex.Status.Feasible || cplex.getStatus() == IloCplex.Status.Optimal) {
			////////////////////////////////////////////////// EVAL SOL /////////////////////////////////////////////////////////////////
			xsol = new double[data.numCounty][data.numCounty];
			for (int k = 0; k < countVars; k++) {
				int origin = pairs.get(k)[0];
				int destination = pairs.get(k)[1];
				xsol[origin][destination] = cplex.getValue(x[k]);
			}
			getStatsNoMovement(data, xsol); 

			////////////////////////////////////////////////// PRINT RESULTS /////////////////////////////////////////////////////////////////

			System.out.println("******************************  RESULTS ******************************************");
			double tFlow = 0;
			double tDist = 0;
			for (int k = 0; k < countVars; k++) {
				tFlow += cplex.getValue(x[k]);
				int origin = pairs.get(k)[0];
				int destination = pairs.get(k)[1];
				tDist += data.distMatrix[origin][destination]*cplex.getValue(x[k]);
			}

			System.out.println("TOTAL FLOW:,"+ tFlow);
			System.out.println("TOTAL DISTANCE:,"+ tDist);
			System.out.println("TOTAL R WINS:,"+ cplex.getObjValue());


			// Export output
			java.io.PrintStream ps = new java.io.PrintStream( new java.io.FileOutputStream("2024_"+limit+"_"+lambda+".txt"));
			ps.println("TOTAL FLOW:,"+ tFlow);
			ps.println("TOTAL DISTANCE:,"+ tDist);

			/*
			ps.println("ELECTORAL_VOTES_R:,"+ cplex.getValue(elecR));
			ps.println("ELECTORAL_VOTES_D:,"+ cplex.getValue(elecD));

			int totalR = 0;
			int totalD = 0;
			for (int s = 0; s < data.numStates; s++) {
				totalR += cplex.getValue(TR[s]);
				totalD += cplex.getValue(TD[s]);
			}

			ps.println("TOTAL_VOTES_R:,"+ totalR);
			ps.println("TOTAL_VOTES_D:,"+ totalD);

			ps.println("State, ElectoralVotes, TotalRVotes, TotalDVotes, Winner");
			for (int s = 0; s < data.numStates; s++) {
				if(cplex.getValue(vR[s])>0.01) {
					ps.println(data.states[s]+","+data.eVotes[s]+","+cplex.getValue(TR[s])+","+cplex.getValue(TD[s])+",R");
				}
				else {
					ps.println(data.states[s]+","+data.eVotes[s]+","+cplex.getValue(TR[s])+","+cplex.getValue(TD[s])+",D");
				}		
			}
			if(cplex.getValue(vRE[0])>0.01) {
				ps.println("MAINE_SD1"+","+1+","+cplex.getValue(TRE[0])+","+cplex.getValue(TDE[0])+",R");
			}
			else {
				ps.println("MAINE_SD1"+","+1+","+cplex.getValue(TRE[0])+","+cplex.getValue(TDE[0])+",D");
			}
			if(cplex.getValue(vRE[1])>0.01) {
				ps.println("MAINE_SD2"+","+1+","+cplex.getValue(TRE[1])+","+cplex.getValue(TDE[1])+",R");
			}
			else {
				ps.println("MAINE_SD2"+","+1+","+cplex.getValue(TRE[1])+","+cplex.getValue(TDE[1])+",D");
			}
			if(cplex.getValue(vRE[2])>0.01) {
				ps.println("NEBRASKA_SD1"+","+1+","+cplex.getValue(TRE[2])+","+cplex.getValue(TDE[2])+",R");
			}
			else {
				ps.println("NEBRASKA_SD1"+","+1+","+cplex.getValue(TRE[2])+","+cplex.getValue(TDE[2])+",D");
			}
			if(cplex.getValue(vRE[3])>0.01) {
				ps.println("NEBRASKA_SD2"+","+1+","+cplex.getValue(TRE[3])+","+cplex.getValue(TDE[3])+",R");
			}
			else {
				ps.println("NEBRASKA_SD2"+","+1+","+cplex.getValue(TRE[3])+","+cplex.getValue(TDE[3])+",D");
			}
			if(cplex.getValue(vRE[4])>0.01) {
				ps.println("NEBRASKA_SD3"+","+1+","+cplex.getValue(TRE[4])+","+cplex.getValue(TDE[4])+",R");
			}
			else {
				ps.println("NEBRASKA_SD3"+","+1+","+cplex.getValue(TRE[4])+","+cplex.getValue(TDE[4])+",D");
			}
			 */
			ps.println("FLOW VARS: Origin name, Origin ID, Origin State, Origin Pop, Destination name, Destination ID, Destination State, Destination Pop, Distance, Flow");
			for (int k = 0; k < countVars; k++) {
				int origin = pairs.get(k)[0];
				int destination = pairs.get(k)[1];
				if(cplex.getValue(x[k])>0.0001) {
					ps.println(data.countyList.get(origin).name+","+data.countyList.get(origin).ID+","+data.countyList.get(origin).state+","+data.countyList.get(origin).pop+","+data.countyList.get(destination).name+","+data.countyList.get(destination).ID+","+data.countyList.get(destination).state+","+data.countyList.get(destination).pop+","+data.distMatrix[origin][destination]+","+cplex.getValue(x[k]));
				}
			}
			/*	
			ps.println("County name, County ID, County total votes R, County total votes D");
			for (int i = 0; i < data.numCounty; i++) {
				ps.println(data.countyList.get(i).name+","+data.countyList.get(i).ID+","+cplex.getValue(zR[i])+","+cplex.getValue(zD[i]));
			}
			/* 
		}

		cplex.clearModel();
		cplex.end();	

	}

	void modelTargetSAA(DataHandler data, double lambda, double limit, int peopleLimit) throws IloException, FileNotFoundException {
		int M = 100000000;
		int numSim = 100; //data.countyList.get(0).rVotesSim.length;

		cplex = new IloCplex();
		cplex.setParam(IloCplex.Param.TimeLimit, 3600);
		cplex.setParam(IloCplex.IntParam.MIPEmphasis, 2);
		//cplex.setOut(null);

		//Compute admissible pairs
		int countVars = 0;
		ArrayList<int[]> pairs = new ArrayList<int[]>();
		for (int i = 0; i < data.numCounty; i++) {
			for (int j = 0; j < data.numCounty; j++) {
				if(!data.countyList.get(i).state.equals(data.countyList.get(j).state) && data.distMatrix[i][j] <= limit
						&& data.countyList.get(i).pop >= 50000 && data.countyList.get(j).pop >= 50000 && (data.countyList.get(i).avgRfrac >= 0.5 || data.countyList.get(i).avgDfrac >= 0.5) ) { //Don't move inside the same state!
					int[] aux = new int[2];
					aux[0] = i;
					aux[1] = j;
					pairs.add(aux);
					countVars++;
				}
			}
		}
		System.out.println("TOTAL NUMBER OF FLOW VARS: "+countVars);

		// Define decision variables
		IloNumVar[] xR = new IloNumVar[countVars];
		IloNumVar[] xD = new IloNumVar[countVars];
		IloNumVar[][] zR = new IloNumVar[data.numCounty][numSim];				// Votes per county
		IloNumVar[][] zD = new IloNumVar[data.numCounty][numSim];
		IloNumVar[][] TR = new IloNumVar[data.numStates][numSim];				// Votes per state
		IloNumVar[][] TD = new IloNumVar[data.numStates][numSim];
		IloNumVar elecR[] = new IloNumVar[numSim];
		IloNumVar elecD[] = new IloNumVar[numSim];
		IloNumVar[][] vD = new IloNumVar[data.numStates][numSim];				//Binary aux
		IloNumVar[][] vR = new IloNumVar[data.numStates][numSim];
		// Special districts 0=Maine1, 1=Maine2, ... 4=Nebraska3
		IloNumVar[][] TRE = new IloNumVar[5][numSim];				// Votes per special district
		IloNumVar[][] TDE = new IloNumVar[5][numSim];
		IloNumVar[][] vDE = new IloNumVar[5][numSim];				//Binary aux
		IloNumVar[][] vRE = new IloNumVar[5][numSim];
		// Prefered party wins per simulation
		IloNumVar partyWins[] = new IloNumVar[numSim];

		for (int k = 0; k < countVars; k++) {
			int i = pairs.get(k)[0];
			xR[k] = cplex.numVar(0, data.countyList.get(i).pop, IloNumVarType.Float); //No negativity of the county count takes care of the pool of people to be moved!
			xD[k] = cplex.numVar(0, data.countyList.get(i).pop, IloNumVarType.Float);
		}

		for (int q = 0; q < numSim; q++) {

			for (int i = 0; i < data.numCounty; i++) {
				zR[i][q] = cplex.numVar(0, Double.MAX_VALUE, IloNumVarType.Float);
				zD[i][q] = cplex.numVar(0, Double.MAX_VALUE, IloNumVarType.Float);	
			}
			for (int s = 0; s < data.numStates; s++) {
				TR[s][q] = cplex.numVar(0, Double.MAX_VALUE, IloNumVarType.Float);
				TD[s][q] = cplex.numVar(0, Double.MAX_VALUE, IloNumVarType.Float);
				vD[s][q] = cplex.numVar(0, 1, IloNumVarType.Bool);
				vR[s][q] = cplex.numVar(0, 1, IloNumVarType.Bool);
			}
			for (int s = 0; s < 5; s++) {
				TRE[s][q] = cplex.numVar(0, Double.MAX_VALUE, IloNumVarType.Float);
				TDE[s][q] = cplex.numVar(0, Double.MAX_VALUE, IloNumVarType.Float);
				vDE[s][q] = cplex.numVar(0, 1, IloNumVarType.Bool);
				vRE[s][q] = cplex.numVar(0, 1, IloNumVarType.Bool);
			}
			elecR[q] = cplex.numVar(0, 9999, IloNumVarType.Int); 
			elecD[q] = cplex.numVar(0, 9999, IloNumVarType.Int); 
			partyWins[q] = cplex.numVar(0, 1, IloNumVarType.Bool); 
		}
		// Total people moved out of a county
		for (int i = 0; i < data.numCounty; i++) {
			IloLinearNumExpr expr = cplex.linearNumExpr();
			for (int k = 0; k < countVars; k++) {
				int origin = pairs.get(k)[0];
				if(origin == i) {
					expr.addTerm(1, xR[k]);
					expr.addTerm(1, xD[k]);
				}
			}
			cplex.addLe(expr, lambda*data.countyList.get(i).pop);
		}
		// Total people moved into a county
		for (int i = 0; i < data.numCounty; i++) {
			IloLinearNumExpr expr = cplex.linearNumExpr();
			for (int k = 0; k < countVars; k++) {
				int dest = pairs.get(k)[1];
				if(dest == i) {
					expr.addTerm(1, xR[k]);
					expr.addTerm(1, xD[k]);
				}
			}
			cplex.addLe(expr, lambda*data.countyList.get(i).pop);
		}
		// Total votes for a county
		for (int q = 0; q < numSim; q++) {
			IloLinearNumExpr[] exprR = new IloLinearNumExpr[data.numCounty];
			IloLinearNumExpr[] exprD = new IloLinearNumExpr[data.numCounty];
			// Initialize constraints
			for (int i = 0; i < data.numCounty; i++) {
				exprR[i] = cplex.linearNumExpr();
				exprD[i] = cplex.linearNumExpr();
			}
			// Consider movement
			for (int k = 0; k < countVars; k++) {
				int origin = pairs.get(k)[0];
				int destination = pairs.get(k)[1];
				exprR[origin].addTerm(1, xR[k]);
				exprR[destination].addTerm(-1, xR[k]);
				exprD[origin].addTerm(1, xD[k]);
				exprD[destination].addTerm(-1, xD[k]);
			}

			// Add coefficients and rhs for all counties
			for (int i = 0; i < data.numCounty; i++) {
				County c = data.countyList.get(i);
				exprR[i].addTerm(1, zR[i][q]);
				exprD[i].addTerm(1, zD[i][q]);
				cplex.addEq(exprR[i], c.rVotesSim[q]);
				cplex.addEq(exprD[i], c.dVotesSim[q]);

			}
		}

		// Total votes for a state
		for (int q = 0; q < numSim; q++) {
			for (int s = 0; s < data.numStates; s++) {
				IloLinearNumExpr exprR = cplex.linearNumExpr();
				exprR.addTerm(1, TR[s][q]);
				for (int i = 0; i < data.numCounty; i++) {
					if(data.countyList.get(i).state.equals(data.states[s])) {
						exprR.addTerm(-1, zR[i][q]);
						//System.out.println("MATCH COUNTY "+i+" "+data.states[s]);
					}
				}
				cplex.addEq(exprR, 0);

				IloLinearNumExpr exprD = cplex.linearNumExpr();
				exprD.addTerm(1, TD[s][q]);
				for (int i = 0; i < data.numCounty; i++) {
					if(data.countyList.get(i).state.equals(data.states[s])) {
						exprD.addTerm(-1, zD[i][q]);
						//System.out.println("MATCH COUNTY "+i+" "+data.states[s]);
					}
				}
				cplex.addEq(exprD, 0);

			}
			// Total votes for special districts
			for (int s = 0; s < 5; s++) {
				IloLinearNumExpr exprR = cplex.linearNumExpr();
				exprR.addTerm(1, TRE[s][q]);
				if(s==0) {
					for (int i = 0; i < data.maine1.size(); i++) {
						exprR.addTerm(-1, zR[data.maine1.get(i)][q]);
						//System.out.println("MAINE 1 COUNTY "+data.countyList.get(data.maine1.get(i)).name);
					}
					//Only half of the votes for shared county 23011
					exprR.addTerm(0.5, zR[data.countyIndex.get(23011)][q]);
				}
				else if(s==1) {
					for (int i = 0; i < data.maine2.size(); i++) {
						exprR.addTerm(-1, zR[data.maine2.get(i)][q]);
						//System.out.println("MAINE 2 COUNTY "+data.countyList.get(data.maine2.get(i)).name);
					}
					//Only half of the votes for shared county 23011
					exprR.addTerm(0.5, zR[data.countyIndex.get(23011)][q]);
				}
				else if(s==2) {
					for (int i = 0; i < data.nebraska1.size(); i++) {
						exprR.addTerm(-1, zR[data.nebraska1.get(i)][q]);
					}
					//Only half of the votes for shared county 31153
					exprR.addTerm(0.5, zR[data.countyIndex.get(31153)][q]);
				}
				else if(s==3) {
					for (int i = 0; i < data.nebraska2.size(); i++) {
						exprR.addTerm(-1, zR[data.nebraska2.get(i)][q]);
					}
					//Only half of the votes for shared county 31153
					exprR.addTerm(0.5, zR[data.countyIndex.get(31153)][q]);
				}
				else if(s==4) {
					for (int i = 0; i < data.nebraska3.size(); i++) {
						exprR.addTerm(-1, zR[data.nebraska3.get(i)][q]);
					}
				}

				cplex.addEq(exprR, 0);


				IloLinearNumExpr exprD = cplex.linearNumExpr();
				exprD.addTerm(1, TDE[s][q]);
				if(s==0) {	
					for (int i = 0; i < data.maine1.size(); i++) {
						exprD.addTerm(-1, zD[data.maine1.get(i)][q]);
						//System.out.println("MAINE 1 COUNTY "+data.countyList.get(data.maine1.get(i)).name);
					}
					//Only half of the votes for shared county 23011
					exprD.addTerm(0.5, zD[data.countyIndex.get(23011)][q]);
				}
				else if(s==1) {	
					for (int i = 0; i < data.maine2.size(); i++) {
						exprD.addTerm(-1, zD[data.maine2.get(i)][q]);
						//System.out.println("MAINE 2 COUNTY "+data.countyList.get(data.maine2.get(i)).name);
					}
					//Only half of the votes for shared county 23011
					exprD.addTerm(0.5, zD[data.countyIndex.get(23011)][q]);
				}
				else if(s==2) {	
					for (int i = 0; i < data.nebraska1.size(); i++) {
						exprD.addTerm(-1, zD[data.nebraska1.get(i)][q]);
					}
					//Only half of the votes for shared county 31153
					exprD.addTerm(0.5, zD[data.countyIndex.get(31153)][q]);
				}
				else if(s==3) {	
					for (int i = 0; i < data.nebraska2.size(); i++) {
						exprD.addTerm(-1, zD[data.nebraska2.get(i)][q]);
					}
					//Only half of the votes for shared county 31153
					exprD.addTerm(0.5, zD[data.countyIndex.get(31153)][q]);
				}
				else if(s==4) {	
					for (int i = 0; i < data.nebraska3.size(); i++) {
						exprD.addTerm(-1, zD[data.nebraska3.get(i)][q]);
					}
				}

				cplex.addEq(exprD, 0);

			}
		}
		// Elec votes count
		for (int q = 0; q < numSim; q++) {

			// Every state must be won
			for (int s = 0; s < data.numStates; s++) {
				IloLinearNumExpr expr = cplex.linearNumExpr();
				expr.addTerm(1, vR[s][q]);
				expr.addTerm(1, vD[s][q]);
				cplex.addEq(expr, 1);
			}
			// Every SD must be won
			for (int s = 0; s < 5; s++) {
				IloLinearNumExpr expr = cplex.linearNumExpr();
				expr.addTerm(1, vRE[s][q]);
				expr.addTerm(1, vDE[s][q]);
				cplex.addEq(expr, 1);
			}

			// Binary activations R
			for (int s = 0; s < data.numStates; s++) {
				IloLinearNumExpr expr = cplex.linearNumExpr();
				expr.addTerm(1, TR[s][q]);
				expr.addTerm(-1, TD[s][q]);
				expr.addTerm(-1, vR[s][q]);
				expr.addTerm(-M, vR[s][q]);
				cplex.addGe(expr, -M);
			}
			for (int s = 0; s < 5; s++) {
				IloLinearNumExpr expr = cplex.linearNumExpr();
				expr.addTerm(1, TRE[s][q]);
				expr.addTerm(-1, TDE[s][q]);
				expr.addTerm(-1, vRE[s][q]);
				expr.addTerm(-M, vRE[s][q]);
				cplex.addGe(expr, -M);
			}

			// Binary activations D
			for (int s = 0; s < data.numStates; s++) {
				IloLinearNumExpr expr = cplex.linearNumExpr();
				expr.addTerm(1, TD[s][q]);
				expr.addTerm(-1, TR[s][q]);
				expr.addTerm(-1, vD[s][q]);
				expr.addTerm(-M, vD[s][q]);
				cplex.addGe(expr, -M);
			}
			for (int s = 0; s < 5; s++) {
				IloLinearNumExpr expr = cplex.linearNumExpr();
				expr.addTerm(1, TDE[s][q]);
				expr.addTerm(-1, TRE[s][q]);
				expr.addTerm(-1, vDE[s][q]);
				expr.addTerm(-M, vDE[s][q]);
				cplex.addGe(expr, -M);
			}

			//Total Electoral votes for R
			IloLinearNumExpr exprA = cplex.linearNumExpr();
			exprA.addTerm(1, elecR[q]);
			for (int s = 0; s < data.numStates; s++) {
				exprA.addTerm(-data.eVotes[s], vR[s][q]);
			}
			for (int s = 0; s < 5; s++) {
				exprA.addTerm(-1, vRE[s][q]);
			}
			cplex.addEq(exprA, 0);

			//Total Electoral votes for D
			IloLinearNumExpr exprB = cplex.linearNumExpr();
			exprB.addTerm(1, elecD[q]);
			for (int s = 0; s < data.numStates; s++) {
				exprB.addTerm(-data.eVotes[s], vD[s][q]);
			}
			for (int s = 0; s < 5; s++) {
				exprB.addTerm(-1, vDE[s][q]);
			}
			cplex.addEq(exprB, 0);

		}
		// total number of people moved
		IloLinearNumExpr exprP = cplex.linearNumExpr();
		for (int k = 0; k < countVars; k++) {
			exprP.addTerm(1, xR[k]);
			exprP.addTerm(1, xD[k]);
		}
		cplex.addLe(exprP, peopleLimit);

		// Prefered party wins count
		for (int q = 0; q < numSim; q++) {
			IloLinearNumExpr expr = cplex.linearNumExpr();
			expr.addTerm(1, elecR[q]);
			expr.addTerm(-270, partyWins[q]);
			cplex.addGe(expr, 0);
		}

		// Objective	
		IloLinearNumExpr obj = cplex.linearNumExpr();
		for (int q = 0; q < numSim; q++) {
			obj.addTerm(1, partyWins[q]);
		}
		cplex.addMaximize(obj);
		cplex.solve();

		if(cplex.getStatus() == IloCplex.Status.Infeasible || cplex.getStatus() == IloCplex.Status.InfeasibleOrUnbounded){
			System.out.println("		MODEL TARGET SAA is infeasible!");
		}
		else if(cplex.getStatus() == IloCplex.Status.Feasible || cplex.getStatus() == IloCplex.Status.Optimal) {
			////////////////////////////////////////////////// EVAL SOL /////////////////////////////////////////////////////////////////
			double[][] xRsol = new double[data.numCounty][data.numCounty];
			for (int k = 0; k < countVars; k++) {
				int origin = pairs.get(k)[0];
				int destination = pairs.get(k)[1];
				xRsol[origin][destination] = cplex.getValue(xR[k]);
			}
			//evalSolFluid(data, xsol, peopleLimit); 

			////////////////////////////////////////////////// PRINT RESULTS /////////////////////////////////////////////////////////////////

			System.out.println("******************************  RESULTS ******************************************");
			double tFlow = 0;
			double tDist = 0;
			for (int k = 0; k < countVars; k++) {
				tFlow += cplex.getValue(xR[k]);
				tFlow += cplex.getValue(xD[k]);
			}

			System.out.println("TOTAL FLOW:,"+ tFlow);
			System.out.println("TOTAL R WINS:,"+ cplex.getObjValue());

			/*
			System.out.println("ELECTORAL_VOTES_R:,"+ Math.round(cplex.getValue(elecR)));
			System.out.println("ELECTORAL_VOTES_D:,"+ Math.round(cplex.getValue(elecD)));
			// Export output
			java.io.PrintStream ps = new java.io.PrintStream( new java.io.FileOutputStream("2024_"+limit+"_"+lambda+".txt"));
			ps.println("TOTAL FLOW:,"+ tFlow);
			ps.println("TOTAL DISTANCE:,"+ tDist);
			ps.println("ELECTORAL_VOTES_R:,"+ cplex.getValue(elecR));
			ps.println("ELECTORAL_VOTES_D:,"+ cplex.getValue(elecD));


			int totalR = 0;
			int totalD = 0;
			for (int s = 0; s < data.numStates; s++) {
				totalR += cplex.getValue(TR[s]);
				totalD += cplex.getValue(TD[s]);
			}

			ps.println("TOTAL_VOTES_R:,"+ totalR);
			ps.println("TOTAL_VOTES_D:,"+ totalD);

			ps.println("State, ElectoralVotes, TotalRVotes, TotalDVotes, Winner");
			for (int s = 0; s < data.numStates; s++) {
				if(cplex.getValue(vR[s])>0.01) {
					ps.println(data.states[s]+","+data.eVotes[s]+","+cplex.getValue(TR[s])+","+cplex.getValue(TD[s])+",R");
				}
				else {
					ps.println(data.states[s]+","+data.eVotes[s]+","+cplex.getValue(TR[s])+","+cplex.getValue(TD[s])+",D");
				}		
			}
			if(cplex.getValue(vRE[0])>0.01) {
				ps.println("MAINE_SD1"+","+1+","+cplex.getValue(TRE[0])+","+cplex.getValue(TDE[0])+",R");
			}
			else {
				ps.println("MAINE_SD1"+","+1+","+cplex.getValue(TRE[0])+","+cplex.getValue(TDE[0])+",D");
			}
			if(cplex.getValue(vRE[1])>0.01) {
				ps.println("MAINE_SD2"+","+1+","+cplex.getValue(TRE[1])+","+cplex.getValue(TDE[1])+",R");
			}
			else {
				ps.println("MAINE_SD2"+","+1+","+cplex.getValue(TRE[1])+","+cplex.getValue(TDE[1])+",D");
			}
			if(cplex.getValue(vRE[2])>0.01) {
				ps.println("NEBRASKA_SD1"+","+1+","+cplex.getValue(TRE[2])+","+cplex.getValue(TDE[2])+",R");
			}
			else {
				ps.println("NEBRASKA_SD1"+","+1+","+cplex.getValue(TRE[2])+","+cplex.getValue(TDE[2])+",D");
			}
			if(cplex.getValue(vRE[3])>0.01) {
				ps.println("NEBRASKA_SD2"+","+1+","+cplex.getValue(TRE[3])+","+cplex.getValue(TDE[3])+",R");
			}
			else {
				ps.println("NEBRASKA_SD2"+","+1+","+cplex.getValue(TRE[3])+","+cplex.getValue(TDE[3])+",D");
			}
			if(cplex.getValue(vRE[4])>0.01) {
				ps.println("NEBRASKA_SD3"+","+1+","+cplex.getValue(TRE[4])+","+cplex.getValue(TDE[4])+",R");
			}
			else {
				ps.println("NEBRASKA_SD3"+","+1+","+cplex.getValue(TRE[4])+","+cplex.getValue(TDE[4])+",D");
			}

			ps.println("FLOW VARS: Origin name, Origin ID, Origin State, Origin Pop, Destination name, Destination ID, Destination State, Destination Pop, Distance, Flow");
			for (int k = 0; k < countVars; k++) {
				int origin = pairs.get(k)[0];
				int destination = pairs.get(k)[1];
				if(cplex.getValue(x[k])>0.0001) {
					ps.println(data.countyList.get(origin).name+","+data.countyList.get(origin).ID+","+data.countyList.get(origin).state+","+data.countyList.get(origin).pop+","+data.countyList.get(destination).name+","+data.countyList.get(destination).ID+","+data.countyList.get(destination).state+","+data.countyList.get(destination).pop+","+data.distMatrix[origin][destination]+","+cplex.getValue(x[k]));
				}
			}
			ps.println("County name, County ID, County total votes R, County total votes D");
			for (int i = 0; i < data.numCounty; i++) {
				ps.println(data.countyList.get(i).name+","+data.countyList.get(i).ID+","+cplex.getValue(zR[i])+","+cplex.getValue(zD[i]));
			}
			 */
		}

		cplex.clearModel();
		cplex.end();	

	}

	void modelFluidSAAP(DataHandler data, double lambda, double limit, int pLimit, String party) throws IloException, FileNotFoundException {
		int M = 100000000;
		int numSim = data.countyList.get(0).rVotesSim.length;

		cplex = new IloCplex();
		cplex.setParam(IloCplex.Param.TimeLimit, 3600);
		//cplex.setOut(null);

		//Compute admissible pairs
		int countVars = 0;
		ArrayList<int[]> pairs = new ArrayList<int[]>();
		for (int i = 0; i < data.numCounty; i++) {
			for (int j = 0; j < data.numCounty; j++) {
				if(!data.countyList.get(i).state.equals(data.countyList.get(j).state) && data.distMatrix[i][j] <= limit
						&& data.countyList.get(i).pop >= 50000 && data.countyList.get(j).pop >= 50000 && (data.countyList.get(i).avgRfrac >= 0.5 || data.countyList.get(i).avgDfrac >= 0.5) ) { //Don't move inside the same state!
					int[] aux = new int[2];
					aux[0] = i;
					aux[1] = j;
					pairs.add(aux);
					countVars++;
				}
			}
		}
		System.out.println("TOTAL NUMBER OF FLOW VARS: "+countVars);

		// Define decision variables
		IloNumVar[] x = new IloNumVar[countVars];
		IloNumVar[][] zR = new IloNumVar[data.numCounty][numSim];				// Votes per county
		IloNumVar[][] zD = new IloNumVar[data.numCounty][numSim];
		IloNumVar[][] TR = new IloNumVar[data.numStates][numSim];				// Votes per state
		IloNumVar[][] TD = new IloNumVar[data.numStates][numSim];
		IloNumVar elecR[] = new IloNumVar[numSim];
		IloNumVar elecD[] = new IloNumVar[numSim];
		IloNumVar[][] vD = new IloNumVar[data.numStates][numSim];				//Binary aux
		IloNumVar[][] vR = new IloNumVar[data.numStates][numSim];
		// Special districts 0=Maine1, 1=Maine2, ... 4=Nebraska3
		IloNumVar[][] TRE = new IloNumVar[5][numSim];				// Votes per special district
		IloNumVar[][] TDE = new IloNumVar[5][numSim];
		IloNumVar[][] vDE = new IloNumVar[5][numSim];				//Binary aux
		IloNumVar[][] vRE = new IloNumVar[5][numSim];
		// Prefered party wins per simulation
		IloNumVar partyWins[] = new IloNumVar[numSim];

		for (int k = 0; k < countVars; k++) {
			int i = pairs.get(k)[0];
			x[k] = cplex.numVar(0, data.countyList.get(i).pop, IloNumVarType.Float);
		}

		for (int q = 0; q < numSim; q++) {

			for (int i = 0; i < data.numCounty; i++) {
				zR[i][q] = cplex.numVar(0, Double.MAX_VALUE, IloNumVarType.Float);
				zD[i][q] = cplex.numVar(0, Double.MAX_VALUE, IloNumVarType.Float);	
			}
			for (int s = 0; s < data.numStates; s++) {
				TR[s][q] = cplex.numVar(0, Double.MAX_VALUE, IloNumVarType.Float);
				TD[s][q] = cplex.numVar(0, Double.MAX_VALUE, IloNumVarType.Float);
				vD[s][q] = cplex.numVar(0, 1, IloNumVarType.Bool);
				vR[s][q] = cplex.numVar(0, 1, IloNumVarType.Bool);
			}
			for (int s = 0; s < 5; s++) {
				TRE[s][q] = cplex.numVar(0, Double.MAX_VALUE, IloNumVarType.Float);
				TDE[s][q] = cplex.numVar(0, Double.MAX_VALUE, IloNumVarType.Float);
				vDE[s][q] = cplex.numVar(0, 1, IloNumVarType.Bool);
				vRE[s][q] = cplex.numVar(0, 1, IloNumVarType.Bool);
			}
			elecR[q] = cplex.numVar(0, 9999, IloNumVarType.Int); 
			elecD[q] = cplex.numVar(0, 9999, IloNumVarType.Int); 
			partyWins[q] = cplex.numVar(0, 1, IloNumVarType.Bool); 
		}
		// Total people moved out of a county
		for (int i = 0; i < data.numCounty; i++) {
			IloLinearNumExpr expr = cplex.linearNumExpr();
			for (int k = 0; k < countVars; k++) {
				int origin = pairs.get(k)[0];
				if(origin == i) {
					expr.addTerm(1, x[k]);
				}
			}
			cplex.addLe(expr, lambda*data.countyList.get(i).pop);
		}
		// Total people moved into a county
		for (int i = 0; i < data.numCounty; i++) {
			IloLinearNumExpr expr = cplex.linearNumExpr();
			for (int k = 0; k < countVars; k++) {
				int dest = pairs.get(k)[1];
				if(dest == i) {
					expr.addTerm(1, x[k]);
				}
			}
			cplex.addLe(expr, lambda*data.countyList.get(i).pop);
		}
		// Total votes for a county
		for (int q = 0; q < numSim; q++) {
			IloLinearNumExpr[] exprR = new IloLinearNumExpr[data.numCounty];
			IloLinearNumExpr[] exprD = new IloLinearNumExpr[data.numCounty];
			// Initialize constraints
			for (int i = 0; i < data.numCounty; i++) {
				exprR[i] = cplex.linearNumExpr();
				exprD[i] = cplex.linearNumExpr();
			}
			// Consider movement
			for (int k = 0; k < countVars; k++) {
				int origin = pairs.get(k)[0];
				int destination = pairs.get(k)[1];
				double fracRorigin = data.countyList.get(origin).rVotesSim[q]/(data.countyList.get(origin).pop+0.0);
				double fracDorigin = data.countyList.get(origin).dVotesSim[q]/(data.countyList.get(origin).pop+0.0);
				exprR[origin].addTerm(fracRorigin, x[k]);
				exprR[destination].addTerm(-fracRorigin, x[k]);
				exprD[origin].addTerm(fracDorigin, x[k]);
				exprD[destination].addTerm(-fracDorigin, x[k]);
			}

			// Add coefficients and rhs for all counties
			for (int i = 0; i < data.numCounty; i++) {
				County c = data.countyList.get(i);
				exprR[i].addTerm(1, zR[i][q]);
				exprD[i].addTerm(1, zD[i][q]);
				cplex.addEq(exprR[i], c.rVotesSim[q]);
				cplex.addEq(exprD[i], c.dVotesSim[q]);

			}
		}

		// Total votes for a state
		for (int q = 0; q < numSim; q++) {
			for (int s = 0; s < data.numStates; s++) {
				IloLinearNumExpr exprR = cplex.linearNumExpr();
				exprR.addTerm(1, TR[s][q]);
				for (int i = 0; i < data.numCounty; i++) {
					if(data.countyList.get(i).state.equals(data.states[s])) {
						exprR.addTerm(-1, zR[i][q]);
						//System.out.println("MATCH COUNTY "+i+" "+data.states[s]);
					}
				}
				cplex.addEq(exprR, 0);

				IloLinearNumExpr exprD = cplex.linearNumExpr();
				exprD.addTerm(1, TD[s][q]);
				for (int i = 0; i < data.numCounty; i++) {
					if(data.countyList.get(i).state.equals(data.states[s])) {
						exprD.addTerm(-1, zD[i][q]);
						//System.out.println("MATCH COUNTY "+i+" "+data.states[s]);
					}
				}
				cplex.addEq(exprD, 0);

			}
			// Total votes for special districts
			for (int s = 0; s < 5; s++) {
				IloLinearNumExpr exprR = cplex.linearNumExpr();
				exprR.addTerm(1, TRE[s][q]);
				if(s==0) {
					for (int i = 0; i < data.maine1.size(); i++) {
						exprR.addTerm(-1, zR[data.maine1.get(i)][q]);
						//System.out.println("MAINE 1 COUNTY "+data.countyList.get(data.maine1.get(i)).name);
					}
					//Only half of the votes for shared county 23011
					exprR.addTerm(0.5, zR[data.countyIndex.get(23011)][q]);
				}
				else if(s==1) {
					for (int i = 0; i < data.maine2.size(); i++) {
						exprR.addTerm(-1, zR[data.maine2.get(i)][q]);
						//System.out.println("MAINE 2 COUNTY "+data.countyList.get(data.maine2.get(i)).name);
					}
					//Only half of the votes for shared county 23011
					exprR.addTerm(0.5, zR[data.countyIndex.get(23011)][q]);
				}
				else if(s==2) {
					for (int i = 0; i < data.nebraska1.size(); i++) {
						exprR.addTerm(-1, zR[data.nebraska1.get(i)][q]);
					}
					//Only half of the votes for shared county 31153
					exprR.addTerm(0.5, zR[data.countyIndex.get(31153)][q]);
				}
				else if(s==3) {
					for (int i = 0; i < data.nebraska2.size(); i++) {
						exprR.addTerm(-1, zR[data.nebraska2.get(i)][q]);
					}
					//Only half of the votes for shared county 31153
					exprR.addTerm(0.5, zR[data.countyIndex.get(31153)][q]);
				}
				else if(s==4) {
					for (int i = 0; i < data.nebraska3.size(); i++) {
						exprR.addTerm(-1, zR[data.nebraska3.get(i)][q]);
					}
				}

				cplex.addEq(exprR, 0);


				IloLinearNumExpr exprD = cplex.linearNumExpr();
				exprD.addTerm(1, TDE[s][q]);
				if(s==0) {	
					for (int i = 0; i < data.maine1.size(); i++) {
						exprD.addTerm(-1, zD[data.maine1.get(i)][q]);
						//System.out.println("MAINE 1 COUNTY "+data.countyList.get(data.maine1.get(i)).name);
					}
					//Only half of the votes for shared county 23011
					exprD.addTerm(0.5, zD[data.countyIndex.get(23011)][q]);
				}
				else if(s==1) {	
					for (int i = 0; i < data.maine2.size(); i++) {
						exprD.addTerm(-1, zD[data.maine2.get(i)][q]);
						//System.out.println("MAINE 2 COUNTY "+data.countyList.get(data.maine2.get(i)).name);
					}
					//Only half of the votes for shared county 23011
					exprD.addTerm(0.5, zD[data.countyIndex.get(23011)][q]);
				}
				else if(s==2) {	
					for (int i = 0; i < data.nebraska1.size(); i++) {
						exprD.addTerm(-1, zD[data.nebraska1.get(i)][q]);
					}
					//Only half of the votes for shared county 31153
					exprD.addTerm(0.5, zD[data.countyIndex.get(31153)][q]);
				}
				else if(s==3) {	
					for (int i = 0; i < data.nebraska2.size(); i++) {
						exprD.addTerm(-1, zD[data.nebraska2.get(i)][q]);
					}
					//Only half of the votes for shared county 31153
					exprD.addTerm(0.5, zD[data.countyIndex.get(31153)][q]);
				}
				else if(s==4) {	
					for (int i = 0; i < data.nebraska3.size(); i++) {
						exprD.addTerm(-1, zD[data.nebraska3.get(i)][q]);
					}
				}

				cplex.addEq(exprD, 0);

			}
		}
		// Elec votes count
		for (int q = 0; q < numSim; q++) {

			// Every state must be won
			for (int s = 0; s < data.numStates; s++) {
				IloLinearNumExpr expr = cplex.linearNumExpr();
				expr.addTerm(1, vR[s][q]);
				expr.addTerm(1, vD[s][q]);
				cplex.addEq(expr, 1);
			}
			// Every SD must be won
			for (int s = 0; s < 5; s++) {
				IloLinearNumExpr expr = cplex.linearNumExpr();
				expr.addTerm(1, vRE[s][q]);
				expr.addTerm(1, vDE[s][q]);
				cplex.addEq(expr, 1);
			}

			// Binary activations R
			for (int s = 0; s < data.numStates; s++) {
				IloLinearNumExpr expr = cplex.linearNumExpr();
				expr.addTerm(1, TR[s][q]);
				expr.addTerm(-1, TD[s][q]);
				expr.addTerm(-1, vR[s][q]);
				expr.addTerm(-M, vR[s][q]);
				cplex.addGe(expr, -M);
			}
			for (int s = 0; s < 5; s++) {
				IloLinearNumExpr expr = cplex.linearNumExpr();
				expr.addTerm(1, TRE[s][q]);
				expr.addTerm(-1, TDE[s][q]);
				expr.addTerm(-1, vRE[s][q]);
				expr.addTerm(-M, vRE[s][q]);
				cplex.addGe(expr, -M);
			}

			// Binary activations D
			for (int s = 0; s < data.numStates; s++) {
				IloLinearNumExpr expr = cplex.linearNumExpr();
				expr.addTerm(1, TD[s][q]);
				expr.addTerm(-1, TR[s][q]);
				expr.addTerm(-1, vD[s][q]);
				expr.addTerm(-M, vD[s][q]);
				cplex.addGe(expr, -M);
			}
			for (int s = 0; s < 5; s++) {
				IloLinearNumExpr expr = cplex.linearNumExpr();
				expr.addTerm(1, TDE[s][q]);
				expr.addTerm(-1, TRE[s][q]);
				expr.addTerm(-1, vDE[s][q]);
				expr.addTerm(-M, vDE[s][q]);
				cplex.addGe(expr, -M);
			}

			//Total Electoral votes for R
			IloLinearNumExpr exprA = cplex.linearNumExpr();
			exprA.addTerm(1, elecR[q]);
			for (int s = 0; s < data.numStates; s++) {
				exprA.addTerm(-data.eVotes[s], vR[s][q]);
			}
			for (int s = 0; s < 5; s++) {
				exprA.addTerm(-1, vRE[s][q]);
			}
			cplex.addEq(exprA, 0);

			//Total Electoral votes for D
			IloLinearNumExpr exprB = cplex.linearNumExpr();
			exprB.addTerm(1, elecD[q]);
			for (int s = 0; s < data.numStates; s++) {
				exprB.addTerm(-data.eVotes[s], vD[s][q]);
			}
			for (int s = 0; s < 5; s++) {
				exprB.addTerm(-1, vDE[s][q]);
			}
			cplex.addEq(exprB, 0);

		}	
		// Prefered party wins count
		for (int q = 0; q < numSim; q++) {
			IloLinearNumExpr expr = cplex.linearNumExpr();
			if(party.equals("R")) {expr.addTerm(1, elecR[q]);}		//Use for R
			else{ expr.addTerm(1, elecD[q]);}						//Use for D
			expr.addTerm(-270, partyWins[q]);
			cplex.addGe(expr, 0);
		}

		// Objective prob	
		IloLinearNumExpr probProxy = cplex.linearNumExpr();
		for (int q = 0; q < numSim; q++) {
			probProxy.addTerm(1, partyWins[q]);
		}
		cplex.addGe(probProxy, pLimit);

		// Objective people	
		IloLinearNumExpr objPeople = cplex.linearNumExpr();
		for (int k = 0; k < countVars; k++) {
			objPeople.addTerm(1, x[k]);
		}

		cplex.addMinimize(objPeople);
		cplex.solve();

		if(cplex.getStatus() == IloCplex.Status.Infeasible || cplex.getStatus() == IloCplex.Status.InfeasibleOrUnbounded){
			System.out.println("		MODEL FLUID SAA is infeasible!");
		}
		else if(cplex.getStatus() == IloCplex.Status.Feasible || cplex.getStatus() == IloCplex.Status.Optimal) {
			////////////////////////////////////////////////// EVAL SOL /////////////////////////////////////////////////////////////////
			xsol = new double[data.numCounty][data.numCounty];
			for (int k = 0; k < countVars; k++) {
				int origin = pairs.get(k)[0];
				int destination = pairs.get(k)[1];
				xsol[origin][destination] = cplex.getValue(x[k]);
			}
			//evalSolFluid(data, xsol, peopleLimit); 

			////////////////////////////////////////////////// PRINT RESULTS /////////////////////////////////////////////////////////////////

			System.out.println("******************************  RESULTS ******************************************");
			double tFlow = 0;
			double tDist = 0;
			for (int k = 0; k < countVars; k++) {
				tFlow += cplex.getValue(x[k]);
				int origin = pairs.get(k)[0];
				int destination = pairs.get(k)[1];
				tDist += data.distMatrix[origin][destination]*cplex.getValue(x[k]);
			}

			//System.out.println("TOTAL FLOW:,"+ tFlow);
			//System.out.println("TOTAL DISTANCE:,"+ tDist);
			flow = tFlow;

			// Export output
			java.io.PrintStream ps = new java.io.PrintStream( new java.io.FileOutputStream("Fluid_"+party+"_"+pLimit+".txt"));
			ps.println("TOTAL_FLOW:,"+ tFlow);
			ps.println("Origin_name,Origin_ID,Origin_State,Origin_Pop,Destination_name,Destination_ID,Destination_State,Destination_Pop,Distance,Flow");
			for (int k = 0; k < countVars; k++) {
				int origin = pairs.get(k)[0];
				int destination = pairs.get(k)[1];
				if(cplex.getValue(x[k])>0.0001) {
					ps.println(data.countyList.get(origin).name+","+data.countyList.get(origin).ID+","+data.countyList.get(origin).state+","+data.countyList.get(origin).pop+","+data.countyList.get(destination).name+","+data.countyList.get(destination).ID+","+data.countyList.get(destination).state+","+data.countyList.get(destination).pop+","+data.distMatrix[origin][destination]+","+cplex.getValue(x[k]));
				}
			}
		}

		cplex.clearModel();
		cplex.end();	

	}

	//Use simulations to evaluate a flow solution
	public void getStatsNoMovement(DataHandler data, double[][] xsol) throws FileNotFoundException {
		int numSim = data.countyList.get(0).rVotesSim.length;
		double[] totalVotes = new double[numSim];
		double[] totalRvotes = new double[numSim];
		double[] totalDvotes = new double[numSim];
		double[] elecVotesR = new double[numSim];
		double[] elecVotesD = new double[numSim];
		double probRwin = 0;
		double probDwin = 0;

		java.io.PrintStream ps = new java.io.PrintStream( new java.io.FileOutputStream("SimulationResults_KAMALA.txt"));
		ps.println("SimID, TotalVotes, TotalRVotes, TotalDVotes, ElecVotesR, ElecVotesD, RWins, Dwins");
		for (int k = 0; k < numSim; k++) {
			ArrayList<Double> aux =  calcVotesFluid(data, xsol, k); 
			totalRvotes[k] = aux.get(0);
			totalDvotes[k] = aux.get(1);
			elecVotesR[k] = aux.get(2);
			elecVotesD[k] = aux.get(3);
			totalVotes[k] = totalRvotes[k]+totalDvotes[k]; 
			int Rwins = 0; int Dwins = 0;
			if(elecVotesR[k] >= 270) {Rwins = 1; probRwin++;}
			if(elecVotesD[k] >= 270) {Dwins = 1; probDwin++;}
			ps.println(k+","+totalVotes[k]+","+totalRvotes[k]+","+totalDvotes[k]+","+elecVotesR[k]+","+elecVotesD[k]+","+Rwins+","+Dwins);
		}
		ps.println("PROB R WINS: "+probRwin);
		ps.println("PROB D WINS: "+probDwin);

		System.out.println("PROB R WINS: "+probRwin);
		System.out.println("PROB D WINS: "+probDwin);
	}

	private ArrayList<Double> calcVotesFluid(DataHandler data, double[][] xsol, int k) {
		ArrayList<Double> results = new ArrayList<Double>();

		double[] zR = new double[data.numCounty];				// Votes per county
		double[] zD = new double[data.numCounty];
		double[] TR = new double[data.numStates];				// Votes per state
		double[] TD = new double[data.numStates];
		double elecR = 0;
		double elecD = 0;
		// Special districts 0=Maine1, 1=Maine2, ... 4=Nebraska3
		double[] TRE = new double[5];				// Votes per special district
		double[] TDE = new double[5];

		// Count votes per county
		for (int o = 0; o < data.numCounty; o++) {
			zR[o] += data.countyList.get(o).rVotesSim[k];
			zD[o] += data.countyList.get(o).dVotesSim[k];
			double fracRorigin = data.countyList.get(o).rVotesSim[k]/(data.countyList.get(o).pop+0.0);
			double fracDorigin = data.countyList.get(o).dVotesSim[k]/(data.countyList.get(o).pop+0.0);
			for (int d = 0; d < data.numCounty; d++) {
				if(xsol[o][d]>0.001) {
					zR[o]-= fracRorigin*xsol[o][d];
					zR[d]+= fracRorigin*xsol[o][d];
					zD[o]-= fracDorigin*xsol[o][d];
					zD[d]+= fracDorigin*xsol[o][d];
				}
			}
			//System.out.println("COUNTY "+o+" "+zR[o]+" vs "+zD[o]);	
		}

		// Total votes for a state
		double totalRvotes = 0;
		double totalDvotes = 0;
		for (int s = 0; s < data.numStates; s++) {
			for (int i = 0; i < data.numCounty; i++) {
				if(data.countyList.get(i).state.equals(data.states[s])) {
					TR[s] += zR[i];
					TD[s] += zD[i];
				}
			}

			totalRvotes+=TR[s];
			totalDvotes+=TD[s];
		}

		// Total votes for special districts
		for (int s = 0; s < 5; s++) {
			if(s==0) {
				for (int i = 0; i < data.maine1.size(); i++) {
					TRE[s] += zR[data.maine1.get(i)];
					TDE[s] += zD[data.maine1.get(i)];
				}
				//Only half of the votes for shared county 23011
				TRE[s] -= 0.5*zR[data.countyIndex.get(23011)];
				TDE[s] -= 0.5*zD[data.countyIndex.get(23011)];
			}
			else if(s==1) {
				for (int i = 0; i < data.maine2.size(); i++) {
					TRE[s] += zR[data.maine2.get(i)];
					TDE[s] += zD[data.maine2.get(i)];
				}
				//Only half of the votes for shared county 23011
				TRE[s] -= 0.5*zR[data.countyIndex.get(23011)];
				TDE[s] -= 0.5*zD[data.countyIndex.get(23011)];
			}
			else if(s==2) {
				for (int i = 0; i < data.nebraska1.size(); i++) {
					TRE[s] += zR[data.nebraska1.get(i)];
					TDE[s] += zD[data.nebraska1.get(i)];
				}
				//Only half of the votes for shared county 31153
				TRE[s] -= 0.5*zR[data.countyIndex.get(31153)];
				TDE[s] -= 0.5*zD[data.countyIndex.get(31153)];
			}
			else if(s==3) {
				for (int i = 0; i < data.nebraska2.size(); i++) {
					TRE[s] += zR[data.nebraska2.get(i)];
					TDE[s] += zD[data.nebraska2.get(i)];
				}
				//Only half of the votes for shared county 31153
				TRE[s] -= 0.5*zR[data.countyIndex.get(31153)];
				TDE[s] -= 0.5*zD[data.countyIndex.get(31153)];
			}
			else if(s==4) {
				for (int i = 0; i < data.nebraska3.size(); i++) {
					TRE[s] += zR[data.nebraska3.get(i)];
					TDE[s] += zD[data.nebraska3.get(i)];
				}
			}

		}

		// Total Electoral Votes
		for (int s = 0; s < data.numStates; s++) {
			if(TR[s] > TD[s]+0.1) {
				elecR+=data.eVotes[s];
			}
			else {
				elecD+=data.eVotes[s];
			}
		}
		for (int s = 0; s < 5; s++) {
			if(TRE[s] > TDE[s]+0.1) {
				elecR+=1;
			}
			else {
				elecD+=1;
			}
		}

		results.add(totalRvotes);
		results.add(totalDvotes);
		results.add(elecR);
		results.add(elecD);



		return results;
	}

	public void updateSim(DataHandler data) throws IloException, FileNotFoundException {
		int numSim = data.countyList.get(0).rVotesSim.length;
		// Bayesian update for each state using the polling data
		for (int s = 0; s < data.numStates; s++) {
			checkStateSim(data, s, numSim);
		}


		// Get avg rates and min pops, and state pops

		for (int i = 0; i < data.numCounty; i++) {
			data.statePop[data.stateIndex.get(data.countyList.get(i).state)] += data.countyList.get(i).pop;
			// Get min
			for (int q = 0; q < numSim; q++) {
				if(data.countyList.get(i).rVotesSim[q] < data.countyList.get(i).minRvotes) {
					data.countyList.get(i).minRvotes = data.countyList.get(i).rVotesSim[q]; 
				}
				if(data.countyList.get(i).dVotesSim[q] < data.countyList.get(i).minDvotes) {
					data.countyList.get(i).minDvotes = data.countyList.get(i).dVotesSim[q]; 
				}
			}
			//System.out.println("COUNTY "+i+" MIN VOTES: "+data.countyList.get(i).minRvotes+" and "+data.countyList.get(i).minDvotes);

			//Get avg
			if(data.countyList.get(i).pop >= 1) {
				for (int q = 0; q < numSim; q++) {
					data.countyList.get(i).avgRfrac += (data.countyList.get(i).rVotesSim[q]/(data.countyList.get(i).pop+0.0));
					data.countyList.get(i).avgDfrac += (data.countyList.get(i).dVotesSim[q]/(data.countyList.get(i).pop+0.0));
				}

				data.countyList.get(i).avgRfrac = data.countyList.get(i).avgRfrac/(numSim+0.0);
				data.countyList.get(i).avgDfrac = data.countyList.get(i).avgDfrac/(numSim+0.0);

				//System.out.println("AVG FRAC "+data.countyList.get(i).avgRfrac+" and "+data.countyList.get(i).avgDfrac);

			}

		}

//Output Post Calibrated Sims
/*		java.io.PrintStream ps = new java.io.PrintStream( new java.io.FileOutputStream("OptimizationSetRace.csv"));
		ps.print("fip,party");
		for (int k = 0; k < numSim; k++) {
			ps.print(",sim_"+k);	
		}
		ps.println();
		for (int i = 0; i < data.numCounty; i++) {
			ps.print(data.countyList.get(i).ID+",DEMOCRAT");
			for (int k = 0; k < numSim; k++) {
				ps.print(","+data.countyList.get(i).dVotesSim[k]);	
			}
			ps.println();

			ps.print(data.countyList.get(i).ID+",REPUBLICAN");
			for (int k = 0; k < numSim; k++) {
				ps.print(","+data.countyList.get(i).rVotesSim[k]);	
			}
			ps.println();

		}
		 
*/

	}

	private void checkStateSim(DataHandler data, int s, int numSim) throws IloException {
		if(s==19) {//Especial model for Maine
			adjustMaine(data, numSim);
		}
		else if(s==27) {//Especial model for Nebraska
			adjustNebraska(data, numSim);
		}
		else {
			// Compute votes per state for each simulation
			double[] TR = new double[numSim];				// Votes per state
			double[] TD = new double[numSim];
			double probRwins = 0;
			double probDwins = 0;
			int rWins = 0;
			int dWins = 0;
			for (int k = 0; k < numSim; k++) {
				// Count votes per county
				double[] zR = new double[data.numCounty];		// Votes per county
				double[] zD = new double[data.numCounty];
				for (int o = 0; o < data.numCounty; o++) {
					zR[o] = data.countyList.get(o).rVotesSim[k];
					zD[o] = data.countyList.get(o).dVotesSim[k];
				}

				// Total votes for a state
				for (int i = 0; i < data.numCounty; i++) {
					if(data.countyList.get(i).state.equals(data.states[s])) {
						TR[k] += zR[i];
						TD[k] += zD[i];
					}
				}
				if(TR[k] > TD[k]) {rWins++;}
				else {dWins++;}

				if(TR[k] < 0 || TD[k]<0) {System.out.println("ERROR "+data.states[s]+" TotalRvotes "+TR[k]+" TotalDvotes "+TD[k]+" sim_"+k);}
			}
			probRwins = rWins/(numSim+0.0);
			probDwins = dWins/(numSim+0.0);

			System.out.print(data.states[s]+",Initial_RWINS_Prob:,"+probRwins+",Initial_DWINS_Prob:,"+probDwins+",PredR:,"+ data.PredR[s]+",PredD:,"+ data.PredD[s]+",");
			if(probRwins < data.PredR[s]-0.01) {
				adjustStateR(data, s, numSim, TR, TD);			
			}
			else if(probDwins < data.PredD[s]-0.01) {
				adjustStateD(data, s, numSim, TR, TD);
			}
			System.out.println();

		}


	}

	private void adjustStateR(DataHandler data, int s, int numSim, double[] R, double[] D) throws IloException {
		int M = 10000000;
		cplex = new IloCplex();
		cplex.setOut(null);
		// Define decision variables
		IloNumVar[] newR = new IloNumVar[numSim];
		IloNumVar[] newD = new IloNumVar[numSim];
		IloNumVar[] vR = new IloNumVar[numSim];
		IloNumVar alpha = cplex.numVar(0, Double.MAX_VALUE, IloNumVarType.Float);

		for (int k = 0; k < numSim; k++) {
			newR[k] = cplex.numVar(0, Double.MAX_VALUE, IloNumVarType.Float);
			newD[k] = cplex.numVar(0, Double.MAX_VALUE, IloNumVarType.Float);
			vR[k] = cplex.numVar(0, 1, IloNumVarType.Bool);
		}
		// New votes
		for (int k = 0; k < numSim; k++) {
			IloLinearNumExpr expr = cplex.linearNumExpr();
			expr.addTerm(1, newR[k]);
			expr.addTerm(-D[k], alpha);
			cplex.addEq(expr, R[k]);

			IloLinearNumExpr exprD = cplex.linearNumExpr();
			exprD.addTerm(1, newD[k]);
			exprD.addTerm(D[k], alpha);
			cplex.addEq(exprD, D[k]);
		}

		// Binary activations R
		for (int k = 0; k < numSim; k++) {
			IloLinearNumExpr expr = cplex.linearNumExpr();
			expr.addTerm(1, newR[k]);
			expr.addTerm(-1, newD[k]);
			expr.addTerm(-1, vR[k]);
			expr.addTerm(-M, vR[k]);
			cplex.addGe(expr, -M);
		}
		// Match prediction
		IloLinearNumExpr expr = cplex.linearNumExpr();
		for (int k = 0; k < numSim; k++) {
			expr.addTerm(1, vR[k]);
		}
		cplex.addGe(expr, data.PredR[s]*numSim);

		// Objective	
		IloLinearNumExpr obj = cplex.linearNumExpr();
		obj.addTerm(1, alpha);
		cplex.addMinimize(obj);
		cplex.solve();

		if(cplex.getStatus() == IloCplex.Status.Infeasible || cplex.getStatus() == IloCplex.Status.InfeasibleOrUnbounded){
			System.out.println("		Adjust State is infeasible!");
		}
		else if(cplex.getStatus() == IloCplex.Status.Optimal) {
			int newRwins = 0;
			for (int k = 0; k < numSim; k++) {
				newRwins+=(int)Math.round(cplex.getValue(vR[k]));
			}


			System.out.print("alpha:,"+cplex.getValue(alpha)+",newProbRWins:,"+newRwins/(numSim+0.0)+",newProbDWins:,"+(1-newRwins/(numSim+0.0)));
		}
		// Adjust all counties
		double alphaVal = cplex.getValue(alpha);
		for (int i = 0; i < data.numCounty; i++) {
			if(data.countyList.get(i).state.equals(data.states[s])) {
				for (int k = 0; k < numSim; k++) {
					data.countyList.get(i).rVotesSim[k] = data.countyList.get(i).rVotesSim[k] + alphaVal*data.countyList.get(i).dVotesSim[k];
					data.countyList.get(i).dVotesSim[k] = (1-alphaVal)*data.countyList.get(i).dVotesSim[k];
				}	
			}
		}


		cplex.clearModel();
		cplex.end();	



	}

	private void adjustStateRevaporate(DataHandler data, int s, int numSim, double[] R, double[] D) throws IloException {
		int M = 10000000;
		cplex = new IloCplex();
		cplex.setOut(null);
		// Define decision variables
		IloNumVar[] newR = new IloNumVar[numSim];
		IloNumVar[] newD = new IloNumVar[numSim];
		IloNumVar[] vR = new IloNumVar[numSim];
		IloNumVar alpha = cplex.numVar(0, Double.MAX_VALUE, IloNumVarType.Float);
	
		for (int k = 0; k < numSim; k++) {
			newR[k] = cplex.numVar(0, Double.MAX_VALUE, IloNumVarType.Float);
			newD[k] = cplex.numVar(0, Double.MAX_VALUE, IloNumVarType.Float);
			vR[k] = cplex.numVar(0, 1, IloNumVarType.Bool);
		}
		// New votes
		for (int k = 0; k < numSim; k++) {
			IloLinearNumExpr expr = cplex.linearNumExpr();
			expr.addTerm(1, newR[k]);
			cplex.addEq(expr, R[k]);
	
			IloLinearNumExpr exprD = cplex.linearNumExpr();
			exprD.addTerm(1, newD[k]);
			exprD.addTerm(D[k], alpha);
			cplex.addEq(exprD, D[k]);
		}
	
		// Binary activations R
		for (int k = 0; k < numSim; k++) {
			IloLinearNumExpr expr = cplex.linearNumExpr();
			expr.addTerm(1, newR[k]);
			expr.addTerm(-1, newD[k]);
			expr.addTerm(-1, vR[k]);
			expr.addTerm(-M, vR[k]);
			cplex.addGe(expr, -M);
		}
		// Match prediction
		IloLinearNumExpr expr = cplex.linearNumExpr();
		for (int k = 0; k < numSim; k++) {
			expr.addTerm(1, vR[k]);
		}
		cplex.addGe(expr, data.PredR[s]*numSim);
	
		// Objective	
		IloLinearNumExpr obj = cplex.linearNumExpr();
		obj.addTerm(1, alpha);
		cplex.addMinimize(obj);
		cplex.solve();
	
		if(cplex.getStatus() == IloCplex.Status.Infeasible || cplex.getStatus() == IloCplex.Status.InfeasibleOrUnbounded){
			System.out.println("		Adjust State is infeasible!");
		}
		else if(cplex.getStatus() == IloCplex.Status.Optimal) {
			int newRwins = 0;
			for (int k = 0; k < numSim; k++) {
				newRwins+=(int)Math.round(cplex.getValue(vR[k]));
			}
	
	
			System.out.print("alpha:,"+cplex.getValue(alpha)+",newProbRWins:,"+newRwins/(numSim+0.0)+",newProbDWins:,"+(1-newRwins/(numSim+0.0)));
		}
		// Adjust all counties
		double alphaVal = cplex.getValue(alpha);
		for (int i = 0; i < data.numCounty; i++) {
			if(data.countyList.get(i).state.equals(data.states[s])) {
				for (int k = 0; k < numSim; k++) {
					data.countyList.get(i).dVotesSim[k] = (1-alphaVal)*data.countyList.get(i).dVotesSim[k];
				}	
			}
		}
	
	
		cplex.clearModel();
		cplex.end();	
	
	
	
	}

	private void adjustStateD(DataHandler data, int s, int numSim, double[] R, double[] D) throws IloException {
		int M = 10000000;
		cplex = new IloCplex();
		cplex.setOut(null);
		// Define decision variables
		IloNumVar[] newR = new IloNumVar[numSim];
		IloNumVar[] newD = new IloNumVar[numSim];
		IloNumVar[] vD = new IloNumVar[numSim];
		IloNumVar alpha = cplex.numVar(0, Double.MAX_VALUE, IloNumVarType.Float);

		for (int k = 0; k < numSim; k++) {
			newR[k] = cplex.numVar(0, Double.MAX_VALUE, IloNumVarType.Float);
			newD[k] = cplex.numVar(0, Double.MAX_VALUE, IloNumVarType.Float);
			vD[k] = cplex.numVar(0, 1, IloNumVarType.Bool);
		}
		// New votes
		for (int k = 0; k < numSim; k++) {
			IloLinearNumExpr expr = cplex.linearNumExpr();
			expr.addTerm(1, newD[k]);
			expr.addTerm(-R[k], alpha);
			cplex.addEq(expr, D[k]);

			IloLinearNumExpr exprR = cplex.linearNumExpr();
			exprR.addTerm(1, newR[k]);
			exprR.addTerm(R[k], alpha);
			cplex.addEq(exprR, R[k]);
		}

		// Binary activations R
		for (int k = 0; k < numSim; k++) {
			IloLinearNumExpr expr = cplex.linearNumExpr();
			expr.addTerm(1, newD[k]);
			expr.addTerm(-1, newR[k]);
			expr.addTerm(-1, vD[k]);
			expr.addTerm(-M, vD[k]);
			cplex.addGe(expr, -M);
		}
		// Match prediction
		IloLinearNumExpr expr = cplex.linearNumExpr();
		for (int k = 0; k < numSim; k++) {
			expr.addTerm(1, vD[k]);
		}
		cplex.addGe(expr, data.PredD[s]*numSim);

		// Objective	
		IloLinearNumExpr obj = cplex.linearNumExpr();
		obj.addTerm(1, alpha);
		cplex.addMinimize(obj);
		cplex.solve();

		if(cplex.getStatus() == IloCplex.Status.Infeasible || cplex.getStatus() == IloCplex.Status.InfeasibleOrUnbounded){
			System.out.println("		Adjust State is infeasible!");
		}
		else if(cplex.getStatus() == IloCplex.Status.Optimal) {
			int newDwins = 0;
			for (int k = 0; k < numSim; k++) {
				newDwins+=(int)Math.round(cplex.getValue(vD[k]));
			}


			System.out.print("alpha:,"+cplex.getValue(alpha)+",newProbRWins:,"+(1-newDwins/(numSim+0.0))+",newProbDWins:,"+newDwins/(numSim+0.0));

		}
		// Adjust all counties
		double alphaVal = cplex.getValue(alpha);
		for (int i = 0; i < data.numCounty; i++) {
			if(data.countyList.get(i).state.equals(data.states[s])) {
				for (int k = 0; k < numSim; k++) {
					data.countyList.get(i).dVotesSim[k] = data.countyList.get(i).dVotesSim[k] + alphaVal*data.countyList.get(i).rVotesSim[k];
					data.countyList.get(i).rVotesSim[k] = (1-alphaVal)*data.countyList.get(i).rVotesSim[k];
				}	
			}
		}


		cplex.clearModel();
		cplex.end();	



	}

	private void adjustStateDevaporate(DataHandler data, int s, int numSim, double[] R, double[] D) throws IloException {
		int M = 10000000;
		cplex = new IloCplex();
		cplex.setOut(null);
		// Define decision variables
		IloNumVar[] newR = new IloNumVar[numSim];
		IloNumVar[] newD = new IloNumVar[numSim];
		IloNumVar[] vD = new IloNumVar[numSim];
		IloNumVar alpha = cplex.numVar(0, Double.MAX_VALUE, IloNumVarType.Float);
	
		for (int k = 0; k < numSim; k++) {
			newR[k] = cplex.numVar(0, Double.MAX_VALUE, IloNumVarType.Float);
			newD[k] = cplex.numVar(0, Double.MAX_VALUE, IloNumVarType.Float);
			vD[k] = cplex.numVar(0, 1, IloNumVarType.Bool);
		}
		// New votes
		for (int k = 0; k < numSim; k++) {
			IloLinearNumExpr expr = cplex.linearNumExpr();
			expr.addTerm(1, newD[k]);
			cplex.addEq(expr, D[k]);
	
			IloLinearNumExpr exprR = cplex.linearNumExpr();
			exprR.addTerm(1, newR[k]);
			exprR.addTerm(R[k], alpha);
			cplex.addEq(exprR, R[k]);
		}
	
		// Binary activations R
		for (int k = 0; k < numSim; k++) {
			IloLinearNumExpr expr = cplex.linearNumExpr();
			expr.addTerm(1, newD[k]);
			expr.addTerm(-1, newR[k]);
			expr.addTerm(-1, vD[k]);
			expr.addTerm(-M, vD[k]);
			cplex.addGe(expr, -M);
		}
		// Match prediction
		IloLinearNumExpr expr = cplex.linearNumExpr();
		for (int k = 0; k < numSim; k++) {
			expr.addTerm(1, vD[k]);
		}
		cplex.addGe(expr, data.PredD[s]*numSim);
	
		// Objective	
		IloLinearNumExpr obj = cplex.linearNumExpr();
		obj.addTerm(1, alpha);
		cplex.addMinimize(obj);
		cplex.solve();
	
		if(cplex.getStatus() == IloCplex.Status.Infeasible || cplex.getStatus() == IloCplex.Status.InfeasibleOrUnbounded){
			System.out.println("		Adjust State is infeasible!");
		}
		else if(cplex.getStatus() == IloCplex.Status.Optimal) {
			int newDwins = 0;
			for (int k = 0; k < numSim; k++) {
				newDwins+=(int)Math.round(cplex.getValue(vD[k]));
			}
	
	
			System.out.print("alpha:,"+cplex.getValue(alpha)+",newProbRWins:,"+(1-newDwins/(numSim+0.0))+",newProbDWins:,"+newDwins/(numSim+0.0));
	
		}
		// Adjust all counties
		double alphaVal = cplex.getValue(alpha);
		for (int i = 0; i < data.numCounty; i++) {
			if(data.countyList.get(i).state.equals(data.states[s])) {
				for (int k = 0; k < numSim; k++) {
					data.countyList.get(i).rVotesSim[k] = (1-alphaVal)*data.countyList.get(i).rVotesSim[k];
				}	
			}
		}
	
	
		cplex.clearModel();
		cplex.end();	
	
	
	
	}

	private void adjustMaine(DataHandler data, int numSim) throws IloException {
		// Get the total number of votes for the especial districts
		double[][] TRE = new double[2][numSim];
		double[][] TDE = new double[2][numSim];
		double probRwins0 = 0;
		double probDwins0 = 0;
		int rWins0 = 0;
		int dWins0 = 0;
		double probRwins1 = 0;
		double probDwins1 = 0;
		int rWins1 = 0;
		int dWins1 = 0;
		int pop0 = 0;
		int pop1 = 0;
		int rWins = 0;
		int dWins = 0;
		double probRwins = 0;
		double probDwins = 0;

		for (int i = 0; i < data.maine1.size(); i++) {
			pop0+=data.countyList.get(data.maine1.get(i)).pop;
		}
		for (int i = 0; i < data.maine2.size(); i++) {
			pop1 += data.countyList.get(data.maine2.get(i)).pop;
		}

		for (int q = 0; q < numSim; q++) {
			for (int s = 0; s <= 1; s++) {
				if(s==0) {
					for (int i = 0; i < data.maine1.size(); i++) {
						TRE[s][q] += data.countyList.get(data.maine1.get(i)).rVotesSim[q];
					}
					//Only half of the votes for shared county 23011
					TRE[s][q] -= 0.5*data.countyList.get(data.countyIndex.get(23011)).rVotesSim[q];
				}
				else if(s==1) {
					for (int i = 0; i < data.maine2.size(); i++) {
						TRE[s][q] += data.countyList.get(data.maine2.get(i)).rVotesSim[q];
					}
					//Only half of the votes for shared county 23011
					TRE[s][q] -= 0.5*data.countyList.get(data.countyIndex.get(23011)).rVotesSim[q];
				}

				if(s==0) {
					for (int i = 0; i < data.maine1.size(); i++) {
						TDE[s][q] += data.countyList.get(data.maine1.get(i)).dVotesSim[q];
					}
					//Only half of the votes for shared county 23011
					TDE[s][q]-=0.5*data.countyList.get(data.countyIndex.get(23011)).dVotesSim[q];
				}
				else if(s==1) {
					for (int i = 0; i < data.maine2.size(); i++) {
						TDE[s][q] +=data.countyList.get(data.maine2.get(i)).dVotesSim[q];
					}
					//Only half of the votes for shared county 23011
					TDE[s][q]-=0.5*data.countyList.get(data.countyIndex.get(23011)).dVotesSim[q];
				}
			}

			if(TRE[0][q] > TDE[0][q]) {rWins0++;}
			else {dWins0++;}
			if(TRE[1][q] > TDE[1][q]) {rWins1++;}
			else {dWins1++;}
			if(TRE[0][q]+TRE[1][q] > TDE[0][q]+TDE[1][q]) {rWins++;}
			else {dWins++;}

		}
		probRwins = rWins/(numSim+0.0);
		probDwins = dWins/(numSim+0.0);
		probRwins0 = rWins0/(numSim+0.0);
		probDwins0 = dWins0/(numSim+0.0);
		probRwins1 = rWins1/(numSim+0.0);
		probDwins1 = dWins1/(numSim+0.0);

		int M = 10000000;
		cplex = new IloCplex();
		cplex.setParam(IloCplex.Param.TimeLimit, 300);
		cplex.setOut(null);

		// Define decision variables
		IloNumVar[][] newR = new IloNumVar[2][numSim];
		IloNumVar[][] newD = new IloNumVar[2][numSim];
		IloNumVar[][] vD = new IloNumVar[3][numSim];		//Index 2 is the complete state
		IloNumVar[][] vR = new IloNumVar[3][numSim];
		IloNumVar[] alphaR = new IloNumVar[2];
		IloNumVar [] alphaD = new IloNumVar[2];

		for (int s = 0; s <= 1; s++) {
			alphaR[s] = cplex.numVar(0, Double.MAX_VALUE, IloNumVarType.Float); //Flip from D to R
			alphaD[s] = cplex.numVar(0, Double.MAX_VALUE, IloNumVarType.Float); //Flip from R to D
		}

		for (int k = 0; k < numSim; k++) {
			for (int s = 0; s <= 1; s++) {
				newR[s][k] = cplex.numVar(0, Double.MAX_VALUE, IloNumVarType.Float);
				newD[s][k] = cplex.numVar(0, Double.MAX_VALUE, IloNumVarType.Float);
				vD[s][k] = cplex.numVar(0, 1, IloNumVarType.Bool);
				vR[s][k] = cplex.numVar(0, 1, IloNumVarType.Bool);
			}
			vD[2][k] = cplex.numVar(0, 1, IloNumVarType.Bool);
			vR[2][k] = cplex.numVar(0, 1, IloNumVarType.Bool);
		}
		// New votes
		for (int k = 0; k < numSim; k++) {
			for (int s = 0; s <= 1; s++) {
				IloLinearNumExpr expr = cplex.linearNumExpr();
				expr.addTerm(1, newD[s][k]);
				expr.addTerm(-1, alphaD[s]);
				expr.addTerm(1, alphaR[s]);
				cplex.addEq(expr, TDE[s][k]);

				IloLinearNumExpr exprR = cplex.linearNumExpr();
				exprR.addTerm(1, newR[s][k]);
				exprR.addTerm(-1, alphaR[s]);
				exprR.addTerm(1, alphaD[s]);
				cplex.addEq(exprR, TRE[s][k]);
			}
		}

		// Binary activations D
		for (int k = 0; k < numSim; k++) {
			for (int s = 0; s <= 1; s++) {
				IloLinearNumExpr expr = cplex.linearNumExpr();
				expr.addTerm(1, newD[s][k]);
				expr.addTerm(-1, newR[s][k]);
				expr.addTerm(-1, vD[s][k]);
				expr.addTerm(-M, vD[s][k]);
				cplex.addGe(expr, -M);
			}
			// For the complete state
			IloLinearNumExpr exprSD = cplex.linearNumExpr();
			exprSD.addTerm(1, newD[0][k]);
			exprSD.addTerm(1, newD[1][k]);
			exprSD.addTerm(-1, newR[0][k]);
			exprSD.addTerm(-1, newR[1][k]);
			exprSD.addTerm(-1, vD[2][k]);
			exprSD.addTerm(-M, vD[2][k]);
			cplex.addGe(exprSD, -M);
		}
		// Binary activations R
		for (int k = 0; k < numSim; k++) {
			for (int s = 0; s <= 1; s++) {
				IloLinearNumExpr expr = cplex.linearNumExpr();
				expr.addTerm(1, newR[s][k]);
				expr.addTerm(-1, newD[s][k]);
				expr.addTerm(-1, vR[s][k]);
				expr.addTerm(-M, vR[s][k]);
				cplex.addGe(expr, -M);
			}
			// For the complete state
			IloLinearNumExpr exprSR = cplex.linearNumExpr();
			exprSR.addTerm(1, newR[0][k]);
			exprSR.addTerm(1, newR[1][k]);
			exprSR.addTerm(-1, newD[0][k]);
			exprSR.addTerm(-1, newD[1][k]);
			exprSR.addTerm(-1, vR[2][k]);
			exprSR.addTerm(-M, vR[2][k]);
			cplex.addGe(exprSR, -M);
		}
		//1 winner constraints
		for (int k = 0; k < numSim; k++) {
			for (int s = 0; s <= 2; s++) {
				IloLinearNumExpr expr = cplex.linearNumExpr();
				expr.addTerm(1, vD[s][k]);
				expr.addTerm(1, vR[s][k]);
				cplex.addEq(expr, 1);
			}
		}

		// Match prediction 0
		IloLinearNumExpr expr = cplex.linearNumExpr();
		for (int k = 0; k < numSim; k++) {
			expr.addTerm(1, vD[0][k]);
		}
		cplex.addGe(expr, data.NSPredDE[0]*numSim-100);
		cplex.addLe(expr, data.NSPredDE[0]*numSim+100);

		IloLinearNumExpr exprR = cplex.linearNumExpr();
		for (int k = 0; k < numSim; k++) {
			exprR.addTerm(1, vR[0][k]);
		}
		cplex.addGe(exprR, data.NSPredRE[0]*numSim-100);
		cplex.addLe(exprR, data.NSPredRE[0]*numSim+100);

		// Match prediction 1
		IloLinearNumExpr expr1 = cplex.linearNumExpr();
		for (int k = 0; k < numSim; k++) {
			expr1.addTerm(1, vD[1][k]);
		}
		cplex.addGe(expr1, data.NSPredDE[1]*numSim-100);
		cplex.addLe(expr1, data.NSPredDE[1]*numSim+100);

		IloLinearNumExpr exprR1 = cplex.linearNumExpr();
		for (int k = 0; k < numSim; k++) {
			exprR1.addTerm(1, vR[1][k]);
		}
		cplex.addGe(exprR1, data.NSPredRE[1]*numSim-100);
		cplex.addLe(exprR1, data.NSPredRE[1]*numSim+100);

		// Match prediction STATE
		IloLinearNumExpr exprS = cplex.linearNumExpr();
		for (int k = 0; k < numSim; k++) {
			exprS.addTerm(1, vD[2][k]);
		}
		cplex.addGe(exprS, data.PredD[19]*numSim-25);
		cplex.addLe(exprS, data.PredD[19]*numSim+25);

		IloLinearNumExpr exprSR1 = cplex.linearNumExpr();
		for (int k = 0; k < numSim; k++) {
			exprSR1.addTerm(1, vR[2][k]);
		}
		cplex.addGe(exprSR1, data.PredR[19]*numSim-25);
		cplex.addLe(exprSR1, data.PredR[19]*numSim+25);


		// Objective	
		IloLinearNumExpr obj = cplex.linearNumExpr();
		obj.addTerm(1, alphaR[0]);
		obj.addTerm(1, alphaR[1]);
		obj.addTerm(1, alphaD[0]);
		obj.addTerm(1, alphaD[1]);
		cplex.addMinimize(obj);
		cplex.solve();

		if(cplex.getStatus() == IloCplex.Status.Infeasible || cplex.getStatus() == IloCplex.Status.InfeasibleOrUnbounded){
			System.out.println("		Adjust MAINE is infeasible!");
		}
		else {
			int newDwins0 = 0;
			int newRwins0 = 0;
			int newDwins1 = 0;
			int newRwins1 = 0;
			int newDS = 0;
			int newRS = 0;


			for (int k = 0; k < numSim; k++) {
				newDwins0+=(int)Math.round(cplex.getValue(vD[0][k]));
				newRwins0+=(int)Math.round(cplex.getValue(vR[0][k]));
				newDwins1+=(int)Math.round(cplex.getValue(vD[1][k]));
				newRwins1+=(int)Math.round(cplex.getValue(vR[1][k]));
				newDS+=(int)Math.round(cplex.getValue(vD[2][k]));
				newRS+=(int)Math.round(cplex.getValue(vR[2][k]));
			}

			System.out.println("MAINE"+",Initial_RWINS_Prob:,"+probRwins+",Initial_DWINS_Prob:,"+probDwins+",NSPredR:,"+ data.PredR[19]+",NSPredD:,"+ data.PredD[19]+",alpha:,0,newProbRwins:,"+newRS/(numSim+0.0)+",newDprob:,"+newDS/(numSim+0.0));
			System.out.println("MAINE_D1"+",Initial_RWINS_Prob:,"+probRwins0+",Initial_DWINS_Prob:,"+probDwins0+",NSPredR:,"+ data.NSPredRE[0]+",NSPredD:,"+ data.NSPredDE[0]+",alpha:,"+(cplex.getValue(alphaD[0])/(pop0+0.0)+cplex.getValue(alphaR[0])/(pop0+0.0)));
			System.out.println("MAINE_D2"+",Initial_RWINS_Prob:,"+probRwins1+",Initial_DWINS_Prob:,"+probDwins1+",NSPredR:,"+ data.NSPredRE[1]+",NSPredD:,"+ data.NSPredDE[1]+",alpha:,"+(cplex.getValue(alphaD[1])/(pop1+0.0)+cplex.getValue(alphaR[1])/(pop1+0.0)));

		}

		// Adjust all counties
		double alphaValD0 = cplex.getValue(alphaD[0])/(pop0+0.0);
		double alphaValD1 = cplex.getValue(alphaD[1])/(pop1+0.0);
		double alphaValR0 = cplex.getValue(alphaR[0])/(pop0+0.0);
		double alphaValR1 = cplex.getValue(alphaR[1])/(pop1+0.0);

		for (int k = 0; k < numSim; k++) {
			for (int i = 0; i < data.maine1.size(); i++) {
				data.countyList.get(data.maine1.get(i)).dVotesSim[k] = data.countyList.get(data.maine1.get(i)).dVotesSim[k] + alphaValD0*data.countyList.get(data.maine1.get(i)).rVotesSim[k];
				data.countyList.get(data.maine1.get(i)).rVotesSim[k] = (1-alphaValD0)*data.countyList.get(data.maine1.get(i)).rVotesSim[k];
				data.countyList.get(data.maine1.get(i)).rVotesSim[k] = data.countyList.get(data.maine1.get(i)).rVotesSim[k] + alphaValR0*data.countyList.get(data.maine1.get(i)).dVotesSim[k];
				data.countyList.get(data.maine1.get(i)).dVotesSim[k] = (1-alphaValR0)*data.countyList.get(data.maine1.get(i)).dVotesSim[k];
			}
			for (int i = 0; i < data.maine2.size(); i++) {
				data.countyList.get(data.maine2.get(i)).dVotesSim[k] = data.countyList.get(data.maine2.get(i)).dVotesSim[k] + alphaValD1*data.countyList.get(data.maine2.get(i)).rVotesSim[k];
				data.countyList.get(data.maine2.get(i)).rVotesSim[k] = (1-alphaValD1)*data.countyList.get(data.maine2.get(i)).rVotesSim[k];
				data.countyList.get(data.maine2.get(i)).rVotesSim[k] = data.countyList.get(data.maine2.get(i)).rVotesSim[k] + alphaValR1*data.countyList.get(data.maine2.get(i)).dVotesSim[k];
				data.countyList.get(data.maine2.get(i)).dVotesSim[k] = (1-alphaValR1)*data.countyList.get(data.maine2.get(i)).dVotesSim[k];
			}
		}


		cplex.clearModel();
		cplex.end();	



	}

	private void adjustMaineEvaporate(DataHandler data, int numSim) throws IloException {
		// Get the total number of votes for the especial districts
		double[][] TRE = new double[2][numSim];
		double[][] TDE = new double[2][numSim];
		double probRwins0 = 0;
		double probDwins0 = 0;
		int rWins0 = 0;
		int dWins0 = 0;
		double probRwins1 = 0;
		double probDwins1 = 0;
		int rWins1 = 0;
		int dWins1 = 0;
		int pop0 = 0;
		int pop1 = 0;
		int rWins = 0;
		int dWins = 0;
		double probRwins = 0;
		double probDwins = 0;
	
		for (int i = 0; i < data.maine1.size(); i++) {
			pop0+=data.countyList.get(data.maine1.get(i)).pop;
		}
		for (int i = 0; i < data.maine2.size(); i++) {
			pop1 += data.countyList.get(data.maine2.get(i)).pop;
		}
	
		for (int q = 0; q < numSim; q++) {
			for (int s = 0; s <= 1; s++) {
				if(s==0) {
					for (int i = 0; i < data.maine1.size(); i++) {
						TRE[s][q] += data.countyList.get(data.maine1.get(i)).rVotesSim[q];
					}
					//Only half of the votes for shared county 23011
					TRE[s][q] -= 0.5*data.countyList.get(data.countyIndex.get(23011)).rVotesSim[q];
				}
				else if(s==1) {
					for (int i = 0; i < data.maine2.size(); i++) {
						TRE[s][q] += data.countyList.get(data.maine2.get(i)).rVotesSim[q];
					}
					//Only half of the votes for shared county 23011
					TRE[s][q] -= 0.5*data.countyList.get(data.countyIndex.get(23011)).rVotesSim[q];
				}
	
				if(s==0) {
					for (int i = 0; i < data.maine1.size(); i++) {
						TDE[s][q] += data.countyList.get(data.maine1.get(i)).dVotesSim[q];
					}
					//Only half of the votes for shared county 23011
					TDE[s][q]-=0.5*data.countyList.get(data.countyIndex.get(23011)).dVotesSim[q];
				}
				else if(s==1) {
					for (int i = 0; i < data.maine2.size(); i++) {
						TDE[s][q] +=data.countyList.get(data.maine2.get(i)).dVotesSim[q];
					}
					//Only half of the votes for shared county 23011
					TDE[s][q]-=0.5*data.countyList.get(data.countyIndex.get(23011)).dVotesSim[q];
				}
			}
	
			if(TRE[0][q] > TDE[0][q]) {rWins0++;}
			else {dWins0++;}
			if(TRE[1][q] > TDE[1][q]) {rWins1++;}
			else {dWins1++;}
			if(TRE[0][q]+TRE[1][q] > TDE[0][q]+TDE[1][q]) {rWins++;}
			else {dWins++;}
	
		}
		probRwins = rWins/(numSim+0.0);
		probDwins = dWins/(numSim+0.0);
		probRwins0 = rWins0/(numSim+0.0);
		probDwins0 = dWins0/(numSim+0.0);
		probRwins1 = rWins1/(numSim+0.0);
		probDwins1 = dWins1/(numSim+0.0);
	
		int M = 10000000;
		cplex = new IloCplex();
		cplex.setParam(IloCplex.Param.TimeLimit, 300);
		cplex.setOut(null);
	
		// Define decision variables
		IloNumVar[][] newR = new IloNumVar[2][numSim];
		IloNumVar[][] newD = new IloNumVar[2][numSim];
		IloNumVar[][] vD = new IloNumVar[3][numSim];		//Index 2 is the complete state
		IloNumVar[][] vR = new IloNumVar[3][numSim];
		IloNumVar[] alphaR = new IloNumVar[2];
		IloNumVar [] alphaD = new IloNumVar[2];
	
		for (int s = 0; s <= 1; s++) {
			alphaR[s] = cplex.numVar(0, Double.MAX_VALUE, IloNumVarType.Float); //Flip from D to R
			alphaD[s] = cplex.numVar(0, Double.MAX_VALUE, IloNumVarType.Float); //Flip from R to D
		}
	
		for (int k = 0; k < numSim; k++) {
			for (int s = 0; s <= 1; s++) {
				newR[s][k] = cplex.numVar(0, Double.MAX_VALUE, IloNumVarType.Float);
				newD[s][k] = cplex.numVar(0, Double.MAX_VALUE, IloNumVarType.Float);
				vD[s][k] = cplex.numVar(0, 1, IloNumVarType.Bool);
				vR[s][k] = cplex.numVar(0, 1, IloNumVarType.Bool);
			}
			vD[2][k] = cplex.numVar(0, 1, IloNumVarType.Bool);
			vR[2][k] = cplex.numVar(0, 1, IloNumVarType.Bool);
		}
		// New votes
		for (int k = 0; k < numSim; k++) {
			for (int s = 0; s <= 1; s++) {
				IloLinearNumExpr expr = cplex.linearNumExpr();
				expr.addTerm(1, newD[s][k]);
				expr.addTerm(1, alphaR[s]);
				cplex.addEq(expr, TDE[s][k]);
	
				IloLinearNumExpr exprR = cplex.linearNumExpr();
				exprR.addTerm(1, newR[s][k]);
				exprR.addTerm(1, alphaD[s]);
				cplex.addEq(exprR, TRE[s][k]);
			}
		}
	
		// Binary activations D
		for (int k = 0; k < numSim; k++) {
			for (int s = 0; s <= 1; s++) {
				IloLinearNumExpr expr = cplex.linearNumExpr();
				expr.addTerm(1, newD[s][k]);
				expr.addTerm(-1, newR[s][k]);
				expr.addTerm(-1, vD[s][k]);
				expr.addTerm(-M, vD[s][k]);
				cplex.addGe(expr, -M);
			}
			// For the complete state
			IloLinearNumExpr exprSD = cplex.linearNumExpr();
			exprSD.addTerm(1, newD[0][k]);
			exprSD.addTerm(1, newD[1][k]);
			exprSD.addTerm(-1, newR[0][k]);
			exprSD.addTerm(-1, newR[1][k]);
			exprSD.addTerm(-1, vD[2][k]);
			exprSD.addTerm(-M, vD[2][k]);
			cplex.addGe(exprSD, -M);
		}
		// Binary activations R
		for (int k = 0; k < numSim; k++) {
			for (int s = 0; s <= 1; s++) {
				IloLinearNumExpr expr = cplex.linearNumExpr();
				expr.addTerm(1, newR[s][k]);
				expr.addTerm(-1, newD[s][k]);
				expr.addTerm(-1, vR[s][k]);
				expr.addTerm(-M, vR[s][k]);
				cplex.addGe(expr, -M);
			}
			// For the complete state
			IloLinearNumExpr exprSR = cplex.linearNumExpr();
			exprSR.addTerm(1, newR[0][k]);
			exprSR.addTerm(1, newR[1][k]);
			exprSR.addTerm(-1, newD[0][k]);
			exprSR.addTerm(-1, newD[1][k]);
			exprSR.addTerm(-1, vR[2][k]);
			exprSR.addTerm(-M, vR[2][k]);
			cplex.addGe(exprSR, -M);
		}
		//1 winner constraints
		for (int k = 0; k < numSim; k++) {
			for (int s = 0; s <= 2; s++) {
				IloLinearNumExpr expr = cplex.linearNumExpr();
				expr.addTerm(1, vD[s][k]);
				expr.addTerm(1, vR[s][k]);
				cplex.addEq(expr, 1);
			}
		}
	
		// Match prediction 0
		IloLinearNumExpr expr = cplex.linearNumExpr();
		for (int k = 0; k < numSim; k++) {
			expr.addTerm(1, vD[0][k]);
		}
		cplex.addGe(expr, data.NSPredDE[0]*numSim-100);
		cplex.addLe(expr, data.NSPredDE[0]*numSim+100);
	
		IloLinearNumExpr exprR = cplex.linearNumExpr();
		for (int k = 0; k < numSim; k++) {
			exprR.addTerm(1, vR[0][k]);
		}
		cplex.addGe(exprR, data.NSPredRE[0]*numSim-100);
		cplex.addLe(exprR, data.NSPredRE[0]*numSim+100);
	
		// Match prediction 1
		IloLinearNumExpr expr1 = cplex.linearNumExpr();
		for (int k = 0; k < numSim; k++) {
			expr1.addTerm(1, vD[1][k]);
		}
		cplex.addGe(expr1, data.NSPredDE[1]*numSim-100);
		cplex.addLe(expr1, data.NSPredDE[1]*numSim+100);
	
		IloLinearNumExpr exprR1 = cplex.linearNumExpr();
		for (int k = 0; k < numSim; k++) {
			exprR1.addTerm(1, vR[1][k]);
		}
		cplex.addGe(exprR1, data.NSPredRE[1]*numSim-100);
		cplex.addLe(exprR1, data.NSPredRE[1]*numSim+100);
	
		// Match prediction STATE
		IloLinearNumExpr exprS = cplex.linearNumExpr();
		for (int k = 0; k < numSim; k++) {
			exprS.addTerm(1, vD[2][k]);
		}
		cplex.addGe(exprS, data.PredD[19]*numSim-25);
		cplex.addLe(exprS, data.PredD[19]*numSim+25);
	
		IloLinearNumExpr exprSR1 = cplex.linearNumExpr();
		for (int k = 0; k < numSim; k++) {
			exprSR1.addTerm(1, vR[2][k]);
		}
		cplex.addGe(exprSR1, data.PredR[19]*numSim-25);
		cplex.addLe(exprSR1, data.PredR[19]*numSim+25);
	
	
		// Objective	
		IloLinearNumExpr obj = cplex.linearNumExpr();
		obj.addTerm(1, alphaR[0]);
		obj.addTerm(1, alphaR[1]);
		obj.addTerm(1, alphaD[0]);
		obj.addTerm(1, alphaD[1]);
		cplex.addMinimize(obj);
		cplex.solve();
	
		if(cplex.getStatus() == IloCplex.Status.Infeasible || cplex.getStatus() == IloCplex.Status.InfeasibleOrUnbounded){
			System.out.println("		Adjust MAINE is infeasible!");
		}
		else {
			int newDwins0 = 0;
			int newRwins0 = 0;
			int newDwins1 = 0;
			int newRwins1 = 0;
			int newDS = 0;
			int newRS = 0;
	
	
			for (int k = 0; k < numSim; k++) {
				newDwins0+=(int)Math.round(cplex.getValue(vD[0][k]));
				newRwins0+=(int)Math.round(cplex.getValue(vR[0][k]));
				newDwins1+=(int)Math.round(cplex.getValue(vD[1][k]));
				newRwins1+=(int)Math.round(cplex.getValue(vR[1][k]));
				newDS+=(int)Math.round(cplex.getValue(vD[2][k]));
				newRS+=(int)Math.round(cplex.getValue(vR[2][k]));
			}
	
			System.out.println("MAINE"+",Initial_RWINS_Prob:,"+probRwins+",Initial_DWINS_Prob:,"+probDwins+",NSPredR:,"+ data.PredR[19]+",NSPredD:,"+ data.PredD[19]+",alpha:,0,newProbRwins:,"+newRS/(numSim+0.0)+",newDprob:,"+newDS/(numSim+0.0));
			System.out.println("MAINE_D1"+",Initial_RWINS_Prob:,"+probRwins0+",Initial_DWINS_Prob:,"+probDwins0+",NSPredR:,"+ data.NSPredRE[0]+",NSPredD:,"+ data.NSPredDE[0]+",alpha:,"+(cplex.getValue(alphaD[0])/(pop0+0.0)+cplex.getValue(alphaR[0])/(pop0+0.0)));
			System.out.println("MAINE_D2"+",Initial_RWINS_Prob:,"+probRwins1+",Initial_DWINS_Prob:,"+probDwins1+",NSPredR:,"+ data.NSPredRE[1]+",NSPredD:,"+ data.NSPredDE[1]+",alpha:,"+(cplex.getValue(alphaD[1])/(pop1+0.0)+cplex.getValue(alphaR[1])/(pop1+0.0)));
	
		}
	
		// Adjust all counties
		double alphaValD0 = cplex.getValue(alphaD[0])/(pop0+0.0);
		double alphaValD1 = cplex.getValue(alphaD[1])/(pop1+0.0);
		double alphaValR0 = cplex.getValue(alphaR[0])/(pop0+0.0);
		double alphaValR1 = cplex.getValue(alphaR[1])/(pop1+0.0);
	
		for (int k = 0; k < numSim; k++) {
			for (int i = 0; i < data.maine1.size(); i++) {
				data.countyList.get(data.maine1.get(i)).rVotesSim[k] = (1-alphaValD0)*data.countyList.get(data.maine1.get(i)).rVotesSim[k];
				data.countyList.get(data.maine1.get(i)).dVotesSim[k] = (1-alphaValR0)*data.countyList.get(data.maine1.get(i)).dVotesSim[k];
			}
			for (int i = 0; i < data.maine2.size(); i++) {
				data.countyList.get(data.maine2.get(i)).rVotesSim[k] = (1-alphaValD1)*data.countyList.get(data.maine2.get(i)).rVotesSim[k];
				data.countyList.get(data.maine2.get(i)).dVotesSim[k] = (1-alphaValR1)*data.countyList.get(data.maine2.get(i)).dVotesSim[k];
			}
		}
	
	
		cplex.clearModel();
		cplex.end();	
	
	
	
	}

	private void adjustNebraska(DataHandler data, int numSim) throws IloException {
		// Get the total number of votes for the especial districts
		double[][] TRE = new double[3][numSim];
		double[][] TDE = new double[3][numSim];
		double probRwins0 = 0;
		double probDwins0 = 0;
		int rWins0 = 0;
		int dWins0 = 0;
		double probRwins1 = 0;
		double probDwins1 = 0;
		int rWins1 = 0;
		int dWins1 = 0;
		double probRwins2 = 0;
		double probDwins2 = 0;
		int rWins2 = 0;
		int dWins2 = 0;
		int pop0 = 0;
		int pop1 = 0;
		int pop2 = 0;
		int rWins = 0;
		int dWins = 0;
		double probRwins = 0;
		double probDwins = 0;

		for (int i = 0; i < data.nebraska1.size(); i++) {
			pop0+=data.countyList.get(data.nebraska1.get(i)).pop;
		}
		for (int i = 0; i < data.nebraska2.size(); i++) {
			pop1 += data.countyList.get(data.nebraska2.get(i)).pop;
		}
		for (int i = 0; i < data.nebraska3.size(); i++) {
			pop2 += data.countyList.get(data.nebraska3.get(i)).pop;
		}

		for (int q = 0; q < numSim; q++) {
			for (int s = 0; s <= 2; s++) {
				if(s==0) {
					for (int i = 0; i < data.nebraska1.size(); i++) {
						TRE[s][q] += data.countyList.get(data.nebraska1.get(i)).rVotesSim[q];
					}
					//Only half of the votes for shared county 31153
					TRE[s][q] -= 0.5*data.countyList.get(data.countyIndex.get(31153)).rVotesSim[q];
				}
				else if(s==1) {
					for (int i = 0; i < data.nebraska2.size(); i++) {
						TRE[s][q] += data.countyList.get(data.nebraska2.get(i)).rVotesSim[q];
					}
					//Only half of the votes for shared county 31153
					TRE[s][q] -= 0.5*data.countyList.get(data.countyIndex.get(31153)).rVotesSim[q];
				}
				else if(s==2) {
					for (int i = 0; i < data.nebraska3.size(); i++) {
						TRE[s][q] += data.countyList.get(data.nebraska3.get(i)).rVotesSim[q];
					}
				}

				if(s==0) {
					for (int i = 0; i < data.nebraska1.size(); i++) {
						TDE[s][q] += data.countyList.get(data.nebraska1.get(i)).dVotesSim[q];
					}
					//Only half of the votes for shared county 31153
					TDE[s][q]-=0.5*data.countyList.get(data.countyIndex.get(31153)).dVotesSim[q];
				}
				else if(s==1) {
					for (int i = 0; i < data.nebraska2.size(); i++) {
						TDE[s][q] +=data.countyList.get(data.nebraska2.get(i)).dVotesSim[q];
					}
					//Only half of the votes for shared county 31153
					TDE[s][q]-=0.5*data.countyList.get(data.countyIndex.get(31153)).dVotesSim[q];
				}
				else if(s==2) {
					for (int i = 0; i < data.nebraska3.size(); i++) {
						TDE[s][q] +=data.countyList.get(data.nebraska3.get(i)).dVotesSim[q];
					}
				}
			}

			if(TRE[0][q] > TDE[0][q]) {rWins0++;}
			else {dWins0++;}
			if(TRE[1][q] > TDE[1][q]) {rWins1++;}
			else {dWins1++;}
			if(TRE[2][q] > TDE[2][q]) {rWins2++;}
			else {dWins2++;}
			if(TRE[0][q]+TRE[1][q]+TRE[2][q] > TDE[0][q]+TDE[1][q]+TDE[2][q]) {rWins++;}
			else {dWins++;}
		}


		probRwins = rWins/(numSim+0.0);
		probDwins = dWins/(numSim+0.0);
		probRwins0 = rWins0/(numSim+0.0);
		probDwins0 = dWins0/(numSim+0.0);
		probRwins1 = rWins1/(numSim+0.0);
		probDwins1 = dWins1/(numSim+0.0);
		probRwins2 = rWins2/(numSim+0.0);
		probDwins2 = dWins2/(numSim+0.0);

		int M = 10000000;
		cplex = new IloCplex();
		cplex.setParam(IloCplex.Param.TimeLimit, 300);
		cplex.setOut(null);

		// Define decision variables
		IloNumVar[][] newR = new IloNumVar[3][numSim];
		IloNumVar[][] newD = new IloNumVar[3][numSim];
		IloNumVar[][] vD = new IloNumVar[4][numSim];		//Index 3 is the complete state
		IloNumVar[][] vR = new IloNumVar[4][numSim];
		IloNumVar[] alphaR = new IloNumVar[3];
		IloNumVar [] alphaD = new IloNumVar[3];

		for (int s = 0; s <= 2; s++) {
			alphaR[s] = cplex.numVar(0, Double.MAX_VALUE, IloNumVarType.Float); //Flip from D to R
			alphaD[s] = cplex.numVar(0, Double.MAX_VALUE, IloNumVarType.Float); //Flip from R to D
		}

		for (int k = 0; k < numSim; k++) {
			for (int s = 0; s <= 2; s++) {
				newR[s][k] = cplex.numVar(0, Double.MAX_VALUE, IloNumVarType.Float);
				newD[s][k] = cplex.numVar(0, Double.MAX_VALUE, IloNumVarType.Float);
				vD[s][k] = cplex.numVar(0, 1, IloNumVarType.Bool);
				vR[s][k] = cplex.numVar(0, 1, IloNumVarType.Bool);
			}
			vD[3][k] = cplex.numVar(0, 1, IloNumVarType.Bool);
			vR[3][k] = cplex.numVar(0, 1, IloNumVarType.Bool);
		}
		// New votes
		for (int k = 0; k < numSim; k++) {
			for (int s = 0; s <= 2; s++) {
				IloLinearNumExpr expr = cplex.linearNumExpr();
				expr.addTerm(1, newD[s][k]);
				expr.addTerm(-1, alphaD[s]);
				expr.addTerm(1, alphaR[s]);
				cplex.addEq(expr, TDE[s][k]);

				IloLinearNumExpr exprR = cplex.linearNumExpr();
				exprR.addTerm(1, newR[s][k]);
				exprR.addTerm(-1, alphaR[s]);
				exprR.addTerm(1, alphaD[s]);
				cplex.addEq(exprR, TRE[s][k]);
			}
		}

		// Binary activations D
		for (int k = 0; k < numSim; k++) {
			for (int s = 0; s <= 2; s++) {
				IloLinearNumExpr expr = cplex.linearNumExpr();
				expr.addTerm(1, newD[s][k]);
				expr.addTerm(-1, newR[s][k]);
				expr.addTerm(-1, vD[s][k]);
				expr.addTerm(-M, vD[s][k]);
				cplex.addGe(expr, -M);
			}
			// For the complete state
			IloLinearNumExpr exprSD = cplex.linearNumExpr();
			exprSD.addTerm(1, newD[0][k]);
			exprSD.addTerm(1, newD[1][k]);
			exprSD.addTerm(1, newD[2][k]);
			exprSD.addTerm(-1, newR[0][k]);
			exprSD.addTerm(-1, newR[1][k]);
			exprSD.addTerm(-1, newR[2][k]);
			exprSD.addTerm(-1, vD[3][k]);
			exprSD.addTerm(-M, vD[3][k]);
			cplex.addGe(exprSD, -M);
		}
		// Binary activations R
		for (int k = 0; k < numSim; k++) {
			for (int s = 0; s <= 2; s++) {
				IloLinearNumExpr expr = cplex.linearNumExpr();
				expr.addTerm(1, newR[s][k]);
				expr.addTerm(-1, newD[s][k]);
				expr.addTerm(-1, vR[s][k]);
				expr.addTerm(-M, vR[s][k]);
				cplex.addGe(expr, -M);
			}
			// For the complete state
			IloLinearNumExpr exprSR = cplex.linearNumExpr();
			exprSR.addTerm(1, newR[0][k]);
			exprSR.addTerm(1, newR[1][k]);
			exprSR.addTerm(1, newR[2][k]);
			exprSR.addTerm(-1, newD[0][k]);
			exprSR.addTerm(-1, newD[1][k]);
			exprSR.addTerm(-1, newD[2][k]);
			exprSR.addTerm(-1, vR[3][k]);
			exprSR.addTerm(-M, vR[3][k]);
			cplex.addGe(exprSR, -M);
		}
		//1 winner constraints
		for (int k = 0; k < numSim; k++) {
			for (int s = 0; s <= 3; s++) {
				IloLinearNumExpr expr = cplex.linearNumExpr();
				expr.addTerm(1, vD[s][k]);
				expr.addTerm(1, vR[s][k]);
				cplex.addEq(expr, 1);
			}
		}

		// Match prediction 0
		IloLinearNumExpr expr = cplex.linearNumExpr();
		for (int k = 0; k < numSim; k++) {
			expr.addTerm(1, vD[0][k]);
		}
		cplex.addGe(expr, data.NSPredDE[2]*numSim-100);
		cplex.addLe(expr, data.NSPredDE[2]*numSim+100);

		IloLinearNumExpr exprR = cplex.linearNumExpr();
		for (int k = 0; k < numSim; k++) {
			exprR.addTerm(1, vR[0][k]);
		}
		cplex.addGe(exprR, data.NSPredRE[2]*numSim-100);
		cplex.addLe(exprR, data.NSPredRE[2]*numSim+100);

		// Match prediction 1
		IloLinearNumExpr expr1 = cplex.linearNumExpr();
		for (int k = 0; k < numSim; k++) {
			expr1.addTerm(1, vD[1][k]);
		}
		cplex.addGe(expr1, data.NSPredDE[3]*numSim-100);
		cplex.addLe(expr1, data.NSPredDE[3]*numSim+100);

		IloLinearNumExpr exprR1 = cplex.linearNumExpr();
		for (int k = 0; k < numSim; k++) {
			exprR1.addTerm(1, vR[1][k]);
		}
		cplex.addGe(exprR1, data.NSPredRE[3]*numSim-100);
		cplex.addLe(exprR1, data.NSPredRE[3]*numSim+100);

		// Match prediction 2
		IloLinearNumExpr expr2 = cplex.linearNumExpr();
		for (int k = 0; k < numSim; k++) {
			expr2.addTerm(1, vD[2][k]);
		}
		cplex.addGe(expr2, data.NSPredDE[4]*numSim-100);
		cplex.addLe(expr2, data.NSPredDE[4]*numSim+100);

		IloLinearNumExpr exprR2 = cplex.linearNumExpr();
		for (int k = 0; k < numSim; k++) {
			exprR2.addTerm(1, vR[2][k]);
		}
		cplex.addGe(exprR2, data.NSPredRE[4]*numSim-100);
		cplex.addLe(exprR2, data.NSPredRE[4]*numSim+100);

		// Match prediction STATE
		IloLinearNumExpr exprS = cplex.linearNumExpr();
		for (int k = 0; k < numSim; k++) {
			exprS.addTerm(1, vD[3][k]);
		}
		cplex.addGe(exprS, data.PredD[27]*numSim-25);
		cplex.addLe(exprS, data.PredD[27]*numSim+25);

		IloLinearNumExpr exprSR1 = cplex.linearNumExpr();
		for (int k = 0; k < numSim; k++) {
			exprSR1.addTerm(1, vR[3][k]);
		}
		cplex.addGe(exprSR1, data.PredR[27]*numSim-25);
		cplex.addLe(exprSR1, data.PredR[27]*numSim+25);

		// Objective	
		IloLinearNumExpr obj = cplex.linearNumExpr();
		obj.addTerm(1, alphaR[0]);
		obj.addTerm(1, alphaR[1]);
		obj.addTerm(1, alphaR[2]);
		obj.addTerm(1, alphaD[0]);
		obj.addTerm(1, alphaD[1]);
		obj.addTerm(1, alphaD[2]);
		cplex.addMinimize(obj);
		cplex.solve();

		if(cplex.getStatus() == IloCplex.Status.Infeasible || cplex.getStatus() == IloCplex.Status.InfeasibleOrUnbounded){
			System.out.println("		Adjust NEBRASKA is infeasible!");
		}
		else {
			int newDwins0 = 0;
			int newRwins0 = 0;
			int newDwins1 = 0;
			int newRwins1 = 0;
			int newDwins2 = 0;
			int newRwins2 = 0;
			int newDS = 0;
			int newRS = 0;


			for (int k = 0; k < numSim; k++) {
				newDwins0+=(int)Math.round(cplex.getValue(vD[0][k]));
				newRwins0+=(int)Math.round(cplex.getValue(vR[0][k]));
				newDwins1+=(int)Math.round(cplex.getValue(vD[1][k]));
				newRwins1+=(int)Math.round(cplex.getValue(vR[1][k]));
				newDwins2+=(int)Math.round(cplex.getValue(vD[2][k]));
				newRwins2+=(int)Math.round(cplex.getValue(vR[2][k]));
				newDS+=(int)Math.round(cplex.getValue(vD[3][k]));
				newRS+=(int)Math.round(cplex.getValue(vR[3][k]));
			}

			System.out.println("NEBRASKA"+",Initial_RWINS_Prob:,"+probRwins+",Initial_DWINS_Prob:,"+probDwins+",NSPredR:,"+ data.PredR[27]+",NSPredD:,"+ data.PredD[27]+",alpha:,0,newProbRwins:,"+newRS/(numSim+0.0)+",newDprob:,"+newDS/(numSim+0.0));
			System.out.println("NEBRASKA_D1"+",Initial_RWINS_Prob:,"+probRwins0+",Initial_DWINS_Prob:,"+probDwins0+",NSPredR:,"+ data.NSPredRE[2]+",NSPredD:,"+ data.NSPredDE[2]+",alpha:,"+(cplex.getValue(alphaD[0])/(pop0+0.0)+cplex.getValue(alphaR[0])/(pop0+0.0)));
			System.out.println("NEBRASKA_D2"+",Initial_RWINS_Prob:,"+probRwins1+",Initial_DWINS_Prob:,"+probDwins1+",NSPredR:,"+ data.NSPredRE[3]+",NSPredD:,"+ data.NSPredDE[3]+",alpha:,"+(cplex.getValue(alphaD[1])/(pop1+0.0)+cplex.getValue(alphaR[1])/(pop1+0.0)));
			System.out.println("NEBRASKA_D3"+",Initial_RWINS_Prob:,"+probRwins1+",Initial_DWINS_Prob:,"+probDwins1+",NSPredR:,"+ data.NSPredRE[4]+",NSPredD:,"+ data.NSPredDE[4]+",alpha:,"+(cplex.getValue(alphaD[2])/(pop1+0.0)+cplex.getValue(alphaR[2])/(pop1+0.0)));

		}



		// Adjust all counties
		double alphaValD0 = cplex.getValue(alphaD[0])/(pop0+0.0);
		double alphaValD1 = cplex.getValue(alphaD[1])/(pop1+0.0);
		double alphaValD2 = cplex.getValue(alphaD[2])/(pop1+0.0);
		double alphaValR0 = cplex.getValue(alphaR[0])/(pop0+0.0);
		double alphaValR1 = cplex.getValue(alphaR[1])/(pop1+0.0);
		double alphaValR2 = cplex.getValue(alphaR[2])/(pop1+0.0);

		for (int k = 0; k < numSim; k++) {
			for (int i = 0; i < data.nebraska1.size(); i++) {
				data.countyList.get(data.nebraska1.get(i)).dVotesSim[k] = data.countyList.get(data.nebraska1.get(i)).dVotesSim[k] + alphaValD0*data.countyList.get(data.nebraska1.get(i)).rVotesSim[k];
				data.countyList.get(data.nebraska1.get(i)).rVotesSim[k] = (1-alphaValD0)*data.countyList.get(data.nebraska1.get(i)).rVotesSim[k];
				data.countyList.get(data.nebraska1.get(i)).rVotesSim[k] = data.countyList.get(data.nebraska1.get(i)).rVotesSim[k] + alphaValR0*data.countyList.get(data.nebraska1.get(i)).dVotesSim[k];
				data.countyList.get(data.nebraska1.get(i)).dVotesSim[k] = (1-alphaValR0)*data.countyList.get(data.nebraska1.get(i)).dVotesSim[k];
			}
			for (int i = 0; i < data.nebraska2.size(); i++) {
				data.countyList.get(data.nebraska2.get(i)).dVotesSim[k] = data.countyList.get(data.nebraska2.get(i)).dVotesSim[k] + alphaValD1*data.countyList.get(data.nebraska2.get(i)).rVotesSim[k];
				data.countyList.get(data.nebraska2.get(i)).rVotesSim[k] = (1-alphaValD1)*data.countyList.get(data.nebraska2.get(i)).rVotesSim[k];
				data.countyList.get(data.nebraska2.get(i)).rVotesSim[k] = data.countyList.get(data.nebraska2.get(i)).rVotesSim[k] + alphaValR1*data.countyList.get(data.nebraska2.get(i)).dVotesSim[k];
				data.countyList.get(data.nebraska2.get(i)).dVotesSim[k] = (1-alphaValR1)*data.countyList.get(data.nebraska2.get(i)).dVotesSim[k];
			}
			for (int i = 0; i < data.nebraska3.size(); i++) {
				data.countyList.get(data.nebraska3.get(i)).dVotesSim[k] = data.countyList.get(data.nebraska3.get(i)).dVotesSim[k] + alphaValD2*data.countyList.get(data.nebraska3.get(i)).rVotesSim[k];
				data.countyList.get(data.nebraska3.get(i)).rVotesSim[k] = (1-alphaValD2)*data.countyList.get(data.nebraska3.get(i)).rVotesSim[k];
				data.countyList.get(data.nebraska3.get(i)).rVotesSim[k] = data.countyList.get(data.nebraska3.get(i)).rVotesSim[k] + alphaValR2*data.countyList.get(data.nebraska3.get(i)).dVotesSim[k];
				data.countyList.get(data.nebraska3.get(i)).dVotesSim[k] = (1-alphaValR2)*data.countyList.get(data.nebraska3.get(i)).dVotesSim[k];
			}
		}

		cplex.clearModel();
		cplex.end();	


	}

	private void adjustNebraskaEvaporate(DataHandler data, int numSim) throws IloException {
		// Get the total number of votes for the especial districts
		double[][] TRE = new double[3][numSim];
		double[][] TDE = new double[3][numSim];
		double probRwins0 = 0;
		double probDwins0 = 0;
		int rWins0 = 0;
		int dWins0 = 0;
		double probRwins1 = 0;
		double probDwins1 = 0;
		int rWins1 = 0;
		int dWins1 = 0;
		double probRwins2 = 0;
		double probDwins2 = 0;
		int rWins2 = 0;
		int dWins2 = 0;
		int pop0 = 0;
		int pop1 = 0;
		int pop2 = 0;
		int rWins = 0;
		int dWins = 0;
		double probRwins = 0;
		double probDwins = 0;
	
		for (int i = 0; i < data.nebraska1.size(); i++) {
			pop0+=data.countyList.get(data.nebraska1.get(i)).pop;
		}
		for (int i = 0; i < data.nebraska2.size(); i++) {
			pop1 += data.countyList.get(data.nebraska2.get(i)).pop;
		}
		for (int i = 0; i < data.nebraska3.size(); i++) {
			pop2 += data.countyList.get(data.nebraska3.get(i)).pop;
		}
	
		for (int q = 0; q < numSim; q++) {
			for (int s = 0; s <= 2; s++) {
				if(s==0) {
					for (int i = 0; i < data.nebraska1.size(); i++) {
						TRE[s][q] += data.countyList.get(data.nebraska1.get(i)).rVotesSim[q];
					}
					//Only half of the votes for shared county 31153
					TRE[s][q] -= 0.5*data.countyList.get(data.countyIndex.get(31153)).rVotesSim[q];
				}
				else if(s==1) {
					for (int i = 0; i < data.nebraska2.size(); i++) {
						TRE[s][q] += data.countyList.get(data.nebraska2.get(i)).rVotesSim[q];
					}
					//Only half of the votes for shared county 31153
					TRE[s][q] -= 0.5*data.countyList.get(data.countyIndex.get(31153)).rVotesSim[q];
				}
				else if(s==2) {
					for (int i = 0; i < data.nebraska3.size(); i++) {
						TRE[s][q] += data.countyList.get(data.nebraska3.get(i)).rVotesSim[q];
					}
				}
	
				if(s==0) {
					for (int i = 0; i < data.nebraska1.size(); i++) {
						TDE[s][q] += data.countyList.get(data.nebraska1.get(i)).dVotesSim[q];
					}
					//Only half of the votes for shared county 31153
					TDE[s][q]-=0.5*data.countyList.get(data.countyIndex.get(31153)).dVotesSim[q];
				}
				else if(s==1) {
					for (int i = 0; i < data.nebraska2.size(); i++) {
						TDE[s][q] +=data.countyList.get(data.nebraska2.get(i)).dVotesSim[q];
					}
					//Only half of the votes for shared county 31153
					TDE[s][q]-=0.5*data.countyList.get(data.countyIndex.get(31153)).dVotesSim[q];
				}
				else if(s==2) {
					for (int i = 0; i < data.nebraska3.size(); i++) {
						TDE[s][q] +=data.countyList.get(data.nebraska3.get(i)).dVotesSim[q];
					}
				}
			}
	
			if(TRE[0][q] > TDE[0][q]) {rWins0++;}
			else {dWins0++;}
			if(TRE[1][q] > TDE[1][q]) {rWins1++;}
			else {dWins1++;}
			if(TRE[2][q] > TDE[2][q]) {rWins2++;}
			else {dWins2++;}
			if(TRE[0][q]+TRE[1][q]+TRE[2][q] > TDE[0][q]+TDE[1][q]+TDE[2][q]) {rWins++;}
			else {dWins++;}
		}
	
	
		probRwins = rWins/(numSim+0.0);
		probDwins = dWins/(numSim+0.0);
		probRwins0 = rWins0/(numSim+0.0);
		probDwins0 = dWins0/(numSim+0.0);
		probRwins1 = rWins1/(numSim+0.0);
		probDwins1 = dWins1/(numSim+0.0);
		probRwins2 = rWins2/(numSim+0.0);
		probDwins2 = dWins2/(numSim+0.0);
	
		int M = 10000000;
		cplex = new IloCplex();
		cplex.setParam(IloCplex.Param.TimeLimit, 300);
		cplex.setOut(null);
	
		// Define decision variables
		IloNumVar[][] newR = new IloNumVar[3][numSim];
		IloNumVar[][] newD = new IloNumVar[3][numSim];
		IloNumVar[][] vD = new IloNumVar[4][numSim];		//Index 3 is the complete state
		IloNumVar[][] vR = new IloNumVar[4][numSim];
		IloNumVar[] alphaR = new IloNumVar[3];
		IloNumVar [] alphaD = new IloNumVar[3];
	
		for (int s = 0; s <= 2; s++) {
			alphaR[s] = cplex.numVar(0, Double.MAX_VALUE, IloNumVarType.Float); //Flip from D to R
			alphaD[s] = cplex.numVar(0, Double.MAX_VALUE, IloNumVarType.Float); //Flip from R to D
		}
	
		for (int k = 0; k < numSim; k++) {
			for (int s = 0; s <= 2; s++) {
				newR[s][k] = cplex.numVar(0, Double.MAX_VALUE, IloNumVarType.Float);
				newD[s][k] = cplex.numVar(0, Double.MAX_VALUE, IloNumVarType.Float);
				vD[s][k] = cplex.numVar(0, 1, IloNumVarType.Bool);
				vR[s][k] = cplex.numVar(0, 1, IloNumVarType.Bool);
			}
			vD[3][k] = cplex.numVar(0, 1, IloNumVarType.Bool);
			vR[3][k] = cplex.numVar(0, 1, IloNumVarType.Bool);
		}
		// New votes
		for (int k = 0; k < numSim; k++) {
			for (int s = 0; s <= 2; s++) {
				IloLinearNumExpr expr = cplex.linearNumExpr();
				expr.addTerm(1, newD[s][k]);
				expr.addTerm(1, alphaR[s]);
				cplex.addEq(expr, TDE[s][k]);
	
				IloLinearNumExpr exprR = cplex.linearNumExpr();
				exprR.addTerm(1, newR[s][k]);
				exprR.addTerm(1, alphaD[s]);
				cplex.addEq(exprR, TRE[s][k]);
			}
		}
	
		// Binary activations D
		for (int k = 0; k < numSim; k++) {
			for (int s = 0; s <= 2; s++) {
				IloLinearNumExpr expr = cplex.linearNumExpr();
				expr.addTerm(1, newD[s][k]);
				expr.addTerm(-1, newR[s][k]);
				expr.addTerm(-1, vD[s][k]);
				expr.addTerm(-M, vD[s][k]);
				cplex.addGe(expr, -M);
			}
			// For the complete state
			IloLinearNumExpr exprSD = cplex.linearNumExpr();
			exprSD.addTerm(1, newD[0][k]);
			exprSD.addTerm(1, newD[1][k]);
			exprSD.addTerm(1, newD[2][k]);
			exprSD.addTerm(-1, newR[0][k]);
			exprSD.addTerm(-1, newR[1][k]);
			exprSD.addTerm(-1, newR[2][k]);
			exprSD.addTerm(-1, vD[3][k]);
			exprSD.addTerm(-M, vD[3][k]);
			cplex.addGe(exprSD, -M);
		}
		// Binary activations R
		for (int k = 0; k < numSim; k++) {
			for (int s = 0; s <= 2; s++) {
				IloLinearNumExpr expr = cplex.linearNumExpr();
				expr.addTerm(1, newR[s][k]);
				expr.addTerm(-1, newD[s][k]);
				expr.addTerm(-1, vR[s][k]);
				expr.addTerm(-M, vR[s][k]);
				cplex.addGe(expr, -M);
			}
			// For the complete state
			IloLinearNumExpr exprSR = cplex.linearNumExpr();
			exprSR.addTerm(1, newR[0][k]);
			exprSR.addTerm(1, newR[1][k]);
			exprSR.addTerm(1, newR[2][k]);
			exprSR.addTerm(-1, newD[0][k]);
			exprSR.addTerm(-1, newD[1][k]);
			exprSR.addTerm(-1, newD[2][k]);
			exprSR.addTerm(-1, vR[3][k]);
			exprSR.addTerm(-M, vR[3][k]);
			cplex.addGe(exprSR, -M);
		}
		//1 winner constraints
		for (int k = 0; k < numSim; k++) {
			for (int s = 0; s <= 3; s++) {
				IloLinearNumExpr expr = cplex.linearNumExpr();
				expr.addTerm(1, vD[s][k]);
				expr.addTerm(1, vR[s][k]);
				cplex.addEq(expr, 1);
			}
		}
	
		// Match prediction 0
		IloLinearNumExpr expr = cplex.linearNumExpr();
		for (int k = 0; k < numSim; k++) {
			expr.addTerm(1, vD[0][k]);
		}
		cplex.addGe(expr, data.NSPredDE[2]*numSim-100);
		cplex.addLe(expr, data.NSPredDE[2]*numSim+100);
	
		IloLinearNumExpr exprR = cplex.linearNumExpr();
		for (int k = 0; k < numSim; k++) {
			exprR.addTerm(1, vR[0][k]);
		}
		cplex.addGe(exprR, data.NSPredRE[2]*numSim-100);
		cplex.addLe(exprR, data.NSPredRE[2]*numSim+100);
	
		// Match prediction 1
		IloLinearNumExpr expr1 = cplex.linearNumExpr();
		for (int k = 0; k < numSim; k++) {
			expr1.addTerm(1, vD[1][k]);
		}
		cplex.addGe(expr1, data.NSPredDE[3]*numSim-100);
		cplex.addLe(expr1, data.NSPredDE[3]*numSim+100);
	
		IloLinearNumExpr exprR1 = cplex.linearNumExpr();
		for (int k = 0; k < numSim; k++) {
			exprR1.addTerm(1, vR[1][k]);
		}
		cplex.addGe(exprR1, data.NSPredRE[3]*numSim-100);
		cplex.addLe(exprR1, data.NSPredRE[3]*numSim+100);
	
		// Match prediction 2
		IloLinearNumExpr expr2 = cplex.linearNumExpr();
		for (int k = 0; k < numSim; k++) {
			expr2.addTerm(1, vD[2][k]);
		}
		cplex.addGe(expr2, data.NSPredDE[4]*numSim-100);
		cplex.addLe(expr2, data.NSPredDE[4]*numSim+100);
	
		IloLinearNumExpr exprR2 = cplex.linearNumExpr();
		for (int k = 0; k < numSim; k++) {
			exprR2.addTerm(1, vR[2][k]);
		}
		cplex.addGe(exprR2, data.NSPredRE[4]*numSim-100);
		cplex.addLe(exprR2, data.NSPredRE[4]*numSim+100);
	
		// Match prediction STATE
		IloLinearNumExpr exprS = cplex.linearNumExpr();
		for (int k = 0; k < numSim; k++) {
			exprS.addTerm(1, vD[3][k]);
		}
		cplex.addGe(exprS, data.PredD[27]*numSim-25);
		cplex.addLe(exprS, data.PredD[27]*numSim+25);
	
		IloLinearNumExpr exprSR1 = cplex.linearNumExpr();
		for (int k = 0; k < numSim; k++) {
			exprSR1.addTerm(1, vR[3][k]);
		}
		cplex.addGe(exprSR1, data.PredR[27]*numSim-25);
		cplex.addLe(exprSR1, data.PredR[27]*numSim+25);
	
		// Objective	
		IloLinearNumExpr obj = cplex.linearNumExpr();
		obj.addTerm(1, alphaR[0]);
		obj.addTerm(1, alphaR[1]);
		obj.addTerm(1, alphaR[2]);
		obj.addTerm(1, alphaD[0]);
		obj.addTerm(1, alphaD[1]);
		obj.addTerm(1, alphaD[2]);
		cplex.addMinimize(obj);
		cplex.solve();
	
		if(cplex.getStatus() == IloCplex.Status.Infeasible || cplex.getStatus() == IloCplex.Status.InfeasibleOrUnbounded){
			System.out.println("		Adjust NEBRASKA is infeasible!");
		}
		else {
			int newDwins0 = 0;
			int newRwins0 = 0;
			int newDwins1 = 0;
			int newRwins1 = 0;
			int newDwins2 = 0;
			int newRwins2 = 0;
			int newDS = 0;
			int newRS = 0;
	
	
			for (int k = 0; k < numSim; k++) {
				newDwins0+=(int)Math.round(cplex.getValue(vD[0][k]));
				newRwins0+=(int)Math.round(cplex.getValue(vR[0][k]));
				newDwins1+=(int)Math.round(cplex.getValue(vD[1][k]));
				newRwins1+=(int)Math.round(cplex.getValue(vR[1][k]));
				newDwins2+=(int)Math.round(cplex.getValue(vD[2][k]));
				newRwins2+=(int)Math.round(cplex.getValue(vR[2][k]));
				newDS+=(int)Math.round(cplex.getValue(vD[3][k]));
				newRS+=(int)Math.round(cplex.getValue(vR[3][k]));
			}
	
			System.out.println("NEBRASKA"+",Initial_RWINS_Prob:,"+probRwins+",Initial_DWINS_Prob:,"+probDwins+",NSPredR:,"+ data.PredR[27]+",NSPredD:,"+ data.PredD[27]+",alpha:,0,newProbRwins:,"+newRS/(numSim+0.0)+",newDprob:,"+newDS/(numSim+0.0));
			System.out.println("NEBRASKA_D1"+",Initial_RWINS_Prob:,"+probRwins0+",Initial_DWINS_Prob:,"+probDwins0+",NSPredR:,"+ data.NSPredRE[2]+",NSPredD:,"+ data.NSPredDE[2]+",alpha:,"+(cplex.getValue(alphaD[0])/(pop0+0.0)+cplex.getValue(alphaR[0])/(pop0+0.0)));
			System.out.println("NEBRASKA_D2"+",Initial_RWINS_Prob:,"+probRwins1+",Initial_DWINS_Prob:,"+probDwins1+",NSPredR:,"+ data.NSPredRE[3]+",NSPredD:,"+ data.NSPredDE[3]+",alpha:,"+(cplex.getValue(alphaD[1])/(pop1+0.0)+cplex.getValue(alphaR[1])/(pop1+0.0)));
			System.out.println("NEBRASKA_D3"+",Initial_RWINS_Prob:,"+probRwins1+",Initial_DWINS_Prob:,"+probDwins1+",NSPredR:,"+ data.NSPredRE[4]+",NSPredD:,"+ data.NSPredDE[4]+",alpha:,"+(cplex.getValue(alphaD[2])/(pop1+0.0)+cplex.getValue(alphaR[2])/(pop1+0.0)));
	
		}
	
	
	
		// Adjust all counties
		double alphaValD0 = cplex.getValue(alphaD[0])/(pop0+0.0);
		double alphaValD1 = cplex.getValue(alphaD[1])/(pop1+0.0);
		double alphaValD2 = cplex.getValue(alphaD[2])/(pop1+0.0);
		double alphaValR0 = cplex.getValue(alphaR[0])/(pop0+0.0);
		double alphaValR1 = cplex.getValue(alphaR[1])/(pop1+0.0);
		double alphaValR2 = cplex.getValue(alphaR[2])/(pop1+0.0);
	
		for (int k = 0; k < numSim; k++) {
			for (int i = 0; i < data.nebraska1.size(); i++) {
				data.countyList.get(data.nebraska1.get(i)).rVotesSim[k] = (1-alphaValD0)*data.countyList.get(data.nebraska1.get(i)).rVotesSim[k];
				data.countyList.get(data.nebraska1.get(i)).dVotesSim[k] = (1-alphaValR0)*data.countyList.get(data.nebraska1.get(i)).dVotesSim[k];
			}
			for (int i = 0; i < data.nebraska2.size(); i++) {
				data.countyList.get(data.nebraska2.get(i)).rVotesSim[k] = (1-alphaValD1)*data.countyList.get(data.nebraska2.get(i)).rVotesSim[k];
				data.countyList.get(data.nebraska2.get(i)).dVotesSim[k] = (1-alphaValR1)*data.countyList.get(data.nebraska2.get(i)).dVotesSim[k];
			}
			for (int i = 0; i < data.nebraska3.size(); i++) {
				data.countyList.get(data.nebraska3.get(i)).rVotesSim[k] = (1-alphaValD2)*data.countyList.get(data.nebraska3.get(i)).rVotesSim[k];
				data.countyList.get(data.nebraska3.get(i)).dVotesSim[k] = (1-alphaValR2)*data.countyList.get(data.nebraska3.get(i)).dVotesSim[k];
			}
		}
	
		cplex.clearModel();
		cplex.end();	
	
	
	}

	void modelTargetSAAProb(DataHandler data, double lambda, double limit, int pLimit, String party, String poll) throws IloException, FileNotFoundException {
		int M = 100000000;
		int numSim = data.countyList.get(0).rVotesSim.length;
		cplex = new IloCplex();
		cplex.setParam(IloCplex.Param.TimeLimit, 3600);
		//cplex.setOut(null);
		//Compute admissible pairs
		int countVars = 0;
		ArrayList<int[]> pairs = new ArrayList<int[]>();
		for (int i = 0; i < data.numCounty; i++) {
			for (int j = 0; j < data.numCounty; j++) {
				if(!data.countyList.get(i).state.equals("MAINE") && !data.countyList.get(j).state.equals("MAINE") && !data.countyList.get(i).state.equals("NEBRASKA") && !data.countyList.get(j).state.equals("NEBRASKA") ) { //Don't move from MAINE or NEBRASKA
					if(!data.countyList.get(i).state.equals(data.countyList.get(j).state) && data.distMatrix[i][j] <= limit && data.countyList.get(i).pop >= 10000 && data.countyList.get(j).pop >= 10000) { //Don't move inside the same state!
						if( (party.equals("R") && Math.abs(0.5-data.PredR[data.stateIndex.get(data.countyList.get(j).state)]) <= 0.3) || (party.equals("D") && Math.abs(0.5-data.PredD[data.stateIndex.get(data.countyList.get(j).state)]) <= 0.3) ) {// Go to swing states		
							if((party.equals("R") && Math.abs(0.5-data.PredR[data.stateIndex.get(data.countyList.get(i).state)]) > 0.3) || (party.equals("D") && Math.abs(0.5-data.PredD[data.stateIndex.get(data.countyList.get(i).state)]) > 0.3)) {// From polarized states
								//System.out.println("WEPA "+data.countyList.get(i).state+" to "+data.countyList.get(j).state);	
								int[] aux = new int[2];
								aux[0] = i;
								aux[1] = j;
								pairs.add(aux);
								countVars++;
							}
						}
					}
				}
			}
		}
		System.out.println("TOTAL NUMBER OF FLOW VARS: "+countVars);
		// Define decision variables
		IloNumVar[] xR = new IloNumVar[countVars];
		IloNumVar[] xD = new IloNumVar[countVars];
		IloNumVar[][] zR = new IloNumVar[data.numCounty][numSim]; // Votes per county
		IloNumVar[][] zD = new IloNumVar[data.numCounty][numSim];
		IloNumVar[][] TR = new IloNumVar[data.numStates][numSim]; // Votes per state
		IloNumVar[][] TD = new IloNumVar[data.numStates][numSim];
		IloNumVar elecR[] = new IloNumVar[numSim];
		IloNumVar elecD[] = new IloNumVar[numSim];
		IloNumVar[][] vD = new IloNumVar[data.numStates][numSim]; //Binary aux
		IloNumVar[][] vR = new IloNumVar[data.numStates][numSim];
		// Special districts 0=Maine1, 1=Maine2, ... 4=Nebraska3
		IloNumVar[][] TRE = new IloNumVar[5][numSim]; // Votes per special district
		IloNumVar[][] TDE = new IloNumVar[5][numSim];
		IloNumVar[][] vDE = new IloNumVar[5][numSim]; //Binary aux
		IloNumVar[][] vRE = new IloNumVar[5][numSim];
		// Prefered party wins per simulation
		IloNumVar partyWins[] = new IloNumVar[numSim];
		for (int k = 0; k < countVars; k++) {
			int i = pairs.get(k)[0];
			if(party.equals("R")) {
				xR[k] = cplex.numVar(0, data.countyList.get(i).minRvotes, IloNumVarType.Float); 
				xD[k] = cplex.numVar(0, 0, IloNumVarType.Float);
			}
			else {
				xR[k] = cplex.numVar(0, 0, IloNumVarType.Float); 
				xD[k] = cplex.numVar(0, data.countyList.get(i).minDvotes, IloNumVarType.Float);
			}
		}
		for (int q = 0; q < numSim; q++) {
			for (int i = 0; i < data.numCounty; i++) {
				zR[i][q] = cplex.numVar(0, Double.MAX_VALUE, IloNumVarType.Float);
				zD[i][q] = cplex.numVar(0, Double.MAX_VALUE, IloNumVarType.Float);
			}
			for (int s = 0; s < data.numStates; s++) {
				TR[s][q] = cplex.numVar(0, Double.MAX_VALUE, IloNumVarType.Float);
				TD[s][q] = cplex.numVar(0, Double.MAX_VALUE, IloNumVarType.Float);
				vD[s][q] = cplex.numVar(0, 1, IloNumVarType.Bool);
				vR[s][q] = cplex.numVar(0, 1, IloNumVarType.Bool);
			}
			for (int s = 0; s < 5; s++) {
				TRE[s][q] = cplex.numVar(0, Double.MAX_VALUE, IloNumVarType.Float);
				TDE[s][q] = cplex.numVar(0, Double.MAX_VALUE, IloNumVarType.Float);
				vDE[s][q] = cplex.numVar(0, 1, IloNumVarType.Bool);
				vRE[s][q] = cplex.numVar(0, 1, IloNumVarType.Bool);
			}
			elecR[q] = cplex.numVar(0, 9999, IloNumVarType.Int); 
			elecD[q] = cplex.numVar(0, 9999, IloNumVarType.Int); 
			partyWins[q] = cplex.numVar(0, 1, IloNumVarType.Bool); 
		}
		// Total people moved out of a county
		for (int i = 0; i < data.numCounty; i++) {
			IloLinearNumExpr expr = cplex.linearNumExpr();
			for (int k = 0; k < countVars; k++) {
				int origin = pairs.get(k)[0];
				if(origin == i) {
					expr.addTerm(1, xR[k]);
					expr.addTerm(1, xD[k]);
				}
			}
			cplex.addLe(expr, lambda*data.countyList.get(i).pop);
		}
		// Total people moved into a county
		for (int i = 0; i < data.numCounty; i++) {
			IloLinearNumExpr expr = cplex.linearNumExpr();
			for (int k = 0; k < countVars; k++) {
				int dest = pairs.get(k)[1];
				if(dest == i) {
					expr.addTerm(1, xR[k]);
					expr.addTerm(1, xD[k]);
				}
			}
			cplex.addLe(expr, lambda*data.countyList.get(i).pop);
		}

		// Total votes for a county
		for (int q = 0; q < numSim; q++) {
			IloLinearNumExpr[] exprR = new IloLinearNumExpr[data.numCounty];
			IloLinearNumExpr[] exprD = new IloLinearNumExpr[data.numCounty];
			// Initialize constraints
			for (int i = 0; i < data.numCounty; i++) {
				exprR[i] = cplex.linearNumExpr();
				exprD[i] = cplex.linearNumExpr();
			}
			// Consider movement
			for (int k = 0; k < countVars; k++) {
				int origin = pairs.get(k)[0];
				int destination = pairs.get(k)[1];
				exprR[origin].addTerm(1, xR[k]);
				exprR[destination].addTerm(-1, xR[k]);
				exprD[origin].addTerm(1, xD[k]);
				exprD[destination].addTerm(-1, xD[k]);
			}
			// Add coefficients and rhs for all counties
			for (int i = 0; i < data.numCounty; i++) {
				County c = data.countyList.get(i);
				exprR[i].addTerm(1, zR[i][q]);
				exprD[i].addTerm(1, zD[i][q]);
				cplex.addEq(exprR[i], c.rVotesSim[q]);
				cplex.addEq(exprD[i], c.dVotesSim[q]);
			}
		}
		// Total votes for a state
		for (int q = 0; q < numSim; q++) {
			for (int s = 0; s < data.numStates; s++) {
				IloLinearNumExpr exprR = cplex.linearNumExpr();
				exprR.addTerm(1, TR[s][q]);
				for (int i = 0; i < data.numCounty; i++) {
					if(data.countyList.get(i).state.equals(data.states[s])) {
						exprR.addTerm(-1, zR[i][q]);
						//System.out.println("MATCH COUNTY "+i+" "+data.states[s]);
					}
				}
				cplex.addEq(exprR, 0);
				IloLinearNumExpr exprD = cplex.linearNumExpr();
				exprD.addTerm(1, TD[s][q]);
				for (int i = 0; i < data.numCounty; i++) {
					if(data.countyList.get(i).state.equals(data.states[s])) {
						exprD.addTerm(-1, zD[i][q]);
						//System.out.println("MATCH COUNTY "+i+" "+data.states[s]);
					}
				}
				cplex.addEq(exprD, 0);
			}
			// Total votes for special districts
			for (int s = 0; s < 5; s++) {
				IloLinearNumExpr exprR = cplex.linearNumExpr();
				exprR.addTerm(1, TRE[s][q]);
				if(s==0) {
					for (int i = 0; i < data.maine1.size(); i++) {
						exprR.addTerm(-1, zR[data.maine1.get(i)][q]);
						//System.out.println("MAINE 1 COUNTY "+data.countyList.get(data.maine1.get(i)).name);
					}
					//Only half of the votes for shared county 23011
					exprR.addTerm(0.5, zR[data.countyIndex.get(23011)][q]);
				}
				else if(s==1) {
					for (int i = 0; i < data.maine2.size(); i++) {
						exprR.addTerm(-1, zR[data.maine2.get(i)][q]);
						//System.out.println("MAINE 2 COUNTY "+data.countyList.get(data.maine2.get(i)).name);
					}
					//Only half of the votes for shared county 23011
					exprR.addTerm(0.5, zR[data.countyIndex.get(23011)][q]);
				}
				else if(s==2) {
					for (int i = 0; i < data.nebraska1.size(); i++) {
						exprR.addTerm(-1, zR[data.nebraska1.get(i)][q]);
					}
					//Only half of the votes for shared county 31153
					exprR.addTerm(0.5, zR[data.countyIndex.get(31153)][q]);
				}
				else if(s==3) {
					for (int i = 0; i < data.nebraska2.size(); i++) {
						exprR.addTerm(-1, zR[data.nebraska2.get(i)][q]);
					}
					//Only half of the votes for shared county 31153
					exprR.addTerm(0.5, zR[data.countyIndex.get(31153)][q]);
				}
				else if(s==4) {
					for (int i = 0; i < data.nebraska3.size(); i++) {
						exprR.addTerm(-1, zR[data.nebraska3.get(i)][q]);
					}
				}
				cplex.addEq(exprR, 0);
				IloLinearNumExpr exprD = cplex.linearNumExpr();
				exprD.addTerm(1, TDE[s][q]);
				if(s==0) {
					for (int i = 0; i < data.maine1.size(); i++) {
						exprD.addTerm(-1, zD[data.maine1.get(i)][q]);
						//System.out.println("MAINE 1 COUNTY "+data.countyList.get(data.maine1.get(i)).name);
					}
					//Only half of the votes for shared county 23011
					exprD.addTerm(0.5, zD[data.countyIndex.get(23011)][q]);
				}
				else if(s==1) {
					for (int i = 0; i < data.maine2.size(); i++) {
						exprD.addTerm(-1, zD[data.maine2.get(i)][q]);
						//System.out.println("MAINE 2 COUNTY "+data.countyList.get(data.maine2.get(i)).name);
					}
					//Only half of the votes for shared county 23011
					exprD.addTerm(0.5, zD[data.countyIndex.get(23011)][q]);
				}
				else if(s==2) {
					for (int i = 0; i < data.nebraska1.size(); i++) {
						exprD.addTerm(-1, zD[data.nebraska1.get(i)][q]);
					}
					//Only half of the votes for shared county 31153
					exprD.addTerm(0.5, zD[data.countyIndex.get(31153)][q]);
				}
				else if(s==3) {
					for (int i = 0; i < data.nebraska2.size(); i++) {
						exprD.addTerm(-1, zD[data.nebraska2.get(i)][q]);
					}
					//Only half of the votes for shared county 31153
					exprD.addTerm(0.5, zD[data.countyIndex.get(31153)][q]);
				}
				else if(s==4) {
					for (int i = 0; i < data.nebraska3.size(); i++) {
						exprD.addTerm(-1, zD[data.nebraska3.get(i)][q]);
					}
				}
				cplex.addEq(exprD, 0);
			}
		}
		// Elec votes count
		for (int q = 0; q < numSim; q++) {
			// Binary activations R
			for (int s = 0; s < data.numStates; s++) {
				IloLinearNumExpr expr = cplex.linearNumExpr();
				expr.addTerm(1, TR[s][q]);
				expr.addTerm(-1, TD[s][q]);
				expr.addTerm(-1, vR[s][q]);
				expr.addTerm(-M, vR[s][q]);
				cplex.addGe(expr, -M);
			}
			for (int s = 0; s < 5; s++) {
				IloLinearNumExpr expr = cplex.linearNumExpr();
				expr.addTerm(1, TRE[s][q]);
				expr.addTerm(-1, TDE[s][q]);
				expr.addTerm(-1, vRE[s][q]);
				expr.addTerm(-M, vRE[s][q]);
				cplex.addGe(expr, -M);
			}
			// Binary activations D
			for (int s = 0; s < data.numStates; s++) {
				IloLinearNumExpr expr = cplex.linearNumExpr();
				expr.addTerm(1, TD[s][q]);
				expr.addTerm(-1, TR[s][q]);
				expr.addTerm(-1, vD[s][q]);
				expr.addTerm(-M, vD[s][q]);
				cplex.addGe(expr, -M);
			}
			for (int s = 0; s < 5; s++) {
				IloLinearNumExpr expr = cplex.linearNumExpr();
				expr.addTerm(1, TDE[s][q]);
				expr.addTerm(-1, TRE[s][q]);
				expr.addTerm(-1, vDE[s][q]);
				expr.addTerm(-M, vDE[s][q]);
				cplex.addGe(expr, -M);
			}
			//Total Electoral votes for R
			IloLinearNumExpr exprA = cplex.linearNumExpr();
			exprA.addTerm(1, elecR[q]);
			for (int s = 0; s < data.numStates; s++) {
				exprA.addTerm(-data.eVotes[s], vR[s][q]);
			}
			for (int s = 0; s < 5; s++) {
				exprA.addTerm(-1, vRE[s][q]);
			}
			cplex.addEq(exprA, 0);
			//Total Electoral votes for D
			IloLinearNumExpr exprB = cplex.linearNumExpr();
			exprB.addTerm(1, elecD[q]);
			for (int s = 0; s < data.numStates; s++) {
				exprB.addTerm(-data.eVotes[s], vD[s][q]);
			}
			for (int s = 0; s < 5; s++) {
				exprB.addTerm(-1, vDE[s][q]);
			}
			cplex.addEq(exprB, 0);
		}
		// Prefered party wins count
		for (int q = 0; q < numSim; q++) {
			IloLinearNumExpr expr = cplex.linearNumExpr();
			if(party.equals("R")) {expr.addTerm(1, elecR[q]);}		//Use for R
			else{ expr.addTerm(1, elecD[q]);}						//Use for D
			expr.addTerm(-270, partyWins[q]);
			cplex.addGe(expr, 0);
		}

		// total number of people moved
		IloLinearNumExpr exprP = cplex.linearNumExpr();
		for (int k = 0; k < countVars; k++) {
			exprP.addTerm(1, xR[k]);
			exprP.addTerm(1, xD[k]);
		}
		cplex.addLe(exprP, pLimit);

		// Prob proxy
		IloLinearNumExpr objProb = cplex.linearNumExpr();
		for (int q = 0; q < numSim; q++) {
			objProb.addTerm(1, partyWins[q]);
		}

		cplex.addMaximize(objProb);
		cplex.solve();

		if(cplex.getStatus() == IloCplex.Status.Infeasible || cplex.getStatus() == IloCplex.Status.InfeasibleOrUnbounded){
			System.out.println(" MODEL TARGET SAAPROB is infeasible!");
		}
		else if(cplex.getStatus() == IloCplex.Status.Feasible || cplex.getStatus() == IloCplex.Status.Optimal) {
			double probProxy = cplex.getObjValue();
			xRsol = new double[data.numCounty][data.numCounty];
			xDsol = new double[data.numCounty][data.numCounty];
			xsol = new double[data.numCounty][data.numCounty];
			ArrayList<int[]> activeFlows = new ArrayList<int[]>();
			double tDist = 0;

			for (int k = 0; k < countVars; k++) {
				int origin = pairs.get(k)[0];
				int destination = pairs.get(k)[1];
				xRsol[origin][destination] = cplex.getValue(xR[k]);
				xDsol[origin][destination] = cplex.getValue(xD[k]);
				tDist+=data.distMatrix[origin][destination]*(xRsol[origin][destination]+xDsol[origin][destination]);
				if(xRsol[origin][destination]+xDsol[origin][destination]>0.0001) {
					activeFlows.add(pairs.get(k));
				}
			}
			if(party.equals("R")) {
				probRin = (probProxy/numSim+0.0);
			}
			else {
				probDin = (probProxy/numSim+0.0);
			}
			System.out.println("WEPAAAAAAAAAAAAAAAAAAAAAAAAAAAAA: "+tDist);
			cplex.clearModel();
			cplex.end();

			//postOpttune(data, lambda, activeFlows, party, pLimit); // Turn off when lambda < 1
			postMinDist(data, activeFlows, lambda, party);

			////////////////////////////////////////////////// PRINT RESULTS /////////////////////////////////////////////////////////////////
			java.io.PrintStream ps = new java.io.PrintStream( new java.io.FileOutputStream("Sol_"+party+"_"+pLimit+"_"+lambda+"_"+poll+".txt"));
			System.out.println("******************************  RESULTS ******************************************");
			double tFlow = 0;
			for (int i = 0; i < data.numCounty; i++) {
				for (int j = 0; j < data.numCounty; j++) {
					tFlow += xRsol[i][j];
					tFlow += xDsol[i][j];
				}
			}
			System.out.println("TOTAL_FLOW:,"+ tFlow);


			flow = tFlow;
			ps.println("TOTAL_FLOW:,"+ tFlow);
			ps.println("Origin_name,Origin_ID,Origin_State,Origin_Pop,Destination_name,Destination_ID,Destination_State,Destination_Pop,Distance,FlowR,FlowD");
			for (int i = 0; i < data.numCounty; i++) {
				for (int j = 0; j < data.numCounty; j++) {
					if(xRsol[i][j]+xDsol[i][j]>0.0001) {
						ps.println(data.countyList.get(i).name+","+data.countyList.get(i).ID+","+data.countyList.get(i).state+","+data.countyList.get(i).pop+","+data.countyList.get(j).name+","+data.countyList.get(j).ID+","+data.countyList.get(j).state+","+data.countyList.get(j).pop+","+data.distMatrix[i][j]+","+xRsol[i][j]+","+xDsol[i][j]);
					}
				}
			}

		}
	}

	private void postMinDist(DataHandler data, ArrayList<int[]> activeFlows, double lambda, String party) throws IloException {
		cplexB = new IloCplex();
		// Define decision variables
		IloNumVar[][] xR = new IloNumVar[data.numCounty][data.numCounty];
		IloNumVar[][] xD = new IloNumVar[data.numCounty][data.numCounty];
		for (int i = 0; i < data.numCounty; i++) {
			IloLinearNumExpr expr = cplexB.linearNumExpr();
			for (int j = 0; j < data.numCounty; j++) {
				if(data.distMatrix[i][j] <= 100 && i!=j && !data.countyList.get(i).state.equals(data.countyList.get(j).state)) {
					if(party.equals("R")) {
						xD[i][j] = cplexB.numVar(0, 0, IloNumVarType.Float);
						xR[i][j] = cplexB.numVar(0, data.countyList.get(i).minRvotes, IloNumVarType.Float);
					}
					else if(party.equals("D")) {
						xR[i][j] = cplexB.numVar(0, 0, IloNumVarType.Float); 
						xD[i][j]= cplexB.numVar(0, data.countyList.get(i).minDvotes, IloNumVarType.Float);
					}
					expr.addTerm(1, xR[i][j]);
					expr.addTerm(1, xD[i][j]);
				}
			}
			//Total out of county i <= limit
			cplexB.addLe(expr, data.countyList.get(i).minRvotes);
			cplexB.addLe(expr, lambda*data.countyList.get(i).pop);
		}
		
		// Total into any county <= limit
		for (int j = 0; j < data.numCounty; j++) {
			IloLinearNumExpr exprIN = cplexB.linearNumExpr();
			for (int i = 0; i < data.numCounty; i++) {
				if(data.distMatrix[i][j] <= 100 && i!=j && !data.countyList.get(i).state.equals(data.countyList.get(j).state)) {
						exprIN.addTerm(1, xR[i][j]);
						exprIN.addTerm(1, xD[i][j]);					
				}
			}
			cplexB.addLe(exprIN, lambda*data.countyList.get(j).pop);
		}


		//Maintain same IN flow for all states
		double flowStateIn[] = new double[data.numStates];
		double flowStateOut[] = new double[data.numStates];
		for (int k = 0; k < activeFlows.size(); k++) {
			int origin = activeFlows.get(k)[0];
			int destination = activeFlows.get(k)[1];
			int oState = data.stateIndex.get(data.countyList.get(origin).state);
			int dState = data.stateIndex.get(data.countyList.get(destination).state);
			flowStateOut[oState] += (xRsol[origin][destination]+xDsol[origin][destination]);
			flowStateIn[dState] += (xRsol[origin][destination]+xDsol[origin][destination]);
			//System.out.println("FLOW INTO "+data.countyList.get(destination).state+" : "+flowStateIn[dState]);
			//System.out.println("FLOW OUT "+data.countyList.get(origin).state+" : "+flowStateIn[dState]);
		}

		IloLinearNumExpr[] exprIN = new IloLinearNumExpr[data.numStates];
		IloLinearNumExpr[] exprOUT = new IloLinearNumExpr[data.numStates];
		for (int j = 0; j < data.numStates; j++) {
			exprIN[j] = cplexB.linearNumExpr();
			exprOUT[j] = cplexB.linearNumExpr();
		}
		for (int i = 0; i < data.numCounty; i++) {
			for (int j = 0; j < data.numCounty; j++) {
				if(data.distMatrix[i][j] <= 100 && i!=j && !data.countyList.get(i).state.equals(data.countyList.get(j).state)) {
					int oState = data.stateIndex.get(data.countyList.get(i).state);
					int dState = data.stateIndex.get(data.countyList.get(j).state);
					exprIN[dState].addTerm(1, xR[i][j]);
					exprIN[dState].addTerm(1, xD[i][j]);

					exprOUT[oState].addTerm(1, xR[i][j]);
					exprOUT[oState].addTerm(1, xD[i][j]);
				}
			}
		}
		for (int j = 0; j < data.numStates; j++) {
			cplexB.addEq(exprIN[j], flowStateIn[j]);
			cplexB.addEq(exprOUT[j], flowStateOut[j]);
		}

		// dist
		IloLinearNumExpr objDist = cplexB.linearNumExpr();
		for (int i = 0; i < data.numCounty; i++) {
			for (int j = 0; j < data.numCounty; j++) {
				if(data.distMatrix[i][j] <= 100 && i!=j && !data.countyList.get(i).state.equals(data.countyList.get(j).state)) {
					objDist.addTerm(data.distMatrix[i][j], xR[i][j]);
					objDist.addTerm(data.distMatrix[i][j], xD[i][j]);
				}
			}
		}
		
		cplexB.addMinimize(objDist);
		cplexB.solve();

		xRsol = new double[data.numCounty][data.numCounty];
		xDsol = new double[data.numCounty][data.numCounty];

		double tDist = 0;
		for (int i = 0; i < data.numCounty; i++) {
			for (int j = 0; j < data.numCounty; j++) {
				if(data.distMatrix[i][j] <= 100 && i!=j && !data.countyList.get(i).state.equals(data.countyList.get(j).state) ) {
					xRsol[i][j] = cplexB.getValue(xR[i][j]);
					xDsol[i][j] = cplexB.getValue(xD[i][j]);
					tDist+=data.distMatrix[i][j]*(xRsol[i][j]+xDsol[i][j]);
					if(cplexB.getValue(xR[i][j])+0*cplexB.getValue(xD[i][j])>0.1) {
						System.out.println("MOVE FROM "+i+" to "+j+": "+(cplexB.getValue(xR[i][j])+cplexB.getValue(xD[i][j])));
					}
				}
			}
		}
		System.out.println("WEPUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUU: "+tDist);
		System.out.println("WEPOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO: "+cplexB.getObjValue());
		cplexB.clearModel();
		cplexB.end();

	}

	private void postOpttune(DataHandler data, double lambda, ArrayList<int[]> activeFlows, String party, int pLimit) {

		if(party.equals("R")) {
			//To Allocate remaining capacity
			double tFlow = 0;
			for (int k = 0; k < activeFlows.size(); k++) {
				int origin = activeFlows.get(k)[0];
				int destination = activeFlows.get(k)[1];
				tFlow+=xRsol[origin][destination];
				double auxFlow = 0;
				if(xRsol[origin][destination]<1000) {
					auxFlow = xRsol[origin][destination];
					xRsol[origin][destination] = 0;
					System.out.println(auxFlow+" REMOVED FROM "+origin+" to "+destination);
				}
				if(auxFlow>0) {
					// Relocate the flow
					boolean success = false;
					for (int q = 0; q < activeFlows.size(); q++) {
						int originPrime = activeFlows.get(q)[0];
						int destinationPrime = activeFlows.get(q)[1];
						if(xRsol[originPrime][destinationPrime] >= 1000 && xRsol[originPrime][destinationPrime]+auxFlow <= data.countyList.get(originPrime).minRvotes) {
							xRsol[originPrime][destinationPrime] = xRsol[originPrime][destinationPrime]+auxFlow;
							q = 100000;
							success = true;
							System.out.println(auxFlow+" ADDED "+originPrime+" to "+destinationPrime+" new val: "+xRsol[originPrime][destinationPrime]+" limit: "+data.countyList.get(originPrime).minRvotes);
						}
					}
					if(!success) {//Post opt failed!
						xRsol[origin][destination] = auxFlow;
						System.out.println("Failed to adjust flows!!!");
					}
				}
			}
			// Allocate missing cap
			double auxFlow = pLimit - tFlow;
			if(auxFlow>0) {
				// Relocate the flow
				for (int q = 0; q < activeFlows.size(); q++) {
					int originPrime = activeFlows.get(q)[0];
					int destinationPrime = activeFlows.get(q)[1];
					if(xRsol[originPrime][destinationPrime]+auxFlow <= data.countyList.get(originPrime).minRvotes) {
						xRsol[originPrime][destinationPrime] = xRsol[originPrime][destinationPrime]+auxFlow;
						q =100000;
						System.out.println(auxFlow+" ADDED "+originPrime+" to "+destinationPrime+" new val: "+xRsol[originPrime][destinationPrime]+" limit: "+data.countyList.get(originPrime).minRvotes);
					}
				}
			}
		}
		////////////////////////////////////////////////////////// DEMOCRATS /////////////////////////////////////////////
		else {
			//To Allocate remaining capacity
			double tFlow = 0;
			for (int k = 0; k < activeFlows.size(); k++) {
				int origin = activeFlows.get(k)[0];
				int destination = activeFlows.get(k)[1];
				tFlow+=xDsol[origin][destination];
				double auxFlow = 0;
				if(xDsol[origin][destination]<1000) {
					auxFlow = xDsol[origin][destination];
					xDsol[origin][destination] = 0;
					System.out.println(auxFlow+" REMOVED FROM "+origin+" to "+destination);
				}
				if(auxFlow>0) {
					// Relocate the flow
					boolean success = false;
					for (int q = 0; q < activeFlows.size(); q++) {
						int originPrime = activeFlows.get(q)[0];
						int destinationPrime = activeFlows.get(q)[1];
						if(xDsol[originPrime][destinationPrime] >= 1000 && xDsol[originPrime][destinationPrime]+auxFlow <= data.countyList.get(originPrime).minDvotes) {
							xDsol[originPrime][destinationPrime] = xDsol[originPrime][destinationPrime]+auxFlow;
							q =100000;
							success = true;
							System.out.println(auxFlow+" ADDED "+originPrime+" to "+destinationPrime+" new val: "+xDsol[originPrime][destinationPrime]+" limit: "+data.countyList.get(originPrime).minDvotes);
						}
					}
					if(!success) {//Post opt failed!
						xDsol[origin][destination] = auxFlow;
						System.out.println("Failed to adjust flows!!!");
					}
				}
			}
			// Allocate missing cap
			double auxFlow = pLimit - tFlow;
			if(auxFlow>0) {
				// Relocate the flow
				for (int q = 0; q < activeFlows.size(); q++) {
					int originPrime = activeFlows.get(q)[0];
					int destinationPrime = activeFlows.get(q)[1];
					if(xDsol[originPrime][destinationPrime]+auxFlow <= data.countyList.get(originPrime).minDvotes) {
						xDsol[originPrime][destinationPrime] = xDsol[originPrime][destinationPrime]+auxFlow;
						q =100000;
						System.out.println(auxFlow+" ADDED "+originPrime+" to "+destinationPrime+" new val: "+xDsol[originPrime][destinationPrime]+" limit: "+data.countyList.get(originPrime).minDvotes);
					}
				}
			}
		}


	}

	public void evalOutSampleT(DataHandler data, double lambda, int p, String party, String poll, int q) throws FileNotFoundException {
		int numSim = data.countyList.get(0).rVotesSim.length;
		double[] totalVotes = new double[numSim];
		double[] totalRvotes = new double[numSim];
		double[] totalDvotes = new double[numSim];
		double[] elecVotesR = new double[numSim];
		double[] elecVotesD = new double[numSim];
		double probRwin = 0;
		double probDwin = 0;
		probRout = 0;
		probDout = 0;

		java.io.PrintStream ps = new java.io.PrintStream( new java.io.FileOutputStream("Stats_"+party+"_"+p+"_"+lambda+"_"+poll+"_"+q+".txt"));
		ps.println("StateName,ProbRWinsBefore,ProbDWinsBefore,ProbRwinsAfter,ProbDWinsAfter,TotalVotes,TotalRVotes,TotalDVotes,QR25,QR975,QD25,QD975");
		for (int s = 0; s < data.numStates; s++) {
			checkStatesProbsT(data, s, numSim, ps);
		}
		checkEspecialDistricts(data, numSim, ps);

		ps.println("SimID, TotalVotes, TotalRVotes, TotalDVotes, ElecVotesR, ElecVotesD, RWins, Dwins");
		for (int k = 0; k < numSim; k++) {
			ArrayList<Double> aux =  calcVotesTarget(data, xRsol, xDsol, k); 
			totalRvotes[k] = aux.get(0);
			totalDvotes[k] = aux.get(1);
			elecVotesR[k] = aux.get(2);
			elecVotesD[k] = aux.get(3);
			totalVotes[k] = totalRvotes[k]+totalDvotes[k]; 
			int Rwins = 0; int Dwins = 0;
			if(elecVotesR[k] >= 270) {Rwins = 1; probRwin++;}
			if(elecVotesD[k] >= 270) {Dwins = 1; probDwin++;}
			ps.println(k+","+totalVotes[k]+","+totalRvotes[k]+","+totalDvotes[k]+","+elecVotesR[k]+","+elecVotesD[k]+","+Rwins+","+Dwins);
		}
		ps.println("OUT_SAMPLE_PROB_R_WINS:,"+probRwin);
		ps.println("OUT_SAMPLE_PROB_D_WINS:,"+probDwin);

		System.out.println("OUT SAMPLE PROB R WINS: "+probRwin);
		System.out.println("OUT SAMPLE PROB D WINS: "+probDwin);
		probRout = (probRwin/numSim+0.0);
		probDout = (probDwin/numSim+0.0);

	}

	private void checkEspecialDistricts(DataHandler data, int numSim, PrintStream ps) {
		double[][] TRE = new double[5][numSim];				// Votes per state
		double[][] TDE = new double[5][numSim];
		double[] probRwins = new double[5];
		double[] probDwins = new double[5];
		int[] rWins = new int[5];
		int[] dWins = new int[5];

		for (int k = 0; k < numSim; k++) {
			// Count votes per county
			double[] zR = new double[data.numCounty];		// Votes per county
			double[] zD = new double[data.numCounty];
			for (int o = 0; o < data.numCounty; o++) {
				zR[o] = data.countyList.get(o).rVotesSim[k];
				zD[o] = data.countyList.get(o).dVotesSim[k];
			}

			// Total votes for special districts
			for (int s = 0; s < 5; s++) {
				if(s==0) {
					for (int i = 0; i < data.maine1.size(); i++) {
						TRE[s][k] += zR[data.maine1.get(i)];
						TDE[s][k] += zD[data.maine1.get(i)];
					}
					//Only half of the votes for shared county 23011
					TRE[s][k] -= 0.5*zR[data.countyIndex.get(23011)];
					TDE[s][k] -= 0.5*zD[data.countyIndex.get(23011)];
				}
				else if(s==1) {
					for (int i = 0; i < data.maine2.size(); i++) {
						TRE[s][k] += zR[data.maine2.get(i)];
						TDE[s][k] += zD[data.maine2.get(i)];
					}
					//Only half of the votes for shared county 23011
					TRE[s][k] -= 0.5*zR[data.countyIndex.get(23011)];
					TDE[s][k] -= 0.5*zD[data.countyIndex.get(23011)];
				}
				else if(s==2) {
					for (int i = 0; i < data.nebraska1.size(); i++) {
						TRE[s][k] += zR[data.nebraska1.get(i)];
						TDE[s][k] += zD[data.nebraska1.get(i)];
					}
					//Only half of the votes for shared county 31153
					TRE[s][k] -= 0.5*zR[data.countyIndex.get(31153)];
					TDE[s][k] -= 0.5*zD[data.countyIndex.get(31153)];
				}
				else if(s==3) {
					for (int i = 0; i < data.nebraska2.size(); i++) {
						TRE[s][k] += zR[data.nebraska2.get(i)];
						TDE[s][k] += zD[data.nebraska2.get(i)];
					}
					//Only half of the votes for shared county 31153
					TRE[s][k] -= 0.5*zR[data.countyIndex.get(31153)];
					TDE[s][k] -= 0.5*zD[data.countyIndex.get(31153)];
				}
				else if(s==4) {
					for (int i = 0; i < data.nebraska3.size(); i++) {
						TRE[s][k] += zR[data.nebraska3.get(i)];
						TDE[s][k] += zD[data.nebraska3.get(i)];
					}
				}

			}

			for (int s = 0; s < 5; s++) {
				if(TRE[s][k] > TDE[s][k]) {rWins[s]++;}
				else {dWins[s]++;}
			}
			//if(TR[k] < 0 || TD[k]<0) {System.out.println("ERROR "+data.states[s]+" TotalRvotes "+TR[k]+" TotalDvotes "+TD[k]+" sim_"+k);}
		}
		for (int s = 0; s < 5; s++) {
			probRwins[s] = rWins[s]/(numSim+0.0);
			probDwins[s] = dWins[s]/(numSim+0.0);
			//ps.print("BEFORE_DISTRICT,"+s+","+probRwins[s]+","+probDwins[s]+",");
		}


		// Compute votes per state for each simulation AFTER MOVE
		TRE = new double[5][numSim];				// Votes per state
		TDE = new double[5][numSim];
		probRwins = new double[5];
		probDwins = new double[5];
		rWins = new int[5];
		dWins = new int[5];
		double[] avgRvotes = new double[5];
		double[] avgDvotes = new double[5];

		for (int k = 0; k < numSim; k++) {
			double[] zR = new double[data.numCounty];		// Votes per county
			double[] zD = new double[data.numCounty];
			// Count votes per county
			for (int o = 0; o < data.numCounty; o++) {
				zR[o] += data.countyList.get(o).rVotesSim[k];
				zD[o] += data.countyList.get(o).dVotesSim[k];
				for (int d = 0; d < data.numCounty; d++) {
					if(xRsol[o][d]>0.001 || xDsol[o][d]>0.001) {
						zR[o]-= xRsol[o][d];
						zR[d]+= xRsol[o][d];
						zD[o]-= xDsol[o][d];
						zD[d]+= xDsol[o][d];
					}
				}
				//System.out.println("COUNTY "+o+" "+zR[o]+" vs "+zD[o]);	
			}

			// Total votes for special districts
			for (int s = 0; s < 5; s++) {
				if(s==0) {
					for (int i = 0; i < data.maine1.size(); i++) {
						TRE[s][k] += zR[data.maine1.get(i)];
						TDE[s][k] += zD[data.maine1.get(i)];
					}
					//Only half of the votes for shared county 23011
					TRE[s][k] -= 0.5*zR[data.countyIndex.get(23011)];
					TDE[s][k] -= 0.5*zD[data.countyIndex.get(23011)];
				}
				else if(s==1) {
					for (int i = 0; i < data.maine2.size(); i++) {
						TRE[s][k] += zR[data.maine2.get(i)];
						TDE[s][k] += zD[data.maine2.get(i)];
					}
					//Only half of the votes for shared county 23011
					TRE[s][k] -= 0.5*zR[data.countyIndex.get(23011)];
					TDE[s][k] -= 0.5*zD[data.countyIndex.get(23011)];
				}
				else if(s==2) {
					for (int i = 0; i < data.nebraska1.size(); i++) {
						TRE[s][k] += zR[data.nebraska1.get(i)];
						TDE[s][k] += zD[data.nebraska1.get(i)];
					}
					//Only half of the votes for shared county 31153
					TRE[s][k] -= 0.5*zR[data.countyIndex.get(31153)];
					TDE[s][k] -= 0.5*zD[data.countyIndex.get(31153)];
				}
				else if(s==3) {
					for (int i = 0; i < data.nebraska2.size(); i++) {
						TRE[s][k] += zR[data.nebraska2.get(i)];
						TDE[s][k] += zD[data.nebraska2.get(i)];
					}
					//Only half of the votes for shared county 31153
					TRE[s][k] -= 0.5*zR[data.countyIndex.get(31153)];
					TDE[s][k] -= 0.5*zD[data.countyIndex.get(31153)];
				}
				else if(s==4) {
					for (int i = 0; i < data.nebraska3.size(); i++) {
						TRE[s][k] += zR[data.nebraska3.get(i)];
						TDE[s][k] += zD[data.nebraska3.get(i)];
					}
				}

			}

			for (int s = 0; s < 5; s++) {
				if(TRE[s][k] > TDE[s][k]) {rWins[s]++;}
				else {dWins[s]++;}
				avgRvotes[s]+=TRE[s][k];
				avgDvotes[s]+=TDE[s][k];
			}


		}

		for (int s = 0; s < 5; s++) {
			probRwins[s] = rWins[s]/(numSim+0.0);
			probDwins[s] = dWins[s]/(numSim+0.0);
			avgRvotes[s] = avgRvotes[s]/(numSim+0.0);
			avgDvotes[s] = avgDvotes[s]/(numSim+0.0);

			//Get quantiles
			Arrays.sort(TRE[s]);
			double qR25 = (0.975)*TRE[s][124]+(0.025)*TRE[s][125];
			double qR975 = (0.025)*TRE[s][4874]+(0.975)*TRE[s][4875];

			Arrays.sort(TDE[s]);
			double qD25 = (0.975)*TDE[s][124]+(0.025)*TDE[s][125];
			double qD975 = (0.025)*TDE[s][4874]+(0.975)*TDE[s][4875];


			ps.println("AFTER_DISTRICT,"+s+probRwins[s]+","+probDwins[s]+","+avgRvotes[s]+","+avgDvotes[s]+","+qR25+","+qR975+","+qD25+","+qD975);

			System.out.println(data.states[s]+"DONE!");
		}


	}

	private ArrayList<Double> calcVotesTarget(DataHandler data, double[][] xR, double[][] xD, int k) {
		ArrayList<Double> results = new ArrayList<Double>();

		double[] zR = new double[data.numCounty];				// Votes per county
		double[] zD = new double[data.numCounty];
		double[] TR = new double[data.numStates];				// Votes per state
		double[] TD = new double[data.numStates];
		double elecR = 0;
		double elecD = 0;
		// Special districts 0=Maine1, 1=Maine2, ... 4=Nebraska3
		double[] TRE = new double[5];				// Votes per special district
		double[] TDE = new double[5];

		// Count votes per county
		for (int o = 0; o < data.numCounty; o++) {
			zR[o] += data.countyList.get(o).rVotesSim[k];
			zD[o] += data.countyList.get(o).dVotesSim[k];
			for (int d = 0; d < data.numCounty; d++) {
				if(xR[o][d]>0.001 || xD[o][d]>0.001) {
					zR[o]-= xR[o][d];
					zR[d]+= xR[o][d];
					zD[o]-= xD[o][d];
					zD[d]+= xD[o][d];
				}
			}

			//System.out.println("COUNTY "+o+" "+zR[o]+" vs "+zD[o]);	
		}

		// Total votes for a state
		double totalRvotes = 0;
		double totalDvotes = 0;
		for (int s = 0; s < data.numStates; s++) {
			for (int i = 0; i < data.numCounty; i++) {
				if(data.countyList.get(i).state.equals(data.states[s])) {
					TR[s] += zR[i];
					TD[s] += zD[i];
				}
			}

			totalRvotes+=TR[s];
			totalDvotes+=TD[s];
		}

		// Total votes for special districts
		for (int s = 0; s < 5; s++) {
			if(s==0) {
				for (int i = 0; i < data.maine1.size(); i++) {
					TRE[s] += zR[data.maine1.get(i)];
					TDE[s] += zD[data.maine1.get(i)];
				}
				//Only half of the votes for shared county 23011
				TRE[s] -= 0.5*zR[data.countyIndex.get(23011)];
				TDE[s] -= 0.5*zD[data.countyIndex.get(23011)];
			}
			else if(s==1) {
				for (int i = 0; i < data.maine2.size(); i++) {
					TRE[s] += zR[data.maine2.get(i)];
					TDE[s] += zD[data.maine2.get(i)];
				}
				//Only half of the votes for shared county 23011
				TRE[s] -= 0.5*zR[data.countyIndex.get(23011)];
				TDE[s] -= 0.5*zD[data.countyIndex.get(23011)];
			}
			else if(s==2) {
				for (int i = 0; i < data.nebraska1.size(); i++) {
					TRE[s] += zR[data.nebraska1.get(i)];
					TDE[s] += zD[data.nebraska1.get(i)];
				}
				//Only half of the votes for shared county 31153
				TRE[s] -= 0.5*zR[data.countyIndex.get(31153)];
				TDE[s] -= 0.5*zD[data.countyIndex.get(31153)];
			}
			else if(s==3) {
				for (int i = 0; i < data.nebraska2.size(); i++) {
					TRE[s] += zR[data.nebraska2.get(i)];
					TDE[s] += zD[data.nebraska2.get(i)];
				}
				//Only half of the votes for shared county 31153
				TRE[s] -= 0.5*zR[data.countyIndex.get(31153)];
				TDE[s] -= 0.5*zD[data.countyIndex.get(31153)];
			}
			else if(s==4) {
				for (int i = 0; i < data.nebraska3.size(); i++) {
					TRE[s] += zR[data.nebraska3.get(i)];
					TDE[s] += zD[data.nebraska3.get(i)];
				}
			}

		}

		// Total Electoral Votes
		for (int s = 0; s < data.numStates; s++) {
			if(TR[s] > TD[s]+0.1) {
				elecR+=data.eVotes[s];
			}
			else {
				elecD+=data.eVotes[s];
			}
		}
		for (int s = 0; s < 5; s++) {
			if(TRE[s] > TDE[s]+0.1) {
				elecR+=1;
			}
			else {
				elecD+=1;
			}
		}

		results.add(totalRvotes);
		results.add(totalDvotes);
		results.add(elecR);
		results.add(elecD);



		return results;
	}

	private void checkStatesProbsT(DataHandler data, int s, int numSim, PrintStream ps){
		// Compute votes per state for each simulation BEFORE MOVE
		double[] TR = new double[numSim];				// Votes per state
		double[] TD = new double[numSim];
		double probRwins = 0;
		double probDwins = 0;
		int rWins = 0;
		int dWins = 0;

		for (int k = 0; k < numSim; k++) {
			// Count votes per county
			double[] zR = new double[data.numCounty];		// Votes per county
			double[] zD = new double[data.numCounty];
			for (int o = 0; o < data.numCounty; o++) {
				zR[o] = data.countyList.get(o).rVotesSim[k];
				zD[o] = data.countyList.get(o).dVotesSim[k];
			}

			// Total votes for a state
			for (int i = 0; i < data.numCounty; i++) {
				if(data.countyList.get(i).state.equals(data.states[s])) {
					TR[k] += zR[i];
					TD[k] += zD[i];
				}
			}
			if(TR[k] > TD[k]) {rWins++;}
			else {dWins++;}

			//if(TR[k] < 0 || TD[k]<0) {System.out.println("ERROR "+data.states[s]+" TotalRvotes "+TR[k]+" TotalDvotes "+TD[k]+" sim_"+k);}
		}
		probRwins = rWins/(numSim+0.0);
		probDwins = dWins/(numSim+0.0);

		ps.print(data.states[s]+","+probRwins+","+probDwins+",");

		// Compute votes per state for each simulation AFTER MOVE
		TR = new double[numSim];				// Votes per state
		TD = new double[numSim];
		probRwins = 0;
		probDwins = 0;
		rWins = 0;
		dWins = 0;
		double avgRvotes = 0;
		double avgDvotes = 0;

		for (int k = 0; k < numSim; k++) {
			double[] zR = new double[data.numCounty];		// Votes per county
			double[] zD = new double[data.numCounty];
			// Count votes per county
			for (int o = 0; o < data.numCounty; o++) {
				zR[o] += data.countyList.get(o).rVotesSim[k];
				zD[o] += data.countyList.get(o).dVotesSim[k];
				for (int d = 0; d < data.numCounty; d++) {
					if(xRsol[o][d]>0.001 || xDsol[o][d]>0.001) {
						zR[o]-= xRsol[o][d];
						zR[d]+= xRsol[o][d];
						zD[o]-= xDsol[o][d];
						zD[d]+= xDsol[o][d];
					}
				}
				//System.out.println("COUNTY "+o+" "+zR[o]+" vs "+zD[o]);	
			}

			// Total votes for a state
			for (int i = 0; i < data.numCounty; i++) {
				if(data.countyList.get(i).state.equals(data.states[s])) {
					TR[k] += zR[i];
					TD[k] += zD[i];
				}
			}

			if(TR[k] > TD[k]) {rWins++;}
			else {dWins++;}
			//Get avg
			avgRvotes+=TR[k];
			avgDvotes+=TD[k];

			//if(TR[k] < 0 || TD[k]<0) {System.out.println("ERROR "+data.states[s]+" TotalRvotes "+TR[k]+" TotalDvotes "+TD[k]+" sim_"+k);}
		}

		probRwins = rWins/(numSim+0.0);
		probDwins = dWins/(numSim+0.0);
		avgRvotes = avgRvotes/(numSim+0.0);
		avgDvotes = avgDvotes/(numSim+0.0);
		double totalV = avgRvotes+avgDvotes;

		//Get quantiles
		Arrays.sort(TR);
		double qR25 = (0.975)*TR[124]+(0.025)*TR[125];
		double qR975 = (0.025)*TR[4874]+(0.975)*TR[4875];

		Arrays.sort(TD);
		double qD25 = (0.975)*TD[124]+(0.025)*TD[125];
		double qD975 = (0.025)*TD[4874]+(0.975)*TD[4875];


		ps.println(probRwins+","+probDwins+","+totalV+","+avgRvotes+","+avgDvotes+","+qR25+","+qR975+","+qD25+","+qD975);
		System.out.println(data.states[s]+"DONE!");

	}







}

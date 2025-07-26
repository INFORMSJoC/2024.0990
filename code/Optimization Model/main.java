package FluidModels;
import java.io.IOException;

import org.apache.poi.openxml4j.exceptions.InvalidFormatException;

import ilog.concert.IloException;

public class main {

	public static void main(String[] args) throws InvalidFormatException, IOException, IloException {
		// TODO Auto-generated method stub

		long aTime = System.nanoTime();
		java.io.PrintStream ps = new java.io.PrintStream( new java.io.FileOutputStream("Results.txt", true));
		//ps.println("GroundTruth,TotalFlow,lambda,modelProb,outSampleProb_NS,outSampleProb_HILL,outSampleProb_RACE");
		ps.println("GroundTruth,TotalFlow,lambda,modelProb,outSampleProb");

		DataHandler data = new DataHandler("InputData", 2024, 1); // 1 for NS, 2 for Hill, 3 for RACE
		data.readDist("county_dist.csv");

		int[] P = {50000}; //Movement limit
		double[] Deltas = {1};
		// Select the party to optimize and the ground truth 
		String party = "R";
		String GT = "NS";

		for (int l = 0; l < Deltas.length; l++) {
			double delta = Deltas[l];

			for (int p = 0; p < P.length; p++) {
				// Read the initial simulations
				data.readSim("Calibration1000.csv", 1000);
				// Read data from the polls
				data.updatePredictions("InputData", 1); // 1 for NS, 2 for Hill, 3 for RACE
				// Run the bayesian upates and the optimization flow model
				AlgorithmHandler alg = new AlgorithmHandler(data);
				alg.updateSim(data);
				alg.modelTargetSAAProb(data, delta, 100, P[p], party, GT);

				// Evaluation out of sample with new data
				if(party.equals("R")) {ps.print(GT+","+alg.flow+","+delta+","+alg.probRin);}
				else{ps.print(GT+","+alg.flow+","+delta+","+alg.probDin);}
				// Check robustness with each poll 
				for(int q = 1; q <= 3; q++) { // 1 for NS, 2 for Hill, 3 for RACE
					// Read the evaluation sims
					data.readValidation("Calibration5000.csv", 5000);
					// Read data from the polls 
					data.updatePredictions("InputData", q);
					// Run the bayesian update 
					alg.updateSim(data);
					// Evaluate the solution from the opt model with the sims 
					alg.evalOutSampleT(data, delta, P[p], party, GT, q);

					if(party.equals("R")) {ps.print(","+alg.probRout);}
					else{ps.print(","+alg.probDout);}
				}
				ps.println();
				System.out.println(" Execution Time: "+((System.nanoTime()-aTime)/1000000000)); 
			}
		}
	}
}

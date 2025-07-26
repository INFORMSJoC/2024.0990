package FluidModels;

public class County {

	int ID;
	String name;
	String state;
	int pop;
	
	double[] rVotesSim;
	double[] dVotesSim;
	double avgRfrac;
	double avgDfrac;
	double minRvotes;
	double minDvotes;
	
	public County() {
		
		ID = -1;
		name = "";
		state= "";
		pop = 0;
		avgRfrac = 0;
		avgDfrac = 0;
		minRvotes = 1000000000;
		minDvotes = 1000000000;
	}

	public void print() {
		System.out.println(name+","+ID+","+state+","+pop);
		
	}

	public void reset(int numSim) {
		rVotesSim = new double[numSim];
		dVotesSim = new double[numSim];
		avgRfrac = 0;
		avgDfrac = 0;
	}

}

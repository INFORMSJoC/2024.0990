package FluidModels;
import java.io.BufferedReader;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileReader;
import java.io.IOException;
import java.io.InputStream;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.Iterator;
import java.util.Random;

import org.apache.poi.openxml4j.exceptions.InvalidFormatException;
import org.apache.poi.ss.usermodel.Cell;
import org.apache.poi.ss.usermodel.Row;
import org.apache.poi.ss.usermodel.Sheet;
import org.apache.poi.ss.usermodel.Workbook;
import org.apache.poi.ss.usermodel.WorkbookFactory;

public class DataHandler {

	int numCounty;							// Number of counties 
	int numStates = 51;						// Number of States
	String[]  states = new String[numStates];		// State names
	HashMap<String, Integer> stateIndex;		// Map between ID and position in the List
	int[] eVotes = new int[numStates]; 			// Electoral votes for each state
	ArrayList<County> countyList;				// Counties with all their info
	HashMap<Integer, Integer> countyIndex;		// Map between ID and position in the List
	double [][] distMatrix;							// Distance matrix
	double [] PredR = new double[numStates];		// NS prediction of R wining probability
	double [] PredD = new double[numStates];		// NS prediction of D wining probability
	double [] statePop = new double[numStates];

	// Special districts 0=Maine1, 1=Maine2, 2=Nebraska1, 3=Nebraska2, 4=Nebraska3
	ArrayList<Integer> maine1;						// Maine special district 1
	ArrayList<Integer> maine2;						// Maine special district 2
	ArrayList<Integer> nebraska1;				// Maine special district 2
	ArrayList<Integer> nebraska2;				// Maine special district 2
	ArrayList<Integer> nebraska3;				// Maine special district 2
	double [] NSPredRE = new double[5];		// NS prediction of R wining probability
	double [] NSPredDE = new double[5];		// NS prediction of D wining probability


	public DataHandler(String fileName, int year, int poll) throws InvalidFormatException, IOException {
		maine1 = new ArrayList<Integer>();
		maine2 = new ArrayList<Integer>();
		nebraska1 = new ArrayList<Integer>();
		nebraska2 = new ArrayList<Integer>();
		nebraska3 = new ArrayList<Integer>();
		stateIndex = new HashMap<String, Integer>();

		if(year == 2024) {
			InputStream inp = new FileInputStream("data/"+fileName+".xlsx");

			Workbook wb = WorkbookFactory.create(inp);
			// READ STATE data
			Sheet sheet = wb.getSheetAt(2);
			Row row = sheet.getRow(0);
			Cell cell = row.getCell(0);

			for(int i = 0; i< states.length; i++) {
				row = sheet.getRow(i+1);
				cell = row.getCell(0);
				states[i] = cell.getStringCellValue();
				stateIndex.put(states[i], i);

				row = sheet.getRow(i+1);
				cell = row.getCell(1);
				eVotes[i] = (int) cell.getNumericCellValue();

				//CHANGE HERE TO READ PREDICTIONS NS (6 and 7), Hill (8 and 9), RaceToWhiteHouse (10 and 11)
				int index1 = -1;
				int index2 = -1;
				if(poll ==1) {index1 = 6; index2 = 7;}
				if(poll ==2) {index1 = 8; index2 = 9;}
				if(poll ==3) {index1 = 10; index2 = 11;}
				row = sheet.getRow(i+1);		
				cell = row.getCell(index1);
				PredR[i] = cell.getNumericCellValue();

				row = sheet.getRow(i+1);
				cell = row.getCell(index2);
				PredD[i] = cell.getNumericCellValue();

				//System.out.println("Hello I just read "+i+": "+states[i]+" and "+eVotes[i]+" PollR: "+PredR[i]+" PollD: "+PredD[i]+" Pollster: "+poll );
			}

			for(int i = states.length; i< states.length+5; i++) {
				//CHANGE HERE TO READ PREDICTIONS NS (6 and 7), Hill (8 and 9), RaceToWhiteHouse (10 and 11)
				int index1 = -1;
				int index2 = -1;
				if(poll ==1) {index1 = 6; index2 = 7;}
				if(poll ==2) {index1 = 8; index2 = 9;}
				if(poll ==3) {index1 = 10; index2 = 11;}
				row = sheet.getRow(i+1);
				cell = row.getCell(index1);
				NSPredRE[i-states.length] = cell.getNumericCellValue();

				row = sheet.getRow(i+1);
				cell = row.getCell(index2);
				NSPredDE[i-states.length] = cell.getNumericCellValue();
				countyIndex = new HashMap<Integer, Integer>();

				//System.out.println("Hello I just read District"+(i-states.length)+" NSR: "+NSPredRE[i-states.length]+" NSD: "+NSPredDE[i-states.length] );
			}

			// READ COUNTY data
			countyList = new ArrayList<>();
			numCounty = 0;
			countyIndex = new HashMap<Integer, Integer>();
			sheet = wb.getSheetAt(0);

			for(int i = 0; i < 3150; i++) {
				County c = new County();
				row = sheet.getRow(i+1);
				cell = row.getCell(0);
				c.name = cell.getStringCellValue();

				row = sheet.getRow(i+1);
				cell = row.getCell(1);
				c.ID = (int) cell.getNumericCellValue();

				row = sheet.getRow(i+1);
				cell = row.getCell(2);
				c.state = cell.getStringCellValue();

				row = sheet.getRow(i+1);
				cell = row.getCell(4);
				c.pop = (int) cell.getNumericCellValue();

				countyList.add(c);
				countyIndex.put(c.ID, numCounty);
				numCounty++;
				//c.print();

			}
			//Read special districts for Maine
			sheet = wb.getSheetAt(3);
			for(int i = 1; i <= 6; i++) {

				row = sheet.getRow(i);
				cell = row.getCell(2);
				int ID = (int) cell.getNumericCellValue();
				int index = countyIndex.get(ID);
				maine1.add(index);
				//System.out.println("MAINE SPECIAL 1 "+countyList.get(index).name+" "+countyList.get(index).ID);
			}
			for(int i = 1; i <= 11; i++) {

				row = sheet.getRow(i);
				cell = row.getCell(3);
				int ID = (int) cell.getNumericCellValue();
				int index = countyIndex.get(ID);
				maine2.add(index);
				//System.out.println("MAINE SPECIAL 2 "+countyList.get(index).name+" "+countyList.get(index).ID);
			}

			//Read special districts for Nebraska
			sheet = wb.getSheetAt(4);
			for(int i = 1; i <= 16; i++) {

				row = sheet.getRow(i);
				cell = row.getCell(3);
				int ID = (int) cell.getNumericCellValue();
				int index = countyIndex.get(ID);
				nebraska1.add(index);
				//System.out.println("NEBRASKA SPECIAL 1 "+countyList.get(index).name+" "+countyList.get(index).ID);
			}
			for(int i = 1; i <= 2; i++) {

				row = sheet.getRow(i);
				cell = row.getCell(4);
				int ID = (int) cell.getNumericCellValue();
				int index = countyIndex.get(ID);
				nebraska2.add(index);
				//if(nebraska1.contains(index)) {System.out.println(countyList.get(index).ID+" SHARED 1 and 2");}
				//System.out.println("NEBRASKA SPECIAL 2 "+countyList.get(index).name+" "+countyList.get(index).ID);
			}
			for(int i = 1; i <= 72; i++) {

				row = sheet.getRow(i);
				cell = row.getCell(5);
				int ID = (int) cell.getNumericCellValue();
				int index = countyIndex.get(ID);
				nebraska3.add(index);
				//System.out.println("NEBRASKA SPECIAL 3 "+countyList.get(index).name+" "+countyList.get(index).ID);
				//if(nebraska1.contains(index)) {System.out.println(countyList.get(index).ID+" SHARED 1 and 3");}
				//if(nebraska2.contains(index)) {System.out.println(countyList.get(index).ID+" SHARED 2 and 3");}
			}

			//for (int i = 0; i < numCounty; i++) {
			//		countyList.get(i).print();	
			//}


			System.out.println("DONE READING TOTAL NUMBER OF COUNTIES "+numCounty);
			inp.close();
		}

	}

	public void readDist(String name) throws IOException {
		File file = new File("data/"+name);
		BufferedReader bufRdr = new BufferedReader(new FileReader(file));
		distMatrix = new double[numCounty][numCounty];

		// Read first line
		String line = bufRdr.readLine();
		String[] tokens = line.split(" ");

		while (line!=null) {
			line = bufRdr.readLine();
			if(line!=null) {
				//System.out.println(line);
				tokens = line.split(",");
				String[] IDs = tokens[0].split("_");
				int ID1 = Integer.parseInt(IDs[0]);
				int ID2 = Integer.parseInt(IDs[1]);
				//System.out.println(ID1+" and "+ID2);
				if(countyIndex.get(ID1)!=null && countyIndex.get(ID2)!= null) {
					distMatrix[countyIndex.get(ID1)][countyIndex.get(ID2)] = Double.parseDouble(tokens[1]);
					//System.out.println("WEPA "+distMatrix[countyIndex.get(ID1)][countyIndex.get(ID2)]);
				}
			}
		}
		System.out.println("DONE READING DISTANCES");

		int countMissing = 0;
		for (int i = 0; i < numCounty; i++) {
			for (int j = 0; j < numCounty; j++) {
				if(i!=j && distMatrix[i][j]==0) { //Counties with no distance found
					countMissing++;
					distMatrix[i][j]=100000;
				}
			}
		}
		System.out.println("MISSING ENTRIES: "+countMissing);


	}


	public void readSim(String name, int numSim) throws IOException {
		// Reset all info in counties
		for (int i = 0; i < numCounty; i++) {
			countyList.get(i).reset(numSim);
		}

		File file = new File("data/"+name);
		BufferedReader bufRdr = new BufferedReader(new FileReader(file));

		// Read first line
		String line = bufRdr.readLine();
		String[] tokens = line.split(" ");

		while (line!=null) {
			line = bufRdr.readLine();
			if(line!=null) {
				tokens = line.split(",");
				int ID = Integer.parseInt(tokens[1]);
				int index = countyIndex.get(ID);
				String party = tokens[2];
				//System.out.println(ID+" "+index+" "+party);
				for (int k = 3; k < Math.min(tokens.length, numSim+3); k++) {//Read numSim
					double votes = Double.parseDouble(tokens[k]);
					//System.out.println("sim "+(k-3)+" "+votes);
					if(party.equals("DEMOCRAT")) {
						countyList.get(index).dVotesSim[k-3] = votes;
					}
					else if(party.equals("REPUBLICAN")) {
						countyList.get(index).rVotesSim[k-3] = votes;
					}
				}

			}
		}
		//System.out.println("DONE READING SIMULATIONS");

	}

	public void readValidation(String name, int numSim) throws IOException {
		// Reset all info in counties
		for (int i = 0; i < numCounty; i++) {
			countyList.get(i).reset(numSim);
		}

		File file = new File("data/"+name);
		BufferedReader bufRdr = new BufferedReader(new FileReader(file));

		// Read first line
		String line = bufRdr.readLine();
		String[] tokens = line.split(" ");

		while (line!=null) {
			line = bufRdr.readLine();
			if(line!=null) {
				tokens = line.split(",");
				int ID = Integer.parseInt(tokens[1]);
				int index = countyIndex.get(ID);
				String party = tokens[2];
				//System.out.println(ID+" "+index+" "+party);
				for (int k = 3; k < Math.min(tokens.length, numSim+3); k++) {//Read numSim
					double votes = Double.parseDouble(tokens[k]);
					//System.out.println("sim "+(k-3)+" "+votes);
					if(party.equals("DEMOCRAT")) {
						countyList.get(index).dVotesSim[k-3] = votes;
					}
					else if(party.equals("REPUBLICAN")) {
						countyList.get(index).rVotesSim[k-3] = votes;
					}
				}

			}
		}

		//System.out.println("DONE READING VALIDATION");

	}

	public void updatePredictions(String fileName, int poll) throws InvalidFormatException, IOException {
		PredR = new double[numStates];		// NS prediction of R wining probability
		PredD = new double[numStates];		// NS prediction of D wining probability
		NSPredRE = new double[5];		// NS prediction of R wining probability
		NSPredDE = new double[5];		// NS prediction of D wining probability

		InputStream inp = new FileInputStream("data/"+fileName+".xlsx");
		Workbook wb = WorkbookFactory.create(inp);
		// READ STATE data
		Sheet sheet = wb.getSheetAt(2);
		Row row = sheet.getRow(0);
		Cell cell = row.getCell(0);

		for(int i = 0; i< states.length; i++) {
			//CHANGE HERE TO READ PREDICTIONS NS (6 and 7), Hill (8 and 9), RaceToWhiteHouse (10 and 11)
			int index1 = -1;
			int index2 = -1;
			if(poll == 1) {index1 = 6; index2 = 7;}
			if(poll == 2) {index1 = 8; index2 = 9;}
			if(poll == 3) {index1 = 10; index2 = 11;}
			row = sheet.getRow(i+1);		
			cell = row.getCell(index1);
			PredR[i] = cell.getNumericCellValue();

			row = sheet.getRow(i+1);
			cell = row.getCell(index2);
			PredD[i] = cell.getNumericCellValue();
			System.out.println("Hello I just updated "+i+": "+states[i]+" and "+eVotes[i]+" PollR: "+PredR[i]+" PollD: "+PredD[i]+" Pollster: "+poll );
		}
		for(int i = states.length; i< states.length+5; i++) {
			//CHANGE HERE TO READ PREDICTIONS NS (6 and 7), Hill (8 and 9), RaceToWhiteHouse (10 and 11)
			int index1 = -1;
			int index2 = -1;
			if(poll ==1) {index1 = 6; index2 = 7;}
			if(poll ==2) {index1 = 8; index2 = 9;}
			if(poll ==3) {index1 = 10; index2 = 11;}
			row = sheet.getRow(i+1);
			cell = row.getCell(index1);
			NSPredRE[i-states.length] = cell.getNumericCellValue();

			row = sheet.getRow(i+1);
			cell = row.getCell(index2);
			NSPredDE[i-states.length] = cell.getNumericCellValue();

			System.out.println("Hello I just read District"+(i-states.length)+" PollR: "+NSPredRE[i-states.length]+" PollD: "+NSPredDE[i-states.length] );
		}

		inp.close();

	}

	public void readSimSample(String name, int numSim, int seed) throws IOException {
		Random r = new Random(seed);
		// Reset all info in counties
		for (int i = 0; i < numCounty; i++) {
			countyList.get(i).reset(numSim);
		}

		File file = new File("data/"+name);
		BufferedReader bufRdr = new BufferedReader(new FileReader(file));

		// Read first line
		String line = bufRdr.readLine();
		String[] tokens = line.split(" ");
		int[] simSelected = new int[5000];
		int numSelected = 0;
		int[] indexSelected  = new int[numSim];
		// Do the sampling
		while (numSelected < numSim) {//Read numSim
			int k = 3+r.nextInt(1000-3);
			if(simSelected[k] == 0) {	
				indexSelected[numSelected] = k;
				simSelected[k] = 1;
				numSelected++;
			}
		}
		System.out.println("WEPA "+Arrays.toString(indexSelected));

		while (line!=null) {
			line = bufRdr.readLine();
			if(line!=null) {
				tokens = line.split(",");
				int ID = Integer.parseInt(tokens[1]);
				int index = countyIndex.get(ID);
				String party = tokens[2];
				//System.out.println(ID+" "+index+" "+party);

				for (int q=0; q<numSim; q++) {//Read numSim
					int k = indexSelected[q];
						double votes = Double.parseDouble(tokens[k]);
						//System.out.println("sim "+(k-3)+" "+votes);
						if(party.equals("DEMOCRAT")) {
							countyList.get(index).dVotesSim[q] = votes;
						}
						else if(party.equals("REPUBLICAN")) {
							countyList.get(index).rVotesSim[q] = votes;
						}

					
				}

			}
		}
		//System.out.println("DONE READING SIMULATIONS");

	}





}

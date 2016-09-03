package MainProcess;

import java.io.FileNotFoundException;

import protein.proteinSet;

public class simplySwissIntAct {

	public static void main(String[] args) throws FileNotFoundException {
		// TODO Auto-generated method stub
		proteinSet measure = new proteinSet();
		measure.loadSwissIntAct("../InFile/SwissOri/uniprot_sprot201401.dat");
		measure.outputIntActProtein("IntAct201401");
	}
}

package MainProcess;

import java.io.FileNotFoundException;

import protein.proteinSet;

public class parseTrembl {

	public static void main(String[] args) throws FileNotFoundException {
		// TODO Auto-generated method stub
		String data = args[1];
		proteinSet.parseUniprotDatFile(data);
	}

}

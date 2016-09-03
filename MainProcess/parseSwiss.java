package MainProcess;

import java.io.FileNotFoundException;
import protein.proteinSet;

public class parseSwiss {

	public static void main(String[] args) throws FileNotFoundException 
	{
		// TODO Auto-generated method stub
		
		String data = args[1];
		String swissFileName = "../InFile/SwissOri/uniprot_sprot" + data +".dat";
		proteinSet.LoadSwissAccessMap(swissFileName);
		proteinSet.OutputAccess2NameMap("../InFile/Swiss/ac2Name" + data);
		proteinSet.OutputAccess2AccessMap("../InFile/Swiss/ac2ac" + data);
		

		proteinSet Swiss = new proteinSet();
		Swiss.loadSwissAnnotation(swissFileName);
		Swiss.OutputAnnotation("../InFile/Swiss/Ann" + data);
		
		
		Swiss.loadSwissProteinSequence(swissFileName);
		Swiss.OutputFastaAccessSequence("../InFile/Swiss/Swiss" + data + ".fasta");
	}

}

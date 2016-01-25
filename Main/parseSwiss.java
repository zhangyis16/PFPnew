package Main;

import java.io.FileNotFoundException;
import protein.proteinSet;

public class parseSwiss {

	public static void main(String[] args) throws FileNotFoundException 
	{
		// TODO Auto-generated method stub
		
		String data = args[1];
		String swissFileName = "../InFile/Swiss/uniprot_sprot" + data +".dat";
		proteinSet.LoadSwissAccessMap(swissFileName);
		proteinSet.OutputAccess2NameMap("../InOutFile/ac2Name" + data);
		proteinSet.OutputAccess2AccessMap("../InOutFile/ac2ac" + data);
		

		proteinSet Swiss = new proteinSet();
		Swiss.parseSwissAnnotation(swissFileName);
		Swiss.OutputAnnotation("../InOutFile/SwissAnnotation" + data);
		
		
		Swiss.parseSwissProteinSequence(swissFileName);
		Swiss.OutputFastaSequence("../InOutFile/Swiss" + data + ".fasta");
	}

}

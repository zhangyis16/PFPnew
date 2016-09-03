package MainProcess;

import java.io.FileNotFoundException;

import protein.proteinSet;

public class simplyPubMed 
{

	public static void main(String[] args) throws FileNotFoundException 
	{
		// TODO Auto-generated method stub
		proteinSet measure = new proteinSet();
		measure.loadSwissPubMed("../InFile/SwissOri/uniprot_sprot201401.dat");
		measure.outputPubMedID("PubMed201401");
	}

}

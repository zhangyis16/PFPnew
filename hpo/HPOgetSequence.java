package hpo;

import java.io.FileNotFoundException;
import java.util.HashMap;

import common.Load;
import protein.proteinSet;

public class HPOgetSequence {

	public static void main(String[] args) throws FileNotFoundException {
		// TODO Auto-generated method stub
		HPOSet aHPOSet = new HPOSet();
		aHPOSet.Load("../InFile/HPO/Ontology/HPO2013-09-17.obo");
		proteinSet.LoadAccess2NameMap("../InFile/Swiss/ac2Name201307");
		proteinSet Pset = new proteinSet();

		
		

		
		
		//Pset.addProteinFromFile("../InFile/HPO/Measure/CAFA2/hpoType1Name");
		
		Pset.AddHPOAnnotation("../InFile/HPO/Train/CAFA2/Ann");

		//Pset.mapAccess();
		//Pset.OutputProteinAccess("../InFile/HPO/Measure/CAFA2/Ann");

		//Pset.OutputHPOAnnotation("../InFile/HPO/Measure/CAFA2/Ann");
		
		Pset.loadFastaSequence("../InFile/Swiss/Swiss201401.fasta");
		Pset.loadTableSequence("../InFile/Trembl/SeqTrembl201401");
		Pset.OutputFastaAccessSequence("../InFile/HPO/Train/CAFA2/Seq");
	}

}

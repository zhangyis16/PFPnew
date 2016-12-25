package tempor;

import java.io.FileNotFoundException;

import Main.learningOfGO;
import protein.GoSet;
import protein.proteinSet;

public class demo {

	public static void main(String[] args) throws FileNotFoundException 
	{
		// TODO Auto-generated method stub
		learningOfGO.aGoSet = new GoSet("../InFile/GeneOntology/gene_ontology_edit.obo.2013-06-15");
		
		String Tstart = "201301";
		String Tend =   "201401";
		proteinSet proteinSetStart = new proteinSet();
		proteinSet proteinSetEnd = new proteinSet();
		proteinSet AnnType1 = new proteinSet();
		
		
		proteinSet.LoadAccess2NameMap("../InFile/Swiss/ac2Name" + Tend);
		proteinSet.LoadSwissMapAccess2UniAccess("../InFile/Swiss/ac2ac" + Tend);
		
		proteinSetStart.AddAnnotation("../InFile/Swiss/Ann" + Tstart);
		proteinSetEnd.AddAnnotation("../InFile/Swiss/Ann" + Tend);
		//AnnType1 = proteinSet.getNewProtein(proteinSetStart,proteinSetEnd);
		
		AnnType1 = proteinSet.getNewProtein(Tstart, Tend);
		
		AnnType1.filterCAFA2Species();
		
		AnnType1.UpdateIndex();
		AnnType1.OutputProteinAccess("../InFile/temporary/Ann");	
		AnnType1.calMFBPCCsize(learningOfGO.aGoSet);
		
		proteinSet MFSub =	AnnType1.getSubAnnotation('F');
		MFSub.OutputIDAnnSpecies("MFOIDAnnSpecies");
		
		proteinSet BPSub =	AnnType1.getSubAnnotation('P');
		BPSub.OutputIDAnnSpecies("BPOIDAnnSpecies");
		
		proteinSet CCSub =	AnnType1.getSubAnnotation('C');
		CCSub.OutputIDAnnSpecies("CCOIDAnnSpecies");
		
		
		
		AnnType1.loadFastaSequence("../InFile/Swiss/Swiss201401.fasta");
		AnnType1.OutputFastaAccessSequence("seq");
	}
}

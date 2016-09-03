package MainProcess;

import java.io.FileNotFoundException;

import Main.learning;
import protein.GoSet;
import protein.proteinSet;

public class getCAFA2Type1NewProtein {

	
	public static void main(String[] args) throws FileNotFoundException {
		// TODO Auto-generated method stub
		String Tstart = "201401";
		String Tend =   "201409";
		learning.aGoSet = new GoSet("../InFile/GeneOntology/gene_ontology_edit.obo.2013-06-15");
		

		
		proteinSet CAFA2target = new proteinSet();
		proteinSet CAFA2type1benchmark = new proteinSet();
		
		proteinSet AnnType1 = proteinSet.getNewProtein(Tstart, Tend);
		AnnType1.removeGoNotIn(learning.aGoSet);
		
		System.out.println(AnnType1.size());
		
		AnnType1.filterCAFA2Species();
		
		System.out.println(AnnType1.size());
		proteinSet.LoadAccess2NameMap("../InFile/Swiss/ac2Name" + "201401");
		
		
		CAFA2target.addProteinFromFile("../InFile/CAFA2target100816");
		CAFA2target.OutputProteinAccess("../InFile/CAFA2targetAccess");
		AnnType1.getIntersection(CAFA2target);
		System.out.println(AnnType1.size());
		
		AnnType1.OutputAnnotation("myCAFA2measure");
		CAFA2type1benchmark.AddAnnotation("../InFile/CAFA2/type1annotation");
		proteinSet.compare2Set(AnnType1, CAFA2type1benchmark, "./" ,learning.aGoSet);
	}
}

package Main;

import java.io.FileNotFoundException;

import protein.GoSet;
import protein.proteinSet;

public class getType1NewProtein {

	
	public static void main(String[] args) throws FileNotFoundException {
		// TODO Auto-generated method stub
		String t1 = "201409";
		String t0 = "201401";
		learning.aGoSet = new GoSet("../InFile/gene ontology/gene_ontology_edit.obo.2013-07-01");
		

		
		proteinSet CAFA2target = new proteinSet();
		proteinSet CAFA2type1benchmark = new proteinSet();
		
		proteinSet AnnType1 = proteinSet.getNewProtein(t0, t1, learning.aGoSet);
		System.out.println(AnnType1.size());
		AnnType1.filterCAFA2Species();
		System.out.println(AnnType1.size());
		proteinSet.LoadAccess2NameMap("../InFile/Swiss/ac2Name" + t1);
		CAFA2target.AddProtein("../InFile/CAFA2target");
		AnnType1.getIntersection(CAFA2target);
		System.out.println(AnnType1.size());
		
		AnnType1.OutputAnnotation("myCAFA2measure");
		CAFA2type1benchmark.AddAnnotation("../InFile/CAFA2/type1annotation");
		proteinSet.compare2Set(AnnType1, CAFA2type1benchmark, "./" ,learning.aGoSet);
	}

}

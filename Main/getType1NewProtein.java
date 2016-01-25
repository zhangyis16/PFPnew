package Main;

import java.io.FileNotFoundException;

import protein.GoSet;
import protein.proteinSet;

public class getType1NewProtein {

	
	public static void main(String[] args) throws FileNotFoundException {
		// TODO Auto-generated method stub
		learning.aGoSet = new GoSet("../InFile/gene ontology/gene_ontology_edit.obo.2013-07-01");
		
		proteinSet Swiss201401 = new proteinSet();
		proteinSet Swiss201409 = new proteinSet();
		proteinSet AnnType1 = new proteinSet();
		proteinSet CAFA2target = new proteinSet();
		proteinSet CAFA2type1benchmark = new proteinSet();
		
		proteinSet.LoadSwissMapAccess2UniAccess("../InOutFile/ac2ac201401");
		proteinSet.LoadAccess2NameMap("../InFile/Swiss/ac2Name201401");
		
		Swiss201401.AddAnnotationInSwiss("../InFile/Swiss/SwissAnnotation201401");
		Swiss201401.AddAnnotationInSwiss("../InFile/Goa/Annotation201401");
		Swiss201401.AddAnnotationInSwiss("../InFile/GODB/GODB201401");
		Swiss201401.eraserProteinOnly5515();
		Swiss201401.filterCAFA2Species();
		Swiss201401.removeGoNotIn(learning.aGoSet);
		
		
		Swiss201409.AddAnnotationInSwiss("../InFile/Swiss/SwissAnnotation201409");
		Swiss201409.AddAnnotationInSwiss("../InFile/Goa/Annotation201409");
		Swiss201409.AddAnnotationInSwiss("../InFile/GODB/GODB201409");
		Swiss201409.eraserProteinOnly5515();
		Swiss201409.filterCAFA2Species();
		Swiss201409.removeGoNotIn(learning.aGoSet);
		
		AnnType1 = proteinSet.getNewProtein(Swiss201409,Swiss201401);
		
		proteinSet.LoadAccess2NameMap("../InFile/Swiss/ac2Name201409");
		
		CAFA2target.AddProtein("../InFile/CAFA2target");
		AnnType1.getIntersection(CAFA2target);
		
		CAFA2type1benchmark.AddAnnotation("..InFile//CAFA2/type1annotation");
		proteinSet.compareSpecies(AnnType1,CAFA2type1benchmark,"Resulttype1.csv");
		proteinSet.listCHABIE(AnnType1,CAFA2type1benchmark,"chabietype1.csv");
		
		//Swiss201401.removeLabelRareProtein(3);   //去掉标注少于2的蛋白
		//System.out.println(Swiss201401.size());
		//Swiss201401.statSpecies();
		//Swiss201401.statLabelNum();
		//Swiss201401.addFather(Go201307);
		//Swiss201401.statLabelNum();
		//System.out.println(Swiss201401.size());
	}

}

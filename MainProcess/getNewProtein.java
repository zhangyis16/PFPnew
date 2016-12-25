package MainProcess;
import protein.*;

import java.util.*;
import java.io.*;

public class getNewProtein
{
	public static void main(String[] args) throws FileNotFoundException
	{
		proteinSet Ann201401 =   new proteinSet();
		proteinSet Ann201409 = 	 new proteinSet();
		proteinSet AnnType1  = 	 new proteinSet();
		proteinSet Type2MF = 	 new proteinSet();
		proteinSet Type2BP = 	 new proteinSet();
		proteinSet Type2CC = 	 new proteinSet();
		proteinSet CAFA2target = new proteinSet();
		proteinSet CAFA2type1benchmark = new proteinSet();
		proteinSet CAFA2type2MFbenchmark = new proteinSet();
		proteinSet CAFA2type2BPbenchmark = new proteinSet();
		proteinSet CAFA2type2CCbenchmark = new proteinSet();
		GoSet Go201307 = new GoSet("../InFile/gene ontology/gene_ontology_edit.obo.2013-07-01");
		
		Ann201401.AddAnnotation("../InFile/Swiss/Ann201401");
		Ann201401.AddAnnotation("../InFile/Goa/Ann201401");
		Ann201401.AddAnnotation("../InFile/GODB/Ann201401");
		
		Ann201409.AddAnnotation("../InFile/Swiss/Ann201409");
		Ann201409.AddAnnotation("../InFile/Goa/Ann201409");
		Ann201409.AddAnnotation("../InFile/GODB/Ann201409");
		
		Ann201401.filterProteinOnly5515();
		Ann201409.filterProteinOnly5515();
		
		AnnType1 = proteinSet.getNewProtein(Ann201409,Ann201401);
		System.out.println(AnnType1.size());
		
		CAFA2target.addProteinFromFile("../InFile/CAFA2target");
		
		AnnType1.getIntersection(CAFA2target);
		
		System.out.println(AnnType1.size());
		
		AnnType1.OutputAnnotation("myCAFA2measure");
		
		
		
		CAFA2type1benchmark.AddAnnotation("CAFA2/type1annotation");
		
//		Type2MF = proteinSet.getLimitKnowProtein(Ann201410,Ann201401,'F');
//		Type2BP = proteinSet.getLimitKnowProtein(Ann201410,Ann201401,'P');
//  	Type2CC = proteinSet.getLimitKnowProtein(Ann201410,Ann201401,'C');
//		Type2MF.getIntersection(CAFA2target);
//		Type2BP.getIntersection(CAFA2target);
//		Type2CC.getIntersection(CAFA2target);
		
//		Type2BP.OutputAnnotationList("BPmylist");
		
//		CAFA2type2MFbenchmark.AddAnnotation("CAFA2/MFOtype2Annotation");
//		CAFA2type2BPbenchmark.AddAnnotation("CAFA2/BPOtype2Annotation");
//		CAFA2type2CCbenchmark.AddAnnotation("CAFA2/CCOtype2Annotation");
		
		proteinSet.compareSpecies(AnnType1,CAFA2type1benchmark,"Resulttype1.csv");
//		proteinSet.compareSpecies(Type2MF,CAFA2type2MFbenchmark,"ResultMF.csv");
//		proteinSet.compareSpecies(Type2BP,CAFA2type2BPbenchmark,"ResultBP.csv");
//		proteinSet.compareSpecies(Type2CC,CAFA2type2CCbenchmark,"ResultCC.csv");
		
//		proteinSet.listCHABIE(AnnType1,CAFA2type1benchmark,"chabietype1.csv");
//		proteinSet.listCHABIE(Type2MF,CAFA2type2MFbenchmark,"chabieMF.csv");
//		proteinSet.listCHABIE(Type2BP,CAFA2type2BPbenchmark,"chabieBP.csv");
//		proteinSet.listCHABIE(Type2CC,CAFA2type2CCbenchmark,"chabieCC.csv");
		
//		AnnType1.getIntersection(CAFA2type1benchmark);
//		CAFA2type1benchmark.getIntersection(AnnType1);   //取交叉蛋白
//		AnnType1.eraserAnnotation(5515);
//		CAFA2type1benchmark.eraserAnnotation(5515);		
//		proteinSet.compare(AnnType1,CAFA2type1benchmark,"type1");   
		
/*		Type2MF.getIntersection(CAFA2type2MFbenchmark);
		CAFA2type2MFbenchmark.getIntersection(Type2MF);   //取交叉蛋白
		Type2MF.eraserAnnotation(5515);
		CAFA2type2MFbenchmark.eraserAnnotation(5515);		
		proteinSet.compare(Type2MF,CAFA2type2MFbenchmark,"MFtype2");
		
		Type2BP.getIntersection(CAFA2type2BPbenchmark);
		CAFA2type2BPbenchmark.getIntersection(Type2BP);   //取交叉蛋白

		proteinSet.compare(Type2BP,CAFA2type2BPbenchmark,"BPtype2");
		
		Type2CC.getIntersection(CAFA2type2CCbenchmark);
		CAFA2type2CCbenchmark.getIntersection(Type2CC);   //取交叉蛋白

		proteinSet.compare(Type2CC,CAFA2type2CCbenchmark,"CCtype2");   */
		/*AnnType1.OutputAnnotation("../OutFile/type1");
		Type2MF.OutputAnnotation("../OutFile/type2MF");
		Type2BP.OutputAnnotation("../OutFile/type2BP");
		Type2CC.OutputAnnotation("../OutFile/type2CC");   */
	}
}
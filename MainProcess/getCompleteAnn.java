package MainProcess;

import java.io.FileNotFoundException;

import Main.learningOfGO;
import protein.GoSet;
import protein.proteinSet;

public class getCompleteAnn {

	public static void main(String[] args) throws FileNotFoundException 
	{
		// TODO Auto-generated method stub
		learningOfGO.aGoSet = new GoSet("../InFile/GeneOntology/gene_ontology_edit.obo.2013-06-15");
		
		proteinSet measure = new proteinSet();
		measure.AddAnnotation("../InFile/Measure/CAFA2/Ann");
		System.out.println(measure.size());
		System.out.println(measure.sizeOfAnn());
		
		measure.calMFBPCCsize(learningOfGO.aGoSet);
		
		measure.removeAllAnn();
		
		measure.AddAnnotationInProteinSet("../InFile/Swiss/Ann201509");
		measure.AddAnnotationInProteinSet("../InFile/Goa/Ann201509");
		measure.AddAnnotationInProteinSet("../InFile/GODB/Ann201509");
		measure.removeProteinHaveNoGOAnn();
		
		System.out.println(measure.size());
		System.out.println(measure.sizeOfAnn());
		measure.OutputType2Annotation("CAFA2AnnType2",learningOfGO.aGoSet);
		measure.OutputType3Annotation("CAFA2AnnType3",learningOfGO.aGoSet);
	}

}

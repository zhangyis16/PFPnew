package MainProcess;

import java.io.FileNotFoundException;

import Main.learning;
import protein.GoSet;
import protein.proteinSet;

public class getCompleteAnn {

	public static void main(String[] args) throws FileNotFoundException 
	{
		// TODO Auto-generated method stub
		learning.aGoSet = new GoSet("../InFile/GeneOntology/gene_ontology_edit.obo.2013-06-15");
		
		proteinSet measure = new proteinSet();
		measure.AddAnnotation("../InFile/Measure/CAFA2/Ann");
		System.out.println(measure.size());
		System.out.println(measure.sizeOfAnn());
		
		measure.calMFBPCCsize(learning.aGoSet);
		
		measure.removeAllAnn();
		
		measure.AddAnnotationInProteinSet("../InFile/Swiss/Ann201509");
		measure.AddAnnotationInProteinSet("../InFile/Goa/Ann201509");
		measure.AddAnnotationInProteinSet("../InFile/GODB/Ann201509");
		measure.removeProteinHaveNoAnn();
		
		System.out.println(measure.size());
		System.out.println(measure.sizeOfAnn());
		measure.OutputType2Annotation("CAFA2AnnType2",learning.aGoSet);
		measure.OutputType3Annotation("CAFA2AnnType3",learning.aGoSet);
	}

}

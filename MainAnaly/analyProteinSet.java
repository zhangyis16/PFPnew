package MainAnaly;

import java.io.FileNotFoundException;

import Main.learning;
import protein.proteinCommon;
import protein.proteinSet;

public class analyProteinSet {

	public static void main(String[] args) throws FileNotFoundException {
		// TODO Auto-generated method stub
		learning.aGoSet.Load("../InFile/GeneOntology/gene_ontology_edit.obo.2013-06-15");
		proteinSet measure = new proteinSet();
		proteinSet.LoadAccess2NameMap("../InFile/Swiss/ac2Name201401");
		measure.AddAnnotation("../InFile/Train/201401/Ann");
		measure.removeGoNotIn(learning.aGoSet);
		measure.addFather(learning.aGoSet);
		//measure.getAnnotationsize();
		measure.calMFBPCCsize(learning.aGoSet);
		
	}
}

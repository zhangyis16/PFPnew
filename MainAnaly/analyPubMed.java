package MainAnaly;

import java.io.FileNotFoundException;

import Main.learning;
import protein.proteinSet;

public class analyPubMed {

	public static void main(String[] args) throws FileNotFoundException {
		// TODO Auto-generated method stub
		learning.aGoSet.Load("../InFile/GeneOntology/gene_ontology_edit.obo.2013-06-15");
		proteinSet measure = new proteinSet();
		proteinSet all201401 = new proteinSet();
		all201401.loadPubMed("../InFile/Swiss/PubMed201401");
		measure.AddAnnotation("../InFile/Train/201401/Ann");
		all201401.AddAnnotation("../InFile/Train/201401/Ann");
		all201401.getIntersection(measure);
		all201401.calMFBPCCsize(learning.aGoSet);
		
		all201401.analyPubmed("../InFile/Train/201401/Pubmed.csv");
		all201401.outputPubMedID("../InFile/Train/201401/PubMedList");
		
	}

}

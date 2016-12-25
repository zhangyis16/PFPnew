package MainAnaly;

import java.io.FileNotFoundException;

import Main.learningOfGO;
import protein.proteinCommon;
import protein.proteinSet;

public class analyProteinSet {

	public static void main(String[] args) throws FileNotFoundException {
		// TODO Auto-generated method stub
		learningOfGO.aGoSet.Load("../InFile/GeneOntology/gene_ontology_edit.obo.2013-06-15");
		proteinSet measure = new proteinSet();
		proteinSet.LoadAccess2NameMap("../InFile/Swiss/ac2Name201401");
		measure.AddAnnotation("../InFile/Train/201401/Ann");
		
		System.out.println(measure.size());
		System.out.println(measure.sizeOfAnn());
		
		
		//measure.removeGoNotIn(learning.aGoSet);
		//measure.addFather(learning.aGoSet);
		//measure.getAnnotationsize();
		
		measure.removeGoNotIn(learningOfGO.aGoSet);
		measure.addGOFather(learningOfGO.aGoSet);
		measure.removeAnnotation(8150,3674,5575);
		measure.statSpeciesGoNum(proteinCommon.SetMySpecies, "201401SpeciesGoNum");
		
		measure.calMFBPCCsize(learningOfGO.aGoSet);
		measure.statGoFrequency("201401");
		
		
		
		/*proteinSet MFSub =	measure.getSubAnnotation('F');
		MFSub.OutputIDAnnSpecies("MFOIDAnnSpecies");
		proteinSet BPSub =	measure.getSubAnnotation('P');
		BPSub.OutputIDAnnSpecies("BPOIDAnnSpecies");
		proteinSet CCSub =	measure.getSubAnnotation('C');
		CCSub.OutputIDAnnSpecies("CCOIDAnnSpecies");    */
		
	}
}

package MainProcess;

import java.io.FileNotFoundException;

import Main.learningOfGO;
import protein.GoSet;
import protein.proteinSet;

public class evalutionCAFA2IEA {

	public static void main(String[] args) throws FileNotFoundException 
	{
		// TODO Auto-generated method stub
		learningOfGO.aGoSet = new GoSet("../InFile/gene ontology/gene_ontology_edit.obo.2013-07-01");
		proteinSet.LoadAccess2NameMap("../InFile/Swiss/ac201307");
		proteinSet measure = new proteinSet();
		
		measure.AddAnnotation("../InFile/CAFA2/type1annotation");
		measure.AddIEA_Annotation("../InFile/Swiss/SwissIEA_Ann201401");
		measure.AddIEA_Annotation("../InFile/Goa/GoaIEA_Ann201401");
		
		measure.addGOFather(learningOfGO.aGoSet);
		measure.addIEAFather(learningOfGO.aGoSet);
		measure.removeAnnotation(8150,3674,5575);
		measure.removeIEA_Annotation(8150,3674,5575);
		measure.evalutionIEAPR('F');
		measure.evalutionIEAPR('P');
		measure.evalutionIEAPR('C');
	}

}

package Main;

import java.io.FileNotFoundException;

import protein.GoSet;
import protein.proteinSet;

public class divide3SubSet {

	public static void main(String[] args) throws FileNotFoundException {
		// TODO Auto-generated method stub
		String GoVersion = "2013-06-15";
		learning.aGoSet = new GoSet("../InFile/gene ontology/gene_ontology_edit.obo." + GoVersion);
		proteinSet train =   new proteinSet();
		train.AddAnnotation("../InFile/Goa/Ann201409");
		train.removeGoNotIn(learning.aGoSet);
		proteinSet MFO =   train.getSubAnnotation('F');
		proteinSet BPO =   train.getSubAnnotation('P');
		proteinSet CCO =   train.getSubAnnotation('C');
		MFO.OutputAnnotation("goa.uniprot-goa.MFO");
		BPO.OutputAnnotation("goa.uniprot-goa.BPO");
		CCO.OutputAnnotation("goa.uniprot-goa.CCO");
	}

}

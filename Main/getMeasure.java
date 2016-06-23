package Main;

import java.io.FileNotFoundException;

import protein.GoSet;
import protein.proteinSet;

public class getMeasure {
	
	public static void main(String[] args) throws FileNotFoundException
	{
		learning.aGoSet = new GoSet("../InFile/gene ontology/gene_ontology_edit.obo.2013-07-01");
		proteinSet.LoadAccess2NameMap("../InFile/Swiss/ac2Name201307");
		proteinSet measure = new proteinSet();
		measure.AddAnnotation("../InFile/Measure/type1Ann");
		measure.loadFastaSequence("../InFile/Swiss/Sequence201401.fasta");
		measure.OutputAnnotation("../InFile/Measure/CAFA2Mea");
		measure.OutputFastaSequence("../InFile/Measure/201401seq");
	}
}

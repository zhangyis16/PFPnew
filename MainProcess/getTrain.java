package MainProcess;

import java.io.FileNotFoundException;

import Main.learning;
import protein.GoSet;
import protein.proteinCommon;
import protein.proteinSet;

public class getTrain {

	public static void main(String[] args) throws FileNotFoundException {
		learning.aGoSet = new GoSet("../InFile/GeneOntology/gene_ontology_edit.obo.2013-06-15");

		String t = "201301";
		
		
		proteinSet.LoadSwissMapAccess2UniAccess("../InFile/Swiss/ac2ac" + t);
		proteinSet.LoadAccess2NameMap("../InFile/Swiss/ac2Name" + t);
		proteinSet Train = new proteinSet();
		
		Train.AddAnnotation("../InFile/Swiss/Ann" + t);
		Train.AddAnnotationInSwiss("../InFile/Goa/Ann" + t);
		Train.AddAnnotationInSwiss("../InFile/GODB/Ann" + t);
		
		Train.eraserProteinOnly5515();   
		
		//Train.removeGoNotIn(learning.aGoSet);
		Train.loadFastaSequence("../InFile/Swiss/Swiss" + t + ".fasta");
		
		
		
		Train.OutputAnnotation("../InFile/temporary/Ann");
		Train.OutputFastaAccessSequence("../InFile/temporary/Seq");
		System.out.println(Train.size());

	}

}

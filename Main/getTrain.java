package Main;

import java.io.FileNotFoundException;

import protein.GoSet;
import protein.proteinCommon;
import protein.proteinSet;

public class getTrain {

	public static void main(String[] args) throws FileNotFoundException {
		//GoSet Go201307 = new GoSet("../InFile/gene ontology/gene_ontology_edit.obo.2013-07-01");

		learning.aGoSet = new GoSet("../InFile/gene ontology/gene_ontology_edit.obo.2014-01-01");
		
		
		
		proteinSet.LoadSwissMapAccess2UniAccess("../InFile/Swiss/ac2ac201401");
		proteinSet.LoadAccess2NameMap("../InFile/Swiss/ac2Name201401");
		proteinSet Swiss201401 = new proteinSet();
		Swiss201401.AddAnnotation("../InFile/Swiss/SwissAnnotation201401");
		Swiss201401.AddAnnotationInSwiss("../InFile/Goa/Annotation201401");
		Swiss201401.eraserProteinOnly5515();   
		
		//Swiss201401.statSpecies();
		//Swiss201401.loadFastaSequence("../InFile/Swiss/Sequence201401.fasta");
		//Swiss201401.statAminoAcid();
		
		//Swiss201401.statSpeciesStringLength(proteinCommon.SetMySpecies);
		
		//Swiss201401.statSpecies();
		
		learning.aGoSet.AddSon();
		learning.aGoSet.addDepth();
		
		
/*		learning.aGoSet.stat();
		
		learning.aGoSet.outputDepth("Depth.csv");
		learning.aGoSet.statMinDepth();
		learning.aGoSet.statMaxDepth();   */
		//Swiss201401.statLabelNum();

		
		
		
		Swiss201401.removeGoNotIn(learning.aGoSet);
		learning.aGoSet.AddAllFather();
		Swiss201401.removeRedunAnn(learning.aGoSet);
		System.out.println(Swiss201401.AnnotationSize());
		Swiss201401.addFather(learning.aGoSet);
		
		//Swiss201401.statSpeciesGoNum(proteinCommon.SetMySpecies,"spGOdisWithoutFather.csv");
		
		Swiss201401.statGoFrequency();
		
		//Swiss201401.statLabelNum();
		
		//Swiss201401.statLabelDepth();
	}

}
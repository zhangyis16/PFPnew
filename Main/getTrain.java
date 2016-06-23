package Main;

import java.io.FileNotFoundException;

import protein.GoSet;
import protein.proteinCommon;
import protein.proteinSet;

public class getTrain {

	public static void main(String[] args) throws FileNotFoundException {
		learning.aGoSet = new GoSet("../InFile/gene ontology/gene_ontology_edit.obo.2013-07-01");

		
		
		
		proteinSet.LoadSwissMapAccess2UniAccess("../InFile/Swiss/ac2ac201401");
		proteinSet.LoadAccess2NameMap("../InFile/Swiss/ac2Name201401");
		proteinSet Train201401 = new proteinSet();
		
		Train201401.AddAnnotation("../InFile/Swiss/Ann201401");
		Train201401.AddAnnotationInSwiss("../InFile/Goa/Ann201401");
		Train201401.AddAnnotationInSwiss("../InFile/GODB/Ann201401");
		
		Train201401.eraserProteinOnly5515();   
		
		//Swiss201401.statSpecies();
		//Swiss201401.loadFastaSequence("../InFile/Swiss/Sequence201401.fasta");
		//Swiss201401.statAminoAcid();
		
		//Swiss201401.statSpeciesStringLength(proteinCommon.SetMySpecies);
		
		//Swiss201401.statSpecies();
		
		//learning.aGoSet.AddSon();
		//learning.aGoSet.addDepth();
		
		
/*		learning.aGoSet.stat();
		
		learning.aGoSet.outputDepth("Depth.csv");
		learning.aGoSet.statMinDepth();
		learning.aGoSet.statMaxDepth();   */
		//Swiss201401.statLabelNum();

		
		
		
		Train201401.removeGoNotIn(learning.aGoSet);
		Train201401.loadFastaSequence("../InFile/Swiss/Swiss201401.fasta");
		
		Train201401.OutputAnnotation("../InFile/train/201401");
		Train201401.OutputFastaSequence("../InFile/train/201401seq");
		System.out.println(Train201401.size());
		//learning.aGoSet.AddAllFather();
		//Swiss201401.removeRedunAnn(learning.aGoSet);
		//System.out.println(Swiss201401.AnnotationSize());
		//Train201401.addFather(learning.aGoSet);
		//Train201401.removeAnnotation(8150,3674,5575);
		//System.out.println(Train201401.AnnotationSize());
		//Swiss201401.statSpeciesGoNum(proteinCommon.SetMySpecies,"spGOdisWithoutFather.csv");
		
		//Train201401.statGoFrequency();
		
		//Swiss201401.statLabelNum();
		
		//Swiss201401.statLabelDepth();
	}

}

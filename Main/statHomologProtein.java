package Main;

import java.io.IOException;
import java.util.Arrays;
import java.util.TreeSet;

import protein.GoSet;
import protein.proteinSet;

public class statHomologProtein {

	public static void main(String[] args) throws IOException, InterruptedException 
	{
		String species1 = "_HUMAN";
		String species2 = "_MOUSE";
		// TODO Auto-generated method stub
		String TrainDirectory = "../InFile/Train/201401/";
		proteinSet.LoadAccess2NameMap("../InFile/Swiss/ac2Name201401");
		learningOfGO.aGoSet = new GoSet("../InFile/GeneOntology/gene_ontology_edit.obo.2013-06-15");
		proteinSet train =   new proteinSet();
		proteinSet trainHuman =   new proteinSet();
		proteinSet trainMouse =   new proteinSet();
		train.AddAnnotation(TrainDirectory + "Ann");
		train.removeGoNotIn(learningOfGO.aGoSet);
		train.sortProteinBaseName();
		
		TreeSet<String> spset = new TreeSet<String>(Arrays.asList
				(species1,species2));
		train.filterSpecies(spset);
		train.filterHomologProtein();
		train.OutputProteinName("201401HomoProtein");
		
		train.loadFastaSequence("../InFile/Train/201401/Seq");
		trainHuman = train.getSpeciesSubSet(species1);
		trainMouse = train.getSpeciesSubSet(species2);
		
		trainHuman.addGOFather(learningOfGO.aGoSet);
		trainMouse.addGOFather(learningOfGO.aGoSet);
		
		proteinSet.compareHomoProtein(trainHuman, trainMouse ,learningOfGO.aGoSet);
		//trainHuman.OutputFastaSequence("HUMAN");
		//trainMouse.OutputFastaSequence("MOUSE");
		//trainHuman.OutputAnnotationName("HumanAnn");
		//trainMouse.OutputAnnotationName("MouseAnn");
	}
}

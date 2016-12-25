package MainProcess;

import java.io.FileNotFoundException;

import Main.learningOfGO;
import protein.GoSet;
import protein.proteinSet;

public class extractSpecies {

	public static void main(String[] args) throws FileNotFoundException {
		// TODO Auto-generated method stub
		learningOfGO.aGoSet = new GoSet("../InFile/gene ontology/gene_ontology_edit.obo." + "2013-06-15");
		
		String AccessVersion = "201401";
		String MeaDirectory = "../InFile/Measure/CAFA2/";
		proteinSet.LoadAccess2NameMap("../InFile/Swiss/ac2Name" + AccessVersion);
		
		proteinSet measure = new proteinSet();
		measure.AddAnnotation(MeaDirectory + "Ann");
		System.out.println("Measure size = " + measure.size());
		measure.loadFastaSequence(MeaDirectory + "Seq");
		measure.UpdateIndex();
		measure.LoadPredScore("../OutFile/KnnResult");
		
		proteinSet newSet = measure.getSpeciesSubSet("_ECOLI");
		
		System.out.println(newSet.size());
		newSet.OutputAnnotation("ECOLI_Ann");
		newSet.OutputFastaNameSequence("ECOLI_seq");
		newSet.OutputGOPredScore("Ecoli_pred");
	}

}

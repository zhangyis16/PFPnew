package MainProcess;

import java.io.FileNotFoundException;

import Main.learningOfGO;
import protein.GoSet;
import protein.proteinSet;

public class Simply {

	public static void main(String[] args) throws FileNotFoundException {
		// TODO Auto-generated method stub
		String TrainDirectory = "../InFile/Train/201401/";
		String GoVersion = "2013-06-15";
		String MeaDirectory   = "../InFile/Measure/201409-201401/";
		
		
		learningOfGO.aGoSet = new GoSet("../InFile/gene ontology/gene_ontology_edit.obo." + GoVersion);
		proteinSet train =   new proteinSet();
		
		train.AddAnnotation("Ann");
		System.out.println(train.size());
		train.removeGoNotIn(learningOfGO.aGoSet);
		System.out.println(train.size());
		train.OutputAnnotation("AnnWithoutFather");
		train.addGOFather(learningOfGO.aGoSet);
		
		
		
	/*	proteinSet measure = new proteinSet();
		measure.AddAnnotation(MeaDirectory+"Ann");
		
		measure.removeGoNotIn(learning.aGoSet);
		measure.OutputAnnotation("AnnWithoutFather");
		measure.addFather(learning.aGoSet);
		measure.OutputAnnotation("AnnWithFather");
		measure.removeAnnotation(8150,3674,5575);   */
	}

}

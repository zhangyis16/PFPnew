package MainProcess;

import java.io.FileNotFoundException;

import Main.learningOfGO;
import protein.GoSet;
import protein.proteinSet;

public class DivideTrainSet {

	public static void main(String[] args) throws FileNotFoundException {
		// TODO Auto-generated method stub
		String TrainDirectory = "../InFile/Train/201401new/";
		String GoVersion = "2013-06-15";
		learningOfGO.aGoSet = new GoSet("../InFile/gene ontology/gene_ontology_edit.obo." + GoVersion);
		
		proteinSet train =   new proteinSet();
		proteinSet trainNew =   new proteinSet();
		proteinSet leaveOneTrain =   new proteinSet();
		train.AddAnnotation(TrainDirectory + "Ann");
		train.loadFastaSequence(TrainDirectory + "Seq");
		
		System.out.println(train.size());
		leaveOneTrain = train.getRandomSubProteinSet(2000);
		
		leaveOneTrain.OutputAnnotation("Ann2000");
		leaveOneTrain.OutputFastaAccessSequence("Seq2000");
		System.out.println(leaveOneTrain.size());
		
		trainNew = proteinSet.getNewProtein(train, leaveOneTrain);
		trainNew.OutputAnnotation("AnnNew");
		trainNew.OutputFastaAccessSequence("SeqNew");
		
		System.out.println(trainNew.size());
	}

}

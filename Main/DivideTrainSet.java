package Main;

import java.io.FileNotFoundException;

import protein.GoSet;
import protein.proteinSet;

public class DivideTrainSet {

	public static void main(String[] args) throws FileNotFoundException {
		// TODO Auto-generated method stub
		String TrainDirectory = "../InFile/Train/201401new/";
		String GoVersion = "2013-06-15";
		learning.aGoSet = new GoSet("../InFile/gene ontology/gene_ontology_edit.obo." + GoVersion);
		
		proteinSet train =   new proteinSet();
		proteinSet trainNew =   new proteinSet();
		proteinSet leaveOneTrain =   new proteinSet();
		train.AddAnnotation(TrainDirectory + "Ann");
		train.loadFastaSequence(TrainDirectory + "Seq");
		
		System.out.println(train.size());
		leaveOneTrain = train.getRandomSubProteinSet(2000);
		
		leaveOneTrain.OutputAnnotation("Ann2000");
		leaveOneTrain.OutputFastaSequence("Seq2000");
		System.out.println(leaveOneTrain.size());
		
		trainNew = proteinSet.getNewProtein(train, leaveOneTrain);
		trainNew.OutputAnnotation("AnnNew");
		trainNew.OutputFastaSequence("SeqNew");
		
		System.out.println(trainNew.size());
	}

}

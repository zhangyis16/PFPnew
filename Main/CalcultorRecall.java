package Main;

import java.io.FileNotFoundException;
import java.io.IOException;

import protein.GoSet;
import protein.proteinSet;

public class CalcultorRecall {

	public static void main(String[] args) throws IOException {
		// TODO Auto-generated method stub
		String GoVersion = "2013-07-01";
		String AccessVersion = "201307";
		String TrainDirectory = "../InFile/Train/201401/";
		String MeaDirectory   = "../InFile/Measure/CAFA2/";
		String blastResult = "../InFile/blastResult/CAFA2|201401.outfmt6";
		
		learningOfGO.aGoSet = new GoSet("../InFile/gene ontology/gene_ontology_edit.obo." + GoVersion);
		proteinSet.LoadAccess2NameMap("../InFile/Swiss/ac2Name" + AccessVersion);
		
		proteinSet measure = new proteinSet();
		measure.AddAnnotation(MeaDirectory + "Ann");

		proteinSet train =   new proteinSet();
		train.AddAnnotation(TrainDirectory + "Ann");
		
		
		train.removeGoNotIn(learningOfGO.aGoSet);
		measure.removeGoNotIn(learningOfGO.aGoSet);
		
		train.addGOFather(learningOfGO.aGoSet);
		train.removeAnnotation(8150,3674,5575);
		
		measure.addGOFather(learningOfGO.aGoSet);
		measure.removeAnnotation(8150,3674,5575);
		
		measure.LoadNoSparseFeatureVector(MeaDirectory + "Feature" ,567);
		measure.setLiblinearFeatureFromNoSparseFeature();        
		train.setPredList(train);
		measure.setPredList(train);
		train.sortPredListBaseLabelOrder();
		measure.sortPredListBaseLabelOrder();

		
			
/*		measure.libLinearPredict("../OutFile/TrainModel/");
		measure.addTopK_L2RCandidate(50);
		measure.clearPredResult();
		System.out.println("liblinear finish");    */
		
		measure.ScoreNaiveBaseline(train);	
		measure.addTopK_L2RCandidate(50,'F');
		measure.clearPredResult();
		System.out.println("naive finish");
		
		measure.addBlastResultBitScore(blastResult);
		measure.BlastKnnBaseline(train);
		measure.addTopK_L2RCandidate(50,'F');
		
		System.out.println("KNN finish");
		
		double recall = measure.CalCandidateRecall('F');
		System.out.println(recall);
	}

}

package Main;

import java.io.FileOutputStream;
import java.io.IOException;
import java.io.PrintWriter;

import liblinear.InvalidInputDataException;
import protein.GoSet;
import protein.proteinSet;

public class learning {
	
	public static GoSet aGoSet;
	public static void main(String[] args) throws IOException, InvalidInputDataException 
	{

		GoSet Go201307 = new GoSet("../InFile/gene ontology/gene_ontology_edit.obo.2013-07-01");
		proteinSet.LoadAccess2NameMap("../InFile/Swiss/ac201307");
		
		proteinSet train =   new proteinSet();
		proteinSet measure = new proteinSet();
		
		train.AddAnnotation("../InFile/Train/201401");
		measure.AddAnnotation("../InFile/Measure/Intersect/type1");
		
	    measure.addFather(Go201307);
		train.addFather(Go201307);
		train.removeAnnotation(8150,3674,5575);
		measure.removeAnnotation(8150,3674,5575);
		
		measure.loadFastaSequence("../InFile/Measure/Intersect/type1measure.fasta");
		train.loadFastaSequence("../InFile/Train/201401train.fasta");
			
		measure.setPredList(train);    
		
		//train.tranSequence2Vector();
		//measure.tranSequence2Vector();
		
		//train.setPredList(train);
		//train.setLiblinearFeature();
		//train.libLinearTrain("../OutFile/TrainModel/");   
		
		
		//measure.setLiblinearFeature();
		
		//measure.libLinearPredict("../OutFile/TrainModel/");
		//measure.OutputPred("../OutFile/predResult");
		//measure.addPredResult("../OutFile/predResult");
		
		//measure.addBlastResultBitScore("../InFile/blastResult/psiBlastResult.out");
		measure.addBlastResultSimility("../InFile/blastResult/psiBlastResult.out");
		//measure.readMSAlignResult("../InFile/blastResult/psiBlastResult.outfmt4");
		//measure.outputMSA("error.out");
		measure.GOtchaBaseline(train);
		
		measure.evalution();
	}
}

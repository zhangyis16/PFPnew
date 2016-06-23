package Main;

import java.io.FileInputStream;

import java.io.IOException;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Properties;

import common.FeatureTran;
import liblinear.InvalidInputDataException;
import protein.GoSet;
import protein.proteinCommon;
import protein.proteinSet;

public class learning {
	
	public static GoSet aGoSet = new GoSet();
	public static String GoDirectory;
	public static String Access2NameDirectory;
	public static String TrainDirectory;
	public static String MeaDirectory;
	public static String blastResult;
	public static String Task;
	public static int ThreadNum;
	public static double liblinearC;
	public static String modelDir;

	public static String predResultFile;
	public static String resultOutPattern;
	public static String resultOutFile;
	public static liblinear.Parameter liblinearParam   = new liblinear.Parameter(liblinear.SolverType.L2R_LR, 1, 1000, 0.05);
	public static int FeatureSize = 8000+400+4;
	
	
	public static void ReadConfig() throws IOException
	{
		FileInputStream in = new FileInputStream("../bin/config.in");
		Properties config = new Properties();
		config.load(in);
		GoDirectory = config.getProperty("GoDirectory");
		Access2NameDirectory = config.getProperty("Access2NameDirectory");
		TrainDirectory = config.getProperty("TrainDirectory");
		MeaDirectory   = config.getProperty("MeaDirectory");
		blastResult = config.getProperty("blastResult");
		Task = config.getProperty("Task");
		ThreadNum = Integer.parseInt(config.getProperty("ThreadNum"));
		liblinearC = Double.parseDouble(config.getProperty("liblinearC"));
		modelDir = config.getProperty("modelDir");
		predResultFile = config.getProperty("predResultFile");
		resultOutPattern = config.getProperty("resultOutPattern");
		resultOutFile = config.getProperty("resultOutFile");
	}
	public static void main(String[] args) throws IOException, InvalidInputDataException 
	{
		learning.ReadConfig();
		
		
		learning.aGoSet.Load(GoDirectory);
		System.out.println("TrainSet = " + TrainDirectory);
		System.out.println("MeasureSet = " + MeaDirectory);
		proteinSet.LoadAccess2NameMap(Access2NameDirectory);
		
		proteinSet measure = new proteinSet();
		measure.AddAnnotation(MeaDirectory + "Ann");
		System.out.println("Measure size = " + measure.size());
		
		proteinSet train =   new proteinSet();
		train.AddAnnotation(TrainDirectory + "Ann");
		train.removeGoNotIn(learning.aGoSet);
		
		System.out.println("Train size = " + train.size());
		measure.removeGoNotIn(learning.aGoSet);
		
		train.addFather(learning.aGoSet);
		train.removeAnnotation(8150,3674,5575);
		
		measure.addFather(learning.aGoSet);
		measure.removeAnnotation(8150,3674,5575);
		
		if (Task.equals("NaiveSpecies"))
		{
			ArrayList<String> SpeciesList = new ArrayList<String>(Arrays.asList
					("_ARATH","_BACSU","_DANRE","_DICDI","_DROME","_ECOLI",
					 "_HUMAN","_MOUSE","_PSEAE","_RAT"  ,"_SCHPO","_XENLA","_YEAST"));
			measure.naiveSpeciesBaseline(train,SpeciesList); 
		}
		if (Task.equals("Naive"))
		{
			measure.naiveBaseline(train);
		}
		if (Task.equals("LiblinearTrain"))
		{
			liblinearParam   = new liblinear.Parameter(liblinear.SolverType.L2R_LR, liblinearC, 1000, 0.05);
			train.loadFastaSequence(TrainDirectory + "Seq");
			System.out.println("Train size = " + train.size());

			
			train.OutputFastaSequence("Seq");
			
			train.tranSequence2TriSparseFeature();
			train.outputLibTrainFile("TrainSparseFeature", 1);
			train.setLiblinearFeatureFromSparseFeature();
			
			train.setPredList(train);
			train.sortPredListBaseFrequency();
			train.libLinearTrain(modelDir,ThreadNum); 
		}
		if (Task.equals("LiblinearPred"))
		{

			measure.loadFastaSequence(MeaDirectory + "Seq");
			measure.tranSequence2TriSparseFeature();
			measure.outputLibTrainFile("MeasureSparseFeature", 1);
			measure.setLiblinearFeatureFromSparseFeature();   
			
			
			measure.setPredList(train);
			measure.sortPredListBaseFrequency();
			measure.libLinearPredict(modelDir);
			
			measure.removeLowPred(2000);
			measure.OutputPredScore(predResultFile);
		}
		
		if (Task.equals("LiblinearPredCheck"))
		{
			proteinCommon.CalTriAcidSimiliar(0.7);
			measure.loadFastaSequence(MeaDirectory + "Seq");
			measure.tranSequence2TriSparseFeature();
			
			measure.setLiblinearFeatureFromSparseFeature();   
			
			measure.libLinearPredict(modelDir,3006);
			measure.calOneLabelAUC(3006);
		}
		
		if (Task.equals("L2R"))
		{
			measure.LoadPredScore("../OutFile/KnnResult");
			measure.evalEveryLabelAUC("everyLabelResult");
			measure.clearPredResult();
		}
		if (Task.equals("Blast"))
		{
			measure.addBlastResultBitScore(blastResult);
			measure.blastBaseline(train);
		}
		
		if (Task.equals("BlastKnn"))
		{
			System.out.println("Begin Blast Weight Knn Score");
			measure.addBlastResultBitScore(blastResult);
			//measure.removeBlastOtherSp();
			measure.GOtchaBaseline(train);
			measure.OutputPredScore(predResultFile);
		}
		
		
		
		//measure.addBlastResultSimility("../InFile/blastResult/psiBlastResult.out");
		//measure.readMSAlignResult("../InFile/blastResult/psiBlastResult.outfmt4");
		//measure.outputMSA("error.out");
		
		
		System.out.println(Task);
		measure.evalutionEveryProtein("everyProteinResult","MFO");
		measure.evalutionFmaxAndAUPR(resultOutFile,resultOutPattern,"all");
		//measure.evalutionSpecies(resultOutFile,"_HUMAN","_MOUSE","_ARATH","_ECOLI","_RAT","_PSEAE","_YEAST");
		
		/*
		measure.setPredList(train);
		measure.sortPredListBaseFrequency();
		measure.evalEveryLabelAUC("everyLabelResult");
		measure.CalTopkRecall(10,50,100);
		*/   
	}
}

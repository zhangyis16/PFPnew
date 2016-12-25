package hpo;

import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.PrintWriter;
import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.Properties;

import common.Output;
import common.Pair;
import multiLabel.MLlearning;
import protein.proteinSet;

public class learningHPO {
	
	public static HPOSet aHPOSet = new HPOSet();
	
	public static String Task = "Naive";
	public static liblinear.Parameter liblinearParam   = new liblinear.Parameter(liblinear.SolverType.L2R_LR, 1, 1000, 0.05);
	public static double liblinearC = 1;
/*	public static String HPODirectory = "../InFile/HPO/Ontology/HPO2013-09-17.obo";
	public static String MeaDirectory = "../InFile/HPO/Measure/CAFA2/";
	public static String TrainDirectory = "../InFile/HPO/Train/CAFA2/";
	public static String Access2NameDirectory = "../InFile/Swiss/ac2Name201307";
	public static String blastResult = "../InFile/blastResult/HPOblast/CAFA2.outfmt6";   */
	
/*	public static String HPODirectory = "../InFile/HPO/Ontology/HPO2016-01-13.obo";
	public static String Access2NameDirectory = "../InFile/Swiss/ac2Name201605";
	public static String MeaDirectory = "../InFile/HPO/Measure/HPO20160125CV1/";
	public static String TrainDirectory = "../InFile/HPO/Train/HPO20160125CV1/";
	public static String blastResult = "../InFile/blastResult/HPOblast/CV1_20160125.outfmt6";   */
	
	public static String HPODirectory;
	public static String Access2NameDirectory;
	public static String MeaDirectory;
	public static String TrainDirectory;
	public static String blastResult;   
	
	public static String modelDir = new String();  
	public static DecimalFormat doubleFormat = new DecimalFormat( "0.00000"); 
	public static Properties config;
	public static int ThreadNum;
	public static String predResultFile;
	
	public static void ReadConfig() throws IOException
	{
		FileInputStream in = new FileInputStream("../bin/HPOconfig.in");
		config = new Properties();
		config.load(in);
		doubleFormat = new DecimalFormat( "0.00000"); 
		HPODirectory = config.getProperty("HPODirectory");
		Access2NameDirectory = config.getProperty("Access2NameDirectory");
		TrainDirectory = config.getProperty("TrainDirectory");
		MeaDirectory   = config.getProperty("MeaDirectory");
		blastResult = config.getProperty("blastResult");
		Task = config.getProperty("Task");
		ThreadNum = Integer.parseInt(config.getProperty("ThreadNum"));
		liblinearC = Double.parseDouble(config.getProperty("liblinearC"));
		modelDir = config.getProperty("modelDir");
		predResultFile = config.getProperty("predResultFile");
	}
	
	public static void main(String[] args) throws IOException 
	{
		// TODO Auto-generated method stub
		learningHPO.ReadConfig();
		aHPOSet.Load(HPODirectory);
		proteinSet.LoadAccess2NameMap(Access2NameDirectory);
		
		System.out.println("HPO size = " + aHPOSet.sizeOfHPO());
		proteinSet measure = new proteinSet();
		measure.AddHPOAnnotation(MeaDirectory + "Ann");
		
		//proteinSet CAFA2Target = new proteinSet();
		//CAFA2Target.addProteinFromFile("../InFile/CAFA2target100816");
		//measure.getIntersection(CAFA2Target);
		//measure.loadFastaSequence("../InFile/Swiss/Swiss201401.fasta");
		//measure.loadFastaSequence("../InFile/Trembl/SeqALL.fasta");
		//measure.OutputFastaAccessSequence(MeaDirectory + "Seq");
		
		
		measure.removeHPONotIn(aHPOSet);				//移除不在这个HPO本体中的标注
		System.out.println("Measure Ann size = " + measure.sizeOfHPOAnn());
		measure.addHPOFather(aHPOSet);					//根据功能的传递关系加入父亲节点
		measure.removeHPOAnnotation(1);					//去掉根节点
		
		proteinSet train =   new proteinSet();
		train.AddHPOAnnotation(TrainDirectory + "Ann");
		//train.loadFastaSequence("../InFile/Swiss/Swiss201401.fasta");
		//train.loadFastaSequence("../InFile/Trembl/SeqALL.fasta");
		//train.OutputFastaAccessSequence(TrainDirectory + "Seq");
		
		train.removeHPONotIn(aHPOSet);
		System.out.println("Train Ann size = " + train.sizeOfHPOAnn());
		train.addHPOFather(aHPOSet);
		train.removeHPOAnnotation(1);
		
		System.out.println("Measure Ann size = " + measure.sizeOfHPOAnn());
		System.out.println("Measure size = " + measure.size());
		System.out.println("Train Ann size = " + train.sizeOfHPOAnn());
		System.out.println("Train size = " + train.size());
		
		if (Task.equals("Naive"))
		{
			measure.ScoreHPONaiveBaseline(train);
			ArrayList<Pair<Integer,Double>> predList = train.getHPONaiveList();
			Output.OutputPairList("HPONaiveResult", predList);
		}
		if (Task.equals("Blast"))
		{
			measure.addBlastResultBitScore(blastResult);
			measure.blastHPOBaseline(train);
		}
		System.out.println("*************");
		if (Task.equals("BlastKnn"))
		{
			measure.addBlastResultBitScore(blastResult);	//载入Blast的预测结果
			measure.BlastKnnHPOBaseline(train);	//进行带权重的KNN预测.
			
			ArrayList<ArrayList<Pair<Integer,Double>>> predScoreList = measure.getHPOPredScoreList();
			MLlearning.addNoScoreLabelNaiveResult(predScoreList, train.getHPONaiveList());
			measure.addHPOPredScoreList(predScoreList);
		}
		
		if (Task.equals("Liblinear3merTrain"))
		{
			liblinearParam   = new liblinear.Parameter(liblinear.SolverType.L2R_LR, liblinearC, 1000, 0.05);
			train.setHPOPredList(train);
			train.loadFastaSequence(TrainDirectory + "Seq");
			ArrayList<Integer> LabelList = train.getHPOLabelList();
			ArrayList<HashSet<Integer>> InstanceLabel = train.getInstanceLabel();
			ArrayList<ArrayList<Pair<Integer, Double>>> FeatureList = train.get3merFeatureList();
			int FeatureSize = 8000;
			MLlearning.MultiLabelLiblinearTrain(LabelList, InstanceLabel, FeatureList, modelDir, FeatureSize, ThreadNum); 
		}
		
		if (Task.equals("Liblinear2merTrain"))
		{
			liblinearParam   = new liblinear.Parameter(liblinear.SolverType.L2R_LR, liblinearC, 1000, 0.05);
			train.setHPOPredList(train);
			train.loadFastaSequence(TrainDirectory + "Seq");
			ArrayList<Integer> LabelList = train.getHPOLabelList();
			ArrayList<HashSet<Integer>> InstanceLabel = train.getInstanceLabel();
			ArrayList<ArrayList<Pair<Integer, Double>>> FeatureList = train.get2merFeatureList();
			int FeatureSize = 400;
			MLlearning.MultiLabelLiblinearTrain(LabelList, InstanceLabel, FeatureList, modelDir, FeatureSize, ThreadNum); 
		}
		
		if (Task.equals("Liblinear2merPred"))
		{
			measure.loadFastaSequence(MeaDirectory + "Seq");
			ArrayList<Integer> LabelList = train.getHPOLabelList();
			//Output.OutputPairList("HPONaiveList", HPONaiveList);
			ArrayList<ArrayList<Pair<Integer, Double>>> SparseFeatureList = measure.get2merFeatureList();
			int FeatureSize = 400;
			ArrayList<ArrayList<Pair<Integer,Double>>> predScoreList = MLlearning.MultiLabelLiblinearPrediction(LabelList, SparseFeatureList, modelDir, FeatureSize);
			measure.addHPOPredScoreList(predScoreList);
			measure.OutputHPOPredScore(predResultFile);
		}
		
		if (Task.equals("Liblinear3merPred"))
		{
			measure.loadFastaSequence(MeaDirectory + "Seq");
			ArrayList<Integer> LabelList = train.getHPOLabelList();
			//Output.OutputPairList("HPONaiveList", HPONaiveList);
			ArrayList<ArrayList<Pair<Integer, Double>>> SparseFeatureList = measure.get3merFeatureList();
			int FeatureSize = 8000;
			ArrayList<ArrayList<Pair<Integer,Double>>> predScoreList = MLlearning.MultiLabelLiblinearPrediction(LabelList, SparseFeatureList, modelDir, FeatureSize);
			measure.addHPOPredScoreList(predScoreList);
			measure.OutputHPOPredScore(predResultFile);
		}
		
		measure.evalutionHPOFmaxAndAUPR();
		measure.setHPOPredList(train);
		measure.evalEveryHPOLabelAUC("hpoAUC");
		
	}

}

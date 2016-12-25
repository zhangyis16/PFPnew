package multiLabel;

import java.io.File;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.List;

import Main.learningOfGO;
import common.Pair;
import liblinear.Feature;
import liblinear.FeatureNode;
import liblinear.Linear;
import protein.oneProtein;
import protein.proteinSet;
import protein.proteinSet.MultiThreadTrain;

public class MLlearning 
{
	public static liblinear.Parameter liblinearParam   = new liblinear.Parameter(liblinear.SolverType.L2R_LR, 1, 1000, 0.05);
	public static ArrayList<Integer> LabelList;
	public static ArrayList<HashSet<Integer>> InstanceLabel;
	public static ArrayList<ArrayList<Pair<Integer, Double>>> FeatureList;
	public static String modelDir = new String();
	public static int  FeatureSize;
	public static int ThreadNum;
	
	public static void MultiLabelLiblinearTrain(ArrayList<Integer> LabelList, ArrayList<HashSet<Integer>> InstanceLabel,
			ArrayList<ArrayList<Pair<Integer, Double>>> FeatureList, String modelDir,int FeatureSize,
			int ThreadNum) throws IOException
	{
		
		MLlearning.LabelList = LabelList;
		MLlearning.InstanceLabel = InstanceLabel;
		MLlearning.FeatureList = FeatureList;
		MLlearning.modelDir = modelDir;
		MLlearning.FeatureSize = FeatureSize;
		MLlearning.ThreadNum = ThreadNum;
		MLlearning.libLinearTrain();
	}
	
    public static void libLinearTrain() throws IOException
    {
    	MultiThreadTrain my = new MultiThreadTrain(0, modelDir);
    	ArrayList<Thread> t = new ArrayList<Thread>();
    	for (int i = 1;i <= ThreadNum ;i++)
    	{
    		t.add(new Thread(my,"thread " + i));
    	}
    	for (int i = 0;i < ThreadNum ;i++)
    	{
    		t.get(i).start();
    	}
    	boolean bo = true;
    	while (bo)
    	{
    		bo = false;
    		for (int i = 0;i < ThreadNum ;i++)
    			if (t.get(i).getState() != Thread.State.TERMINATED) bo = true;
    	}   	
    }
    
	public static class MultiThreadTrain implements Runnable 
	{
		private int count;
		private String Direction  = new String();
		public MultiThreadTrain(int countNum,String aDirection)
		{
			count = countNum;
			Direction = aDirection;
		}
		public void run() {
			// TODO Auto-generated method stub
			while(count<=MLlearning.LabelList.size() - 1)
			{
				try {
					int num = this.get();
					System.out.println(Thread.currentThread().getName()+ "now train label" + LabelList.get(num));
					libLinearTrain(LabelList.get(num),Direction);
					
				} catch (IOException e) 
				{
					// TODO Auto-generated catch block
					e.printStackTrace();
				}
			}
		}
		synchronized public int get(){
			return this.count++;
		}
		
	}
    
	public static void libLinearTrain(int labelNum, String Direction) throws IOException
    {
    	int Ann = labelNum;
        int max_index = MLlearning.FeatureSize;   //表示特征长度
        double bias = -1;
        
        liblinear.Problem   prob    = new liblinear.Problem();
        
        prob.bias = bias;
        prob.l = FeatureList.size(); //表示训练集的大小
        
        prob.n = max_index;            //表示最大特征维数
        if (bias >= 0) {
            prob.n++;
        }
        prob.x = new liblinear.Feature[prob.l][];
        prob.y = new double[prob.l];
        
        for (int i = 0; i < FeatureList.size(); i++)
        {
        		List<Feature> y = new ArrayList<Feature>();
        		for (Pair<Integer, Double> pair : FeatureList.get(i)) 
        		{
        			Feature node = new FeatureNode(pair.getFirst(), pair.getSecond());
        			y.add(node);
        		}
         		liblinear.Feature[] x = new Feature[y.size()];
        		x = y.toArray(x);
        		
        		prob.x[i] = x;
        }

        for (int i = 0; i < InstanceLabel.size(); i++)
        {
        	if (InstanceLabel.get(i).contains(Ann))	
        		prob.y[i] = 1; 
        	
        	else 
        		prob.y[i] = -1;
        	
        }
        if (Ann == 118)
        {
        	PrintWriter Fout = new PrintWriter(new FileOutputStream("TrainFeatureOutputFile"));
        	for (int i = 0;i<FeatureList.size();i++)
        	{
        		if (prob.y[i] == 1) Fout.print("+1 "); else Fout.print("-1 ");
        		for (int j=0;j< prob.x[i].length;j++)
        		{
        			Fout.print(prob.x[i][j].getIndex()+":"+ prob.x[i][j].getValue() + " ");
        		}
        		Fout.println();
        	}
        	Fout.close();
        }
        
        
    	liblinear.Model model = liblinear.Linear.train(prob, MLlearning.liblinearParam);
    	liblinear.Linear.saveModel(new File(Direction + "Label" + labelNum + ".model"), model);
    }

	public static ArrayList<ArrayList<Pair<Integer,Double>>> MultiLabelLiblinearPrediction(ArrayList<Integer> LabelList,
			ArrayList<ArrayList<Pair<Integer, Double>>> SparseFeatureList, String modelDir,int FeatureSize) throws IOException
	{
		ArrayList<ArrayList<Pair<Integer,Double>>> predScoreList = new ArrayList<ArrayList<Pair<Integer,Double>>>();
		for (int i = 0;i<SparseFeatureList.size();i++)
			predScoreList.add(new ArrayList<Pair<Integer,Double>>());
		
		
		ArrayList<Feature[]> liblinearFeatureList = new ArrayList<Feature[]>();
		for (int i = 0;i<SparseFeatureList.size();i++)
		{
			Feature[] liblinearFeature;
			List<Feature> x = new ArrayList<Feature>();
			for (Pair<Integer, Double> pair : SparseFeatureList.get(i)) 
			{
				Feature node = new FeatureNode(pair.getFirst(), pair.getSecond());
				x.add(node);
			}
			liblinearFeature = new Feature[x.size()];
			liblinearFeature = x.toArray(liblinearFeature);
			liblinearFeatureList.add(liblinearFeature);
		}
		
		
		
    	for (int count = 0;count<LabelList.size();count++)
    	{
    		int labelNum = LabelList.get(count);
    		try
    		{
    			liblinear.Model model = Linear.loadModel(new File(modelDir + "Label" + labelNum + ".model"));
    			if (count % 100 ==0)
    				System.out.println("Predict " + count + "th Label" + labelNum);
    			int index = 0;
    			for(Feature[] e:liblinearFeatureList)
    			{
    				double[] prob_estimates = new double[2];
    				double predict_label = Linear.predictProbability(model, e, prob_estimates);
    				predScoreList.get(index).add(new Pair<Integer, Double>(labelNum, prob_estimates[0]));
    				index++;
    			}
    		}
    		catch (IOException exp) 
    		{
    			System.out.println("Not find Lable " + labelNum + " Model File");
    		}
    	}
    	
		return predScoreList;
	}

	public static void addNoScoreLabelNaiveResult(ArrayList<ArrayList<Pair<Integer,Double>>> PredictionList
			,ArrayList<Pair<Integer,Double>> naiveResult)
	{
		double miniScore = 1.0;
		for (int i=0;i<PredictionList.size();i++)
		{
			for (int j=0;j<PredictionList.get(i).size();j++)
			{
				double score = PredictionList.get(i).get(j).getSecond();
				if (score<miniScore) miniScore = score;
			}
		}
		for (int i=0;i<PredictionList.size();i++)
		{
			HashSet<Integer> LabelSet = new HashSet<Integer>();
			for (int j=0;j<PredictionList.get(i).size();j++)
				LabelSet.add(PredictionList.get(i).get(j).getFirst());
			for (int j=0;j<naiveResult.size();j++)
			{
				int Label = naiveResult.get(j).getFirst();
				double score = naiveResult.get(j).getSecond();
				if (!LabelSet.contains(Label)) PredictionList.get(i).add(new Pair<Integer,Double>(Label,score*miniScore));
			}
		}
	}

}

package protein;
import common.*;
import java.util.*;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.PrintWriter;

public class evalution {
	public static Pair<Double,Double> GetFMeasureMax(ArrayList<HashSet<Integer>> answer,
			ArrayList<ArrayList<Pair<Integer,Double>>> predictor,
			ArrayList<Pair<Double,Double>> PreReCallPair) throws FileNotFoundException
	{
		
		ArrayList<labelResult> Result = new ArrayList<labelResult>();
		PreReCallPair.clear();
		
		for (int i = 0; i< answer.size(); i++)
		{
			while (answer.get(i).size() == 0)
			{
				answer.remove(i);
				predictor.remove(i);
			}
		}
		for (int i = 0; i< answer.size(); i++)
		{
            for (Pair<Integer,Double> node : predictor.get(i))
            {
            	int first = node.getFirst();
            	if (answer.get(i).contains(first))
            		Result.add(new labelResult(i,true,node.getSecond()));
            	else
            		Result.add(new labelResult(i,false,node.getSecond()));
            } 
		}
		Collections.sort(Result);
		
		int cutIndex = 0;
		double threshhold = 0;
		double fmax = 0;
		int [] predNum = new int[answer.size()];
		int [] TpNum = new int[answer.size()];
		double recall,presion;
	    int prenum,recallnum;
	    double presum,recallsum,fmeasure;
	    
		for (int i=0;i<predNum.length;i++) 
		{
			predNum[i] = 0;
			TpNum[i] = 0;
		}
		
		while (cutIndex < Result.size())
		{
			int beginIndex = cutIndex;
			int endIndex = cutIndex + 1;
			while (endIndex < Result.size() && 
				(Result.get(beginIndex).getValue() - Result.get(endIndex).getValue()) < Parameter.resolution)
				endIndex++;
			cutIndex = endIndex;
			for (int i = beginIndex;i < endIndex; i++)
			{
				labelResult node = Result.get(i); 
				predNum[node.getIndex()]++;
				if (Result.get(i).getResult()) TpNum[node.getIndex()]++;
			}
	        prenum = 0; recallnum = predNum.length;
	        presum = 0; recallsum = 0;
	        for (int i=0;i<predNum.length;i++) 
	        {
	        	recall = (double)TpNum[i] / answer.get(i).size();
	        	presion = (double)TpNum[i] / predNum[i];
	        	recallsum += recall;
	        	if (predNum[i] > 0) 
	        	{
	        		prenum++;
	        		presum += presion;
	        	}
	        }
	        recallsum /= recallnum;
	        if (prenum > 0) presum /= prenum; else presum = 0;
	        if (presum + recallsum == 0) fmeasure = 0; else fmeasure = 2 * presum * recallsum / (presum+recallsum);
	        //System.out.println("Fmax = " + fmeasure + "Cut" + Result.get(endIndex - 1).getValue());
	        
	        PreReCallPair.add(new Pair<Double,Double>(presum,recallsum));
	        if (fmeasure > fmax)
	        {
	            fmax = fmeasure;
	            threshhold = Result.get(endIndex - 1).getValue();
	        }
		}
		return new Pair<Double,Double>(fmax,threshhold);
	}
	
	
	
	public static Pair<Double,Double> GetFMeasureMaxSimpleVersion(ArrayList<HashSet<Integer>> answer,
									ArrayList<ArrayList<Pair<Integer,Double>>> predictor)
	{
	    double pre,recall,fmax = 0;
	    double threshhold = 0;
	    int prenum,recallnum;
	    double presum,recallsum,fmeasure;
	    for (int o=1;o<=99;o++)
	    {
	        double cut = (double)o/100;
	        prenum = 0; recallnum = 0;
	        presum = 0; recallsum = 0;
	        
	        HashSet<Integer> golist = new HashSet<Integer>();
	        HashSet<Integer> predlist = new HashSet<Integer>();
	        ArrayList<Pair<Integer,Double>> predpair = new ArrayList<Pair<Integer,Double>>();
	        int jiaocha;
	        for (int i = 0,j = 0; i< answer.size(); i++,j++)
	        {
	            predpair = predictor.get(j);
	            golist = answer.get(i);
	            for(Pair<Integer,Double> iter:predpair)
	            {
	            	if (iter.getSecond() >= cut) predlist.add(iter.getFirst());
	            }
	            jiaocha = 0;
	            for (Integer iter:golist)
	            {
	            	if (predlist.contains(iter)) jiaocha++;
	            }
	            if (predlist.size()>0)
	            {
	                pre = (double)jiaocha/predlist.size();
	                prenum++;
	                presum += pre;
	            }
	            recallnum++;
	            recall = (double)jiaocha/golist.size();
	            recallsum += recall;
	        }
	        if (prenum > 0) presum /= prenum; else presum = 0;
	        recallsum /=recallnum;
	        if (presum + recallsum == 0) fmeasure = 0; else fmeasure = 2 * presum * recallsum / (presum+recallsum);
	        //System.out.println(cut + "\t" + fmeasure);
	        if (fmeasure > fmax)
	        {
	            fmax = fmeasure;
	            threshhold = cut;
	        }
	        predlist.clear();
	    }
	    
	    Pair<Double,Double> return_pair = new Pair<Double,Double>(fmax,threshhold);  
		return return_pair;
	}
	

}

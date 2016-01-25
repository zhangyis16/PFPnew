package protein;

import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.Scanner;
import java.util.Set;
import java.util.TreeMap;

import Main.learning;
import common.Pair;
import liblinear.Feature;
import liblinear.FeatureNode;
import liblinear.Linear;

public class oneProtein implements Comparable<oneProtein>
{
	private String sequence = new String();
    private String access = new String();
    private HashSet<Integer> Annotation = new HashSet<Integer>();
    private ArrayList<Pair<Integer,Double>> Prediction = new ArrayList<Pair<Integer,Double>>();
    
    private ArrayList<Pair<Integer,Double>> SequenceVector = new ArrayList<Pair<Integer,Double>>();
    private Feature[] liblinearFeature;
    
    
    private ArrayList<Pair<String,Double>>  blastResult = new ArrayList<Pair<String,Double>>();
    private HashMap<String,Integer> blastAccessIndex = new HashMap<String,Integer>();
    private String OriAlignSequence = new String();
    private ArrayList<String> MSASequence = new ArrayList<String>();
    
    private HashSet<Integer> blastCandidate = new HashSet<Integer>();
    
    public oneProtein(String access)
    {
    	if (access.contains("_")) access = proteinSet.MapName2UniAccess.get(access);
    	this.access = access;
    }
    public oneProtein(String access,HashSet<Integer> aAnnotation)
    {
    	this.access = access;
    	this.Annotation = aAnnotation;
    }
    public oneProtein(String access,String asequence)
    {
    	this.access = access;
    	sequence = asequence;
    }
    
    public void AddAnnotation(int aAnnotation)
    {
    	Annotation.add(aAnnotation);
    }
    public void addFather(GoSet aGoSet)
    {
    	ArrayList<Integer> tempor = new ArrayList<Integer>(Annotation);
    	for (int i=0;i<tempor.size();i++)
    	{
    		int gonum = tempor.get(i);
    		Set<Integer> fatherList = aGoSet.getFatherList(gonum);
    		for (Integer father:fatherList)
    		{
    			if (!Annotation.contains(father))
    			{
    				Annotation.add(father);
    				tempor.add(father);
    			}
    		}
    	}
    }

    public void addPred(Pair<Integer,Double> apred)
    {
    	Prediction.add(apred);
    }
    public void addSequence(String asequence)
    {
    	sequence = asequence;
    }
    public void addBlastResult(String Access,double bitscore)
    {
    	this.blastResult.add(new Pair<String,Double>(Access,bitscore));
    	int index = this.blastResult.size() - 1;
    	this.blastAccessIndex.put(Access, index);
    	this.MSASequence.add("");
    }
    public void addMSASequenceChar(char ch,String Acc)
    {
    	int index = this.blastAccessIndex.get(Acc);
    	String string = this.MSASequence.get(index) + ch;
    	this.MSASequence.set(index, string);
    }
    public void addOriAlignSequenceChar(char ch)
    {
    	this.OriAlignSequence += ch;
    }
    public void addMSAgap()
    {
    	for (int count = 0; count < this.MSASequence.size(); count++)
    	{
    		String e = this.MSASequence.get(count);
    		if (e.length()<OriAlignSequence.length())
    		{
    			for(int cou = 0;cou<OriAlignSequence.length() -this.MSASequence.get(count).length();cou++)
    				e = e + "-";
    		}
    		this.MSASequence.set(count, e);
    	}
    }
    
    public void addBlastCandidate(proteinSet train)
    {
    	for (int count=0;count<blastResult.size();count++)
    	{
    		String Acc = blastResult.get(count).getFirst();
    		for (Integer label:train.getAnnotation(Acc))
    		{
    			blastCandidate.add(label);
    		}
    	}
    }
    public void GoFDR(proteinSet train)
    {
    	this.addBlastCandidate(train);
    	
    }
    
    public void outputBlastResult()
    {
    	for (int count=0;count<blastResult.size();count++)
    		System.out.println(blastResult.get(count).getFirst() + "   " + blastResult.get(count).getSecond());
    }
    public void outputMSA(PrintWriter Fout)
    {
    	Fout.println(this.access + "  " + OriAlignSequence);
    	for (int count=0;count<MSASequence.size();count++)
    		Fout.println(this.blastResult.get(count).getFirst() + "  " + this.MSASequence.get(count));
    }
    public boolean compare(oneProtein other)
    {
    	if (Annotation.size() != other.getAnnotationSize()) return false;
    	for (Integer i:Annotation)
    		if (!other.containAnnotation(i)) return false;
    	return true;
    }
    
    
    public int compareTo(oneProtein o) 
    {
		return this.getName().compareTo(o.getName());
	}
    public boolean containAnnotation(int index)
    {
    	return Annotation.contains(index);
    }

    public void eraserAnnotation(int gonum)
    {
    	Annotation.remove(gonum);
    }
    public String getAccess()
    {
    	return this.access;
    }
    public String getSequence()
    {
    	return this.sequence;
    }
    public Feature[] getliblinearFeature()
    {
    	return this.liblinearFeature;
    }
    public HashSet<Integer> getAnnotation()
    {
    	HashSet<Integer> clone = (HashSet<Integer>)Annotation.clone();
		return clone;
    }
    public int getAnnotationSize()
    {
    	return this.Annotation.size();
    }
    public int getAnnotationSize(char sp)
    {
    	int num = 0;
    	for (Integer Ann:this.Annotation)
    	{
    		if (learning.aGoSet.getSpace(Ann) == sp)
    			num++;
    	}
    	return num;
    }
    public ArrayList<Pair<Integer,Double>> getFeature()
    {
    	ArrayList<Pair<Integer,Double>> clone = (ArrayList<Pair<Integer,Double>>)this.SequenceVector.clone();
		return clone;
    }
    public String getName()
    {
    	return proteinSet.MapUniAccess2Name.get(this.access);
    }
    public ArrayList<Pair<Integer,Double>> getPredList()
    {		
    	return this.Prediction;
    }
    public String getSpecies()
    {
    	return proteinCommon.getSpecies(this.access);
    }
    public HashSet<Integer> getSubAnnotation(char space)
    {
    	HashSet<Integer> newAnnotation = new HashSet<Integer>();
    	for (int node:Annotation)
    	{
    		if (learning.aGoSet.getSpace(node) == space)
    			newAnnotation.add(node);
    	}
    	return newAnnotation;
    }
    public oneProtein getSubProtein(char space)
    {
    	oneProtein aProtein = new oneProtein(this.access,getSubAnnotation(space));
    	for (Pair<Integer,Double> entry : Prediction)
    	{
    		int Annotation = entry.getFirst();
    		if (learning.aGoSet.getSpace(Annotation) == space) aProtein.addPred(entry);
    	}
    	return aProtein;
    }
    public boolean onlyGO5515()
    {
    	if ((Annotation.size()==1) && (Annotation.contains(5515))) 
    		return true; 
    	else 
    		return false;
    }
    public void OutputAnnotation(PrintWriter Fout)
    {
    	for (Integer i:Annotation)
    	Fout.println(this.access + "\t" + proteinCommon.GOInt2Str(i));
    }
    public void OutputAnnotationList(PrintWriter Fout,char c)
    {
    	for (Integer i:Annotation)
    		Fout.print(proteinCommon.GOInt2Str(i) + c);
    }
    public void OutputAnnotationNum(PrintWriter Fout,char c)
    {
    	int i = 0;
    	for (Integer e:Annotation)
    	{
    		i++;
    		if (i <= 1) Fout.print(e);
    		else {
    			Fout.print(c);
    			Fout.print(c + e);
    		}
    	}	
    }
    public void OutputAnnotationSpace(PrintWriter Fout)
    {
    	for (Integer i:Annotation)
    	Fout.println(this.access + "\t" + proteinCommon.GOInt2Str(i) + "\t" + learning.aGoSet.getSpace(i));
    }
    public void OutputFastaSequence(PrintWriter Fout)
    {
    	if (!sequence.equals(""))
    	{
    		Fout.println(">" + this.access);
    		int count = (sequence.length()-1)/60;
    		for (int i=0;i<=count;i++)
    			Fout.println(sequence.substring(i*60, Math.min((i+1)*60,sequence.length())));
    		Fout.println();
    	}
    	else System.out.println(this.access);
    }
    public void OutputPred(PrintWriter Fout)
    {
    	Fout.print(this.access + "\t");
    	for (Pair<Integer,Double> entry : Prediction)
    	{
    		Fout.print(entry.getFirst() + ":"+entry.getSecond()+" ");
    	}
    	Fout.println();
    }
    public void addPredResult(Scanner In)
    {
    	String line = In.nextLine();
    	
    }
    public void setAccess(String Access)
    {
    	this.access = Access;
    }
    public void setLiblinearFeature()
    {
		List<Feature> x = new ArrayList<Feature>();
		for (Pair<Integer,Double> pair:this.getFeature())
		{
    		Feature node = new FeatureNode(pair.getFirst(), pair.getSecond());
            x.add(node);
		}
		this.liblinearFeature = new Feature[x.size()];
		this.liblinearFeature = x.toArray(this.liblinearFeature);
    }
    public void setPredList(ArrayList<Pair<Integer,Double>> PredList)
    {
    	for (Pair<Integer,Double> pair:PredList)
    	{
    		this.Prediction.add(pair);
    	}
    }
    public void setBlastPred(proteinSet train)
    {
    	ArrayList<Integer> ann = new ArrayList<Integer>();
    	HashMap<Integer,Double> pred = new HashMap<Integer,Double>();
    	for (Pair<String,Double> pair:blastResult)
    	{
    		double bitscore = pair.getSecond();
    		String access = pair.getFirst();
    		ann = train.getAnnotation(access);
    		for (int node:ann)
    		{
    			if (!pred.containsKey(node)) 
    			{
    				pred.put(node, bitscore);
    				Prediction.add(new Pair<Integer,Double>(node, bitscore));
    			}
    		}
    	}
    }
    public void setGOtcha(proteinSet train)
    {
    	ArrayList<Integer> ann = new ArrayList<Integer>();
    	HashMap<Integer,Double> pred = new HashMap<Integer,Double>();
    	double SumBitScore = 0;
    	double bitscore;
    	for (Pair<String,Double> pair:blastResult)
    	{
    		bitscore = pair.getSecond();
    		SumBitScore = SumBitScore + bitscore;
    		String access = pair.getFirst();
    		ann = train.getAnnotation(access);
    		for (int node:ann)
    		{
    			if (!pred.containsKey(node)) 
    			{
    				pred.put(node, bitscore);
    			}else
    			{
    				pred.put(node, pred.get(node) + bitscore);
    			}
    		}
    	}
    	for (Map.Entry<Integer, Double> entry : pred.entrySet()) 
    	{
    		this.Prediction.add(new Pair<Integer,Double>(entry.getKey(),(double)entry.getValue()/SumBitScore));
    	}
    	
    }
    public void outputSequenceVector(PrintWriter Fout)
    {
    	for (Pair<Integer,Double> e:this.SequenceVector)
    	{
    		Fout.print(" " + e.getFirst() + ":" + e.getSecond());
    	}
    }
    public void tranSequence2Vector()
    {
    	Map<Integer,Double> map = new TreeMap<Integer,Double>();
		for (int j = 0; j<this.sequence.length()-2; j++)
		{
			String subSeq = this.sequence.substring(j, j+3);
			int dimension = proteinCommon.get3AcidIndex(subSeq);
			double num;
			if (dimension>0)
			{
				if (map.containsKey(dimension))
					num = map.get(dimension);
				else
					num = 0.0;
				map.put(dimension, num + 1.0);
			}
		}
		
		for (Map.Entry<Integer,Double> entry : map.entrySet()) 
    	{
    		this.SequenceVector.add(new Pair<Integer,Double>(entry.getKey(),entry.getValue()));
    	}
    }
    public void filterGoNotIn(GoSet aGo)
    {
    	HashSet<Integer> aAnnotation = new HashSet<Integer>();
    	for (Integer i:this.Annotation)
    	{
    		if (aGo.containNode(i))  aAnnotation.add(i);
    	}
    	this.Annotation = aAnnotation;
    }
    public void removeAnnotation(Set<Integer> args)
    {
    	for (int node : args)
    	{
    		Annotation.remove(node);
    	}
    	Iterator<Pair<Integer, Double>> iter = Prediction.iterator();
    	while (iter.hasNext())
    	{
    		if (args.contains(iter.next().getFirst())) iter.remove();
    	}
    }
    public void libLinearPredict(int labelNum,liblinear.Model model)
    {
    	double[] prob_estimates = new double[2];            
		double predict_label = Linear.predictProbability(model, this.liblinearFeature, prob_estimates);
		this.addPred(new Pair<Integer,Double>(labelNum,prob_estimates[0]));
    }
    public void removeAnnotation(int... args)
    {
    	Set<Integer> SetAnn = new HashSet<Integer>();
    	for (int node : args)
    	{
    		SetAnn.add(node);
    	}
    	this.removeAnnotation(SetAnn);
    }

    public void removeRedunAnn(GoSet aGoSet)
    {
    	ArrayList<Integer> Ann = new ArrayList<Integer>();
    	for (int ann:this.Annotation)
    		Ann.add(ann);
    	for (int i = 0;i<Ann.size()-1;i++)
    		for (int j=0;j<Ann.size();j++)
    		{
    			int ann1 = Ann.get(i);
    			int ann2 = Ann.get(j);
    			if (aGoSet.checkFather(ann1, ann2) != 0)
    			{
    				int father = aGoSet.checkFather(ann1, ann2);
    				this.Annotation.remove(father);
    			}
    		}
    	
    }
}

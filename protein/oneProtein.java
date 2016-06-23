package protein;

import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.Set;

import Main.learning;
import common.FeatureTran;
import common.Pair;
import liblinear.Feature;
import liblinear.FeatureNode;
import liblinear.Linear;

public class oneProtein implements Comparable<oneProtein> {
	private String sequence = new String();
	private String access = new String();

	private HashSet<Integer> Annotation = new HashSet<Integer>();
	private int MFOsize = 0;
	private int BPOsize = 0;
	private int CCOsize = 0;
	private HashSet<Integer> IEA_Annotation = new HashSet<Integer>();
	private ArrayList<Pair<Integer, Double>> PredictionScore = new ArrayList<Pair<Integer, Double>>();

	private ArrayList<Pair<Integer, Double>> SparseFeature = new ArrayList<Pair<Integer, Double>>();
	private ArrayList<Double> NoSparseFeature = new ArrayList<Double>();
	private Feature[] liblinearFeature;

	private ArrayList<Pair<String, Double>> blastResult = new ArrayList<Pair<String, Double>>();
	private HashMap<String, Integer> blastAccessIndex = new HashMap<String, Integer>();
	
	private HashSet<Integer> L2RCandidate  = new HashSet<Integer>();
		
	private String OriAlignSequence = new String();
	private HashSet<Integer> blastCandidate = new HashSet<Integer>();
	private ArrayList<String> MSASequence = new ArrayList<String>();

	public int getMFOsize()
	{
		return this.MFOsize;
	}
	public int getBPOsize()
	{
		return this.BPOsize;
	}
	public int getCCOsize()
	{
		return this.CCOsize;
	}
	public void calMFBPCCsize(GoSet aGoSet)
	{
		this.MFOsize = 0;
		this.BPOsize = 0;
		this.CCOsize = 0;
		for (int node:this.Annotation)
		{
			if (aGoSet.getSpace(node) == 'F') this.MFOsize++;
			if (aGoSet.getSpace(node) == 'P') this.BPOsize++;
			if (aGoSet.getSpace(node) == 'C') this.CCOsize++;
		}
		
	}
	public int getSubAnnSize(char sp)
	{
		if (sp == 'F') return this.MFOsize;
		if (sp == 'P') return this.BPOsize;
		if (sp == 'C') return this.CCOsize;
		return 0;
	}
	public double CalCandidateRecall(char space)
	{
		double recall = 0;
		for (Integer e:L2RCandidate)
		{
			if (Annotation.contains(e)) recall+=1;
		}
		System.out.println(recall);
		recall = (double)recall/Annotation.size();
		return recall;
	}
	
	public void clearPredResult()
	{
		this.PredictionScore.clear();
	}
	public oneProtein(String access) {
		if (access.contains("_"))
			access = proteinSet.MapName2UniAccess.get(access);
		this.access = access;
	}

	public void addTopK_L2RCandidate(int k,char space)
	{
		Comparator<Pair<Integer,Double>> compar = new Comparator<Pair<Integer,Double>>()
		{

			@Override
			public int compare(Pair<Integer, Double> o1, Pair<Integer, Double> o2) {
				// TODO Auto-generated method stub
				if (o1.getSecond()>o2.getSecond()) return -1; else return 1;
			}
		};
		Collections.sort(this.PredictionScore, compar);
		for (int i = 0; i<this.PredictionScore.size();i++)
		{
			int gonum = this.PredictionScore.get(i).getFirst();
			while (learning.aGoSet.getSpace(gonum) == space)
			{
				
			}
			this.L2RCandidate.add(this.PredictionScore.get(i).getFirst());
		}
		
	}
	public oneProtein(String access, HashSet<Integer> aAnnotation) {
		this.access = access;
		this.Annotation = aAnnotation;
	}

	public oneProtein(String access, String asequence) {
		this.access = access;
		sequence = asequence;
	}

	public void AddAnnotation(int aAnnotation) {
		Annotation.add(aAnnotation);
	}

	public void addBlastCandidate(proteinSet train) 
	{
		for (int count = 0; count < blastResult.size(); count++) {
			String Acc = blastResult.get(count).getFirst();
			for (Integer label : train.getAnnotation(Acc)) {
				blastCandidate.add(label);
			}
		}
	}

	public void addBlastResult(String Access, double bitscore) 
	{
		if (!this.blastAccessIndex.containsKey(Access))
		{	
			this.blastResult.add(new Pair<String, Double>(Access, bitscore));
			int index = this.blastResult.size() - 1;
			this.blastAccessIndex.put(Access, index);
		}else
		{
			int index = this.blastAccessIndex.get(Access);
			double score = this.blastResult.get(index).getSecond();
			if (score<bitscore)
				this.blastResult.get(index).setSecond(bitscore);
		}
		
	}
	
	public void removeBlastOtherSp()
	{
		ArrayList<Pair<String, Double>> Result = new ArrayList<Pair<String, Double>>();
		for (int i = 0;i<this.blastResult.size(); i++)
		{
			String Acc = this.blastResult.get(i).getFirst();
			if (proteinCommon.getSpecies(Acc).equals(proteinCommon.getSpecies(this.access))) Result.add(blastResult.get(i));
		}
		blastResult = Result;
	}

	public void addFather(GoSet aGoSet) {
		ArrayList<Integer> tempor = new ArrayList<Integer>(Annotation);
		for (int i = 0; i < tempor.size(); i++) {
			int gonum = tempor.get(i);
			Set<Integer> fatherList = aGoSet.getFatherList(gonum);
			for (Integer father : fatherList) {
				if (!Annotation.contains(father)) {
					Annotation.add(father);
					tempor.add(father);
				}
			}
		}
	}

	public void AddIEA_Annotation(int aAnnotation) {
		IEA_Annotation.add(aAnnotation);
	}

	public void addIEAFather(GoSet aGoSet) {
		ArrayList<Integer> tempor = new ArrayList<Integer>(IEA_Annotation);
		for (int i = 0; i < tempor.size(); i++) {
			int gonum = tempor.get(i);
			Set<Integer> fatherList = aGoSet.getFatherList(gonum);
			for (Integer father : fatherList) {
				if (!IEA_Annotation.contains(father)) {
					IEA_Annotation.add(father);
					tempor.add(father);
				}
			}
		}
	}

	public void addMSAgap() {
		for (int count = 0; count < this.MSASequence.size(); count++) {
			String e = this.MSASequence.get(count);
			if (e.length() < OriAlignSequence.length()) {
				for (int cou = 0; cou < OriAlignSequence.length() - this.MSASequence.get(count).length(); cou++)
					e = e + "-";
			}
			this.MSASequence.set(count, e);
		}
	}

	public void addMSASequenceChar(char ch, String Acc) {
		int index = this.blastAccessIndex.get(Acc);
		String string = this.MSASequence.get(index) + ch;
		this.MSASequence.set(index, string);
	}

	public void SetNoSparseFeature(ArrayList<Double> Fea) 
	{
		this.NoSparseFeature.clear();
		for (Double e:Fea)
		{
			this.NoSparseFeature.add(e);
		}
	}

	public void addOriAlignSequenceChar(char ch) {
		this.OriAlignSequence += ch;
	}

	public void addPred(Pair<Integer, Double> apred) {
		PredictionScore.add(apred);
	}

	public void addSequence(String asequence) {
		sequence = asequence;
	}

	public boolean compare(oneProtein other) {
		if (Annotation.size() != other.getAnnotationSize())
			return false;
		for (Integer i : Annotation)
			if (!other.containAnnotation(i))
				return false;
		return true;
	}

	public int compareTo(oneProtein o) {
		return this.getName().compareTo(o.getName());
	}

	public boolean containAnnotation(int index) {
		return Annotation.contains(index);
	}

	public void eraserAnnotation(int gonum) {
		Annotation.remove(gonum);
	}

	public Pair<Double, Double> evalutionPR(char space) {
		double pre = 0;
		double recall = 0;
		HashSet<Integer> Annotation = this.getSubAnnotation(space);

		HashSet<Integer> IEA_Annotation = this.getSubIEA_Annotation(space);

		for (int ann : Annotation) {
			if (IEA_Annotation.contains(ann))
				recall += 1;
		}
		for (int ann : IEA_Annotation) {
			if (Annotation.contains(ann))
				pre += 1;
		}
		if (Annotation.size() == 0)
			return new Pair<Double, Double>(-1.0, -1.0);
		else
			recall = recall / Annotation.size();

		if (IEA_Annotation.size() > 0)
			pre = pre / IEA_Annotation.size();
		else
			pre = -1;
		return new Pair<Double, Double>(pre, recall);
	}

	public void filterGoNotIn(GoSet aGo) {
		HashSet<Integer> aAnnotation = new HashSet<Integer>();
		for (Integer i : this.Annotation) {
			if (aGo.containNode(i))
				aAnnotation.add(i);
			else
				System.out.println("remove" + this.access + " : " + i);
		}
		this.Annotation = aAnnotation;
	}

	public String getAccess() {
		return this.access;
	}

	public HashSet<Integer> getAnnotation() 
	{
		HashSet<Integer> clone = new HashSet<Integer>();
		for (Integer e:this.Annotation)
			clone.add(e);
		return clone;
	}

	public int getFeatureSize()
	{
		return this.liblinearFeature.length;
	}
	public int getBlastNum()
	{
		return this.blastAccessIndex.size();
	}
	
	public int getAnnotationSize() {
		return this.Annotation.size();
	}

	public int getAnnotationSize(char sp) {
		int num = 0;
		for (Integer Ann : this.Annotation) {
			if (learning.aGoSet.getSpace(Ann) == sp)
				num++;
		}
		return num;
	}
	public int getAnnotationSize(String sp)
	{
		if (sp.equals("MFO")) return this.getAnnotationSize('F');
		if (sp.equals("BPO")) return this.getAnnotationSize('P');
		if (sp.equals("CCO")) return this.getAnnotationSize('C');
		return 0;
	}

	public Feature[] getliblinearFeature() {
		return this.liblinearFeature;
	}

	public String getName() {
		return proteinSet.MapUniAccess2Name.get(this.access);
	}

	public ArrayList<Pair<Integer, Double>> getPredList() {
		return this.PredictionScore;
	}

	public String getSequence() {
		return this.sequence;
	}

	public ArrayList<Pair<Integer, Double>> getSparseFeature() 
	{
		ArrayList<Pair<Integer, Double>> clone = new ArrayList<Pair<Integer, Double>>();
		for (Pair<Integer, Double> e:this.SparseFeature)
			clone.add(e);
		return clone;
	}

	public String getSpecies() {
		return proteinCommon.getSpecies(this.access);
	}

	public HashSet<Integer> getSubAnnotation(char space) {
		HashSet<Integer> newAnnotation = new HashSet<Integer>();
		for (int node : Annotation) {
			if (learning.aGoSet.getSpace(node) == space)
				newAnnotation.add(node);
		}
		return newAnnotation;
	}

	public HashSet<Integer> getSubIEA_Annotation(char space) {
		HashSet<Integer> newAnnotation = new HashSet<Integer>();
		for (int node : this.IEA_Annotation) {
			if (learning.aGoSet.getSpace(node) == space)
				newAnnotation.add(node);
		}
		return newAnnotation;
	}

	public oneProtein getSubProtein(char space) {
		oneProtein aProtein = new oneProtein(this.access, getSubAnnotation(space));
		for (Pair<Integer, Double> entry : PredictionScore) {
			int Annotation = entry.getFirst();
			if (learning.aGoSet.getSpace(Annotation) == space)
				aProtein.addPred(entry);
		}
		return aProtein;
	}

	public void GoFDR(proteinSet train) {
		this.addBlastCandidate(train);

	}

	public void libLinearPredict(int labelNum, liblinear.Model model) 
	{
		double[] prob_estimates = new double[2];
		double predict_label = Linear.predictProbability(model, this.liblinearFeature, prob_estimates);
		this.addPred(new Pair<Integer, Double>(labelNum, prob_estimates[0]));
	}

	public boolean onlyGO5515() {
		if ((Annotation.size() == 1) && (Annotation.contains(5515)))
			return true;
		else
			return false;
	}

	public void OutputAnnotation(PrintWriter Fout) 
	{
		ArrayList<Integer> Ann = new ArrayList<Integer>();
		for (Integer i : Annotation)
		{
			Ann.add(i);
		}
		Collections.sort(Ann);
		for (Integer i : Ann)
		{
			Fout.println(this.access + "\t" + proteinCommon.GOInt2Str(i));
		}
	}

	public void OutputAnnotationList(PrintWriter Fout, char c) {
		Fout.print(this.access);
		ArrayList<Integer> Ann = new ArrayList<Integer>();
		for (Integer i : Annotation) Ann.add(i);
		Collections.sort(Ann);
		for (Integer i : Ann)
			Fout.print(proteinCommon.GOInt2Str(i) + c);
	}

	public void OutputAnnotationNum(PrintWriter Fout, char c) {
		int i = 0;
		for (Integer e : Annotation) {
			i++;
			if (i <= 1)
				Fout.print(e);
			else {
				Fout.print(c);
				Fout.print(c + e);
			}
		}
	}

	public void OutputAnnotationSpace(PrintWriter Fout) {
		for (Integer i : Annotation)
			Fout.println(this.access + "\t" + proteinCommon.GOInt2Str(i) + "\t" + learning.aGoSet.getSpace(i));
	}

	public void outputBlastResult() {
		for (int count = 0; count < blastResult.size(); count++)
			System.out.println(blastResult.get(count).getFirst() + "   " + blastResult.get(count).getSecond());
	}

	public void OutputFastaSequence(PrintWriter Fout) {
		if (!sequence.equals("")) {
			Fout.println(">" + this.access);
			int count = (sequence.length() - 1) / 60;
			for (int i = 0; i <= count; i++)
				Fout.println(sequence.substring(i * 60, Math.min((i + 1) * 60, sequence.length())));
			Fout.println();
		} else
			System.out.println(this.access + "hava no seq");
	}

	public void outputMSA(PrintWriter Fout) {
		Fout.println(this.access + "  " + OriAlignSequence);
		for (int count = 0; count < MSASequence.size(); count++)
			Fout.println(this.blastResult.get(count).getFirst() + "  " + this.MSASequence.get(count));
	}

	public void OutputNoSparseFeature(PrintWriter Fout) {
		Fout.println(this.access);
		for (Double e : this.NoSparseFeature) {
			Fout.print(e + " ");
		}
		Fout.println();
	}

	public void OutputPredScore(PrintWriter Fout) {
		Fout.println(this.access);
		Fout.print(this.PredictionScore.size() + " ");
		for (Pair<Integer, Double> entry : this.PredictionScore) {
			Fout.printf("%d %.4f ",entry.getFirst(),entry.getSecond());
		}
		Fout.println();
	}

	public void OutputSparseFeature(PrintWriter Fout) 
	{
		for (Pair<Integer, Double> e : this.SparseFeature) {
			Fout.print(" " + e.getFirst() + ":" + e.getSecond());
		}
	}

	public void removeAnnotation(int... args) {
		Set<Integer> SetAnn = new HashSet<Integer>();
		for (int node : args) 
		{
			SetAnn.add(node);
		}
		this.removeAnnotation(SetAnn);
	}

	public void removeAnnotation(Set<Integer> args) {
		for (int node : args) {
			Annotation.remove(node);
		}
		Iterator<Pair<Integer, Double>> iter = PredictionScore.iterator();
		while (iter.hasNext()) {
			if (args.contains(iter.next().getFirst()))
				iter.remove();
		}
	}

	public void removeIEA_Annotation(int... args) {
		Set<Integer> SetAnn = new HashSet<Integer>();
		for (int node : args) {
			SetAnn.add(node);
		}
		this.removeIEA_Annotation(SetAnn);
	}

	public void removeIEA_Annotation(Set<Integer> args) {
		for (int node : args) {
			IEA_Annotation.remove(node);
		}
		Iterator<Pair<Integer, Double>> iter = PredictionScore.iterator();
		while (iter.hasNext()) {
			if (args.contains(iter.next().getFirst()))
				iter.remove();
		}
	}

	public void removeRedunAnn(GoSet aGoSet) {
		ArrayList<Integer> Ann = new ArrayList<Integer>();
		for (int ann : this.Annotation)
			Ann.add(ann);
		for (int i = 0; i < Ann.size() - 1; i++)
			for (int j = 0; j < Ann.size(); j++) {
				int ann1 = Ann.get(i);
				int ann2 = Ann.get(j);
				if (aGoSet.checkFather(ann1, ann2) != 0) {
					int father = aGoSet.checkFather(ann1, ann2);
					this.Annotation.remove(father);
				}
			}

	}
	
	public void removeLowPred(int count)
	{
		ArrayList<Pair<Double,Integer>> ArrSort = new ArrayList<Pair<Double,Integer>>();
		for (Pair<Integer,Double> pair:this.PredictionScore)
		{
			ArrSort.add(new Pair<Double,Integer>(pair.getSecond(),pair.getFirst()));
		}
		Collections.sort(ArrSort);
		this.PredictionScore.clear();
		int index = Math.max(0, ArrSort.size() - count);
		for (int i = ArrSort.size() - 1; i >= index; i--)
		{
			this.PredictionScore.add(new Pair<Integer,Double>(ArrSort.get(i).getSecond(),ArrSort.get(i).getFirst()));
		}
	}

	public void setAccess(String Access) {
		this.access = Access;
	}

	public void setBlastPred(proteinSet train) {
		ArrayList<Integer> ann = new ArrayList<Integer>();
		HashMap<Integer, Double> pred = new HashMap<Integer, Double>();
		for (Pair<String, Double> pair : blastResult) {
			double bitscore = pair.getSecond();
			String access = pair.getFirst();
			ann = train.getAnnotation(access);
			for (int node : ann) {
				if (!pred.containsKey(node)) {
					pred.put(node, bitscore);
					PredictionScore.add(new Pair<Integer, Double>(node, bitscore));
				}
			}
		}
	}

	public void setGOtcha(proteinSet train) {
		if (this.blastResult.size() == 0)
		{
			System.out.println(access + " have no blast result");
		}
		ArrayList<Integer> ann = new ArrayList<Integer>();
		HashMap<Integer, Double> pred = new HashMap<Integer, Double>();
		double SumBitScore = 0;
		double bitscore;
		for (Pair<String, Double> pair : blastResult) {
			bitscore = pair.getSecond();
			SumBitScore = SumBitScore + bitscore;
			String access = pair.getFirst();
			if (train.containProtein(access))
			{
				ann = train.getAnnotation(access);
				for (int node : ann) {
					if (!pred.containsKey(node)) {
						pred.put(node, bitscore);
					} else {
						pred.put(node, pred.get(node) + bitscore);
					}
				}
			}
		}
		for (Map.Entry<Integer, Double> entry : pred.entrySet()) {
			this.PredictionScore
					.add(new Pair<Integer, Double>(entry.getKey(), (double) entry.getValue() / SumBitScore));
		}

	}

	public void setLiblinearFeatureFromSparseFeature() 
	{
		List<Feature> x = new ArrayList<Feature>();
		for (Pair<Integer, Double> pair : this.SparseFeature) {
			Feature node = new FeatureNode(pair.getFirst(), pair.getSecond());
			x.add(node);
		}
		this.liblinearFeature = new Feature[x.size()];
		this.liblinearFeature = x.toArray(this.liblinearFeature);
		this.SparseFeature.clear();
	}

	public void setLiblinearFeatureFromNoSparseFeature() 
	{
		List<Feature> x = new ArrayList<Feature>();
		int count = 1;
		for (Double num : this.NoSparseFeature) 
		{
			Feature node = new FeatureNode(count++ , num);
			x.add(node);
		}
		this.liblinearFeature = new Feature[x.size()];
		this.liblinearFeature = x.toArray(this.liblinearFeature);
	}
	
	public void setPredResult(ArrayList<Pair<Integer, Double>> PredList) 
	{
		for (Pair<Integer, Double> pair : PredList) 
		{
			this.PredictionScore.add(pair);
		}
	}

	public void tranSequence2TriSparseFeature() 
	{
		//transfer the sequence to 3 amino acid feature
		
		ArrayList<Double> Fea = new ArrayList<Double>();
		Fea.clear();
		for (int i = 0;i <= 8000 + 400;i++)
			Fea.add(0.0);
		
		for (int j = 0; j < this.sequence.length() - 2; j++) 
		{
			String subSeq = this.sequence.substring(j, j + 3);
			int dimension = proteinCommon.get3AcidIndex(subSeq);
			if (dimension<=8000)
				Fea.set(dimension, Fea.get(dimension)+ (double)500/sequence.length());		
		}
		/*if (learning.isFeaTran)
		{
			System.out.println("Feature Transform Now ......." + this.access);
			Fea = FeatureTran.VecMuitiMatrix(Fea, proteinCommon.TriAcidSimiliar);
		}    */
		for (int j = 0; j < this.sequence.length() - 1; j++) 
		{
			String subSeq = this.sequence.substring(j, j + 2);
			int dimension = proteinCommon.get2AcidIndex(subSeq) + 8000;
			if (dimension<=8000+400)
				Fea.set(dimension, Fea.get(dimension)+ (double)500/sequence.length());		
		}
		for (int i=1;i<=8000+400;i++)
		{
			if (Fea.get(i)>0.001) this.SparseFeature.add(new Pair<Integer, Double>(i, Fea.get(i)));
		}
	}
}

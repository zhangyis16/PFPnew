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

import Main.learningOfGO;
import common.Pair;
import common.Parameter;
import hpo.HPOSet;
import liblinear.Feature;
import liblinear.FeatureNode;
import liblinear.Linear;

public class oneProtein implements Comparable<oneProtein> {
	private String sequence = new String();
	private String access = new String();

	private HashSet<Integer> Annotation = new HashSet<Integer>();
	private HashSet<Integer> HPOAnnotation = new HashSet<Integer>();
	
	private int MFOsize = 0;
	private int BPOsize = 0;
	private int CCOsize = 0;
	private HashSet<Integer> IEA_Annotation = new HashSet<Integer>();
	private ArrayList<Pair<Integer, Double>> PredictionScore = new ArrayList<Pair<Integer, Double>>();
	private ArrayList<Pair<Integer, Double>> HPOPredictionScore = new ArrayList<Pair<Integer, Double>>();
	
	private ArrayList<Pair<Integer, Double>> SparseFeature = new ArrayList<Pair<Integer, Double>>();
	private ArrayList<Double> NoSparseFeature = new ArrayList<Double>();
	private Feature[] liblinearFeature;

	private ArrayList<Pair<String, Double>> blastResult = new ArrayList<Pair<String, Double>>();
	private ArrayList<Pair<String, Double>> Similiar = new ArrayList<Pair<String, Double>>();
	private HashMap<String, Integer> blastAccessIndex = new HashMap<String, Integer>();
	
	
	
	private HashSet<Integer> MFOL2RCandidate  = new HashSet<Integer>();
	private HashSet<Integer> BPOL2RCandidate  = new HashSet<Integer>();
	private HashSet<Integer> CCOL2RCandidate  = new HashSet<Integer>();
	
	private HashMap<Integer,Double> liblinearScore = new HashMap<Integer,Double>();
	private HashMap<Integer,Double> blastKnnScore = new HashMap<Integer,Double>();
	private HashMap<Integer,Double> blastScore = new HashMap<Integer,Double>();
	
	
	private ArrayList<Integer> PubMedID = new ArrayList<Integer>();
	
	private ArrayList<String> IntActProtein = new ArrayList<String>();
	
	private String OriAlignSequence = new String();
	private HashSet<Integer> blastCandidate = new HashSet<Integer>();
	private ArrayList<String> MSASequence = new ArrayList<String>();

	private ArrayList<Pair<String, String>> DatabaseReference = new ArrayList<Pair<String, String>>();
	
	private int integragedYear = 0;
	
	public double calL2RcandidateRecall(char space)
	{
		double recall = 0;
		if (space == 'F') 
		{
			for (Integer ann:MFOL2RCandidate)
			{
				if (this.Annotation.contains(ann))
					recall += 1.0;
			}
			recall /= this.MFOsize;
		}
		if (space == 'P') 
		{
			for (Integer ann:BPOL2RCandidate)
			{
				if (this.Annotation.contains(ann))
					recall += 1.0;
			}
			recall /= this.BPOsize;
		}
		if (space == 'C') 
		{
			for (Integer ann:CCOL2RCandidate)
			{
				if (this.Annotation.contains(ann))
					recall += 1.0;
			}
			recall /= this.CCOsize;
		}
		return recall;
	}
	
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
	public int getStoreAnnotationSize(char sp)
	{
		if (sp == 'F') return this.getMFOsize();
		if (sp == 'P') return this.getBPOsize();
		if (sp == 'C') return this.getCCOsize();
		else return 0;
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
	public int getPubMedsize()
	{
		return this.PubMedID.size();
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
		return 0.0;
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
	public void sortPredScoreFromHighToLow()
	{
		Comparator<Pair<Integer,Double>> compar = new Comparator<Pair<Integer,Double>>()
		{

			@Override
			public int compare(Pair<Integer, Double> o1, Pair<Integer, Double> o2) {
				// TODO Auto-generated method stub
				if (o1.getSecond()>o2.getSecond()) return -1; else 
					if (o1.getSecond()<o2.getSecond())	return 1; else return 0;
			}
		};
		Collections.sort(this.PredictionScore, compar);
	}
	
	public void addTopK_L2RCandidate(int k,char space)
	{
		int count = 0;
		for (int i = 0; i<this.PredictionScore.size();i++)
		{
			int gonum = this.PredictionScore.get(i).getFirst();
			
			if (i>0) 
				if (this.PredictionScore.get(i-1).getSecond() < this.PredictionScore.get(i).getSecond())
					System.out.printf("Prediction Score Sort Error %.4f %.4f",this.PredictionScore.get(i-1).getSecond()
							,this.PredictionScore.get(i).getSecond());
			
			if (learningOfGO.aGoSet.getSpace(gonum) == space)
			{
				count++;
				if (space == 'F') this.MFOL2RCandidate.add(this.PredictionScore.get(i).getFirst());
				if (space == 'P') this.BPOL2RCandidate.add(this.PredictionScore.get(i).getFirst());
				if (space == 'C') this.CCOL2RCandidate.add(this.PredictionScore.get(i).getFirst());	
			}
			if (count>=k) break;

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
	
	public void AddHPOAnnotation(int aAnnotation) {
		HPOAnnotation.add(aAnnotation);
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
	
	public void addSimiliar(String Access, double sim) 
	{
		this.Similiar.add(new Pair<String,Double>(Access,sim));
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

	public void addGOFather(GoSet aGoSet) 
	{
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
	
	public void addHPOFather(HPOSet aHPOSet) 
	{
		ArrayList<Integer> tempor = new ArrayList<Integer>(HPOAnnotation);
		for (int i = 0; i < tempor.size(); i++) 
		{
			int HPOnum = tempor.get(i);
			Set<Integer> fatherList = aHPOSet.getFatherList(HPOnum);
			for (Integer father : fatherList) 
			{
				if (!HPOAnnotation.contains(father)) 
				{
					HPOAnnotation.add(father);
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
	public void addNoScoreLabelRandom(ArrayList<Pair<Integer,Double>> predList,double miniscore)
	{
		HashSet<Integer> predLabelSet = new HashSet<Integer>();
		for (int i = 0;i<this.PredictionScore.size();i++)
		{
			int label = this.PredictionScore.get(i).getFirst();
			predLabelSet.add(label);
		}
		
		for (int i=0;i<predList.size();i++)
		{
			int label = predList.get(i).getFirst();
			double fre = predList.get(i).getSecond();
			if (!predLabelSet.contains(label))
				this.PredictionScore.add(new Pair<Integer,Double>(label,fre*miniscore));
		}
	}

	public void addOriAlignSequenceChar(char ch) {
		this.OriAlignSequence += ch;
	}

	public void addPred(Pair<Integer, Double> apred) 
	{
		PredictionScore.add(apred);
	}

	public void addPredList(ArrayList<Pair<Integer,Double>> aPredList) 
	{
		PredictionScore = aPredList;
	}
	
	public void addHPOPredList(ArrayList<Pair<Integer,Double>> aPredList) 
	{
		HPOPredictionScore = aPredList;
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

	public void filterGoNotIn(GoSet aGo) 
	{
		HashSet<Integer> aAnnotation = new HashSet<Integer>();
		for (Integer i : this.Annotation) {
			if (aGo.containNode(i))
				aAnnotation.add(i);
			else
				System.out.println("remove" + this.access + " : " + i);
		}
		this.Annotation = aAnnotation;
	}
	
	public void filterHPONotIn(HPOSet aHPOSet) 
	{
		HashSet<Integer> aAnnotation = new HashSet<Integer>();
		for (Integer i : this.HPOAnnotation) {
			if (aHPOSet.containNode(i))
				aAnnotation.add(i);
			else
				System.out.println("remove the HPO of " + this.access + " : " + i);
		}
		this.HPOAnnotation = aAnnotation;
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
	
	public HashSet<Integer> getHPOAnnotation() 
	{
		HashSet<Integer> clone = new HashSet<Integer>();
		for (Integer e:this.HPOAnnotation)
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
	
	public int getAnnotationSize() 
	{
		return this.Annotation.size();
	}
	
	public int getHPOAnnotationSize() 
	{
		return this.HPOAnnotation.size();
	}
	
	public int getAnnotationSize(char sp) 
	{
		int num = 0;
		for (Integer Ann : this.Annotation) {
			if (learningOfGO.aGoSet.getSpace(Ann) == sp)
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

	
	public ArrayList<Pair<Integer, Double>> getHPOPredList() {
		return this.HPOPredictionScore;
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
			if (learningOfGO.aGoSet.getSpace(node) == space)
				newAnnotation.add(node);
		}
		return newAnnotation;
	}

	public HashSet<Integer> getSubIEA_Annotation(char space) {
		HashSet<Integer> newAnnotation = new HashSet<Integer>();
		for (int node : this.IEA_Annotation) {
			if (learningOfGO.aGoSet.getSpace(node) == space)
				newAnnotation.add(node);
		}
		return newAnnotation;
	}

	public double miniscore()
	{
		double miniscore = 1.0;
		for (Pair<Integer, Double> entry : PredictionScore)
		{
			if (miniscore>entry.getSecond()) miniscore = entry.getSecond();
		}
		return miniscore;
	}
	public oneProtein getSubProtein(char space) 
	{
		oneProtein aProtein = new oneProtein(this.access, getSubAnnotation(space));
		for (Pair<Integer, Double> entry : PredictionScore) {
			int Annotation = entry.getFirst();
			if (learningOfGO.aGoSet.getSpace(Annotation) == space)
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
	
	public void OutputHPOAnnotation(PrintWriter Fout) 
	{
		ArrayList<Integer> Ann = new ArrayList<Integer>();
		for (Integer i : HPOAnnotation)
		{
			Ann.add(i);
		}
		Collections.sort(Ann);
		for (Integer i : Ann)
		{
			Fout.println(this.access + "\t" + proteinCommon.HPInt2Str(i));
		}
	}
	
	public void OutputAnnotationName(PrintWriter Fout) 
	{
		ArrayList<Integer> Ann = new ArrayList<Integer>();
		for (Integer i : Annotation)
		{
			Ann.add(i);
		}
		Collections.sort(Ann);
		for (Integer i : Ann)
		{
			Fout.println(this.getName() + "\t" + proteinCommon.GOInt2Str(i));
		}
	}
	
	public void OutputType2Annotation(PrintWriter Fout,GoSet aGoSet) 
	{
		ArrayList<Integer> Ann = new ArrayList<Integer>();
		for (Integer i : Annotation)
		{
			Ann.add(i);
		}
		Collections.sort(Ann);
		for (Integer i : Ann)
		{
			if (this.getStoreAnnotationSize(aGoSet.getSpace(i))>0)
				Fout.println(this.access + "\t" + proteinCommon.GOInt2Str(i));
		}
	}
	
	public void OutputType3Annotation(PrintWriter Fout,GoSet aGoSet) 
	{
		ArrayList<Integer> Ann = new ArrayList<Integer>();
		for (Integer i : Annotation)
		{
			Ann.add(i);
		}
		Collections.sort(Ann);
		for (Integer i : Ann)
		{
			if (this.getStoreAnnotationSize(aGoSet.getSpace(i))==0)
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
			Fout.println(this.access + "\t" + proteinCommon.GOInt2Str(i) + "\t" + learningOfGO.aGoSet.getSpace(i));
	}

	public void outputBlastResult() {
		for (int count = 0; count < blastResult.size(); count++)
			System.out.println(blastResult.get(count).getFirst() + "   " + blastResult.get(count).getSecond());
	}

	public void OutputFastaAccessSequence(PrintWriter Fout) 
	{
		if (!sequence.equals(""))
		{
			String name = this.getAccess();
			String seq  = this.sequence;
			proteinCommon.outputFasta(Fout, name, seq);
		} 
		else
			System.out.println(this.access + "  hava no seq");
	}
	
	public void OutputFastaNameSequence(PrintWriter Fout) 
	{
		if (!sequence.equals(""))
		{
			String name = this.getName();
			String seq  = this.sequence;
			proteinCommon.outputFasta(Fout, name, seq);
		} 
		else
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
	
	public void outputPubMedID(PrintWriter Fout) 
	{
		Fout.println(this.access);
		Fout.print(this.PubMedID.size() + " ");
		for (Integer e : this.PubMedID) 
		{
			Fout.print(e + " ");
		}
		Fout.println();
	}

	public void outputIntActProtein(PrintWriter Fout) 
	{
		Fout.println(this.access);
		Fout.print(this.IntActProtein.size() + " ");
		for (String e : this.IntActProtein) 
		{
			Fout.print(e + " ");
		}
		Fout.println();
	}
	
	public void OutputGOPredScore(PrintWriter Fout) {
		Fout.println(this.access);
		Fout.print(this.PredictionScore.size());
		for (Pair<Integer, Double> entry : this.PredictionScore) {
			Fout.printf(" %d:%.4f",entry.getFirst(),entry.getSecond());
		}
		Fout.println();
	}

	public void OutputHPOPredScore(PrintWriter Fout) 
	{
		Fout.println(this.access);
		Fout.print(this.HPOPredictionScore.size());
		for (Pair<Integer, Double> entry : this.HPOPredictionScore) {
			Fout.printf(" %d:%.4f",entry.getFirst(),entry.getSecond());
		}
		Fout.println();
	}
	
	public void OutputSparseFeature(PrintWriter Fout) 
	{
		for (Pair<Integer, Double> e : this.SparseFeature) 
		{
			Fout.print(" " + e.getFirst() + ":" + e.getSecond());
			Fout.printf(" %d:%.4f", e.getFirst(),e.getSecond());
		}
	}
	
	public void outputIntegratedYear(PrintWriter Fout)
	{
		Fout.println(this.access + "\t" + this.integragedYear);
	}
	
	public void OutputRanklibFile(PrintWriter Fout,char space,int index)
	{
		HashSet<Integer> L2RCandidate  = new HashSet<Integer>();
		if (space == 'F') L2RCandidate = this.MFOL2RCandidate;
		if (space == 'P') L2RCandidate = this.BPOL2RCandidate;
		if (space == 'C') L2RCandidate = this.CCOL2RCandidate;
		double score = 0;
		for (Integer e:L2RCandidate)
		{
			if (this.Annotation.contains(e))
				Fout.print("1 qid:");
			else
				Fout.print("0 qid:");
			Fout.print(index);
			Fout.printf(" 1:%.4f", this.liblinearScore.get(e));
			
			if (this.blastKnnScore.containsKey(e)) score = this.blastKnnScore.get(e); else score = 0;
			Fout.printf(" 2:%.4f", score);
			
			if (this.blastScore.containsKey(e)) score = this.blastScore.get(e); else score = 0;
			Fout.printf(" 3:%.4f", score);
			
			Fout.printf(" 4:%.4f", proteinSet.NaiveIndex.get(e));
			Fout.println("# " + access + " ann = " + e);
		}
		
		
		
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
	
	public void removeHPOAnnotation(int... args) 
	{
		Set<Integer> SetAnn = new HashSet<Integer>();
		for (int node : args) 
		{
			SetAnn.add(node);
		}
		this.removeHPOAnnotation(SetAnn);
	}
	
	public void removeAllAnn()
	{
		this.Annotation.clear();
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
	
	public void removeHPOAnnotation(Set<Integer> args) 
	{
		for (int node : args) {
			HPOAnnotation.remove(node);
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

	public void recordliblinearScore()
	{
		for (Pair<Integer, Double> pair : this.PredictionScore)
		{
			int label = pair.getFirst();
			if ((this.MFOL2RCandidate.contains(label))
				|| (this.BPOL2RCandidate.contains(label))
				|| (this.CCOL2RCandidate.contains(label)))
			{
				this.liblinearScore.put(pair.getFirst(), pair.getSecond());
			}
		}
	}
	
	public void recordblastKnnScore()
	{
		for (Pair<Integer, Double> pair : this.PredictionScore)
		{
			int label = pair.getFirst();
			if ((this.MFOL2RCandidate.contains(label))
				|| (this.BPOL2RCandidate.contains(label))
				|| (this.CCOL2RCandidate.contains(label)))
			{
				this.blastKnnScore.put(pair.getFirst(), pair.getSecond());
			}
		}
	}
	
	public void recordblastScore()
	{
		for (Pair<Integer, Double> pair : this.PredictionScore)
		{
			int label = pair.getFirst();
			if ((this.MFOL2RCandidate.contains(label))
				|| (this.BPOL2RCandidate.contains(label))
				|| (this.CCOL2RCandidate.contains(label)))
			{
				this.blastScore.put(pair.getFirst(), pair.getSecond());
			}
		}
	}
	
	
	public void setAccess(String Access) {
		this.access = Access;
	}

	public void setBlastPred(proteinSet train) 
	{
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
	
	public void setHPOBlastPred(proteinSet train) 
	{
		ArrayList<Integer> ann = new ArrayList<Integer>();
		HashMap<Integer, Double> pred = new HashMap<Integer, Double>();
		for (Pair<String, Double> pair : blastResult) {
			double bitscore = pair.getSecond();
			String access = pair.getFirst();
			ann = train.getHPOAnnotation(access);
			for (int node : ann) {
				if (!pred.containsKey(node)) {
					pred.put(node, bitscore);
					HPOPredictionScore.add(new Pair<Integer, Double>(node, bitscore));
				}
			}
		}
	}

	public void setBlastKnn(proteinSet train) 
	{
		if (this.blastResult.size() == 0)
		{
			System.out.println(access + " have no blast result; We use Naive Score");
			this.PredictionScore = train.getNaiveList();
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
	
	public void setBlastKnnHPO(proteinSet train) 
	{
		if (this.blastResult.size() == 0)
		{
			System.out.println(access + " have no blast result; We use Naive Score");
			this.HPOPredictionScore = train.getNaiveList();
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
				ann = train.getHPOAnnotation(access);
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
			this.HPOPredictionScore
					.add(new Pair<Integer, Double>(entry.getKey(), (double) entry.getValue() / SumBitScore));
		}
	}
	
	public void knn(proteinSet train) 
	{
		ArrayList<Integer> ann = new ArrayList<Integer>();
		HashMap<Integer, Double> pred = new HashMap<Integer, Double>();
		double SumBitScore = 0;
		double sim;
		for (Pair<String, Double> pair : this.Similiar) 
		{
			sim = pair.getSecond();
			SumBitScore = SumBitScore + sim;
			String access = pair.getFirst();
			if (train.containProtein(access))
			{
				ann = train.getAnnotation(access);
				for (int node : ann) {
					if (!pred.containsKey(node)) {
						pred.put(node, sim);
					} else {
						pred.put(node, pred.get(node) + sim);
					}
				}
			}
		}
		for (Map.Entry<Integer, Double> entry : pred.entrySet()) 
		{
			double knnScore = 0.0;
			if (SumBitScore>Parameter.resolution)
				knnScore = (double) entry.getValue() / SumBitScore;
			
			this.PredictionScore.add(new Pair<Integer, Double>(entry.getKey(), knnScore));
		}
	}

	public void setLiblinearFeatureFromSparseFeature() 
	{
		List<Feature> x = new ArrayList<Feature>();
		for (Pair<Integer, Double> pair : this.SparseFeature) 
		{
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
	
	public void setHPOPredResult(ArrayList<Pair<Integer, Double>> PredList) 
	{
		for (Pair<Integer, Double> pair : PredList) 
		{
			this.HPOPredictionScore.add(pair);
		}
	}
	
	public void addPubMedID(int PubMedId)
	{
		this.PubMedID.add(PubMedId);
	}
	
	public void addIntAct(String access)
	{
		this.IntActProtein.add(access);
	}
	
	public int getIntegratedYear()
	{
        return this.integragedYear;
	}
	
	public void setIntegratedYear(int year)
	{
		this.integragedYear = year;
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
	
	
	public void OutputIDAnnSpecies(PrintWriter Fout)
	{
		Fout.print(this.access + " ");
		Fout.print(this.getSpecies() + " ");
		//this.OutputAnnotation(Fout);
		for (Integer k:this.Annotation)
		{
			Fout.print(proteinCommon.GOInt2Str(k) + " ");
		}
		
		Fout.println();
	}
}

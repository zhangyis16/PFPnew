package protein;

import java.io.*;
import java.lang.reflect.InvocationTargetException;
import java.lang.reflect.Method;
import java.util.*;

import Main.learningOfGO;
import common.*;
import demo.SetHandle;
import hpo.HPOSet;
import hpo.learningHPO;
import liblinear.Linear;
import multiLabel.evalution;

public class proteinSet
{
    public static HashMap<String,String> MapUniAccess2Name = new HashMap<String,String>();
    public static HashMap<String,String> MapName2UniAccess = new HashMap<String,String>();
    public static HashMap<String,String> MapAccess2UniAccess = new HashMap<String,String>();
    public static HashMap<Integer,Double> NaiveIndex = new HashMap<Integer,Double>();
    public static int FeatureSize;
    private ArrayList<oneProtein> Cell = new ArrayList<oneProtein>();  
    private HashMap<String,Integer> AccessIndex = new HashMap<String,Integer>();
    private ArrayList<Pair<Integer,Double>> predList = new ArrayList<Pair<Integer,Double>>();
    private ArrayList<Pair<Integer,Double>> hpoPredList = new ArrayList<Pair<Integer,Double>>();
    
    private HashMap<String, Object> attribute = new HashMap<String, Object>();
    private int MFOsize;
    private int BPOsize;
    private int CCOsize;
    

    public static void updateNaiveIndex(proteinSet train)
    {
    	ArrayList<Pair<Integer,Double>> predList = train.getNaiveList();
    	for (Pair<Integer,Double> e:predList)
    	{
    		NaiveIndex.put(e.getFirst(), e.getSecond());
    	}
    }
    
    
    public static void compareAnnDiff(proteinSet mySet,proteinSet benchSet,String OutFile) throws FileNotFoundException
    {
    	PrintWriter FoutmySet = new PrintWriter(new FileOutputStream(OutFile + "mySet"));
    	PrintWriter FoutbenchSet = new PrintWriter(new FileOutputStream(OutFile + "benchSet"));
    	for(int i=0;i<benchSet.size();i++)
    	{
    		oneProtein one = benchSet.getProtein(i);
    		oneProtein two = mySet.getProtein(i);
    		if (!one.compare(two))
    		{
    			one.OutputAnnotationList(FoutbenchSet,'\t');
    			FoutbenchSet.println();
    			two.OutputAnnotationList(FoutmySet,'\t');
    			FoutmySet.println();
    		}
    	}
    	FoutmySet.close();
    	FoutbenchSet.close();
    }
    public static void compareSpecies(proteinSet mySet,proteinSet benchSet,String OutFile) throws FileNotFoundException
    {
    	String Access,Species;
    	PrintWriter Fout = new PrintWriter(new FileOutputStream(OutFile));
    	int [] myresult = new int[27];
    	int [] bench = new int[27];
    	int [] overlap = new int[27];
    	int overSum = 0;
    	for (int i=0;i<mySet.size();i++)
    	{
    		Access = mySet.getProteinAccess(i);
    		Species = proteinCommon.getSpecies(Access);
    		int j;
    		for (j=0;j<proteinCommon.ArrayCAFA2Species_List.size();j++)
    			if (Species.equals(proteinCommon.ArrayCAFA2Species_List.get(j))) break;
    		myresult[j]++;
    		if (benchSet.containProtein(Access))
    		{
    			overlap[j]++;
    			overSum++;
    		}
    	}
    	for (int i=0;i<benchSet.size();i++)
    	{
    		Access = benchSet.getProteinAccess(i);
    		Species = proteinCommon.getSpecies(Access);
    		int j;
    		for (j=0;j<proteinCommon.ArrayCAFA2Species_List.size();j++)
    			if (Species.equals(proteinCommon.ArrayCAFA2Species_List.get(j))) break;
    		bench[j]++;
    	}
    	
    	Fout.println(",myresult,bench,overlap");
    	for (int i=0;i<27;i++)
    	{
    		Fout.println(proteinCommon.ArrayCAFA2Species_List.get(i)+","+myresult[i]+","+
    				bench[i]+","+overlap[i]);
    	}
    	Fout.println("," + mySet.size() + ","+benchSet.size() + ","+overSum);
    	Fout.close();
    }
    
    public static proteinSet getLimitKnowProtein(proteinSet EndSet,proteinSet StartSet,char space)
    {
    	proteinSet result = new proteinSet();
    	for(int i=0;i<StartSet.size();i++)
    	{
    		String Name = StartSet.getProteinAccess(i);
    		oneProtein one = StartSet.getProtein(i);
    		if (one.getSubAnnotation(space).size() == 0)
    		{
    			if (EndSet.containProtein(Name))
    			{
    				oneProtein other = EndSet.getProtein(Name);
    				if (other.getSubAnnotation(space).size() > 0)
    					result.addProtein(new oneProtein(Name,other.getSubAnnotation(space)));
    			}
    		}
    	}
    	return result;
    }
    public static proteinSet getNewProtein(proteinSet StartSet,proteinSet EndSet)
    {
    	proteinSet result = new proteinSet();
    	for(int i = 0;i < EndSet.size(); i++)
    	{
    		String Access = EndSet.getProteinAccess(i);
    		oneProtein one = EndSet.getProtein(i);
    		if (!StartSet.containProtein(Access))
    		{
    			result.addProtein(one);			
    		}
    	}
    	return result;
    }
    
    public static void compare2Set(proteinSet mySet,proteinSet benchSet,String Directory,GoSet aGoSet) throws FileNotFoundException
    {
    	compareSpecies(mySet, benchSet,Directory + "Resulttype1.csv");
		listCHABIEprotein(mySet, benchSet,Directory + "chabietype1.csv");
		
		mySet.removeRedunAnn(learningOfGO.aGoSet);
		benchSet.removeRedunAnn(learningOfGO.aGoSet);
		
		mySet.getIntersection(benchSet);
		benchSet.getIntersection(mySet);
		mySet.SortProtein();
		benchSet.SortProtein();
		proteinSet.compareAnnDiff(mySet, benchSet, "CAFA2");
    }
    
    public static proteinSet getCAFA2Protein(String t0 , String t1 ,GoSet aGoSet) throws FileNotFoundException
    {
    	proteinSet proteinSetT0 = new proteinSet();
		proteinSet proteinSetT1 = new proteinSet();
		proteinSet AnnType1 = new proteinSet();
		
		proteinSet.LoadSwissMapAccess2UniAccess("../InFile/Swiss/ac2ac" + t0);
		proteinSet.LoadAccess2NameMap("../InFile/Swiss/ac2Name" + t0);
		
		proteinSetT0.AddAnnotationInSwiss("../InFile/Swiss/Ann" + t0);
		proteinSetT0.AddAnnotationInSwiss("../InFile/Goa/Ann" + t0);
		proteinSetT0.AddAnnotationInSwiss("../InFile/GODB/Ann" + t0);
		proteinSetT0.filterProteinOnly5515();
		proteinSetT0.removeGoNotIn(aGoSet);
		
		
		proteinSetT1.AddAnnotationInSwiss("../InFile/Swiss/Ann" + t1);
		proteinSetT1.AddAnnotationInSwiss("../InFile/Goa/Ann" + t1);
		proteinSetT1.AddAnnotationInSwiss("../InFile/GODB/Ann" + t1);
		proteinSetT1.filterProteinOnly5515();
		proteinSetT1.removeGoNotIn(aGoSet);
		
		AnnType1 = proteinSet.getNewProtein(proteinSetT1,proteinSetT0);
		AnnType1.removeAnnMF_Only5515();
		
		return AnnType1;
    }
    
    public static proteinSet getNewProtein(String Tstart , String Tend) throws FileNotFoundException
    {
    	proteinSet proteinSetStart = new proteinSet();
		proteinSet proteinSetEnd = new proteinSet();
		proteinSet AnnType1 = new proteinSet();
		
		proteinSet.LoadSwissMapAccess2UniAccess("../InFile/Swiss/ac2ac" + Tend);
		proteinSet.LoadAccess2NameMap("../InFile/Swiss/ac2Name" + Tend);
		
		proteinSetStart.AddAnnotation("../InFile/Swiss/Ann" + Tstart);
		proteinSetStart.AddAnnotationInSwiss("../InFile/Goa/Ann" + Tstart);
		proteinSetStart.AddAnnotationInSwiss("../InFile/GODB/Ann" + Tstart);
		proteinSetStart.filterProteinOnly5515();
		System.out.println("proteinSetStart.size = " + proteinSetStart.size());
		
		
		proteinSetEnd.AddAnnotation("../InFile/Swiss/Ann" + Tend);
		proteinSetEnd.AddAnnotationInSwiss("../InFile/Goa/Ann" + Tend);
		proteinSetEnd.AddAnnotationInSwiss("../InFile/GODB/Ann" + Tend);
		proteinSetEnd.filterProteinOnly5515();
		System.out.println("proteinSetEnd.size = " + proteinSetEnd.size());
		
		proteinSetStart.mapAccess();
		AnnType1 = proteinSet.getNewProtein(proteinSetStart,proteinSetEnd);

		return AnnType1;
    }
    
    
    
    public static void listCHABIEprotein(proteinSet mySet,proteinSet benchSet,String OutFile) throws FileNotFoundException
    {
    	PrintWriter Fout = new PrintWriter(new FileOutputStream(OutFile));
    	Fout.println("我多的");
    	for (int i=0;i<mySet.size();i++)
    	{
    		String Name = mySet.getProteinAccess(i);
    		if (!benchSet.containProtein(Name))	Fout.println(mySet.getProteinName(i));
    	}
    	Fout.println("它多的");
    	for (int i=0;i<benchSet.size();i++)
    	{
    		String Name = benchSet.getProteinAccess(i);
    		if (!mySet.containProtein(Name))	Fout.println(Name);
    	}
    	Fout.close();
    }
    public static void LoadAccess2NameMap(String InFile) throws FileNotFoundException
    {
    	Scanner In = new Scanner(new FileInputStream(InFile));
    	MapUniAccess2Name.clear();
    	MapName2UniAccess.clear();
    	String Access,ProteinName;
    	while (In.hasNext())
    	{
    		Access = In.next();
    		ProteinName = In.next();
    		MapUniAccess2Name.put(Access, ProteinName);
    		MapName2UniAccess.put(ProteinName, Access);
    	}
    	In.close();
    }  
    public static void LoadSwissMapAccess2UniAccess(String InFile) throws FileNotFoundException
    {
    	Load.LoadMap(MapAccess2UniAccess, InFile);
    }
    
    public static void LoadSwissAccessMap(String InFile) throws FileNotFoundException
    {
    	System.out.println("Read Swiss" + InFile + " Now");
    	Scanner In = new Scanner(new FileInputStream(InFile));
    	String line = new String();
    	
    	String ProteinName = new String();
    	String Access = new String();
    	String[] strarr = new String[3];
    	int countPro = 0;
    	while (In.hasNext())
    	{
    		line = In.nextLine();
    		strarr = line.split("\\s+",3);
    		if (strarr[0].equals("ID")) 
    		{
    			ProteinName = strarr[1];
    			countPro++;
    			if (countPro%1000 == 0)
    			{
    				System.out.println(ProteinName);
    				System.out.println("We have read" + countPro + "protein");
    			}
    			line = In.nextLine();
    			strarr = line.split("\\s+",3);
    			Access = strarr[1].substring(0, strarr[1].length()-1);
    			MapUniAccess2Name.put(Access, ProteinName);
    			MapAccess2UniAccess.put(Access, Access);
    			if (strarr.length >2)
    			{
    				StringTokenizer st = new StringTokenizer(strarr[2]," ");
    				int m = st.countTokens();
    				for (int count = 0;count < m;count++)
    				{
    					String Acc2 = st.nextToken();
    					Acc2.substring(0, Acc2.length()-1);
    					Acc2 = StringHandle.removeChar(Acc2, ';');
    					MapAccess2UniAccess.put(Acc2, Access);
    				}
    			}
    		}
    		if (strarr[0].equals("AC"))
    		{
    			StringTokenizer st = new StringTokenizer(strarr[1]," ");
    			int m = st.countTokens();
    			for (int count = 0;count < m;count++)
    			{
    				String Acc2 = st.nextToken();
    				Acc2.substring(0, Acc2.length()-1);
    				Acc2 = StringHandle.removeChar(Acc2, ';');
    				MapAccess2UniAccess.put(Acc2, Access);
    			}
    		}
    	}
    	In.close();
    	System.out.println("Read Swiss" + InFile + " Finish");
    }
    
    public static void OutputAccess2NameMap(String OutFile) throws FileNotFoundException
    {
    	PrintWriter Fout = new PrintWriter(new FileOutputStream(OutFile));
    	for (Map.Entry<String, String> entry : MapUniAccess2Name.entrySet()) 
    	{
    		Fout.println(entry.getKey() + "\t" + entry.getValue());
    	}
    	Fout.close();
    }
    public static void OutputAccess2AccessMap(String OutFile) throws FileNotFoundException
    {
    	PrintWriter Fout = new PrintWriter(new FileOutputStream(OutFile));
    	for (Map.Entry<String, String> entry : MapAccess2UniAccess.entrySet()) 
    	{
    		Fout.println(entry.getKey() + "\t" + entry.getValue());
    	}
    	Fout.close();
    }
     
    
    public static void compareHomoProtein(proteinSet proSp1 ,proteinSet proSp2 ,GoSet aGoSet) throws IOException, InterruptedException
    {
    	PrintWriter Fout = new PrintWriter(new FileOutputStream("differentSpecies"));

    	for (int i=0;i<proSp1.size();i++)
    	{
    		String name1 = proSp1.getProteinName(i);
    		String name2 = proSp2.getProteinName(i);
    		int loc = name1.lastIndexOf('_');
    		String name = name1.substring(0, loc);
    		String seq1 = proSp1.getSequence(i);
    		String seq2 = proSp2.getSequence(i);
    		double sim = proteinCommon.blastTwoProteinSeq(seq1, seq2);
    		HashSet<Integer> Set1 =  proSp1.getProtein(i).getSubAnnotation('F');
    		HashSet<Integer> Set2 =  proSp2.getProtein(i).getSubAnnotation('F');
    		double MFOscore = SetHandle.interRatioUnion(Set1, Set2);
    		if ((Set1.size() == 0) || (Set2.size() == 0)) MFOscore = 0.0/0.0;
    		
    		Set1 =  proSp1.getProtein(i).getSubAnnotation('P');
    		Set2 =  proSp2.getProtein(i).getSubAnnotation('P');
    		double BPOscore = SetHandle.interRatioUnion(Set1, Set2);
    		if ((Set1.size() == 0) || (Set2.size() == 0)) BPOscore = 0.0/0.0;
    		
    		Set1 =  proSp1.getProtein(i).getSubAnnotation('C');
    		Set2 =  proSp2.getProtein(i).getSubAnnotation('C');
    		double CCOscore = SetHandle.interRatioUnion(Set1, Set2);
    		if ((Set1.size() == 0) || (Set2.size() == 0)) CCOscore = 0.0/0.0;
    		
    		Fout.printf("%s , %.4f , %.4f , %.4f , %.4f\n", name,sim,MFOscore,BPOscore,CCOscore);
    		System.out.printf("%s , %.4f , %.4f , %.4f , %.4f\n", name,sim,MFOscore,BPOscore,CCOscore);
    	}
    	Fout.close();
    }
    
    public static void parseUniprotDatFile(String data) throws FileNotFoundException
    {
    	String swissFileName = "../InFile/UniprotKB/uniprot_trembl" + data +".dat";
    	
    	Scanner In = new Scanner(new FileInputStream(swissFileName));
    	PrintWriter FoutAC2Name = new PrintWriter("../InFile/UniprotKB/ac2Name" + data);
    	PrintWriter FoutAC2AC =   new PrintWriter("../InFile/UniprotKB/ac2ac" + data);
    	PrintWriter FoutAnn =   new PrintWriter("../InFile/UniprotKB/Ann" + data);
    	PrintWriter FoutSeq =   new PrintWriter("../InFile/UniprotKB/Seq" + data);
    	System.out.println("Read Swiss" + swissFileName + " Now");
    	
    	String line = new String();
    	String ProteinName = new String();
    	String Access = new String();
    	String[] strarr = new String[3];
    	int countPro = 0;
    	while (In.hasNext())
    	{
    		line = In.nextLine();
    		strarr = line.split("\\s+",3);
    		if (strarr[0].equals("ID")) 
    		{
    			ProteinName = strarr[1];
    			countPro++;
    			if (countPro%1000 == 0)
    			{
    				System.out.println(ProteinName);
    				System.out.println("We have read" + countPro + "protein");
    			}
    			line = In.nextLine();
    			strarr = line.split("\\s+",3);
    			Access = strarr[1].substring(0, strarr[1].length()-1);
    			FoutAC2Name.println(Access + '\t' + ProteinName);
    			FoutAC2AC.println(Access + '\t' + Access);
    			if (strarr.length >2)
    			{
    				StringTokenizer st = new StringTokenizer(strarr[2]," ");
    				int m = st.countTokens();
    				for (int count = 0;count < m;count++)
    				{
    					String Acc2 = st.nextToken();
    					Acc2.substring(0, Acc2.length()-1);
    					FoutAC2AC.println(Acc2 + '\t' + Access);
    				}
    			}
    		}
    		if (strarr[0].equals("AC"))
    		{
    			StringTokenizer st = new StringTokenizer(strarr[1]," ");
    			int m = st.countTokens();
    			for (int count = 0;count < m;count++)
    			{
    				String Acc2 = st.nextToken();
    				Acc2.substring(0, Acc2.length()-1);
    				FoutAC2AC.println(Acc2 + '\t' + Access);
    			}
    		}
    		if (strarr[0].equals("DR") && strarr[1].equals("GO;"))
    		{
    			String[] arr = strarr[2].split("; ",3);
    			int gonum = proteinCommon.GOStr2Int(arr[0]);
    			String evidence = arr[2].substring(0,arr[2].indexOf(':'));
    			if (proteinCommon.EvidenceCode.contains(evidence)) 
    				FoutAnn.println(Access + '\t' + proteinCommon.GOInt2Str(gonum));
    		}
    		if (strarr[0].equals("SQ")) 
    		{
    			String str = new String();
    			String sequence = "";
    			do
    			{
    				sequence = sequence + str;
    				str = In.next();
    			}while (!str.equals("//"));
    			FoutSeq.println(Access + '\t' + sequence);
    		}
    	}
    	In.close();
    	System.out.println("Read Swiss" + swissFileName + " Finish");
    	FoutAC2Name.close();
    	FoutAC2AC.close();
    	FoutAnn.close();
    	FoutSeq.close();
    }
    
    


    
    public proteinSet() 
    {
        this.Cell = new ArrayList<oneProtein>();
        this.AccessIndex = new HashMap<String,Integer>();
    }
    public proteinSet(String InFile) throws FileNotFoundException
    {
    	AddAnnotation(InFile);
    }
    
    
    public void mapAccess()
    {
    	for (int i = 0;i<Cell.size();i++)
    	{
    		String acc = Cell.get(i).getAccess();
    		if (proteinSet.MapAccess2UniAccess.containsKey(acc))
    		{
    			String newacc = proteinSet.MapAccess2UniAccess.get(acc);
    			Cell.get(i).setAccess(newacc);
    		}
    	}
    }
    
    public int getSubSetSize(char sp)
    {
    	if (sp == 'F') return this.MFOsize;
    	if (sp == 'P') return this.BPOsize;
    	if (sp == 'C') return this.CCOsize;
    	return 0;
    }
    public oneProtein get(int index)
    {
    	return Cell.get(index);
    }
    public int sizeOfAnn()
    {
    	int Sum = 0;
    	for (oneProtein e:Cell)
    	{
    		Sum += e.getAnnotationSize();
    	}
    	return Sum;
    }
    
    public int sizeOfHPOAnn()
    {
    	int Sum = 0;
    	for (oneProtein e:Cell)
    	{
    		Sum += e.getHPOAnnotationSize();
    	}
    	return Sum;
    }
    public void AddAnnotation(String InFile) throws FileNotFoundException
    {
    	Scanner In = new Scanner(new FileInputStream(InFile));
    	String access,aAnnotation;
    	while (In.hasNext())
    	{
    		access = In.next();
    		if (access.contains("_")) access = proteinSet.MapName2UniAccess.get(access);
    		aAnnotation = In.next();
    		AddAnnotation(access,proteinCommon.GOStr2Int(aAnnotation));
    	}
    	In.close();
    }
    
    public void AddAnnotationInProteinSet(String InFile) throws FileNotFoundException
    {
    	Scanner In = new Scanner(new FileInputStream(InFile));
    	String access,aAnnotation;
    	while (In.hasNext())
    	{
    		access = In.next();
    		if (access.contains("_")) access = proteinSet.MapName2UniAccess.get(access);
    		aAnnotation = In.next();
    		AddAnnotationInProteinSet(access,proteinCommon.GOStr2Int(aAnnotation));
    	}
    	In.close();
    }
    
    public void AddAnnotationInProteinSet(String Access,int Annotation)
    {
    		if (AccessIndex.containsKey(Access)) 
    		{
    			int index = AccessIndex.get(Access);
    			Cell.get(index).AddAnnotation(Annotation);
    		}
    }
    
    
    public void AddIEA_Annotation(String InFile) throws FileNotFoundException
    {
    	Scanner In = new Scanner(new FileInputStream(InFile));
    	String access,aAnnotation;
    	while (In.hasNext())
    	{
    		access = In.next();
    		if (access.contains("_")) access = proteinSet.MapName2UniAccess.get(access);
    		aAnnotation = In.next();
    		AddIEA_Annotation(access,proteinCommon.GOStr2Int(aAnnotation));
    	}
    	In.close();
    }
    public void AddIEA_Annotation(String Access,int Annotation)
    {
    		if (AccessIndex.containsKey(Access)&&(learningOfGO.aGoSet.containNode(Annotation))) 
    		{
    			int index = AccessIndex.get(Access);
    			Cell.get(index).AddIEA_Annotation(Annotation);
    		}
    }
    
    public void AddHPOAnnotation(String InFile) throws FileNotFoundException
    {
    	Scanner In = new Scanner(new FileInputStream(InFile));
    	String access,aAnnotation;
    	while (In.hasNext())
    	{
    		access = In.next();
    		if (access.contains("_")) access = proteinSet.MapName2UniAccess.get(access);
    		aAnnotation = In.next();
    		AddHPOAnnotation(access,proteinCommon.HPStr2Int(aAnnotation));
    	}
    	In.close();
    }
    
    public void AddHPOAnnotation(String Access,int Annotation)
    {
		if (AccessIndex.containsKey(Access)) 
		{
			int index = AccessIndex.get(Access);
			Cell.get(index).AddHPOAnnotation(Annotation);
		}else
		{
			Cell.add(new oneProtein(Access));
			AccessIndex.put(Access, Cell.size()-1);
			int index = AccessIndex.get(Access);
			Cell.get(index).AddHPOAnnotation(Annotation);
		}
    }
    
    
    public void AddAnnotation(String Access,int Annotation)
    {
    		if (AccessIndex.containsKey(Access)) 
    		{
    			int index = AccessIndex.get(Access);
    			Cell.get(index).AddAnnotation(Annotation);
    		}else
    		{
    			Cell.add(new oneProtein(Access));
    			AccessIndex.put(Access, Cell.size()-1);
    			int index = AccessIndex.get(Access);
    			Cell.get(index).AddAnnotation(Annotation);
    		}
    }
    public void AddAnnotationInGOset(String InFile) throws FileNotFoundException
    {
    	Scanner In = new Scanner(new FileInputStream(InFile));
    	String access,aAnnotation;
    	while (In.hasNext())
    	{
    		access = In.next();
    		if (access.contains("_")) access = proteinSet.MapName2UniAccess.get(access);
    		aAnnotation = In.next();
    		AddAnnotationInGOset(access,proteinCommon.GOStr2Int(aAnnotation));
    	}
    	In.close();
    }
    public void AddAnnotationInGOset(String Access,int Annotation)
    {
    	if (learningOfGO.aGoSet.containNode(Annotation))
    	{
    		AddAnnotation(Access,Annotation);
    	}else
    	{
    		
    	}
    }
    public void AddAnnotationInSwiss(String InFile) throws FileNotFoundException
    {
    	Scanner In = new Scanner(new FileInputStream(InFile));
    	String access,aAnnotation;
    	while (In.hasNext())
    	{
    		access = In.next();
    		if (access.contains("_")) access = proteinSet.MapName2UniAccess.get(access);
    		aAnnotation = In.next();
    		if (proteinSet.MapAccess2UniAccess.containsKey(access))
    		{
    			access = proteinSet.MapAccess2UniAccess.get(access);
    			AddAnnotation(access,proteinCommon.GOStr2Int(aAnnotation));
    		}
    	}
    	In.close();
    }
    public void addNoScoreLabelRandom()
    {
    	double miniscore = this.miniscore();
    	System.out.println("miniscore = " + miniscore);
    	for (oneProtein e:Cell)
    	{
    		e.addNoScoreLabelRandom(this.predList,miniscore);
    	}
    }
    public void LoadNoSparseFeatureVector(String InFile,int FeatureNum) throws FileNotFoundException
    {
    	Scanner In = new Scanner(new FileInputStream(InFile));
    	ArrayList<Double> Fea = new ArrayList<Double>();
    	while (In.hasNext())
    	{
    		String Acc = In.nextLine();
    		String Feature = In.nextLine();
    		String[] Featurearr = Feature.split(" ",FeatureNum);
    		int index = this.AccessIndex.get(Acc);
    		Fea.clear();
    		for (int i=0;i<FeatureNum;i++)
    		{
    			double num = Double.parseDouble(Featurearr[i]);
    			Fea.add(num);
    		}
    		this.Cell.get(index).SetNoSparseFeature(Fea);	
    	}
    	In.close();
    }
    public void addGOFather(GoSet aGoSet)
    {
    	for(oneProtein e:Cell)
    		e.addGOFather(aGoSet);
    }
    
    public void addHPOFather(HPOSet aHPOSet)
    {
    	for(oneProtein e:Cell)
    		e.addHPOFather(aHPOSet);
    }
    
    public void addIEAFather(GoSet aGoSet)
    {
    	for(oneProtein e:Cell)
    		e.addIEAFather(aGoSet);
    }
    

    
    public void addGoaAnnotation(String InFile) throws FileNotFoundException
    {
    	Scanner In = new Scanner(new FileInputStream(InFile));
    	String line,access,qualifer;
    	int gonum;
    	while (In.hasNext())
    	{
    		line = In.nextLine();
    		String[] strarr = line.split("\t",6);
    		access = strarr[1];
    		qualifer = strarr[3];
    		gonum = proteinCommon.GOStr2Int(strarr[4]);
    		if (strarr[0].equals("UniProtKB"))
    			if (!qualifer.equals("NOT")) AddAnnotation(access,gonum);
    	}
    	In.close();
    }
    
    
    public void addProtein(String access)
    {
    	this.addProtein(new oneProtein(access));
    }
    public void addProtein(oneProtein one)
    {
    	Cell.add(one);
    }
    public void addProteinFromFile(String InFile) throws FileNotFoundException
    {
    	Scanner In = new Scanner(new FileInputStream(InFile));
    	String Access;
    	while (In.hasNext())
    	{
    		Access = In.next();
    		if (Access.contains("_")) 
    		{
    			if (proteinSet.MapName2UniAccess.containsKey(Access))
    			{
    				Access = proteinSet.MapName2UniAccess.get(Access);
    				Cell.add(new oneProtein(Access));
    			}
    		}
    		else
    		Cell.add(new oneProtein(Access));
    	}
    	In.close();
    	UpdateIndex();
    }
    public void addSequence(String Access, String sequence)
    {
    	if (this.containProtein(Access))
    	{
    		int index = AccessIndex.get(Access);
    		Cell.get(index).addSequence(sequence);
    	}
    }
    public void loadSwissAnnotation(String InFile) throws FileNotFoundException 
    {
    	System.out.println("Read Swiss" + InFile + " Now");
    	Scanner In = new Scanner(new FileInputStream(InFile));
    	String line,ProteinName;
    	String Access = new String();
    	String[] strarr = new String[3];
    	int countPro = 0;
    	while (In.hasNext())
    	{
    		line = In.nextLine();
    		strarr = line.split("\\s+",3);
    		if (strarr[0].equals("ID")) 
    		{
    			ProteinName = strarr[1];
    			countPro++;
    			if (countPro%1000 == 0)
    			{
    				System.out.println(ProteinName);
    				System.out.println("We have read" + countPro + "protein");
    			}
    			line = In.nextLine();
    			strarr = line.split("\\s+",3);
    			Access = strarr[1].substring(0, strarr[1].length()-1);
    		}
    		if (strarr[0].equals("DR") && strarr[1].equals("GO;"))
    		{
    			String[] arr = strarr[2].split("; ",3);
    			int gonum = proteinCommon.GOStr2Int(arr[0]);
    			String evidence = arr[2].substring(0,arr[2].indexOf(':'));
    			if (proteinCommon.EvidenceCode.contains(evidence)) 
    				AddAnnotation(Access,gonum);
    		}
    	}
    	In.close();
    	System.out.println("Read Swiss" + InFile + " Finish");
    }
    public void loadSwissProteinSequence(String InFile) throws FileNotFoundException
    {
    	Cell.clear();
    	System.out.println("Read Swiss" + InFile + " Now");
    	Scanner In = new Scanner(new FileInputStream(InFile));
    	String line = new String();
    	String ProteinName = new String();
    	String Access = new String();
    	String sequence = new String();
    	String str = "";
    	int countPro = 0;
    	while (In.hasNext())
    	{
    		line = In.nextLine();
    		String[] strarr = line.split("\\s+",3);
    		if (strarr[0].equals("ID")) 
    		{
    			ProteinName = strarr[1];
    			countPro++;
    			if (countPro%1000 == 0)
    			{
    				System.out.println(ProteinName);
    				System.out.println("We have read" + countPro + "protein");
    			}
    			line = In.nextLine();
    			strarr = line.split("\\s+",3);
    			Access = strarr[1].substring(0, strarr[1].length()-1);
    		}
    		if (strarr[0].equals("SQ")) 
    		{
    			str = new String();
    			sequence = "";
    			do
    			{
    				sequence = sequence + str;
    				str = In.next();
    			}while (!str.equals("//"));
    			Cell.add(new oneProtein(Access,sequence));
    		}
    	}
    	In.close();
    	System.out.println("Read Swiss" + InFile + " Finish");
    }
    public void addBlastResultBitScore(String InFile) throws FileNotFoundException
    {
    	Scanner In = new Scanner(new FileInputStream(InFile));
    	String line,Access,Access2;
    	double maxbitscore = 0;
    	while (In.hasNext())
    	{
    		line = In.nextLine();
    		if (line.length()>0)
    		if (line.charAt(0) != '#')
    		{
    			String[] strarr = line.split("\\s+",12);
    			Access = strarr[0];
    			if (strarr.length > 11)
    			{
    				Access2 = strarr[1];
    				double bitscore = Double.valueOf(strarr[11]);
    				if (AccessIndex.containsKey(Access))
    				{
    					int index = AccessIndex.get(Access);
    					Cell.get(index).addBlastResult(Access2, bitscore);
    					if (bitscore > maxbitscore && !Access2.equals(Access)) maxbitscore = bitscore;
    				}
    			}
    		}
    	}
    	In.close();
    }
    public void addSimiliar(String InFile) throws FileNotFoundException
    {
    	Scanner In = new Scanner(new FileInputStream(InFile));
    	String line,Access,Access2;
    	while (In.hasNext())
    	{
    		line = In.nextLine();
    		String[] strarr = line.split(" ");
    		Access = strarr[0];
    		Access2 = strarr[1];
    		double sim = Double.parseDouble(strarr[2]);
    		if (this.AccessIndex.containsKey(Access))
    		{
    			int index = AccessIndex.get(Access);
    			Cell.get(index).addSimiliar(Access2,sim);
    		}
    	}
    	In.close();
    }
    public void addBlastResultIdentity(String InFile) throws FileNotFoundException
    {
    	Scanner In = new Scanner(new FileInputStream(InFile));
    	String line,Access,Access2;
    	double maxIdentity = 0;
    	while (In.hasNext())
    	{
    		line = In.nextLine();
    		String[] strarr = line.split("\\s+",12);
    		Access = strarr[0];
    		Access2 = strarr[1];
    		double Identity = Double.valueOf(strarr[2]);
    		int index = AccessIndex.get(Access);
    		Cell.get(index).addBlastResult(Access2, Identity);
    		if (Identity > maxIdentity && !Access2.equals(Access)) maxIdentity = Identity;
    	}
    	In.close();
    }
    public void clear()
    {
    	Cell.clear();
    	AccessIndex.clear();
    }
    public void calMFBPCCsize(GoSet aGoSet)
    {
    	this.MFOsize = 0;
		this.BPOsize = 0;
		this.CCOsize = 0;
    	for (oneProtein e:this.Cell)
    	{
    		e.calMFBPCCsize(aGoSet);
			if (e.getMFOsize() > 0) this.MFOsize++;
			if (e.getBPOsize() > 0) this.BPOsize++;
			if (e.getCCOsize() > 0) this.CCOsize++;
    	}
    	System.out.println("Total size = " + this.size());
    	
		proteinSet MFSub =	this.getSubAnnotation('F');
		System.out.println("MFO  size = " + this.MFOsize + "    MFO Ann : " + MFSub.sizeOfAnn());
		proteinSet BPSub =	this.getSubAnnotation('P');
		System.out.println("BPO  size = " + this.BPOsize + "    BPO Ann : " + BPSub.sizeOfAnn());
		proteinSet CCSub =	this.getSubAnnotation('C');
		System.out.println("CCO  size = " + this.CCOsize + "    CCO Ann : " + CCSub.sizeOfAnn());
    }
    public void readMSAlignResult(String InFile) throws FileNotFoundException  //Read Multi Sequence Alignment Result
    {
    	String line;
    	Scanner In = new Scanner(new FileInputStream(InFile));
    	while (In.hasNext())
    	{
    		line = In.nextLine();
    		
    		int index;
    		if (line.length()>6)
    		if (line.substring(0,6).equals("Query_"))
    		{
    			String[] strarr = line.split("\\s+",4);
    			String oriSeq = strarr[2];
    			index = Integer.parseInt(strarr[0].substring(6)) - 1;
    			
    			System.out.println(index);
    			
    			ArrayList<Pair<String, String>> AccSeq = new ArrayList<Pair<String, String>>();
    			
    			line = In.nextLine();
    			while(!line.equals(""))
    			{
    				String[] strarr2 = line.split("\\s+",2);
    				String acc = strarr2[0];
    				
    				//PrintWriter Fout = new PrintWriter(new FileOutputStream("error.out"));
    				//Fout.println(line);
    				//Fout.close();
    				
    				line = In.nextLine();
    				if (strarr2[1].length()<oriSeq.length()+6) continue;
    				String seq = strarr2[1].substring(6, 6 + oriSeq.length());
    				AccSeq.add(new Pair<String,String>(acc,seq));
    			}
    			for(int i = 0;i < oriSeq.length();i++)
    			{
    				if (Character.isUpperCase(oriSeq.charAt(i)))
    				{
    					Cell.get(index).addOriAlignSequenceChar(oriSeq.charAt(i));
    					for (int j = 0;j < AccSeq.size(); j++)
    					{
    	    				String seq = AccSeq.get(j).getSecond();
    	    				char ch = seq.charAt(i);
    	    				String Acc = AccSeq.get(j).getFirst();
    	    				
    	    				if (Character.isUpperCase(ch))
    	    					Cell.get(index).addMSASequenceChar(ch, Acc);
    	    				else
    	    					Cell.get(index).addMSASequenceChar('-', Acc);
    					}
    				}
    			}
    			Cell.get(index).addMSAgap();
    			//while(true);
    		}
    	}
    	In.close();
    	
    }
    
    public boolean containProtein(String Access)
    {
    	return AccessIndex.containsKey(Access);
    }
    public void eraserAnnotation(int gonum)
    {
    	for(oneProtein e:Cell)
    		e.eraserAnnotation(gonum);
    }
    
    public void filterProteinOnly5515()
    {
    	for(int i=0;i<Cell.size();i++)
    	{
    		while (i<Cell.size() && (Cell.get(i).onlyGO5515()) )
    		{
    			Cell.remove(i);
    		}
    	}
    	UpdateIndex();
    }
        
    public void evalutionIEAPR(char space)
    {
    	int prenum = 0;
    	int recallnum = 0;
    	double presum = 0;
    	double recallsum = 0;
    	for (oneProtein e:Cell)
    	{
    		Pair<Double,Double> result = e.evalutionPR(space);
    		double pre = result.getFirst();
    		double recall = result.getSecond();
    		if (recall>0) {recallnum++; recallsum+=pre;}
    		if (pre>=0) {prenum++; presum+=pre;}
    	}
    	System.out.println(String.format("IEA Predictio number = %d", prenum));
    	System.out.println(String.format("IEA recall number = %d", recallnum));
    	System.out.println(String.format("precsion = %.2f,recall = %.2f ", presum/prenum,recallsum/recallnum));
    }
    
    public void evalutionEveryProtein(String OutFile) throws FileNotFoundException 
    {
    	PrintWriter Fout = new PrintWriter(new FileOutputStream(OutFile));
    	Fout.println("Access,Species,BlastHit,MFO_fmax,MFO_AUPR,BPO_fmax,BPO_AUPR,CCO_fmax,CCO_AUPR");
    	for(oneProtein e:Cell)
    	{
    		
    		proteinSet newSet = new proteinSet();
    		newSet.addProtein(e);
    		Pair<Double,Double> pair = new Pair<Double,Double>();
    		Fout.printf("%s , %s , %d,", e.getAccess(),e.getSpecies(),e.getBlastNum());
    		pair = newSet.getFmaxAndAUPR("MFO");
    		Fout.printf("%.4f , %.4f,", pair.getFirst(),pair.getSecond());
    		pair = newSet.getFmaxAndAUPR("BPO");
    		Fout.printf("%.4f , %.4f,", pair.getFirst(),pair.getSecond());
    		pair = newSet.getFmaxAndAUPR("CCO");
    		Fout.printf("%.4f , %.4f \n", pair.getFirst(),pair.getSecond());
    		
    	}
    	Fout.close();
    }
    public void evalutionEveryProtein(String OutFile,String Space) throws FileNotFoundException
    {
    	PrintWriter Fout = new PrintWriter(new FileOutputStream(OutFile));
    	Fout.println("Access,Species,BlastHit,answer_size,fmax,AUPR");
    	for(oneProtein e:Cell)
    	{
    		
    		proteinSet newSet = new proteinSet();
    		newSet.addProtein(e);
    		Pair<Double,Double> pair = new Pair<Double,Double>();
    		if (e.getAnnotationSize(Space)>0)
    		{
    			Fout.printf("%s , %s , %d, %d ,", e.getAccess(),e.getSpecies(),e.getBlastNum(),e.getAnnotationSize());
    			pair = newSet.getFmaxAndAUPR(Space);
    			Fout.printf("%.4f , %.4f \n", pair.getFirst(),pair.getSecond());
    		}
    		
    	}
    	Fout.close();
    }
    public Pair<Double,Double> getFmaxAndAUPR(String Space) throws FileNotFoundException 
    {
    	ArrayList<HashSet<Integer>> answer = new ArrayList<HashSet<Integer>>();
    	ArrayList<ArrayList<Pair<Integer,Double>>> predictor = new ArrayList<ArrayList<Pair<Integer,Double>>>();
    	proteinSet subSet = new proteinSet();
    	if (Space.equals("MFO"))  subSet = this.getSubAnnotation('F');
    	if (Space.equals("BPO"))  subSet = this.getSubAnnotation('P');
    	if (Space.equals("CCO"))  subSet = this.getSubAnnotation('C');
    	ArrayList<Pair<Double,Double>> ReCallPrePair = new ArrayList<Pair<Double,Double>>();
    	subSet.getAnswerPredictor(answer,predictor);
    	Pair<Double,Double> pair = evalution.GetFMeasureMax(answer, predictor,ReCallPrePair);
    	
    	double aupr = evalution.CalculateAUC(ReCallPrePair);
    	
    	if (this.Cell.get(0).getAccess().equals("P0A6J3"))
    	{
    		Output.OutputPair("MFOAUPRlist",ReCallPrePair,"# Protein  " + "P0A6J3");
    	}
    	return new Pair<Double,Double>(pair.getFirst(),aupr);
    	
    }
    
    public void evalEveryHPOLabelAUC(String OutFile) throws FileNotFoundException
    {
    	PrintWriter Fout = new PrintWriter(new FileOutputStream(OutFile));
    	
    	ArrayList<HashSet<Integer>> answer = new ArrayList<HashSet<Integer>>();
    	ArrayList<ArrayList<Pair<Integer,Double>>> predictor = new ArrayList<ArrayList<Pair<Integer,Double>>>();
    	this.getHPOAnswerPredictor(answer,predictor);
    	
    	ArrayList<Integer> HPOsize = new ArrayList<Integer>();
    	
    	
    	ArrayList<HashMap<Integer, Double>> predictorMap = new ArrayList<HashMap<Integer, Double>>();
    	for (int i = 0; i<predictor.size();i++)
    	{
    		predictorMap.add(new HashMap<Integer,Double>());
    		for (int j = 0; j<predictor.get(i).size();j++)
    		{
    			int lab = predictor.get(i).get(j).getFirst();
    			double score = predictor.get(i).get(j).getSecond();
    			predictorMap.get(i).put(lab, score);
    		}
    	}
    	for (int j=0;j<answer.size();j++)
    	{
    		HPOsize.add(answer.get(j).size());
    	}

		double HPOAUC = 0.0;  int HPOAUCsize = 0;
    	for (Pair<Integer, Double> e:this.hpoPredList)
    	{
    		int label = e.getFirst();
    		double frequency = e.getSecond();
    		ArrayList<Pair<Integer, Double>> predPair = new ArrayList<Pair<Integer, Double>>();
    		predPair.clear();
    		for (int j=0; j<answer.size();j++)
    		{
    			int size = 0;
    			size = HPOsize.get(j);
    			int result = 0;
    			double score = 0;
    			if (answer.get(j).contains(label)) result = 1; else result = 0;
    			if (predictorMap.get(j).containsKey(label)) score = predictorMap.get(j).get(label); else score = 0;
    			if (size>0) predPair.add(new Pair<Integer,Double>(result ,score));
    		}
    		ArrayList<Double> FPRlist = new ArrayList<Double>();
    		ArrayList<Double> TPRlist = new ArrayList<Double>();
    		
    		double AUC = evalution.CalculateROCAUC(predPair,FPRlist,TPRlist);
    		if (AUC>1E-12)
    		{
    			Fout.printf("%d,%.5f,%.4f\n", label ,frequency, AUC);
    			HPOAUC+=AUC;  HPOAUCsize++;
    		}
    	}
    	System.out.println("HPOAverageAUC = " + learningHPO.doubleFormat.format(HPOAUC/HPOAUCsize));
    	Fout.close();
    }
    public void evalEveryLabelAUC(String OutFile) throws FileNotFoundException
    {
    	PrintWriter Fout = new PrintWriter(new FileOutputStream(OutFile));
    	ArrayList<HashSet<Integer>> answer = new ArrayList<HashSet<Integer>>();
    	ArrayList<ArrayList<Pair<Integer,Double>>> predictor = new ArrayList<ArrayList<Pair<Integer,Double>>>();
        
    	ArrayList<Integer> MFsize = new ArrayList<Integer>();
    	ArrayList<Integer> BPsize = new ArrayList<Integer>();
    	ArrayList<Integer> CCsize = new ArrayList<Integer>();
    	
    	this.getAnswerPredictor(answer,predictor);

    	
    	ArrayList<HashMap<Integer, Double>> predictorMap = new ArrayList<HashMap<Integer, Double>>();
    	for (int i = 0; i<predictor.size();i++)
    	{
    		predictorMap.add(new HashMap<Integer,Double>());
    		for (int j = 0; j<predictor.get(i).size();j++)
    		{
    			int lab = predictor.get(i).get(j).getFirst();
    			double score = predictor.get(i).get(j).getSecond();
    			predictorMap.get(i).put(lab, score);
    		}
    	}
    	for (int j=0;j<answer.size();j++)
    	{
			MFsize.add(0);
			BPsize.add(0);
			CCsize.add(0);
    	}
		for (int j=0;j<answer.size();j++)
		{
			for (Integer e:answer.get(j))
			{
				if (learningOfGO.aGoSet.getSpace(e) == 'F') MFsize.set(j, MFsize.get(j) + 1);
				if (learningOfGO.aGoSet.getSpace(e) == 'P') BPsize.set(j, BPsize.get(j) + 1);
				if (learningOfGO.aGoSet.getSpace(e) == 'C') CCsize.set(j, CCsize.get(j) + 1);	
			}
		}
		Fout.println("label,sp,frequency,auc");
    	double MFOAUC = 0.0;  int MFOAUCsize = 0;
    	double BPOAUC = 0.0;  int BPOAUCsize = 0;
    	double CCOAUC = 0.0;  int CCOAUCsize = 0;
    	
    	
    	for (Pair<Integer, Double> e:this.predList)
    	{
    		int label = e.getFirst();
    		char sp = learningOfGO.aGoSet.getSpace(label);
    		double frequency = e.getSecond();
    		ArrayList<Pair<Integer, Double>> predPair = new ArrayList<Pair<Integer, Double>>();
    		predPair.clear();
    		for (int j=0; j<answer.size();j++)
    		{
    			int size = 0;
    			if (sp == 'F') size = MFsize.get(j);
    			if (sp == 'P') size = BPsize.get(j);
    			if (sp == 'C') size = CCsize.get(j);
    			int result = 0;
    			double score = 0;
    			if (answer.get(j).contains(label)) result = 1; else result = 0;
    			if (predictorMap.get(j).containsKey(label)) score = predictorMap.get(j).get(label); else score = 0;
    			if (size>0) predPair.add(new Pair<Integer,Double>(result ,score));
    		}
    		ArrayList<Double> FPRlist = new ArrayList<Double>();
    		ArrayList<Double> TPRlist = new ArrayList<Double>();
    		
    		double AUC = evalution.CalculateROCAUC(predPair,FPRlist,TPRlist);
    		/*debug 
    		if (label == 16740)
    		{
    			PrintWriter fout1 = new PrintWriter(new FileOutputStream("scoreLabel16740"));
    			PrintWriter fout2 = new PrintWriter(new FileOutputStream("ROCpoint16740"));
    			for (Pair<Integer,Double> onePair:predPair)
    			{
    				fout1.println(onePair.getFirst() + "," + onePair.getSecond());
    			}
    			Collections.sort(TPRlist);
    			Collections.sort(FPRlist);
    			for (int i=0;i<FPRlist.size();i++)
    			{
    				fout2.println(FPRlist.get(i) + "," + TPRlist.get(i));
    			}
    			fout1.close();
    			fout2.close();
    		}
    		*/
    		
    		if (AUC>1E-12)
    		{
    			Fout.printf("%d,%c,%.5f,%.4f\n", label , sp , frequency, AUC);
    			if (sp == 'F') {MFOAUC+=AUC;  MFOAUCsize++;}
    			if (sp == 'P') {BPOAUC+=AUC;  BPOAUCsize++;}
    			if (sp == 'C') {CCOAUC+=AUC;  CCOAUCsize++;}
    		}
    	}
    	System.out.println("MFOAverageAUC = " + learningOfGO.doubleFormat.format(MFOAUC/MFOAUCsize));
    	System.out.println("BPOAverageAUC = " + learningOfGO.doubleFormat.format(BPOAUC/BPOAUCsize));
    	System.out.println("CCOAverageAUC = " + learningOfGO.doubleFormat.format(CCOAUC/CCOAUCsize));
    	Fout.close();
    }
    
    
    public void evalutionFmaxAndAUPR(String resultOutFile,String OutPattern,String message) throws FileNotFoundException
    {
    	
    	PrintWriter Fout = new PrintWriter(new FileOutputStream(resultOutFile,true));
    	if ((OutPattern.equals("table")) && (message.equals("all")))
    	{
    		Fout.println("Sp,Onotology,size,AverageAnn,Fmax,Cut,AUPR");
    	}
    	ArrayList<Pair<Double,Double>> ReCallPrePair = new ArrayList<Pair<Double,Double>>();
    	
    	ArrayList<HashSet<Integer>> answer = new ArrayList<HashSet<Integer>>();
    	ArrayList<ArrayList<Pair<Integer,Double>>> predictor = new ArrayList<ArrayList<Pair<Integer,Double>>>();
    	
    	Pair<Double,Double> result;
    	
    	proteinSet MFOsubset = this.getSubAnnotation('F');
    	System.out.println(MFOsubset.sizeOfAnn());
    	
    	MFOsubset.getAnswerPredictor(answer,predictor);
    	
    	if (answer.size()<8000)
    	{
    		result = evalution.GetFMeasureMax(answer, predictor,ReCallPrePair);
    	}else
    	{
    		result = evalution.GetFMeasureMaxSimpleVersion(answer, predictor);
    	}
    	
    	learningOfGO.MFOcut = result.getSecond();
    	
    	Output.OutputPairList("../OutFile/AUPRlist/MFOAUPRlist",ReCallPrePair);
    	if (OutPattern.equals("table"))
    	{
    		double averageAnn = (double)MFOsubset.sizeOfAnn()/MFOsubset.size();
    		
    		Fout.printf("%s , %s , %d , %.3f , %.4f , %.4f , %.4f\n", message,"MFO",MFOsubset.size(),averageAnn,result.getFirst(),result.getSecond()
        			, evalution.CalculateAUC(ReCallPrePair));
    	}
    	
    	else
    	{
    		System.out.printf("SIZE = " + MFOsubset.size() + "  MFO :: Fmax/AUPR = %.4f/%.4f  Cut = %.4f  \n", result.getFirst()
        			, evalution.CalculateAUC(ReCallPrePair),result.getSecond());
    	}
    	
    	proteinSet BPOsubset = this.getSubAnnotation('P');
    	System.out.println(BPOsubset.sizeOfAnn());
    	BPOsubset.getAnswerPredictor(answer,predictor);
    	
    	if (answer.size()<8000)   //测试集过大时没有办法遍历所有的打分？内存不足
    	{
    		result = evalution.GetFMeasureMax(answer, predictor,ReCallPrePair);
    	}else
    	{
    		result = evalution.GetFMeasureMaxSimpleVersion(answer, predictor);
    	}
    	
    	learningOfGO.BPOcut = result.getSecond();
    	
    	Output.OutputPairList("../OutFile/AUPRlist/BPOAUPRlist",ReCallPrePair);
    	if (OutPattern.equals("table"))
    	{
    		double averageAnn = (double)BPOsubset.sizeOfAnn()/BPOsubset.size();
    		Fout.printf("%s , %s , %d , %.3f , %.4f , %.4f , %.4f\n", message,"BPO",BPOsubset.size(),averageAnn,result.getFirst(),result.getSecond()
        			, evalution.CalculateAUC(ReCallPrePair));
    	}else 
    	{
    		System.out.printf("SIZE = " + BPOsubset.size() + "  BPO :: Fmax/AUPR = %.4f/%.4f  Cut = %.4f  \n", result.getFirst()
        			, evalution.CalculateAUC(ReCallPrePair),result.getSecond());
    	}
    	
    	
    	
    	proteinSet CCOsubset = this.getSubAnnotation('C');
    	System.out.println(CCOsubset.sizeOfAnn());
    	CCOsubset.getAnswerPredictor(answer,predictor);
    	
    	
    	if (answer.size()<8000)
    	{
    		result = evalution.GetFMeasureMax(answer, predictor,ReCallPrePair);
    	}else
    	{
    		result = evalution.GetFMeasureMaxSimpleVersion(answer, predictor);
    	}
    	
    	learningOfGO.CCOcut = result.getSecond();
    	
    	Output.OutputPairList("../OutFile/AUPRlist/CCOAUPRlist",ReCallPrePair);

    	if (OutPattern.equals("table"))
    	{
    		double averageAnn = (double)CCOsubset.sizeOfAnn()/CCOsubset.size();
    		Fout.printf("%s , %s , %d , %.3f , %.4f , %.4f , %.4f\n", message,"CCO",CCOsubset.size(),averageAnn,result.getFirst(),result.getSecond()
        			, evalution.CalculateAUC(ReCallPrePair));
    	}else 
    	{
    		System.out.printf("SIZE = " + CCOsubset.size() + "  CCO :: Fmax/AUPR = %.4f/%.4f  Cut = %.4f  \n", result.getFirst()
    			, evalution.CalculateAUC(ReCallPrePair),result.getSecond());
    	}
    	System.out.println();
    	Fout.close();
    }  
    
    public void evalutionFmaxAndAUPR(double MFOcut,double BPOcut,double CCOcut,String resultOutFile,String OutPattern,String message) throws FileNotFoundException
    {
    	
    	PrintWriter Fout = new PrintWriter(new FileOutputStream(resultOutFile,true));
    	if ((OutPattern.equals("table")) && (message.equals("all")))
    	{
    		Fout.println("Sp,Onotology,size,AverageAnn,Fmax,Cut,AUPR");
    	}
    	ArrayList<Pair<Double,Double>> ReCallPrePair = new ArrayList<Pair<Double,Double>>();
    	
    	ArrayList<HashSet<Integer>> answer = new ArrayList<HashSet<Integer>>();
    	ArrayList<ArrayList<Pair<Integer,Double>>> predictor = new ArrayList<ArrayList<Pair<Integer,Double>>>();
    	
    	
    	proteinSet MFOsubset = this.getSubAnnotation('F');
    	System.out.println("MFOAnnsize = " + MFOsubset.sizeOfAnn());
    	
    	MFOsubset.getAnswerPredictor(answer,predictor);
    	double Fmeasure = 0.0;
    	Fmeasure = evalution.GetFMeasureMax(MFOcut,answer, predictor);
    	evalution.GetFMeasureMax(answer, predictor,ReCallPrePair);   //calculator AUPR	
    	Output.OutputPairList("../OutFile/AUPRlist/MFOAUPRlist",ReCallPrePair);
    	
    	if (OutPattern.equals("table"))
    	{
    		double averageAnn = (double)MFOsubset.sizeOfAnn()/MFOsubset.size();
    		Fout.printf("%s , %s , %d , %.3f , %.4f , %.4f , %.4f\n", message,"MFO",MFOsubset.size(),averageAnn,Fmeasure,MFOcut
        		,evalution.CalculateAUC(ReCallPrePair));
    	}
    	
    	else
    	{
    		System.out.printf("SIZE = " + MFOsubset.size() + "  MFO :: Fmeasure = %.4f  Cut = %.4f  \n", Fmeasure,MFOcut);
    	}
    	
    	proteinSet BPOsubset = this.getSubAnnotation('P');
    	System.out.println("BPOAnnsize = " + BPOsubset.sizeOfAnn());
    	BPOsubset.getAnswerPredictor(answer,predictor);
    	Fmeasure = evalution.GetFMeasureMax(BPOcut,answer, predictor);
    	evalution.GetFMeasureMax(answer, predictor,ReCallPrePair);   //calculator AUPR	
    	
    	Output.OutputPairList("../OutFile/AUPRlist/BPOAUPRlist",ReCallPrePair);
    	if (OutPattern.equals("table"))
    	{
    		double averageAnn = (double)BPOsubset.sizeOfAnn()/BPOsubset.size();
    		Fout.printf("%s , %s , %d , %.3f , %.4f , %.4f , %.4f\n", message,"BPO",BPOsubset.size(),averageAnn,Fmeasure,BPOcut
            		,evalution.CalculateAUC(ReCallPrePair));
    	}else 
    	{
    		System.out.printf("SIZE = " + BPOsubset.size() + "  BPO :: Fmeasure = %.4f  Cut = %.4f  \n", Fmeasure,BPOcut);
		}

    	proteinSet CCOsubset = this.getSubAnnotation('C');
    	System.out.println("CCOAnnsize = " + CCOsubset.sizeOfAnn());
    	CCOsubset.getAnswerPredictor(answer,predictor);
    	
    	Fmeasure = evalution.GetFMeasureMax(CCOcut,answer, predictor);
    	evalution.GetFMeasureMax(answer, predictor,ReCallPrePair);   //calculator AUPR	
    	
    	Output.OutputPairList("../OutFile/AUPRlist/CCOAUPRlist",ReCallPrePair);

    	if (OutPattern.equals("table"))
    	{
    		double averageAnn = (double)CCOsubset.sizeOfAnn()/CCOsubset.size();
    		Fout.printf("%s , %s , %d , %.3f , %.4f , %.4f , %.4f\n", message,"CCO",CCOsubset.size(),averageAnn,Fmeasure,CCOcut
            		,evalution.CalculateAUC(ReCallPrePair));
    	}else 
    	{
    		System.out.printf("SIZE = " + CCOsubset.size() + "  CCO :: Fmeasure = %.4f  Cut = %.4f  \n", Fmeasure,CCOcut);
    	}
    	System.out.println();
    	Fout.close();
    }  
    
    
    public void evalutionHPOFmaxAndAUPR() throws FileNotFoundException 
    {
    	ArrayList<Pair<Double,Double>> ReCallPrePair = new ArrayList<Pair<Double,Double>>();
    	ArrayList<HashSet<Integer>> answer = new ArrayList<HashSet<Integer>>();
    	ArrayList<ArrayList<Pair<Integer,Double>>> predictor = new ArrayList<ArrayList<Pair<Integer,Double>>>();
    	Pair<Double,Double> result;
    	

    	
    	this.getHPOAnswerPredictor(answer,predictor);
    	
    	if (answer.size()<8000)
    	{
    		result = evalution.GetFMeasureMax(answer, predictor,ReCallPrePair);
    	}else
    	{
    		result = evalution.GetFMeasureMaxSimpleVersion(answer, predictor);
    	}
    	
    	double HPOcut = result.getSecond();
    	
    	Output.OutputPairList("../OutFile/AUPRlist/MFOAUPRlist",ReCallPrePair);
    	System.out.printf("SIZE = " + this.size() + "  HPO :: Fmax/AUPR = %.4f/%.4f  Cut = %.4f  \n", result.getFirst()
        			, evalution.CalculateAUC(ReCallPrePair),result.getSecond());
    	
    	
    }
    
    
    public void CalTopkRecall(int... args)
    {
    	ArrayList<HashSet<Integer>> answer = new ArrayList<HashSet<Integer>>();
    	ArrayList<ArrayList<Pair<Integer,Double>>> predictor = new ArrayList<ArrayList<Pair<Integer,Double>>>();
    	
    	proteinSet MFOsubset = this.getSubAnnotation('F');
    	MFOsubset.getAnswerPredictor(answer,predictor);
		for (int node : args) 
		{
			evalution.CalTopKrecall(answer, predictor, node);
		}
		
    	proteinSet BPOsubset = this.getSubAnnotation('P');
    	BPOsubset.getAnswerPredictor(answer,predictor);
		for (int node : args) 
		{
			evalution.CalTopKrecall(answer, predictor, node);
		}
		
    	proteinSet CCOsubset = this.getSubAnnotation('C');
    	CCOsubset.getAnswerPredictor(answer,predictor);
		for (int node : args) 
		{
			evalution.CalTopKrecall(answer, predictor, node);
		}
    	
    }
    
    
    public void filterCAFA2Species()
    {
    	for (int i = 0;i<Cell.size();i++)
    	{	
    		oneProtein e = Cell.get(i);
    		String Species = e.getSpecies();
    		if (!proteinCommon.SetCAFA2Species.contains(Species))
    		{
    			Cell.remove(i);
    			i--;
    		}
    	}
    	this.UpdateIndex();
    }
    
    public void filterSpecies(TreeSet<String> speciesSet)
    {
    	for (int i = 0;i<Cell.size();i++)
    	{	
    		oneProtein e = Cell.get(i);
    		String Species = e.getSpecies();
    		if (!speciesSet.contains(Species))
    		{
    			Cell.remove(i);
    			i--;
    		}
    	}
    	this.UpdateIndex();
    }
    
    public void filterNoAccess()
    {
    	for (int i = 0;i<Cell.size();i++)
    	{	
    		oneProtein e = Cell.get(i);
    		String access = e.getAccess();
    		if (!proteinSet.MapUniAccess2Name.containsKey(access))
    		{
    			Cell.remove(i);
    			i--;
    		}
    	}
    }
    public void removeGoNotIn(GoSet aGo)
    {
    	for(oneProtein e:Cell)
    		e.filterGoNotIn(aGo);
    	this.removeProteinHaveNoGOAnn();
    }
    
    public void removeHPONotIn(HPOSet aHPOSet)
    {
    	for(oneProtein e:Cell)
    		e.filterHPONotIn(aHPOSet);
    	this.removeProteinHaveNoHPOAnn();
    }
    
    
    public void removeProteinHaveNoGOAnn()
    {
    	for (int i = 0;i<Cell.size();i++)
    	{	
    		if (Cell.get(i).getAnnotationSize()==0)
    		{
    			System.out.println("Remove Protin" + Cell.get(i).getAccess() + "Because have no GO annotation");
    			Cell.remove(i);
    			i--;
    		}
    	}
    	this.UpdateIndex();
    }
    public void removeProteinHaveNoHPOAnn()
    {
    	for (int i = 0;i<Cell.size();i++)
    	{	
    		if (Cell.get(i).getHPOAnnotationSize()==0)
    		{
    			System.out.println("Remove Protin" + Cell.get(i).getAccess() + "Because have no HPO annotation");
    			Cell.remove(i);
    			i--;
    		}
    	}
    	this.UpdateIndex();
    }
    
    public void getAnswerPredictor(ArrayList<HashSet<Integer>> answer , ArrayList<ArrayList<Pair<Integer,Double>>> predictor)
    {
    	answer.clear();
    	predictor.clear();
    	for(oneProtein e:Cell)
    	{
    		answer.add(new HashSet<Integer>(e.getAnnotation()));
    		predictor.add(e.getPredList());
    	}
    }
    
    public void getHPOAnswerPredictor(ArrayList<HashSet<Integer>> answer , ArrayList<ArrayList<Pair<Integer,Double>>> predictor)
    {
    	answer.clear();
    	predictor.clear();
    	for(oneProtein e:Cell)
    	{
    		answer.add(new HashSet<Integer>(e.getHPOAnnotation()));
    		predictor.add(e.getHPOPredList());
    	}
    }
        

    public void getIntersection(proteinSet other)
    {
    	for(int i=0;i<Cell.size();i++)
    	{
    		while (i<Cell.size() && ( !other.containProtein(Cell.get(i).getAccess()) ))
    		{
    			Cell.remove(i);
    		}
    	}
    	UpdateIndex();
    }
    
    public ArrayList<Integer> getAnnotation(String access)
    {
    	if (AccessIndex.containsKey(access))
    	{
    		int index = AccessIndex.get(access);
    		ArrayList<Integer> ann = new ArrayList<Integer>(Cell.get(index).getAnnotation());
    		return ann;
    	}
    	else
    		return new ArrayList<Integer>();
    }
    
    public ArrayList<Integer> getHPOAnnotation(String access)
    {
    	if (AccessIndex.containsKey(access))
    	{
    		int index = AccessIndex.get(access);
    		ArrayList<Integer> ann = new ArrayList<Integer>(Cell.get(index).getHPOAnnotation());
    		return ann;
    	}
    	else
    		return new ArrayList<Integer>();
    }
    
    public ArrayList<Integer> getAnnotation(int index)
    {
    	if (index<Cell.size())
    	{
    		ArrayList<Integer> ann = new ArrayList<Integer>(Cell.get(index).getAnnotation());
    		return ann;
    	}
    	else
    	{
    		System.out.println("index exceed border");
    		return new ArrayList<Integer>();
    	}
    }
    
    public HashSet<Integer> getAnnotationSet(int index)
    {
    	if (index<Cell.size())
    	{
    		HashSet<Integer> ann = Cell.get(index).getAnnotation();
    		return ann;
    	}
    	else
    	{
    		System.out.println("index exceed border");
    		return new HashSet<Integer>();
    	}
    }
    
    
    public oneProtein getProtein(int index)
    {
    	return Cell.get(index);
    }
    public oneProtein getProtein(String  Access)
    {
    	int index = AccessIndex.get(Access);
    	return Cell.get(index);
    }
    public String getProteinAccess(int index)
    {
    	return Cell.get(index).getAccess();
    }
    public String getProteinName(int index)
    {
    	String Access = Cell.get(index).getAccess();
    	return proteinSet.MapUniAccess2Name.get(Access);
    }
    public String getSequence(int index)
    {
    	return Cell.get(index).getSequence();
    }
    public int getProteinNum()
    {
    	return Cell.size();
    }
    public proteinSet getRandomSubProteinSet(int num)
    {
    	proteinSet leaveOneTrain =   new proteinSet();
    	if (num<=Cell.size())
    	{
    		Set<Integer> set = new HashSet<Integer>();
    		for (int i = 1; i <= num; i++)
    		{
    			Random random = new Random();
    			int s = random.nextInt(Cell.size());
    			while (set.contains(s))
    			{
    				s = random.nextInt(Cell.size());
    			}
    			set.add(s);
    			leaveOneTrain.addProtein(this.getProtein(s));
    		}
    	leaveOneTrain.UpdateIndex();
    	return leaveOneTrain;
    	}
    	else
    	{
    		System.out.println("the SubSet require is too big");
    		return leaveOneTrain;
    	}
    }
    public proteinSet getSubAnnotation(char space)
    {
    	proteinSet result = new proteinSet();
    	for(oneProtein e:Cell)
    	{
    		oneProtein aProtein = e.getSubProtein(space);
    		if (!(aProtein.getAnnotationSize() == 0))
    		result.addProtein(aProtein);
    	}
    	return result;
    }
    public void loadFastaSequence(String InFile) throws FileNotFoundException
    {
    	Scanner In = new Scanner(new FileInputStream(InFile));
    	String Access = new String();
    	String sequence = new String();
    	while (In.hasNext())
    	{
    		String line = In.nextLine();
    		if (line.substring(0,1).equals(">")) 
    		{
    		   Access = line.substring(1);  //System.out.println(Access);
    		   if (Access.contains("_")) Access = proteinSet.MapName2UniAccess.get(Access);
    		}	
    		else
    		{
    			sequence = "";
    			do
    			{
    				sequence = sequence + line;
    				line = In.nextLine();
    			}while (!line.equals(""));
    			addSequence(Access,sequence);
    		}
    	}
    	In.close();
    }
    
    public void loadTableSequence(String InFile) throws FileNotFoundException
    {
    	Scanner In = new Scanner(new FileInputStream(InFile));
    	String Access = new String();
    	String sequence = new String();
    	int count = 0;
    	while (In.hasNext())
    	{
    		count++;
    		if (count % 1000000 == 0) System.out.println(count); 
    		Access = In.next();
    		sequence = In.next();
    		addSequence(Access,sequence);
    	}
    	In.close();
    }

    public double CalCandidateRecall(char space)
    {
    	double recall = 0;
    	for(oneProtein e:Cell)
    	{
    		recall += e.CalCandidateRecall(space);
    	}
    	recall = (double)recall/this.size();
    	return recall;
    }
    
    public void clearPredResult()
    {
    	for(oneProtein e:Cell)
    	{
    		e.clearPredResult();
    	}
    }
    public void addTopK_L2RCandidate(int k,char space)
    {
    	for(oneProtein e:Cell)
    	{
    		e.addTopK_L2RCandidate(k, space);
    	}
    }
    public void naiveBaseline(proteinSet train,double cut) throws FileNotFoundException
    {
    	ArrayList<Pair<Integer,Double>> predList = train.getNaiveList();
    	for (int i = 0;i<predList.size();i++)
    	{
    		if ( (i<predList.size())  && (predList.get(i).getSecond()<cut ))
    		{
    			predList.remove(i);
    		}
    	}
    	for(oneProtein e:Cell)
    	{
    		e.setPredResult(predList);
    	}
    }
    
    public void ScoreNaiveBaseline(proteinSet train) throws FileNotFoundException
    {
    	ArrayList<Pair<Integer,Double>> predList = train.getNaiveList();
    	for(oneProtein e:Cell)
    	{
    		e.setPredResult(predList);
    	}
    }
    
    public void ScoreHPONaiveBaseline(proteinSet train) throws FileNotFoundException
    {
    	ArrayList<Pair<Integer,Double>> predList = train.getHPONaiveList();
    	for(oneProtein e:Cell)
    	{
    		e.setHPOPredResult(predList);
    	}
    }
    
    public void naiveSpeciesBaseline(proteinSet train,ArrayList<String> SpeciesList)
    {
    	HashMap<String,Integer> SpMap = new HashMap<String,Integer>();
    	for (int i = 0;i<SpeciesList.size();i++)
    		SpMap.put(SpeciesList.get(i), i);
    	ArrayList<ArrayList<Pair<Integer,Double>>> SpeciesNaivePredList = 
    			train.getSpeciesNaiveList(SpMap);
    	for (oneProtein e:Cell)
    	{
    		if (SpMap.containsKey(e.getSpecies()))
    		{
    			int index = SpMap.get(e.getSpecies());
    			e.setPredResult(SpeciesNaivePredList.get(index));
    		}else
    		{
    			e.setPredResult(SpeciesNaivePredList.get(SpeciesNaivePredList.size() - 1));
    		}
    	}
    }
    public ArrayList<ArrayList<Pair<Integer,Double>>> getSpeciesNaiveList(HashMap<String,Integer> SpMap)
    {
    	ArrayList<ArrayList<Pair<Integer,Double>>> predList = new ArrayList<ArrayList<Pair<Integer,Double>>>();
    	ArrayList<HashMap<Integer,Integer>> GoIndex = new ArrayList<HashMap<Integer,Integer>>();
    	int[] stat = new int[SpMap.size()+1];
    	System.out.println(SpMap.size());
    	for (int i = 0;i<=SpMap.size();i++)
    	{
    		predList.add(new ArrayList<Pair<Integer,Double>>());
    		GoIndex.add(new HashMap<Integer,Integer>());
    	}
    	for(oneProtein e:Cell)
    	{
    		int index = 0;
    		String sp = e.getSpecies();
    		if (SpMap.containsKey(sp)) index = SpMap.get(sp); else index = SpMap.size();
    		stat[index] += 1;
    		for (int gonum:e.getAnnotation())
    		{
    			if (GoIndex.get(index).containsKey(gonum))
    			{
    				int index1 = GoIndex.get(index).get(gonum);
    				predList.get(index).get(index1).setSecond(predList.get(index).get(index1).getSecond() + 1.0);
    				
    			}else
    			{
    				predList.get(index).add(new Pair<Integer,Double>(gonum,1.0));
    				GoIndex.get(index).put(gonum, predList.get(index).size()-1);
    			}	
    		}
    	}
    	
    	for(int i = 0;i<predList.size();i++)
    	{
    		for (Pair<Integer,Double> pair: predList.get(i))
    		{
    			pair.setSecond(pair.getSecond()/stat[i]);
    		}
    	}
    	return predList;
    }
    
    public void blastBaseline(proteinSet train) throws FileNotFoundException
    {
    	for(oneProtein e:Cell)
    	{
    		e.setBlastPred(train);
    	}
    }
    
    public void blastHPOBaseline(proteinSet train) throws FileNotFoundException
    {
    	for(oneProtein e:Cell)
    	{
    		e.setHPOBlastPred(train);
    	}
    }
    
    public void BlastKnnBaseline(proteinSet train) throws FileNotFoundException
    {
    	for(oneProtein e:Cell)
    	{
    		e.setBlastKnn(train);
    	}
    }
    
    public void BlastKnnHPOBaseline(proteinSet train) throws FileNotFoundException
    {
    	for(oneProtein e:Cell)
    	{
    		e.setBlastKnnHPO(train);
    	}
    }
    
    public void knn(proteinSet train) throws FileNotFoundException
    {
    	for(oneProtein e:Cell)
    	{
    		e.knn(train);
    	}
    }
    public ArrayList<Pair<Integer,Double>> getNaiveList()
    {
    	ArrayList<Pair<Integer,Double>> predList = new ArrayList<Pair<Integer,Double>>();
    	HashMap<Integer,Integer> GoIndex = new HashMap<Integer,Integer>();
    	for(oneProtein e:Cell)
    	{
    		for (int gonum:e.getAnnotation())
    		{
    			if (GoIndex.containsKey(gonum))
    			{
    				int index = GoIndex.get(gonum);
    				predList.get(index).setSecond(predList.get(index).getSecond() + 1.0);
    				
    			}else
    			{
    				predList.add(new Pair<Integer,Double>(gonum,1.0));
    				GoIndex.put(gonum, predList.size()-1);
    			}	
    		}
    	}
    	this.calMFBPCCsize(learningOfGO.aGoSet);
    	for (Pair<Integer,Double> pair:predList)
    	{
    		int size = this.getSubSetSize(learningOfGO.aGoSet.getSpace(pair.getFirst()));
    		pair.setSecond(pair.getSecond()/size);
    	}
    	return predList;
    }
    
    public ArrayList<Pair<Integer,Double>> getHPONaiveList()
    {
    	ArrayList<Pair<Integer,Double>> predList = new ArrayList<Pair<Integer,Double>>();
    	HashMap<Integer,Integer> HPOIndex = new HashMap<Integer,Integer>();
    	for(oneProtein e:Cell)
    	{
    		for (int HPOnum:e.getHPOAnnotation())
    		{
    			if (HPOIndex.containsKey(HPOnum))
    			{
    				int index = HPOIndex.get(HPOnum);
    				predList.get(index).setSecond(predList.get(index).getSecond() + 1.0);
    				
    			}else
    			{
    				predList.add(new Pair<Integer,Double>(HPOnum,1.0));
    				HPOIndex.put(HPOnum, predList.size()-1);
    			}	
    		}
    	}
    	for (Pair<Integer,Double> pair:predList)
    	{
    		pair.setSecond(pair.getSecond()/this.size());
    	}
    	return predList;
    }
    
    

    public void OutputAnnotation(String OutFile) throws FileNotFoundException
    {
    	PrintWriter Fout = new PrintWriter(new FileOutputStream(OutFile));
    	for(oneProtein e:Cell)
    		e.OutputAnnotation(Fout);
    	Fout.close();
    }
    
    public void OutputHPOAnnotation(String OutFile) throws FileNotFoundException
    {
    	PrintWriter Fout = new PrintWriter(new FileOutputStream(OutFile));
    	for(oneProtein e:Cell)
    		e.OutputHPOAnnotation(Fout);
    	Fout.close();
    }
    
    
    public void OutputAnnotationName(String OutFile) throws FileNotFoundException
    {
    	PrintWriter Fout = new PrintWriter(new FileOutputStream(OutFile));
    	for(oneProtein e:Cell)
    		e.OutputAnnotationName(Fout);
    	Fout.close();
    }
    
    public void OutputType2Annotation(String OutFile,GoSet aGoSet) throws FileNotFoundException
    {
    	PrintWriter Fout = new PrintWriter(new FileOutputStream(OutFile));
    	for(oneProtein e:Cell)
    		e.OutputType2Annotation(Fout,aGoSet);
    	Fout.close();
    }

    public void OutputType3Annotation(String OutFile,GoSet aGoSet) throws FileNotFoundException
    {
    	PrintWriter Fout = new PrintWriter(new FileOutputStream(OutFile));
    	for(oneProtein e:Cell)
    		e.OutputType3Annotation(Fout,aGoSet);
    	Fout.close();
    }
    
    public void OutputAnnotationSpace(String OutFile) throws FileNotFoundException
    {
    	PrintWriter Fout = new PrintWriter(new FileOutputStream(OutFile));
    	for(oneProtein e:Cell)
    		e.OutputAnnotationSpace(Fout);
    	Fout.close();
    }
    
    public void OutputNoSparseFeature(String OutFile) throws FileNotFoundException
    {
    	PrintWriter Fout = new PrintWriter(new FileOutputStream(OutFile));
    	for(oneProtein e:Cell)
    		e.OutputNoSparseFeature(Fout);
    	Fout.close();
    }
    
    public void outputMSA(String OutFile) throws FileNotFoundException
    {
    	PrintWriter Fout = new PrintWriter(new FileOutputStream(OutFile));
    	for(oneProtein e:Cell)
    	{
    		e.outputMSA(Fout);
    		Fout.println("***");
    	}
    	Fout.close();
    }
    public void OutputFastaAccessSequence(String OutFile) throws FileNotFoundException
    {
    	PrintWriter Fout = new PrintWriter(new FileOutputStream(OutFile));
    	for(oneProtein e:Cell)
    		e.OutputFastaAccessSequence(Fout);
    	Fout.close();
    }
    
    public void OutputFastaNameSequence(String OutFile) throws FileNotFoundException
    {
    	PrintWriter Fout = new PrintWriter(new FileOutputStream(OutFile));
    	for(oneProtein e:Cell)
    		e.OutputFastaNameSequence(Fout);
    	Fout.close();
    }
    
    public void OutputGOPredScore(String OutFile) throws FileNotFoundException
    {
    	PrintWriter Fout = new PrintWriter(new FileOutputStream(OutFile));
    	for (oneProtein e:Cell)
    	{
    		e.OutputGOPredScore(Fout);
    	}
    	Fout.close();
    }
    
    public void OutputHPOPredScore(String OutFile) throws FileNotFoundException
    {
    	PrintWriter Fout = new PrintWriter(new FileOutputStream(OutFile));
    	for (oneProtein e:Cell)
    	{
    		e.OutputHPOPredScore(Fout);
    	}
    	Fout.close();
    }
    
    public void LoadPredScore(String InFile) throws FileNotFoundException
    {
    	Scanner In = new Scanner(new FileInputStream(InFile));
    	String Acc;
    	int count = 0;
    	int Ann;
    	double score;
    	int index;
    	while (In.hasNext())
    	{
    		Acc = In.nextLine();
    		System.out.println(Acc);
    		index = this.AccessIndex.get(Acc);
    		count = In.nextInt();
    		for (int i=0;i<count;i++)
    		{
    			Ann = In.nextInt();
    			score = In.nextDouble();
    			this.Cell.get(index).addPred(new Pair<Integer,Double>(Ann,score));
    		}
    		
    	}
    	In.close();
    }
    
    public void OutputProteinAccess(String OutFile) throws FileNotFoundException
    {
    	PrintWriter Fout = new PrintWriter(new FileOutputStream(OutFile));
    	for(oneProtein e:Cell)
    		Fout.println(e.getAccess());
    	Fout.close();
    }
    
    public void OutputProteinName(String OutFile) throws FileNotFoundException
    {
    	PrintWriter Fout = new PrintWriter(new FileOutputStream(OutFile));
    	for(oneProtein e:Cell)
    		Fout.println(e.getName());
    	Fout.close();
    }
    
    public void outputLibTrainFile(String OutFile,int Ann) throws FileNotFoundException
    {
    	PrintWriter Fout = new PrintWriter(new FileOutputStream(OutFile));
    	for(oneProtein e:Cell)
    	{
    		if (e.containAnnotation(Ann)) Fout.print(1); else Fout.print(-1);
    		e.OutputSparseFeature(Fout);
    		Fout.println();
    	}
    	Fout.close();
    }
    
    public void OutputPredList(String OutFile) throws FileNotFoundException
    {
    	PrintWriter Fout = new PrintWriter(new FileOutputStream(OutFile));
    	for (Pair<Integer,Double> e:this.predList)
    	{
    		Fout.println(e.getFirst() + " , " + e.getSecond());
    	}
    	Fout.close();
    }
    
    public void tranSequence2TriSparseFeature()
    {
    	for(oneProtein e:Cell)
    		e.tranSequence2TriSparseFeature();
    }
    
    public int size() 
    {
    	return Cell.size();
    }
    public void setPredList(ArrayList<Pair<Integer,Double>> predList)
    {
    	this.predList = predList;
    }
    public void setPredList(proteinSet train)   
    {
    	this.predList = train.getNaiveList();
    	this.sortPredListBaseFrequency();
    }
    
    public void setHPOPredList(proteinSet train)   
    {
    	this.hpoPredList = train.getHPONaiveList();
    	this.sortPredListBaseFrequency();
    }
    
    public ArrayList<Integer> getHPOLabelList()   
    {
    	Comparator<Pair<Integer,Double>> compar = new Comparator<Pair<Integer,Double>>()
		{
			@Override
			public int compare(Pair<Integer, Double> o1, Pair<Integer, Double> o2) {
				// TODO Auto-generated method stub
				return  o2.getSecond().compareTo(o1.getSecond());
			}
		};
    	this.hpoPredList = this.getHPONaiveList();
    	Collections.sort(this.hpoPredList,compar);
    	
    	ArrayList<Integer> LabelSet = new ArrayList<Integer>();
    	for (Pair<Integer,Double> pair:hpoPredList)
    	{
    		LabelSet.add(pair.getFirst());
    	}
    	return LabelSet;
    }
    
    public void sortPredListBaseFrequency()
    {
		Comparator<Pair<Integer,Double>> compar = new Comparator<Pair<Integer,Double>>()
		{
			@Override
			public int compare(Pair<Integer, Double> o1, Pair<Integer, Double> o2) {
				// TODO Auto-generated method stub
				return  o1.getSecond().compareTo(o2.getSecond());
			}
		};
    	
    	Collections.sort(this.predList,compar);
    }
    public void sortPredListBaseLabelOrder()
    {
    	Collections.sort(this.predList);
    }
    public void setLiblinearFeatureFromSparseFeature()
    {
    	for (oneProtein e:Cell)
    		e.setLiblinearFeatureFromSparseFeature();
    }
    
    public void setLiblinearFeatureFromNoSparseFeature()
    {
    	for (oneProtein e:Cell)
    		e.setLiblinearFeatureFromNoSparseFeature();
    }
    
    public void SortProtein()
    {
    	Collections.sort(Cell);
    	UpdateIndex();
    }
    public void removeLowPred(int count)
    {
    	for (oneProtein e:Cell)
    	{
    		e.removeLowPred(count);
    	}
    }
    
    public void removeAnnotation(int... args)
    {
    	for (oneProtein e:Cell)
    		e.removeAnnotation(args);
    }
    
    public void removeHPOAnnotation(int... args)
    {
    	for (oneProtein e:Cell)
    		e.removeHPOAnnotation(args);
    }
    
    public void removeAllAnn()
    {
    	for (oneProtein e:Cell)
    		e.removeAllAnn();
    }
    public void removeIEA_Annotation(int... args)
    {
    	for (oneProtein e:Cell)
    		e.removeIEA_Annotation(args);
    }
    public void evalutionSpecies(String resultOutFile,String... Species) throws FileNotFoundException
    {
    	for (String sp:Species)
    	{
    		proteinSet newSet = this.getSpeciesSubSet(sp);
    		System.out.println(sp + "  evalution" + "  SIZE = " + newSet.size());
    		newSet.evalutionFmaxAndAUPR(resultOutFile,learningOfGO.resultOutPattern,sp);
    	}
    }
    
    public void evalutionSpecies(double MFOcut,double BPOcut,double CCOcut,String resultOutFile,String... Species) throws FileNotFoundException
    {
    	for (String sp:Species)
    	{
    		proteinSet newSet = this.getSpeciesSubSet(sp);
    		System.out.println(sp + "  evalution" + "  SIZE = " + newSet.size());
    		newSet.evalutionFmaxAndAUPR(MFOcut,BPOcut,CCOcut,resultOutFile,learningOfGO.resultOutPattern,sp);
    	}
    }
    
    public proteinSet getSpeciesSubSet(String sp)
    {
    	proteinSet newSet = new proteinSet();
    	for (oneProtein e:Cell)
    	{
    		if (e.getSpecies().equals(sp))
    			newSet.addProtein(e);
    		//System.out.println(e.getSpecies());
    	}
    	return newSet;		
    }
    public void removeAnnotation(Set<Integer> args)
    {
    	for (oneProtein e:Cell)
    		e.removeAnnotation(args);
    }
    public void statAminoAcid()
    {
		Map<Character,Integer> AcidStat = new HashMap<Character,Integer>();
		for (int i = 0;i<this.size();i++)
		{
			String seq = this.getSequence(i);
			for (int j = 0; j<seq.length(); j++)
			{
				char ch = seq.charAt(j);
				if (!AcidStat.containsKey(ch))
					AcidStat.put(ch, 1);
				else
				{
					int num = AcidStat.get(ch);
					AcidStat.put(ch, num + 1);
				}
			}
		}
		int sum = 0;
		for (Map.Entry<Character,Integer> entry : AcidStat.entrySet()) 
			sum += entry.getValue();
		
    	for (Map.Entry<Character,Integer> entry : AcidStat.entrySet()) 
    	{
    		System.out.println(String.format("%c,%d,%.2f%%", entry.getKey(),entry.getValue(),(double)entry.getValue()/sum*100));
    	}
    }
    public void statSpecies()
    {
    	Map<String,Integer> SpeciesStat = new HashMap<String,Integer>();
    	for(oneProtein e:Cell)
    	{
    		String Species = e.getSpecies();
    		if (SpeciesStat.containsKey(Species))
    		{
    			int num = SpeciesStat.get(Species);
    			SpeciesStat.put(Species, num + 1);
    		}
    		else
    			SpeciesStat.put(Species, 1);
    		
    	}
    	for (Map.Entry<String,Integer> entry : SpeciesStat.entrySet()) 
    	{
    		System.out.println(entry.getKey() + "\t" + entry.getValue());
    	}
    }
    public void statLabelNum()
    {
    	Map<Integer,Integer> LabelNumStat = new TreeMap<Integer,Integer>();
    	Map<Integer,Integer> MFLabelNumStat = new TreeMap<Integer,Integer>();
    	Map<Integer,Integer> BPLabelNumStat = new TreeMap<Integer,Integer>();
    	Map<Integer,Integer> CCLabelNumStat = new TreeMap<Integer,Integer>();
    	for(oneProtein e:Cell)
    	{
    		int LabelNum = e.getAnnotationSize();
    		int MFLabelNum = e.getAnnotationSize('F');
    		int BPLabelNum = e.getAnnotationSize('P');
    		int CCLabelNum = e.getAnnotationSize('C');
    		if (LabelNumStat.containsKey(LabelNum))
    		{
    			int num = LabelNumStat.get(LabelNum);
    			LabelNumStat.put(LabelNum, num + 1);
    		}
    		else
    			LabelNumStat.put(LabelNum, 1);
    		
    		if (MFLabelNumStat.containsKey(MFLabelNum))
    		{
    			int num = MFLabelNumStat.get(MFLabelNum);
    			MFLabelNumStat.put(MFLabelNum, num + 1);
    		}
    		else
    			MFLabelNumStat.put(MFLabelNum, 1);
    		
    		if (BPLabelNumStat.containsKey(BPLabelNum))
    		{
    			int num = BPLabelNumStat.get(BPLabelNum);
    			BPLabelNumStat.put(BPLabelNum, num + 1);
    		}
    		else
    			BPLabelNumStat.put(BPLabelNum, 1);
    		
    		if (CCLabelNumStat.containsKey(CCLabelNum))
    		{
    			int num = CCLabelNumStat.get(CCLabelNum);
    			CCLabelNumStat.put(CCLabelNum, num + 1);
    		}
    		else
    			CCLabelNumStat.put(CCLabelNum, 1);
    			
    	}
    	for (Map.Entry<Integer,Integer> entry : LabelNumStat.entrySet()) 
    	{
    		System.out.println(entry.getKey() + "\t" + entry.getValue());
    	}
    	
    	for (Map.Entry<Integer,Integer> entry : MFLabelNumStat.entrySet()) 
    	{
    		System.out.println(entry.getKey() + "\t" + entry.getValue());
    	}
    	
    	for (Map.Entry<Integer,Integer> entry : BPLabelNumStat.entrySet()) 
    	{
    		System.out.println(entry.getKey() + "\t" + entry.getValue());
    	}
    	
    	for (Map.Entry<Integer,Integer> entry : CCLabelNumStat.entrySet()) 
    	{
    		System.out.println(entry.getKey() + "\t" + entry.getValue());
    	}
    	
    }
    public void statLabelDepth()
    {
    	int[] MFdepth = new int[18];
    	int[] BPdepth = new int[18];
    	int[] CCdepth = new int[18];
    	for(oneProtein e:Cell)
    	{
    		Set<Integer> AnnSet = e.getAnnotation();
    		for (Integer ann:AnnSet)
    		{
    			char sp = learningOfGO.aGoSet.getSpace(ann);
    			if (sp == 'F') MFdepth[learningOfGO.aGoSet.getMindepth(ann)]++;
    			if (sp == 'P') BPdepth[learningOfGO.aGoSet.getMindepth(ann)]++;
    			if (sp == 'C') CCdepth[learningOfGO.aGoSet.getMindepth(ann)]++;
    		}
    	}
    	System.out.println("depth,MF,BP,CC");
    	for (int i=1;i<=17;i++)
    	{
    		System.out.println(i + "," + MFdepth[i] + "," + BPdepth[i] + "," + CCdepth[i]);
    	}
    }
    public void statTriAminoAcid(String OutFile) throws FileNotFoundException
    {
    	PrintWriter Fout = new PrintWriter(new FileOutputStream(OutFile));
		Map<String,Integer> AcidStat = new HashMap<String,Integer>();
		for (int i = 0;i<this.size();i++)
		{
			String seq = this.getSequence(i);
			for (int j = 0; j<seq.length()-2; j++)
			{
				String ch = seq.substring(j, j+3);
				if (!AcidStat.containsKey(ch))
					AcidStat.put(ch, 1);
				else
				{
					int num = AcidStat.get(ch);
					AcidStat.put(ch, num + 1);
				}
			}
		}
    	for (Map.Entry<String,Integer> entry : AcidStat.entrySet()) 
    	{
    		Fout.println(entry.getKey() + "\t" + entry.getValue());
    	}
    	Fout.close();
    }
    public void statGoFrequency(String FileNamePrefix) throws FileNotFoundException
    {
    	PrintWriter Fout = new PrintWriter(new FileOutputStream(FileNamePrefix + "GoFrequency.csv"));
    	int[] statGoFre = new int[learningOfGO.aGoSet.size()];
    	int[] statMFOFre = new int[learningOfGO.aGoSet.size()];
    	int[] statBPOFre = new int[learningOfGO.aGoSet.size()];
    	int[] statCCOFre = new int[learningOfGO.aGoSet.size()];
    	 for (oneProtein e:Cell)
    	 {
    		 HashSet<Integer> AnnSet = e.getAnnotation();
    		 for (int ann:AnnSet)
    		 {
    			 statGoFre[learningOfGO.aGoSet.getIndex(ann)]++;
    			 if (learningOfGO.aGoSet.getSpace(ann) == 'F') statMFOFre[learningOfGO.aGoSet.getIndex(ann)]++;
    			 if (learningOfGO.aGoSet.getSpace(ann) == 'P') statBPOFre[learningOfGO.aGoSet.getIndex(ann)]++;
    			 if (learningOfGO.aGoSet.getSpace(ann) == 'C') statCCOFre[learningOfGO.aGoSet.getIndex(ann)]++;
    		 }
    	 }
    	 Fout.println("GO , Frequency , Space");
    	 int SubSize;
    	 for (int i = 0;i < statGoFre.length;i++)
    	 {
    		 SubSize = this.getSubSetSize(learningOfGO.aGoSet.getIndexSpace(i));
    		 
    		 Fout.printf("%s , %.4f , %c\n", learningOfGO.aGoSet.getGOFromIndex(i) ,(double)statGoFre[i]/SubSize,learningOfGO.aGoSet.getIndexSpace(i));
    	 }
    	 common.statiscit.statDistribution(statGoFre, FileNamePrefix + "GOdistribution.csv");
    	 common.statiscit.statDistribution(statMFOFre, FileNamePrefix + "MFOdistribution.csv");
    	 common.statiscit.statDistribution(statBPOFre, FileNamePrefix + "BPOdistribution.csv");
    	 common.statiscit.statDistribution(statCCOFre, FileNamePrefix + "CCOdistribution.csv");
    	 Fout.close();
    }
    public void statSpeciesGoNum(TreeSet<String> spSet, String OutFile) throws FileNotFoundException
    {
    	PrintWriter Fout = new PrintWriter(new FileOutputStream(OutFile));
    	int spNum = spSet.size();
    	double[] GOstat = new double[spNum + 2];
    	double[] MFOstat = new double[spNum + 2];
    	double[] BPOstat = new double[spNum + 2];
    	double[] CCOstat = new double[spNum + 2];
    	double GOsum = 0;
    	double MFOsum = 0;
    	double BPOsum = 0;
    	double CCOsum = 0;
    	
    	
    	int MFOsize = 0;
    	int BPOsize = 0;
    	int CCOsize = 0;
    	
    	int[] GOnum = new int[spNum + 2];
    	int[] MFOnum = new int[spNum + 2];
    	int[] BPOnum = new int[spNum + 2];
    	int[] CCOnum = new int[spNum + 2];
    	
    	TreeMap<String,Integer> MapSpIndex = new TreeMap<String,Integer>();
    	int Index = 1;
    	for (String sp :spSet)
    	{
    		MapSpIndex.put(sp, Index);
    		Index++;
    	}
    	for (oneProtein e:Cell)
    	{
    		String Name = e.getAccess();
    		String sp = proteinCommon.getSpecies(Name);
    		int spIndex = 0;
    		if (spSet.contains(sp)) 
    			spIndex = MapSpIndex.get(sp);
    		else spIndex = spNum + 1;
    		
    		int GoLabelNum = e.getAnnotationSize();
    		int MFOlabel = e.getAnnotationSize('F');
    		int BPOlabel = e.getAnnotationSize('P');
    		int CCOlabel = e.getAnnotationSize('C');
    		GOstat[spIndex]  += GoLabelNum;
    		GOsum += GoLabelNum;
    		MFOstat[spIndex] += MFOlabel;
    		MFOsum += MFOlabel;
    		BPOstat[spIndex] += BPOlabel;
    		BPOsum += BPOlabel;
    		CCOstat[spIndex] += CCOlabel;
    		CCOsum += CCOlabel;
    		GOnum[spIndex]++;
    		if (MFOlabel>0)  {MFOnum[spIndex]++;    MFOsize++;}
    		if (BPOlabel>0)  {BPOnum[spIndex]++;    BPOsize++;}
    		if (CCOlabel>0)  {CCOnum[spIndex]++;    CCOsize++;}
    	}
    	Index = 1;
    	Fout.println("species,GOaverage,GOnum,MFOaverage,MFOnum,BPOaverage,BPOnum,CCOaverage,CCOnum");
		Fout.println(String.format("ALL" + ",%.2f,%d,%.2f,%d,%.2f,%d,%.2f,%d",GOsum/this.size(),this.size(),MFOsum/MFOsize,MFOsize
				,BPOsum/BPOsize,BPOsize,CCOsum/CCOsize,CCOsize));
    	for (String sp :spSet)
    	{
    		Fout.println(String.format(sp + ",%.2f,%d,%.2f,%d,%.2f,%d,%.2f,%d",GOstat[Index]/GOnum[Index],GOnum[Index],MFOstat[Index]/MFOnum[Index],MFOnum[Index]
    				,BPOstat[Index]/BPOnum[Index],BPOnum[Index],CCOstat[Index]/CCOnum[Index],CCOnum[Index]));
    		Index++;
    	}
		Fout.println(String.format("Others,%.2f,%d,%.2f,%d,%.2f,%d,%.2f,%d",GOstat[Index]/GOnum[Index],GOnum[Index],MFOstat[Index]/MFOnum[Index],MFOnum[Index]
				,BPOstat[Index]/BPOnum[Index],BPOnum[Index],CCOstat[Index]/CCOnum[Index],CCOnum[Index]));
		
    	Fout.close();
    }
    public void statSpeciesStringLength(TreeSet<String> spSet)
    {
    	int max = 0;       String maxAccess = new String();
    	int min = 99999;   String minAccess = new String();
    	double lengthSum = 0;
    	int spNum = spSet.size();
    	double[] seqStat = new double[spNum + 2];
    	int[] spDis = new int[spNum + 2];
    	TreeMap<String,Integer> MapSpIndex = new TreeMap<String,Integer>();
    	int Index = 1;
    	for (String sp :spSet)
    	{
    		MapSpIndex.put(sp, Index);
    		Index++;
    	}
    	for (oneProtein e:Cell)
    	{
    		String Name = e.getName();
    		String sp = proteinCommon.getSpecies(Name);
    		int spIndex = 0;
    		if (spSet.contains(sp)) 
    			spIndex = MapSpIndex.get(sp);
    		else spIndex = spNum + 1;
    		int length = e.getSequence().length();
    		seqStat[spIndex] += length;
    		if (length>max) {max = length;  maxAccess = e.getAccess();}
    		if (length<min) {min = length;   minAccess = e.getAccess();}
    		lengthSum += length;
    		spDis[spIndex]++;
    	}
    	Index = 1;
    	System.out.println("Species & average length");
    	for (String sp :spSet)
    	{
    		System.out.println(String.format("%s   &  %.1f", sp,seqStat[Index]/spDis[Index]));
    		Index++;
    	}
    	System.out.println(String.format("Others   &  %.1f", seqStat[Index]/spDis[Index]));
    	
    	System.out.println("max = " + max + "\t" + maxAccess);
    	System.out.println("min = " + min + "\t" + minAccess);
    	System.out.println("average = " + lengthSum/this.size());
    	
    }
    public void libLinearPredict(String Direction) throws IOException
    {
    	for (int count = 0;count<this.predList.size();count++)
    	{
    		int labelNum = predList.get(count).getFirst();
    		try
    		{
    			liblinear.Model model = Linear.loadModel(new File(Direction + "Label" + labelNum + ".model"));
    			if (count % 100 ==0)
    				System.out.println("Predict " + count + "th Label" + labelNum);
    			for(oneProtein e:Cell)	
    				e.libLinearPredict(labelNum,model);
    		}
    		catch (IOException exp) 
    		{
    			System.out.println("Not find Lable " + labelNum + " Model File");
    		}
    	}
    }
    
    public void libLinearPredict(String Direction, int labelNum) throws IOException
    {
    	System.out.println("Predict Label" + labelNum);
    	liblinear.Model model = Linear.loadModel(new File(Direction + "Label" + labelNum + ".model"));
		for(oneProtein e:Cell)
		{
			e.libLinearPredict(labelNum,model);
		}
    }
    
    public void libLinearTrain(String Direction, int ThreadNum) throws IOException
    {
    	this.calMFBPCCsize(learningOfGO.aGoSet);
    	
    	MultiThreadTrain my = new proteinSet.MultiThreadTrain(this.predList.size() - 1,Direction);
    	
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
    
	public class MultiThreadTrain implements Runnable {
		
		private int count;
		private String Direction  = new String();
		public MultiThreadTrain(int countNum,String aDirection)
		{
			count = countNum;
			Direction = aDirection;
		}
		public void run() {
			// TODO Auto-generated method stub
			while(count>=0)
			{
				try {
					int num = this.get();
					System.out.println(Thread.currentThread().getName()+ "now train label" + predList.get(num).getFirst());
					libLinearTrain(predList.get(num).getFirst(),Direction);
					
				} catch (IOException e) {
					// TODO Auto-generated catch block
					e.printStackTrace();
				}
			}
		}
		synchronized public int get(){
			return this.count--;
		}
		
	}
    
	
	public void libLinearTrain(int labelNum, String Direction) throws IOException
    {
    	int Ann = labelNum;
        int max_index = learningOfGO.FeatureSize;
        double bias = -1;
        
        
        liblinear.Problem   prob    = new liblinear.Problem();
        
        prob.bias = bias;
        
        char sp = learningOfGO.aGoSet.getSpace(labelNum);
        prob.l = this.getSubSetSize(sp);
        
        prob.n = max_index;
        if (bias >= 0) {
            prob.n++;
        }
        prob.x = new liblinear.Feature[prob.l][];
        prob.y = new double[prob.l];
        List<liblinear.Feature[]> vx = new ArrayList<liblinear.Feature[]>();
        int count = 0;
        for (int i = 0; i < this.size(); i++)
        {
        	if (this.get(i).getSubAnnSize(sp)>0)
        	{
        		liblinear.Feature[] x = this.Cell.get(i).getliblinearFeature();
        		vx.add(x);
        		prob.x[count] = vx.get(count);
        		count++;
        	}
        }
        count = 0;
        for (int i = 0; i < this.size(); i++)
        	if (this.get(i).getSubAnnSize(sp)>0)
        	{
        		if (this.Cell.get(i).containAnnotation(Ann))	
        			prob.y[count] = 1; else prob.y[count] = -1;
        		count++;
        	}
    	liblinear.Model model = liblinear.Linear.train(prob, learningOfGO.liblinearParam);
    	liblinear.Linear.saveModel(new File(Direction + "Label" + labelNum + ".model"), model);
    }
    
    
	
	public void HPOlibLinearTrain(String Direction, int ThreadNum) throws IOException
    {
    	HPOMultiThreadTrain my = new proteinSet.HPOMultiThreadTrain(this.hpoPredList.size() - 1,Direction);
    	
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
    
	public class HPOMultiThreadTrain implements Runnable {
		
		private int count;
		private String Direction  = new String();
		public HPOMultiThreadTrain(int countNum,String aDirection)
		{
			count = countNum;
			Direction = aDirection;
		}
		public void run() {
			// TODO Auto-generated method stub
			while(count>=0)
			{
				try {
					int num = this.get();
					System.out.println(Thread.currentThread().getName()+ "now train label" + hpoPredList.get(num).getFirst());
					libLinearTrain(hpoPredList.get(num).getFirst(),Direction);
					
				} catch (IOException e) {
					// TODO Auto-generated catch block
					e.printStackTrace();
				}
			}
		}
		synchronized public int get(){
			return this.count--;
		}
		
	}
    
	
	public void HPOlibLinearTrain(int labelNum, String Direction) throws IOException
    {
    	int Ann = labelNum;
		int max_index = proteinSet.FeatureSize;
        double bias = -1;
        
        
        liblinear.Problem   prob    = new liblinear.Problem();
        
        prob.bias = bias;
        
        char sp = learningOfGO.aGoSet.getSpace(labelNum);
        prob.l = this.getSubSetSize(sp);
        
        prob.n = max_index;
        if (bias >= 0) {
            prob.n++;
        }
        prob.x = new liblinear.Feature[prob.l][];
        prob.y = new double[prob.l];
        List<liblinear.Feature[]> vx = new ArrayList<liblinear.Feature[]>();
        int count = 0;
        for (int i = 0; i < this.size(); i++)
        {
        	if (this.get(i).getSubAnnSize(sp)>0)
        	{
        		liblinear.Feature[] x = this.Cell.get(i).getliblinearFeature();
        		vx.add(x);
        		prob.x[count] = vx.get(count);
        		count++;
        	}
        }
        count = 0;
        for (int i = 0; i < this.size(); i++)
        	if (this.get(i).getSubAnnSize(sp)>0)
        	{
        		if (this.Cell.get(i).containAnnotation(Ann))	
        			prob.y[count] = 1; else prob.y[count] = -1;
        		count++;
        	}
    	liblinear.Model model = liblinear.Linear.train(prob, learningOfGO.liblinearParam);
    	liblinear.Linear.saveModel(new File(Direction + "Label" + labelNum + ".model"), model);
    }
	
	public void GoFDR(proteinSet train)
    {
    	for(oneProtein e:Cell)
    	{
    		e.addBlastCandidate(train);
    		
    	}
    }
    public void UpdateIndex()
    {
    	AccessIndex.clear();
    	int index = 0;
    	for(oneProtein e:Cell)
    	{
    		String Access = e.getAccess();
    		AccessIndex.put(Access, index);
    		index++;
    	}
    }
    public void removeLabelRareProtein(int labelnum)
    {
    	for(int i=0;i<Cell.size();i++)
    	{
    		while (i<Cell.size() && (Cell.get(i).getAnnotationSize()<labelnum))
    		{
    			Cell.remove(i);
    		}
    	}
    	UpdateIndex();
    }
    public void removeRedunAnn(GoSet aGoSet)
    {
    	if (!aGoSet.isAddAllFather())
    		aGoSet.AddAllFather();
    	for (oneProtein e:Cell)
    	{
    		e.removeRedunAnn(aGoSet);
    	}
    }
    public void removeNoSeqProtein()
    {
    	for(int i=0;i<Cell.size();i++)
    	{
    		while (i<Cell.size() && (Cell.get(i).getSequence().length()<1))
    		{
    			Cell.remove(i);
    		}
    	}
    	UpdateIndex();
    }
    public void calOneLabelAUC(int labelNum)
    {
    	char space = learningOfGO.aGoSet.getSpace(labelNum);
    	
    }
    public double miniscore()
    {
    	double miniscore = 1.0;
    	if (this.attribute.containsKey("miniscore"))
    		miniscore = (double) this.attribute.get("miniscore");
    	else 
    	{
    		double score = 1.0;
		 	for (oneProtein e:Cell)
		 	{
		 		score = e.miniscore();
		 		if (miniscore>score) miniscore = score;
		 	}
		 	this.attribute.put("miniscore", miniscore);
		}
    	return miniscore;
    }
    
    public void removeBlastOtherSp()
    {
    	for (oneProtein e:Cell)
    	{
    		e.removeBlastOtherSp();
    	}
    }
    
    public void removeAnnMF_Only5515()
    {
    	for (oneProtein e:Cell)
    	{
    		HashSet<Integer> Ann = e.getSubAnnotation('F');
    		if ((Ann.size() == 1) && (Ann.contains(5515))) e.removeAnnotation(5515);
    	}
    }
    
    public void outputPubMedID(String OutFile) throws FileNotFoundException
    {
    	PrintWriter Fout = new PrintWriter(new FileOutputStream(OutFile));
    	for(oneProtein e:Cell)
    		e.outputPubMedID(Fout);
    	Fout.close();
    }
    
    public void outputIntActProtein(String OutFile) throws FileNotFoundException
    {
    	PrintWriter Fout = new PrintWriter(new FileOutputStream(OutFile));
    	for(oneProtein e:Cell)
    		e.outputIntActProtein(Fout);
    	Fout.close();
    }
    
    public void loadPubMed(String InFile) throws FileNotFoundException
    {
    	Scanner In = new Scanner(new FileInputStream(InFile));
    	while(In.hasNext())
    	{
    		String access = In.next();
    		if (!this.AccessIndex.containsKey(access))
    		{
    			this.addProtein(access);
    			this.AccessIndex.put(access, Cell.size() - 1);
    		}
    		int index = this.AccessIndex.get(access);
    		int count = In.nextInt();
    		int PubMedID;
    		for (int i = 1;i<=count;i++)
    		{
    			 PubMedID = In.nextInt();
    			 this.Cell.get(index).addPubMedID(PubMedID);
    		}
    	}
    	In.close();
    }
    public void loadSwissPubMed(String InFile) throws FileNotFoundException
    {
    	Scanner In = new Scanner(new FileInputStream(InFile));
    	String line = new String();
    	String[] strarr = new String[3];
    	int index = 0;
    	while(In.hasNext())
    	{
    		line = In.nextLine();
    		strarr = line.split("\\s+",3);
    		if (strarr[0].equals("ID"))
    		{
    			line = In.nextLine();
    			strarr = line.split("\\s+",3);
    			String access = strarr[1].substring(0, strarr[1].length()-1);
        		if (!this.AccessIndex.containsKey(access))
        		{
        			this.addProtein(access);
        			this.AccessIndex.put(access, Cell.size() - 1);
        		}
        		index = this.AccessIndex.get(access);
    		}
    		if (strarr[0].equals("RX"))
    		{
    			String PubMed = strarr[1].substring(0, strarr[1].length()-1);
    			String PubMedpre = PubMed.substring(0, 6);
    			PubMed = PubMed.substring(7);
    			if (PubMedpre.equals("PubMed"))
    			{
    				int PubMedID = Integer.parseInt(PubMed);
    				this.Cell.get(index).addPubMedID(PubMedID);
    			}
    		}
    	}
    	In.close();
    }
    
    public void loadSwissIntegratedYear(String InFile) throws FileNotFoundException
    {
    	Scanner In = new Scanner(new FileInputStream(InFile));
    	String line = new String();
    	String[] strarr = new String[3];
    	int index = 0;
    	while(In.hasNext())
    	{
    		line = In.nextLine();
    		strarr = line.split("\\s+",3);
    		if (strarr[0].equals("ID"))
    		{
    			line = In.nextLine();
    			strarr = line.split("\\s+",3);
    			String access = strarr[1].substring(0, strarr[1].length()-1);
        		if (!this.AccessIndex.containsKey(access))
        		{
        			this.addProtein(access);
        			this.AccessIndex.put(access, Cell.size() - 1);
        		}
        		index = this.AccessIndex.get(access);
    		}
    		if (strarr[0].equals("DT"))
    		{
    			line = In.nextLine();
    			line = In.nextLine();
    			int len = strarr[1].length();
    			int year = Integer.parseInt(strarr[1].substring(len-5 ,len-1 ));
    			this.Cell.get(index).setIntegratedYear(year);
    			System.out.println(this.Cell.get(index).getAccess() + "\t" + year);
    		}
    	}
    	In.close();
    }
    
    public void loadIntegratedYear(String InFile) throws FileNotFoundException
    {
    	Scanner In = new Scanner(new FileInputStream(InFile));
    	while (In.hasNext())
    	{
    		String access = In.next();
    		int year = In.nextInt();
    		if (this.AccessIndex.containsKey(access))
    		{
    			int index = this.AccessIndex.get(access);
    			this.Cell.get(index).setIntegratedYear(year);
    		}
    	}
    	In.close();
    }
    
    
    public void loadSwissIntAct(String InFile) throws FileNotFoundException
    {
    	Scanner In = new Scanner(new FileInputStream(InFile));
    	String line = new String();
    	String[] strarr = new String[3];
    	int index = 0;
    	while(In.hasNext())
    	{
    		line = In.nextLine();
    		strarr = line.split("\\s+",3);
    		if (strarr[0].equals("ID"))
    		{
    			line = In.nextLine();
    			strarr = line.split("\\s+",3);
    			String access = strarr[1].substring(0, strarr[1].length()-1);
        		if (!this.AccessIndex.containsKey(access))
        		{
        			this.addProtein(access);
        			this.AccessIndex.put(access, Cell.size() - 1);
        		}
        		index = this.AccessIndex.get(access);
    		}
    		if (line.equals("CC   -!- INTERACTION:"))
    		{
    			line = In.nextLine();
    			while ((line.length()>10) &&    line.substring(0,9).equals("CC       "))
    			{
    				strarr = line.split("\\s+",3);
    				String access = new String();
    				int loc = 0;
    				while (  (loc<strarr[1].length()) && (Character.isLetterOrDigit(strarr[1].charAt(loc))))
    				{
    					
    					access = strarr[1].substring(0, loc+1);
    					loc++;
    				}
    				if (!access.equals("Self"))
    					this.Cell.get(index).addIntAct(access);
    				else
    				{
    					access = this.Cell.get(index).getAccess();
    					this.Cell.get(index).addIntAct(access);
    				}
    				line = In.nextLine();
    			}
    		}
    	}
    	In.close();
    }
    
    public void invokeMethodEveryCell(String methodName) throws NoSuchMethodException, SecurityException, IllegalAccessException, IllegalArgumentException, InvocationTargetException
    {
    	Method method = oneProtein.class.getMethod(methodName);
    	for (oneProtein e:Cell)
    		method.invoke(e);
    }
    
    public void invokeMethodEveryCell(String methodName,String OutFile) throws NoSuchMethodException, SecurityException, IllegalAccessException, IllegalArgumentException, InvocationTargetException, FileNotFoundException
    {
    	PrintWriter Fout = new PrintWriter(new FileOutputStream(OutFile));
    	
    	Method method = oneProtein.class.getMethod(methodName , PrintWriter.class);
    	for (oneProtein e:Cell)
    		method.invoke(e,Fout);
    	Fout.close();
    }
    
    
    public void sortProteinBaseName()
    {
		Comparator<oneProtein> compar = new Comparator<oneProtein>()
		{
			@Override
			public int compare(oneProtein o1, oneProtein o2) {
				// TODO Auto-generated method stub
				return  o1.getName().compareTo(o2.getName());
			}
		};
    	
    	Collections.sort(Cell,compar);
    	UpdateIndex();
    }
    
    public void statYear(int beginYear,int endYear,String OutFile) throws FileNotFoundException
    {
    	PrintWriter Fout = new PrintWriter(new FileOutputStream(OutFile));
    	
    	
    	int size = endYear - beginYear + 1;
    	this.removeGoNotIn(learningOfGO.aGoSet);
    	this.calMFBPCCsize(learningOfGO.aGoSet);
    	int[] statMFOFre = new int[size];
    	int[] statBPOFre = new int[size];
    	int[] statCCOFre = new int[size];
    	
    	int[] statMFOnum = new int[size];
    	int[] statBPOnum = new int[size];
    	int[] statCCOnum = new int[size];
    	
    	for (int i = 1;i<Cell.size();i++)
    	{
    		int index = this.Cell.get(i).getIntegratedYear() - beginYear;
    		if (this.Cell.get(i).getMFOsize()>0) 
    		{
    			statMFOFre[index] += this.Cell.get(i).getMFOsize();
    			statMFOnum[index] += 1;
    		}
    		
    		if (this.Cell.get(i).getBPOsize()>0) 
    		{
    			statBPOFre[index] += this.Cell.get(i).getBPOsize();
    			statBPOnum[index] += 1;
    		}
    		
    		if (this.Cell.get(i).getCCOsize()>0) 
    		{
    			statCCOFre[index] += this.Cell.get(i).getCCOsize();
    			statCCOnum[index] += 1;
    		}	
    	}
    	Fout.println("Year,MFOaverage,MFOsize,BPOaverage,BPOsize,CCOaverage,CCOsize");
    	for (int i = 0;i<size;i++)
    	{
    		Fout.printf("%d , %.2f , %d , %.2f , %d , %.2f , %d\n",i+beginYear,(double)statMFOFre[i]/statMFOnum[i],statMFOnum[i],(double)statBPOFre[i]/statBPOnum[i],statBPOnum[i],
    				(double)statCCOFre[i]/statCCOnum[i],statCCOnum[i]);
    	}
    	
    	Fout.close();
    }
    
    public void analyPubmed(String OutFile) throws FileNotFoundException
    {
    	PrintWriter Fout = new PrintWriter(new FileOutputStream(OutFile));
    	for (oneProtein e:Cell)
    	{
    		Fout.println(e.getAccess() + ","+ e.getMFOsize() + "," + e.getBPOsize() + "," + e.getCCOsize() + "," + e.getPubMedsize());
    	}
    	Fout.close();
    }
    
    public void analyPredScore(GoSet aGoSet)
    {
    	for (int o=0;o<=100;o++)
    	{
    		ArrayList<Pair<Integer, Double>> getPredList = Cell.get(o).getPredList();
    		for (int i = 0;i<getPredList.size()-1;i++)
    			for (int j = i+1;j<getPredList.size();j++)
    			{
    				int GO1 = getPredList.get(i).getFirst();
    				int GO2 = getPredList.get(j).getFirst();
    				double score1 = getPredList.get(i).getSecond();
    				double score2 = getPredList.get(j).getSecond();
    				if (aGoSet.checkGO1_isGO2father(GO1, GO2))
    				{
    					if (score1<score2)
    					{
    						System.out.println("father score is smaller than son's");
    						System.out.printf("%d,%.3f   %d,%.3f", GO1,score1,GO2,score2);
    					}
    				}
    			}
    	}
    }
    
    public void recordliblinearScore()
    {
    	for (oneProtein e:Cell)
    	{
    		e.recordliblinearScore();
    	}
    }
    
    public void recordblastKnnScore()
    {
    	for (oneProtein e:Cell)
    	{
    		e.recordblastKnnScore();
    	}
    }
    
    public void recordblastScore()
    {
    	for (oneProtein e:Cell)
    	{
    		e.recordblastScore();
    	}
    }
    
    public void OutputRanklibTrainFile(String OutFile,char space) throws FileNotFoundException
    {
    	PrintWriter Fout = new PrintWriter(new FileOutputStream(OutFile));
    	int index = 0;
    	for(oneProtein e:Cell)
    	{
    		
    		if (e.getSubAnnSize(space)>0)
    		{
    			index++;
    			e.OutputRanklibFile(Fout,space,index);
    		}
    	}
    	Fout.close();
    }
    
    public void OutputRanklibMeasureFile(String OutFile,char space) throws FileNotFoundException
    {
    	PrintWriter Fout = new PrintWriter(new FileOutputStream(OutFile));
    	int index = 0;
    	for(oneProtein e:Cell)
    	{
    		index++;
    		e.OutputRanklibFile(Fout,space,index);
    	}
    	Fout.close();
    }
    
    public void OutputIDAnnSpecies(String OutFile) throws FileNotFoundException
    {
    	PrintWriter Fout = new PrintWriter(new FileOutputStream(OutFile));
    	for(oneProtein e:Cell)
    		e.OutputIDAnnSpecies(Fout);
    	Fout.close();
    }
    
    
    public void filterHomologProtein() throws FileNotFoundException
    {

    	TreeSet<String> nameSet = new TreeSet<String>();
    	for (int i = 1;i<Cell.size();i++)
    	{
    		String name1 = Cell.get(i-1).getName();
    		String name2 = Cell.get(i).getName();
    		int pos1 = name1.lastIndexOf("_");
    		int pos2 = name2.lastIndexOf("_");
    		if (name1.substring(0, pos1).equals(name2.substring(0, pos2)))
    		{
    			nameSet.add(name1);
    			nameSet.add(name2);
    		}
    	}
        for (int i = 0;i<Cell.size();i++)
        {	
        	oneProtein e = Cell.get(i);
        	String access = e.getName();
        	if (!nameSet.contains(access))
        	{
        		Cell.remove(i);
        		i--;
        	}
        }
        this.UpdateIndex();
    	
    }
    public void sortPredScoreFromHighToLow()
    {
    	for(oneProtein e:Cell)
    		e.sortPredScoreFromHighToLow();
    }
    
    public void calL2RcandidateRecall(String message)
    {
    	double MFOrecall = 0;
    	double BPOrecall = 0;
    	double CCOrecall = 0;
    	for(oneProtein e:Cell)
    	{
    		if (e.getMFOsize()>0) 
    			MFOrecall += e.calL2RcandidateRecall('F');
    		if (e.getBPOsize()>0) 
    			BPOrecall += e.calL2RcandidateRecall('P');
    		if (e.getCCOsize()>0) 
    			CCOrecall += e.calL2RcandidateRecall('C');
    	}
    	MFOrecall /= this.MFOsize;
    	BPOrecall /= this.BPOsize;
    	CCOrecall /= this.CCOsize;
    	System.out.println(message);
    	System.out.println("L2RMFO recall = " + MFOrecall);
    	System.out.println("L2RBPO recall = " + BPOrecall);
    	System.out.println("L2RCCO recall = " + CCOrecall);
    }
    
    public void addRankScore(String FeaFile,String ScoreFile) throws FileNotFoundException
    {
    	Scanner InFea = new Scanner(new FileInputStream(FeaFile));
    	Scanner InScore = new Scanner(new FileInputStream(ScoreFile));
    	String access = new String();
    	while (InFea.hasNext())
    	{
    		String line = InFea.nextLine();
    		double score = InScore.nextDouble();
    		int loc1 = line.lastIndexOf('#');
    		int loc2 = line.lastIndexOf('a');
    		int loc3 = line.lastIndexOf('=');
    		
    		int ann = Integer.parseInt(line.substring(loc3 + 2));
    		access = line.substring(loc1 + 2, loc2 - 1);
    		int index = this.AccessIndex.get(access);
    		this.Cell.get(index).addPred(new Pair<Integer,Double>(ann,score));
    		System.out.println(access +" "+ ann +" "+ score);
    	}
    	InFea.close();
    	InScore.close();
    	
    }
    
    public void prepareL2RFileGenerate(String BlastResult,String ProteinSetDirectory,proteinSet train,String message) throws IOException 
    {
    	this.AddAnnotation(ProteinSetDirectory + "Ann");
		
		
    	this.removeGoNotIn(learningOfGO.aGoSet);
    	this.addGOFather(learningOfGO.aGoSet);
    	this.removeAnnotation(8150,3674,5575);
		
		System.out.println("Begin learningOfGO to rank");
		
		this.addBlastResultBitScore(BlastResult);
		
		this.BlastKnnBaseline(train);
		
		
		
		this.sortPredScoreFromHighToLow();
		this.OutputGOPredScore("../OutFile/" + message + "BlastKnnResult");
		this.addTopK_L2RCandidate(50, 'F');
		this.addTopK_L2RCandidate(100, 'P');
		this.addTopK_L2RCandidate(100, 'C');
		
		this.loadFastaSequence(ProteinSetDirectory + "Seq");
		this.tranSequence2TriSparseFeature();
		this.setLiblinearFeatureFromSparseFeature();
		this.setPredList(train);
		
		this.sortPredListBaseFrequency();
		this.clearPredResult();
		this.libLinearPredict(learningOfGO.modelDir);
		
		
		System.out.println("Begin add Top Liblinear L2R Candidate");
		
		this.sortPredScoreFromHighToLow();
		this.addTopK_L2RCandidate(50, 'F');
		this.addTopK_L2RCandidate(100, 'P');
		this.addTopK_L2RCandidate(100, 'C');
		
		
		System.out.println("Begin record liblinear Score");
		this.recordliblinearScore();
		this.removeLowPred(2000);
		this.OutputGOPredScore("../OutFile/" + message + "LiblinearResult");
		
		this.clearPredResult();
		
		System.out.println("Begin blastKnn again");
		
		this.BlastKnnBaseline(train);
		this.recordblastKnnScore();
		this.clearPredResult();
		
		System.out.println("Begin blast again");
		
		this.addBlastResultBitScore(BlastResult);
		this.blastBaseline(train);
		this.recordblastScore();
		this.clearPredResult();
		
		System.out.println("Begin output ranklib file");
		this.calMFBPCCsize(learningOfGO.aGoSet);
    }
    
    public ArrayList<Integer> getEveryProteinHPOAnnSizeList()
    {
    	ArrayList<Integer> sizeList = new ArrayList<Integer>();
    	for(oneProtein e:Cell)
    	{
    		int size = e.getHPOAnnotationSize();
    		sizeList.add(size);
    	}
    	return sizeList;
    }
    public int [] getEveryHPOFrequency(HPOSet aHPOSet)
    {
    	int [] List = new int[aHPOSet.sizeOfHPO()];
    	for(oneProtein e:Cell)
    	{
    		HashSet<Integer> hpo = e.getHPOAnnotation();
    		for (Integer ann : hpo)
    		{
    			int index = aHPOSet.getIndex(ann);
    			List[index]++;
    		}
    	}
		return List;
    }
    
    public ArrayList<HashSet<Integer>> getInstanceLabel()
    {
    	ArrayList<HashSet<Integer>> InstanceLabel = new ArrayList<HashSet<Integer>>();
    	for (oneProtein e:Cell)
    	{
    		InstanceLabel.add(e.getHPOAnnotation());
    	}
    	return InstanceLabel;
    }
    public ArrayList<ArrayList<Pair<Integer, Double>>> get3merFeatureList()
    {
    	ArrayList<ArrayList<Pair<Integer, Double>>> FeatureList = new ArrayList<ArrayList<Pair<Integer, Double>>>();
    	ArrayList<Double> Fea = new ArrayList<Double>();
    	for (oneProtein e:Cell)
    	{
    		String sequence = e.getSequence();
    		Fea.clear();   
    		ArrayList<Pair<Integer, Double>> Feature = new ArrayList<Pair<Integer, Double>>();   //知道乱引用的厉害了吧~
    		for (int i = 0;i <= 8000 ;i++)	Fea.add(0.0);
    		for (int j = 0; j < sequence.length() - 2; j++) 
    		{
    			String subSeq = sequence.substring(j, j + 3);
    			int dimension = proteinCommon.get3AcidIndex(subSeq);
    			if (dimension<=8000)
    				Fea.set(dimension, Fea.get(dimension) + 1.0);		
    		}
    		for (int i=1;i<=8000;i++)
    		{
    			if (Fea.get(i)>0.001) Feature.add(new Pair<Integer, Double>(i, Fea.get(i)));
    		}
    		FeatureList.add(Feature);
    	}
    	return FeatureList;
    }
    
    public ArrayList<ArrayList<Pair<Integer, Double>>> get2merFeatureList()
    {
    	ArrayList<ArrayList<Pair<Integer, Double>>> FeatureList = new ArrayList<ArrayList<Pair<Integer, Double>>>();
    	ArrayList<Double> Fea = new ArrayList<Double>();
		
    	for (oneProtein e:Cell)
    	{
    		String sequence = e.getSequence();
    		//System.out.println(sequence);
    		Fea.clear();   
    		ArrayList<Pair<Integer, Double>> Feature = new ArrayList<Pair<Integer, Double>>();
    		for (int i = 0;i <= 400 ;i++)	Fea.add(0.0);
    		for (int j = 0; j < sequence.length() - 1; j++) 
    		{
    			String subSeq = sequence.substring(j, j + 2);
    			int dimension = proteinCommon.get2AcidIndex(subSeq);
    			if (dimension<=400)
    				Fea.set(dimension, Fea.get(dimension) + 1.0);		
    		}
    		
    		for (int i=1;i<=400;i++)
    		{
    			if (Fea.get(i)>0.001) Feature.add(new Pair<Integer, Double>(i, Fea.get(i)));
    		}
    		FeatureList.add(Feature);
    	}
    	return FeatureList;
    }
    
    public void addPredScoreList(ArrayList<ArrayList<Pair<Integer,Double>>> predScoreList)
    {
    	for (int i = 0;i<Cell.size();i++)
    	{
    		Cell.get(i).addPredList(predScoreList.get(i));
    	}
    }
    public void addHPOPredScoreList(ArrayList<ArrayList<Pair<Integer,Double>>> predScoreList)
    {
    	for (int i = 0;i<Cell.size();i++)
    	{
    		Cell.get(i).addHPOPredList(predScoreList.get(i));
    	}
    }
    public ArrayList<ArrayList<Pair<Integer,Double>>> getHPOPredScoreList()
    {
    	ArrayList<ArrayList<Pair<Integer,Double>>> list = new ArrayList<ArrayList<Pair<Integer,Double>>>();
    	for (int i = 0;i<Cell.size();i++)
    	{
    		list.add(this.Cell.get(i).getHPOPredList());
    	}
    	return list;
    }
}

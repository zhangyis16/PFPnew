package protein;

import java.io.*;
import java.util.*;

import Main.learning;
import common.*;
import liblinear.Linear;

public class proteinSet
{
    public static HashMap<String,String> MapUniAccess2Name = new HashMap<String,String>();
    public static HashMap<String,String> MapName2UniAccess = new HashMap<String,String>();
    public static HashMap<String,String> MapAccess2UniAccess = new HashMap<String,String>();
    

    public static void compare(proteinSet mySet,proteinSet benchSet,String OutFile) throws FileNotFoundException
    {
    	PrintWriter FoutmySet = new PrintWriter(new FileOutputStream(OutFile+"mySet"));
    	PrintWriter FoutbenchSet = new PrintWriter(new FileOutputStream(OutFile+"benchSet"));
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
    					result.AddProtein(new oneProtein(Name,other.getSubAnnotation(space)));
    			}
    		}
    	}
    	return result;
    }
    public static proteinSet getNewProtein(proteinSet EndSet,proteinSet StartSet)
    {
    	proteinSet result = new proteinSet();
    	for(int i = 0;i < EndSet.size(); i++)
    	{
    		String Access = EndSet.getProteinAccess(i);
    		oneProtein one = EndSet.getProtein(i);
    		if (!StartSet.containProtein(Access))
    		{
    			result.AddProtein(one);			
    		}
    	}
    	return result;
    }
    public static void listCHABIE(proteinSet mySet,proteinSet benchSet,String OutFile) throws FileNotFoundException
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
    	while (In.hasNext())
    	{
    		line = In.nextLine();
    		strarr = line.split("\\s+",3);
    		if (strarr[0].equals("ID")) 
    		{
    			ProteinName = strarr[1];
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
    
    private ArrayList<oneProtein> Cell = new ArrayList<oneProtein>();  
    private HashMap<String,Integer> AccessIndex = new HashMap<String,Integer>();
    private ArrayList<Pair<Integer,Double>> predList = new ArrayList<Pair<Integer,Double>>();
    

    
    public proteinSet() 
    {
        this.Cell = new ArrayList<oneProtein>();
        this.AccessIndex = new HashMap<String,Integer>();
    }
    public proteinSet(String InFile) throws FileNotFoundException
    {
    	AddAnnotation(InFile);
    }
    public oneProtein get(int index)
    {
    	return Cell.get(index);
    }
    public int AnnotationSize()
    {
    	int Sum = 0;
    	for (oneProtein e:Cell)
    	{
    		Sum += e.getAnnotationSize();
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
    	if (learning.aGoSet.containNode(Annotation))
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
    
    public void addFather(GoSet aGoSet)
    {
    	for(oneProtein e:Cell)
    		e.addFather(aGoSet);
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
    		if (!qualifer.equals("NOT")) AddAnnotation(access,gonum);
    	}
    	In.close();
    }
    
    public void AddProtein(oneProtein one)
    {
    	Cell.add(one);
    }
    public void AddProtein(String InFile) throws FileNotFoundException
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
    public void parseSwissAnnotation(String InFile) throws FileNotFoundException 
    {
    	System.out.println("Read Swiss" + InFile + " Now");
    	Scanner In = new Scanner(new FileInputStream(InFile));
    	String line,ProteinName;
    	String Access = new String();
    	String[] strarr = new String[3];
    	while (In.hasNext())
    	{
    		line = In.nextLine();
    		strarr = line.split("\\s+",3);
    		if (strarr[0].equals("ID")) 
    		{
    			ProteinName = strarr[1];
    			line = In.nextLine();
    			strarr = line.split("\\s+",3);
    			Access = strarr[1].substring(0, strarr[1].length()-1);
    		}
    		if (strarr[0].equals("DR") && strarr[1].equals("GO;"))
    		{
    			String[] arr = strarr[2].split("; ",3);
    			int gonum = proteinCommon.GOStr2Int(arr[0]);
    			String evidence = arr[2].substring(0,arr[2].indexOf(':'));
    			if (proteinCommon.Evidence8.contains(evidence)) 
    				AddAnnotation(Access,gonum);
    		}
    	}
    	In.close();
    	System.out.println("Read Swiss" + InFile + " Finish");
    }
    public void parseSwissProteinSequence(String InFile) throws FileNotFoundException
    {
    	Cell.clear();
    	System.out.println("Read Swiss" + InFile + " Now");
    	Scanner In = new Scanner(new FileInputStream(InFile));
    	String line = new String();
    	String Name = new String();
    	String Access = new String();
    	String sequence = new String();
    	String str = "";
    	while (In.hasNext())
    	{
    		line = In.nextLine();
    		String[] strarr = line.split("\\s+",3);
    		if (strarr[0].equals("ID")) 
    		{
    			Name = strarr[1];
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
    		String[] strarr = line.split("\\s+",12);
    		Access = strarr[0];
    		Access2 = strarr[1];
    		double bitscore = Double.valueOf(strarr[11]);
    		int index = AccessIndex.get(Access);
    		Cell.get(index).addBlastResult(Access2, bitscore);
    		if (bitscore > maxbitscore && !Access2.equals(Access)) maxbitscore = bitscore;
    	}
    	In.close();
    	System.out.println("MaxBitScore = " + maxbitscore);
    }
    public void addBlastResultSimility(String InFile) throws FileNotFoundException
    {
    	Scanner In = new Scanner(new FileInputStream(InFile));
    	String line,Access,Access2;
    	double maxbitscore = 0;
    	while (In.hasNext())
    	{
    		line = In.nextLine();
    		String[] strarr = line.split("\\s+",12);
    		Access = strarr[0];
    		Access2 = strarr[1];
    		double bitscore = Double.valueOf(strarr[2]);
    		int index = AccessIndex.get(Access);
    		Cell.get(index).addBlastResult(Access2, bitscore);
    		if (bitscore > maxbitscore && !Access2.equals(Access)) maxbitscore = bitscore;
    	}
    	In.close();
    }
    public void clear()
    {
    	Cell.clear();
    	AccessIndex.clear();
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
    
    public void eraserProteinOnly5515()
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
        
    public void evalution() throws FileNotFoundException
    {
    	ArrayList<Pair<Double,Double>> PreReCallPair = new ArrayList<Pair<Double,Double>>();
    	
    	ArrayList<HashSet<Integer>> answer = new ArrayList();
    	ArrayList<ArrayList<Pair<Integer,Double>>> predictor = new ArrayList();
    	
    	this.removeAnnotation(8150,3674,5575);
    	Pair<Double,Double> result;
    	getAnswerPredictor(answer,predictor);
    	result = evalution.GetFMeasureMax(answer, predictor,PreReCallPair);
    	
    	basic.OutputPair("../OutFile/AUPRlist/ALLAUPRlist",PreReCallPair);
    	
    	System.out.println("Fmax = " + result.getFirst() +" Cut = "+  result.getSecond());  
    	
    	proteinSet MFOsubset = this.getSubAnnotation('F');

    	MFOsubset.getAnswerPredictor(answer,predictor);
    	result = evalution.GetFMeasureMax(answer, predictor,PreReCallPair);
    	basic.OutputPair("../OutFile/AUPRlist/MFOAUPRlist",PreReCallPair);
    	System.out.println("MFO :: Fmax = " + result.getFirst() +" Cut = "+  result.getSecond());
    	
    	proteinSet BPOsubset = this.getSubAnnotation('P');
    	BPOsubset.getAnswerPredictor(answer,predictor);
    	result = evalution.GetFMeasureMax(answer, predictor,PreReCallPair);
    	basic.OutputPair("../OutFile/AUPRlist/BPOAUPRlist",PreReCallPair);
    	System.out.println("BPO :: Fmax = " + result.getFirst() +" Cut = "+  result.getSecond());
    	
    	
    	proteinSet CCOsubset = this.getSubAnnotation('C');
    	CCOsubset.getAnswerPredictor(answer,predictor);
    	result = evalution.GetFMeasureMax(answer, predictor,PreReCallPair);
    	basic.OutputPair("../OutFile/AUPRlist/CCOAUPRlist",PreReCallPair);
    	System.out.println("CCO :: Fmax = " + result.getFirst() +" Cut = "+  result.getSecond());
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
    	
    	for (int i = 0;i<Cell.size();i++)
    	{	
    		if (Cell.get(i).getAnnotationSize()==0)
    		{
    			Cell.remove(i);
    			i--;
    		}
    	}
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
    
    public ArrayList<oneProtein> getCell()
    {
    	return Cell;
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
    	int index = AccessIndex.get(access);
    	ArrayList<Integer> ann = new ArrayList<Integer>(Cell.get(index).getAnnotation());
    	return ann;
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
    public proteinSet getSubAnnotation(char space)
    {
    	proteinSet result = new proteinSet();
    	for(oneProtein e:Cell)
    	{
    		oneProtein aProtein = e.getSubProtein(space);
    		if (!(aProtein.getAnnotationSize() == 0))
    		result.AddProtein(aProtein);
    	}
    	return result;
    }
    public void loadFastaSequence(String InFile) throws FileNotFoundException
    {
    	Scanner In = new Scanner(new FileInputStream(InFile));
    	String Name = new String();
    	String Access = new String();
    	String sequence = new String();
    	while (In.hasNext())
    	{
    		String line = In.nextLine();
    		if (line.substring(0,1).equals(">")) 
    		   Access = line.substring(1);  //System.out.println(Access);
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

    public void naiveBaseline(proteinSet train) throws FileNotFoundException
    {
    	ArrayList<Pair<Integer,Double>> predList = train.getNaiveList();
    	for(oneProtein e:Cell)
    	{
    		e.setPredList(predList);
    	}

    }
    public void blastBaseline(proteinSet train) throws FileNotFoundException
    {
    	for(oneProtein e:Cell)
    	{
    		e.setBlastPred(train);
    	}
    }
    public void GOtchaBaseline(proteinSet train) throws FileNotFoundException
    {
    	for(oneProtein e:Cell)
    	{
    		e.setGOtcha(train);
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
    	for (Pair<Integer,Double> pair:predList)
    	{
    		pair.setSecond(pair.getSecond()/Cell.size());
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
    public void OutputAnnotationList(String OutFile) throws FileNotFoundException
    {
    	PrintWriter Fout = new PrintWriter(new FileOutputStream(OutFile));
    	for(oneProtein e:Cell)
    	{
    		e.OutputAnnotationList(Fout,',');
    		e.outputSequenceVector(Fout);
    		Fout.println();
    	}
    	Fout.close();
    }
    public void OutputAnnotationSpace(String OutFile) throws FileNotFoundException
    {
    	PrintWriter Fout = new PrintWriter(new FileOutputStream(OutFile));
    	for(oneProtein e:Cell)
    		e.OutputAnnotationSpace(Fout);
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
    public void OutputFastaSequence(String OutFile) throws FileNotFoundException
    {
    	PrintWriter Fout = new PrintWriter(new FileOutputStream(OutFile));
    	for(oneProtein e:Cell)
    		e.OutputFastaSequence(Fout);
    	Fout.close();
    }
    public void OutputPred(String OutFile) throws FileNotFoundException
    {
    	PrintWriter Fout = new PrintWriter(new FileOutputStream(OutFile));
    	for (oneProtein e:Cell)
    	{
    		e.OutputPred(Fout);
    	}
    	Fout.close();
    }
    
    public void addPredResult(String InFile) throws FileNotFoundException
    {
    	Scanner In = new Scanner(new FileInputStream(InFile));
    	for (oneProtein e:Cell)
    	{
    		e.addPredResult(In);
    	}
    	In.close();
    }
    public void OutputProtein(String OutFile) throws FileNotFoundException
    {
    	PrintWriter Fout = new PrintWriter(new FileOutputStream(OutFile));
    	for(oneProtein e:Cell)
    		Fout.println(e.getAccess());
    	Fout.close();
    }
    
    public void outputLibTrainFile(String OutFile) throws FileNotFoundException
    {
    	this.tranSequence2Vector();
    	PrintWriter Fout = new PrintWriter(new FileOutputStream(OutFile));
    	for(oneProtein e:Cell)
    	{
    		e.OutputAnnotationNum(Fout,',');
    		e.outputSequenceVector(Fout);
    		Fout.println();
    	}
    	Fout.close();
    }
    
    public void outputLibTrainFile(String OutFile,int Ann) throws FileNotFoundException
    {
    	PrintWriter Fout = new PrintWriter(new FileOutputStream(OutFile));
    	for(oneProtein e:Cell)
    	{
    		if (e.containAnnotation(Ann)) Fout.print(1); else Fout.print(-1);
    		e.outputSequenceVector(Fout);
    		Fout.println();
    	}
    	Fout.close();
    }
    
    public void tranSequence2Vector()
    {
    	for(oneProtein e:Cell)
    		e.tranSequence2Vector();
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
    	ArrayList<Pair<Integer,Double>> predList = train.getNaiveList();
		Collections.sort(predList);
		this.predList = predList;
    }
    public void setLiblinearFeature()
    {
    	for (oneProtein e:Cell)
    		e.setLiblinearFeature();
    }
    public void SortProtein()
    {
    	Collections.sort(Cell);
    	UpdateIndex();
    }
    public void removeAnnotation(int... args)
    {
    	for (oneProtein e:Cell)
    		e.removeAnnotation(args);
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
    			char sp = learning.aGoSet.getSpace(ann);
    			if (sp == 'F') MFdepth[learning.aGoSet.getMindepth(ann)]++;
    			if (sp == 'P') BPdepth[learning.aGoSet.getMindepth(ann)]++;
    			if (sp == 'C') CCdepth[learning.aGoSet.getMindepth(ann)]++;
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
    public void statGoFrequency() throws FileNotFoundException
    {
    	int[] statGoFre = new int[learning.aGoSet.size()+1];
    	int[] statMFOFre = new int[learning.aGoSet.size()+1];
    	int[] statBPOFre = new int[learning.aGoSet.size()+1];
    	int[] statCCOFre = new int[learning.aGoSet.size()+1];
    	 for (oneProtein e:Cell)
    	 {
    		 HashSet<Integer> AnnSet = e.getAnnotation();
    		 for (int ann:AnnSet)
    		 {
    			 statGoFre[learning.aGoSet.getIndex(ann)]++;
    			 if (learning.aGoSet.getSpace(ann) == 'F') statMFOFre[learning.aGoSet.getIndex(ann)]++;
    			 if (learning.aGoSet.getSpace(ann) == 'P') statBPOFre[learning.aGoSet.getIndex(ann)]++;
    			 if (learning.aGoSet.getSpace(ann) == 'C') statCCOFre[learning.aGoSet.getIndex(ann)]++;
    		 }
    	 }
    	 common.statiscit.statDistribution(statGoFre, "GOdistribution.csv");
    	 common.statiscit.statDistribution(statMFOFre, "MFOdistribution.csv");
    	 common.statiscit.statDistribution(statBPOFre, "BPOdistribution.csv");
    	 common.statiscit.statDistribution(statCCOFre, "CCOdistribution.csv");
    }
    public void statSpeciesGoNum(TreeSet<String> spSet, String OutFile) throws FileNotFoundException
    {
    	PrintWriter Fout = new PrintWriter(new FileOutputStream(OutFile));
    	int spNum = spSet.size();
    	double[] GOstat = new double[spNum + 2];
    	double[] MFOstat = new double[spNum + 2];
    	double[] BPOstat = new double[spNum + 2];
    	double[] CCOstat = new double[spNum + 2];
    	
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
    		String Name = e.getName();
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
    		MFOstat[spIndex] += MFOlabel;
    		BPOstat[spIndex] += BPOlabel;
    		CCOstat[spIndex] += CCOlabel;
    		GOnum[spIndex]++;
    		if (MFOlabel>0)  MFOnum[spIndex]++;
    		if (BPOlabel>0)  BPOnum[spIndex]++;
    		if (CCOlabel>0)  CCOnum[spIndex]++;
    	}
    	Index = 1;
    	Fout.println("species,GO,MFO,BPO,CCO");
    	for (String sp :spSet)
    	{
    		Fout.println(String.format(sp + ",%.1f,%.1f,%.1f,%.1f",GOstat[Index]/GOnum[Index],MFOstat[Index]/MFOnum[Index]
    				,BPOstat[Index]/BPOnum[Index],CCOstat[Index]/CCOnum[Index]));
    		Index++;
    	}
		Fout.println(String.format("Others,%.1f,%.1f,%.1f,%.1f",GOstat[Index]/GOnum[Index],MFOstat[Index]/MFOnum[Index]
				,BPOstat[Index]/BPOnum[Index],CCOstat[Index]/CCOnum[Index]));
		
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
    		System.out.println("Predict Label" + this.predList.get(count).getFirst());
    		int labelNum = predList.get(count).getFirst();
    		liblinear.Model model = Linear.loadModel(new File(Direction + "Label" + labelNum + ".model"));
    		
    		for(oneProtein e:Cell)
    		{
    			e.libLinearPredict(labelNum,model);
    		}
    	}
    }
    public class MultiThreadLiblinearPred implements Runnable
    {
		private volatile int count;
		private String Direction  = new String();
		public MultiThreadLiblinearPred(int countNum,String aDirection)
		{
			count = countNum;
			Direction = aDirection;
		}
		@Override
		public void run() 
		{
			// TODO Auto-generated method stub
			
		}
    	
    }
    
    public void libLinearTrain(String Direction) throws IOException
    {
    	MultiThreadTrain my = new proteinSet.MultiThreadTrain(this.predList.size()-1,Direction);
    	Thread t1 = new Thread(my,"thread 1");
		Thread t2 = new Thread(my,"thread 2");
		Thread t3 = new Thread(my,"thread 3");
		Thread t4 = new Thread(my,"thread 4");
		t1.start();
		t2.start();
		t3.start();
		t4.start();
		while((t1.getState() != Thread.State.TERMINATED) &&
				(t2.getState() != Thread.State.TERMINATED) &&
				(t3.getState() != Thread.State.TERMINATED) &&
				(t4.getState() != Thread.State.TERMINATED) 
				);    	
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
					libLinearTrain(predList.get(count--).getFirst(),Direction);
				} catch (IOException e) {
					// TODO Auto-generated catch block
					e.printStackTrace();
				}
				System.out.println(Thread.currentThread().getName()+ "now train label" + this.get());
			}
		}
		synchronized public int get(){
			return this.count--;
		}
		
	}
    public void libLinearTrain(int labelNum,String Direction) throws IOException
    {
    	int Ann = labelNum;
        int max_index = 8003;
        double bias = -1;
        
        liblinear.Parameter param   = new liblinear.Parameter(liblinear.SolverType.L2R_LR, 1, 1000, 0.05);
        liblinear.Problem   prob    = new liblinear.Problem();
        
        prob.bias = bias;
        prob.l = this.size();
        prob.n = max_index;
        if (bias >= 0) {
            prob.n++;
        }
        prob.x = new liblinear.Feature[prob.l][];
        prob.y = new double[prob.l];
        List<liblinear.Feature[]> vx = new ArrayList<liblinear.Feature[]>();
        for (int i = 0; i < prob.l; i++)
        {
        	ArrayList<Pair<Integer,Double>> feature = this.Cell.get(i).getFeature();
        	liblinear.Feature[] x = new liblinear.Feature[feature.size()];
        	for (int j = 0;j < feature.size() ;j++)
        		x[j] = new liblinear.FeatureNode(feature.get(j).getFirst(), feature.get(j).getSecond());
        	vx.add(x);
        	prob.x[i] = vx.get(i);
        }
        for (int i = 0; i < prob.l; i++)
        	if (this.Cell.get(i).containAnnotation(Ann))
            prob.y[i] = 1; else prob.y[i] = -1;
        
    	liblinear.Model model = liblinear.Linear.train(prob, param);
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
    	for (oneProtein e:Cell)
    	{
    		e.removeRedunAnn(aGoSet);
    	}
    }
}

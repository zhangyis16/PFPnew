package protein;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.*;

import common.Pair;
public class proteinCommon
{
	public static Map<Integer,Double> NaivePred = new HashMap<Integer,Double>();
	
	public static ArrayList<ArrayList<Pair<Integer,Double>>> TriAcidSimiliar = 
			new ArrayList<ArrayList<Pair<Integer,Double>>>();
	
	public static void outputFasta(PrintWriter Fout,String name,String sequence)
	{
		Fout.println(">" + name);
		int count = (sequence.length() - 1) / 60;
		for (int i = 0; i <= count; i++)
			Fout.println(sequence.substring(i * 60, Math.min((i + 1) * 60, sequence.length())));
		Fout.println();
	}
	public static double blastTwoProteinSeq(String seq1,String seq2) throws IOException, InterruptedException
	{
		PrintWriter Fout1 = new PrintWriter(new FileOutputStream("Seq1.fasta"));
		proteinCommon.outputFasta(Fout1, "name1", seq1);
		Fout1.close();
		
		PrintWriter Fout2 = new PrintWriter(new FileOutputStream("Seq2.fasta"));
		proteinCommon.outputFasta(Fout2, "name2", seq2);
		Fout2.close();
		
		String cmd = "makeblastdb -in Seq1.fasta -parse_seqids -hash_index -dbtype prot";
		Process process1 = Runtime.getRuntime().exec(cmd);
		process1.waitFor();
		
		cmd = "blastp -task blastp -query Seq2.fasta -db Seq1.fasta "
				+ "-out blastResult -outfmt 6";
		Process process2 = Runtime.getRuntime().exec(cmd);
		process2.waitFor();

		Scanner In = new Scanner(new FileInputStream("blastResult"));
		double score = 0.0;
		if (In.hasNext())
		{
			String line = In.nextLine();
			String[] strarr = line.split("\t",4);
			score = Double.valueOf(strarr[2]);
		}else
			score = 0.0;
		
		In.close();
		
		return score;
	}
	public static double CalTriAcidSimiliar(int i,int j ,int SimiliarMatrix[][])
	{
		int i1 = i/400;
		int i2 = (i/20)%20;
		int i3 = i%20;
		int j1 = j/400;
		int j2 = (j/20)%20;
		int j3 = j%20;
		double s12 = (double)(SimiliarMatrix[i1][j1]
					+SimiliarMatrix[i2][j2]
					+SimiliarMatrix[i3][j3]);
		double s1 = (SimiliarMatrix[i1][i1]
					+ SimiliarMatrix[i2][i2]
					+ SimiliarMatrix[i3][i3]);
		double s2 = (SimiliarMatrix[j1][j1]
				+ SimiliarMatrix[j2][j2]
				+ SimiliarMatrix[j3][j3]);		

		double sim = s12/ Math.sqrt(s1*s2);
		return sim;
	}
	public static double CalTriAcidSimiliar(int i,int j)
	{
		return CalTriAcidSimiliar(i,j ,proteinCommon.Blosum62);
	}
	public static void CalTriAcidSimiliar(double cutoff)
	{
		TriAcidSimiliar.add(new ArrayList<Pair<Integer,Double>>());
		int count = 0;
		for (int i=1;i<=8000;i++)
		{
			ArrayList<Pair<Integer,Double>> Ar = new ArrayList<Pair<Integer,Double>>();
			for (int j=1;j<=8000;j++)
			{
				double sim = CalTriAcidSimiliar(i,j);
				if (sim>cutoff)
				{
					Ar.add(new Pair<Integer,Double>(j,sim));
					count++;
				}
			}
			TriAcidSimiliar.add(Ar);
		}
		System.out.println("cutoff = " + cutoff + " we find " + count + " similiar pair");
	}
	
	public final static ArrayList<String> ArrayCAFA2Species_List = new ArrayList<String>(Arrays.asList
	("_ARATH","_BACSU","_DANRE","_DICDI","_DROME","_ECOLI","_HALS3","_HALVD","_HELPY",
     "_HUMAN","_IGNH4","_METJA","_MOUSE","_MYCGE","_NITMS","_PSEAE","_PSEPK","_PSESM",
     "_PYRFU","_RAT",  "_SALCH","_SALTY","_SCHPO","_STRPN","_SULSO","_XENLA","_YEAST")); 
	
	
	public final static TreeSet<String> SetCAFA2Species = new TreeSet<String>(Arrays.asList
	("_ARATH","_BACSU","_DANRE","_DICDI","_DROME","_ECOLI","_HALS3","_HALVD","_HELPY",
     "_HUMAN","_IGNH4","_METJA","_MOUSE","_MYCGE","_NITMS","_PSEAE","_PSEPK","_PSESM",
     "_PYRFU","_RAT",  "_SALCH","_SALTY","_SCHPO","_STRPN","_SULSO","_XENLA","_YEAST"));
	
	public final static TreeSet<String> SetMySpecies = new TreeSet<String>(Arrays.asList
	("_HUMAN","_MOUSE","_ARATH","_YEAST","_SCHPO","_RAT","_ECOLI","_DROME",
     "_DICDI","_DANRE","_XENLA"));
	
	public final static  HashSet<String> EvidenceCode = new HashSet<String>() {/**
		 * 
		 */
		private static final long serialVersionUID = 1L;

	{  
	       add("EXP");  add("IDA");  add("IMP");  add("IPI");
	       add("IGI");  add("IEP");  add("TAS");  add("IC");
	}};  
	
	public final static  Map<Character,Integer> AcidIndex = new HashMap<Character,Integer>() {/**
		 * 
		 */
		private static final long serialVersionUID = 1L;

	{
	       put('A',0);	       put('R',1);  
	       put('N',2);	       put('D',3);
	       put('C',4);	       put('Q',5);
	       put('E',6);	       put('G',7);
	       put('H',8);	       put('I',9);
	       put('L',10);	       put('K',11);
	       put('M',12);	       put('F',13);
	       put('P',14);	       put('S',15);
	       put('T',16);	       put('W',17);
	       put('Y',18);	       put('V',19); 
	       put('B',20);	       put('Z',21);
	       put('X',22);
	}};  
	public final static  Map<Integer,Character> Index2Acid = new HashMap<Integer,Character>() {/**
		 * 
		 */
		private static final long serialVersionUID = 1L;

	{
	       put(0,'A');	       put(1,'R');  
	       put(2,'N');	       put(3,'D');
	       put(4,'C');	       put(5,'Q');
	       put(6,'E');	       put(7,'G');
	       put(8,'H');	       put(9,'I');
	       put(10,'L');	       put(11,'K');
	       put(12,'M');	       put(13,'F');
	       put(14,'P');	       put(15,'S');
	       put(16,'T');	       put(17,'W');
	       put(18,'Y');	       put(19,'V'); 
	       put(20,'B');	       put(21,'Z');
	       put(22,'X');
	}}; 
	
	public final static int Blosum62[][] ={
			 { 4,-1,-2,-2, 0,-1,-1, 0,-2,-1,-1,-1,-1,-2,-1, 1, 0,-3,-2, 0,-2,-1, 0,-4},
			 {-1, 5, 0,-2,-3, 1, 0,-2, 0,-3,-2, 2,-1,-3,-2,-1,-1,-3,-2,-3,-1, 0,-1,-4},
			 {-2, 0, 6, 1,-3, 0, 0, 0, 1,-3,-3, 0,-2,-3,-2, 1, 0,-4,-2,-3, 3, 0,-1,-4},
			 {-2,-2, 1, 6,-3, 0, 2,-1,-1,-3,-4,-1,-3,-3,-1, 0,-1,-4,-3,-3, 4, 1,-1,-4},
			 { 0,-3,-3,-3, 9,-3,-4,-3,-3,-1,-1,-3,-1,-2,-3,-1,-1,-2,-2,-1,-3,-3,-2,-4},
			 {-1, 1, 0, 0,-3, 5, 2,-2, 0,-3,-2, 1, 0,-3,-1, 0,-1,-2,-1,-2, 0, 3,-1,-4},
			 {-1, 0, 0, 2,-4, 2, 5,-2, 0,-3,-3, 1,-2,-3,-1, 0,-1,-3,-2,-2, 1, 4,-1,-4},
			 { 0,-2, 0,-1,-3,-2,-2, 6,-2,-4,-4,-2,-3,-3,-2, 0,-2,-2,-3,-3,-1,-2,-1,-4},
			 {-2, 0, 1,-1,-3, 0, 0,-2, 8,-3,-3,-1,-2,-1,-2,-1,-2,-2, 2,-3, 0, 0,-1,-4},
			 {-1,-3,-3,-3,-1,-3,-3,-4,-3, 4, 2,-3, 1, 0,-3,-2,-1,-3,-1, 3,-3,-3,-1,-4},
			 {-1,-2,-3,-4,-1,-2,-3,-4,-3, 2, 4,-2, 2, 0,-3,-2,-1,-2,-1, 1,-4,-3,-1,-4},
			 {-1, 2, 0,-1,-3, 1, 1,-2,-1,-3,-2, 5,-1,-3,-1, 0,-1,-3,-2,-2, 0, 1,-1,-4},
			 {-1,-1,-2,-3,-1, 0,-2,-3,-2, 1, 2,-1, 5, 0,-2,-1,-1,-1,-1, 1,-3,-1,-1,-4},
			 {-2,-3,-3,-3,-2,-3,-3,-3,-1, 0, 0,-3, 0, 6,-4,-2,-2, 1, 3,-1,-3,-3,-1,-4},
			 {-1,-2,-2,-1,-3,-1,-1,-2,-2,-3,-3,-1,-2,-4, 7,-1,-1,-4,-3,-2,-2,-1,-2,-4},
			 { 1,-1, 1, 0,-1, 0, 0, 0,-1,-2,-2, 0,-1,-2,-1, 4, 1,-3,-2,-2, 0, 0, 0,-4},
			 { 0,-1, 0,-1,-1,-1,-1,-2,-2,-1,-1,-1,-1,-2,-1, 1, 5,-2,-2, 0,-1,-1, 0,-4},
			 {-3,-3,-4,-4,-2,-2,-3,-2,-2,-3,-2,-3,-1, 1,-4,-3,-2,11, 2,-3,-4,-3,-2,-4},
			 {-2,-2,-2,-3,-2,-1,-2,-3, 2,-1,-1,-2,-1, 3,-3,-2,-2, 2, 7,-1,-3,-2,-1,-4},
			 { 0,-3,-3,-3,-1,-2,-2,-3,-3, 3, 1,-2, 1,-1,-2,-2, 0,-3,-1, 4,-3,-2,-1,-4},
			 {-2,-1, 3, 4,-3, 0, 1,-1, 0,-3,-4, 0,-3,-3,-2, 0,-1,-4,-3,-3, 4, 1,-1,-4},
			 {-1, 0, 0, 1,-3, 3, 4,-2, 0,-3,-3, 1,-1,-3,-1, 0,-1,-3,-2,-2, 1, 4,-1,-4},
			 { 0,-1,-1,-1,-2,-1,-1,-1,-1,-1,-1,-1,-1,-1,-2, 0, 0,-2,-1,-1,-1,-1,-1,-4},
			 {-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4, 1}	
	} ;
	public static boolean is20Acid(char ch)
	{
		if (AcidIndex.containsKey(ch))
		{
			int index = AcidIndex.get(ch);
			if ((index>=0) && (index<20))
					return true;
		}
		return false;
			
	}
	public static int get3AcidIndex(String str)  //正常的三个氨基酸序列返回1~8000
	{
		int num = 0;
		char c0 = str.charAt(0);
		char c1 = str.charAt(1);
		char c2 = str.charAt(2);
		if (AcidIndex.containsKey(c0) &&  AcidIndex.containsKey(c1) && AcidIndex.containsKey(c2))
		{
			num = 1 + AcidIndex.get(c0) + AcidIndex.get(c1)*20 + AcidIndex.get(c2)*400;
			if (is20Acid(c0)  &&  is20Acid(c1) && is20Acid(c2) )
			return num;
		}
		if (c0 == 'B') return 8001;
		if (c0 == 'U') return 8002;
		if (c0 == 'Z') return 8003;
		return 0;
		
	}
	
	public static int get2AcidIndex(String str)  //正常的二个氨基酸序列返回1~400
	{
		int num = 0;
		char c0 = str.charAt(0);
		char c1 = str.charAt(1);
		if (AcidIndex.containsKey(c0) &&  AcidIndex.containsKey(c1) )
		{
			num = 1 + AcidIndex.get(c0) + AcidIndex.get(c1)*20;
			if (is20Acid(c0)  &&  is20Acid(c1))
				return num;
		}
		if (c0 == 'B') return 401;
		if (c0 == 'U') return 402;
		if (c0 == 'Z') return 403;
		return 404;
		
	}
	public static String getSpecies(String access)
	{
		String ProteinName;
		if (!access.contains("_"))
		{	
			if (proteinSet.MapUniAccess2Name.containsKey(access))
				ProteinName = proteinSet.MapUniAccess2Name.get(access);
			else 
			{
				ProteinName = "_NotKnow";
				System.out.println(access+"don't have map Species");
			}
		}
		else ProteinName = access;
		return ProteinName.substring(ProteinName.lastIndexOf("_"));
	}
	public static int GOStr2Int(String GOStr)
	{
		if (GOStr.length() == 10)
		{
			String str = GOStr.substring(3,10);
			int num = 0;
			for (int i=0;i<str.length();i++)
			{
				char  item =  str.charAt(i);
				num = num * 10 + item - '0';
			}
			return num;
		}
		else
			return Integer.parseInt(GOStr);
	}
	
	public static int HPStr2Int(String HPStr)
	{
		if (HPStr.length() == 10)
		{
			String str = HPStr.substring(3,10);
			int num = 0;
			for (int i=0;i<str.length();i++)
			{
				char  item =  str.charAt(i);
				num = num * 10 + item - '0';
			}
			return num;
		}
		else
			return Integer.parseInt(HPStr);
	}
	
	public static String GOInt2Str(int GONum)
	{
		String result = "GO:";
		String GO = Integer.toString(GONum);
		for (int i=0;i<7-GO.length();i++) result += "0";
		result += GO;
		return result;
	}
	
	public static String HPInt2Str(int HPNum)
	{
		String result = "HP:";
		String HP = Integer.toString(HPNum);
		for (int i=0;i<7-HP.length();i++) result += "0";
		result += HP;
		return result;
	}
	
	public static boolean isGO_ID(String GOStr)
	{
		String str = GOStr.substring(0,3);
		if ((GOStr.length() == 10) && (str.equals("GO:"))) return true; else return false;
	}
	
	public static boolean isHPO_ID(String HPOStr)
	{
		String str = HPOStr.substring(0,3);
		if ((HPOStr.length() == 10) && (str.equals("HP:"))) return true; else return false;
	}

}

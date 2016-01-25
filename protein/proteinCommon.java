package protein;
import java.io.*;
import java.util.*;
public class proteinCommon
{
	public static Map<Integer,Double> NaivePred = new HashMap<Integer,Double>();
	
	public final static ArrayList<String> ArrayCAFA2Species_List = new ArrayList(Arrays.asList
	("_ARATH","_BACSU","_DANRE","_DICDI","_DROME","_ECOLI","_HALS3","_HALVD","_HELPY",
     "_HUMAN","_IGNH4","_METJA","_MOUSE","_MYCGE","_NITMS","_PSEAE","_PSEPK","_PSESM",
     "_PYRFU","_RAT",  "_SALCH","_SALTY","_SCHPO","_STRPN","_SULSO","_XENLA","_YEAST")); 
	
	
	public final static TreeSet<String> SetCAFA2Species = new TreeSet(Arrays.asList
	("_ARATH","_BACSU","_DANRE","_DICDI","_DROME","_ECOLI","_HALS3","_HALVD","_HELPY",
     "_HUMAN","_IGNH4","_METJA","_MOUSE","_MYCGE","_NITMS","_PSEAE","_PSEPK","_PSESM",
     "_PYRFU","_RAT",  "_SALCH","_SALTY","_SCHPO","_STRPN","_SULSO","_XENLA","_YEAST"));
	
	public final static TreeSet<String> SetMySpecies = new TreeSet(Arrays.asList
	("_HUMAN","_MOUSE","_ARATH","_YEAST","_SCHPO","_RAT","_ECOLI","_DROME","_MYCTU",
     "_CAEEL","_DICDI","_DANRE","_CANAL","_XENLA"));
	
	public final static  Set<String> Evidence8 = new HashSet<String>() {{  
	       add("EXP");  add("IDA");  add("IPI");  add("IMP");  
	       add("IGI");  add("IEP");  add("TAS");  add("IC");
	}};  
	
	final static  Map<Character,Integer> AcidIndex = new HashMap<Character,Integer>() {{
	       put('A',0);	       put('C',1);  
	       put('D',2);	       put('E',3);
	       put('F',4);	       put('G',5);
	       put('H',6);	       put('I',7);
	       put('K',8);	       put('L',9);
	       put('M',10);	       put('N',11);
	       put('P',12);	       put('Q',13);
	       put('R',14);	       put('S',15);
	       put('T',16);	       put('V',17);
	       put('W',18);	       put('Y',19);     
	}};  
	
	public static int get3AcidIndex(String str)
	{
		int num = 0;
		char c0 = str.charAt(0);
		char c1 = str.charAt(1);
		char c2 = str.charAt(2);
		if (AcidIndex.containsKey(c0) &&  AcidIndex.containsKey(c1) && AcidIndex.containsKey(c2))
		{
			num = 1 + AcidIndex.get(c0) + AcidIndex.get(c1)*20 + AcidIndex.get(c2)*400;
			return num;
		}
		else
		{
			if (c0 == 'B') return 8001;
			if (c0 == 'U') return 8002;
			if (c0 == 'Z') return 8003;
			return 0;
		}
	}
	public static String getSpecies(String access)
	{
		String ProteinName;
		if (!access.contains("_"))
		ProteinName = proteinSet.MapUniAccess2Name.get(access);
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
	public static String GOInt2Str(int GONum)
	{
		String result = "GO:";
		String GO = Integer.toString(GONum);
		for (int i=0;i<7-GO.length();i++) result += "0";
		result += GO;
		return result;
	}
	public static boolean isGO_ID(String GOStr)
	{
		String str = GOStr.substring(0,3);
		if ((GOStr.length() == 10) && (str.equals("GO:"))) return true; else return false;
	}

}

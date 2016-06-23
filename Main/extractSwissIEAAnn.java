package Main;

import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.PrintWriter;
import java.util.Scanner;

import protein.proteinCommon;

public class extractSwissIEAAnn {

	public static void main(String[] args) throws FileNotFoundException {
		// TODO Auto-generated method stub
		String InFile = args[0];
    	Scanner In = new Scanner(new FileInputStream(InFile));
    	PrintWriter Fout = new PrintWriter(new FileOutputStream("SwissIEA_Ann"));
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
    			if (!proteinCommon.EvidenceCode.contains(evidence)) 
    				Fout.println(Access + "\t" + gonum);
    		}
    	}
    	In.close();
    	Fout.close();
	}

}

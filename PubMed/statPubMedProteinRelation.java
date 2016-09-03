package PubMed;

import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Scanner;

import common.Pair;

public class statPubMedProteinRelation {

	public static void main(String[] args) throws FileNotFoundException {
		// TODO Auto-generated method stub
		Scanner In = new Scanner(new FileInputStream("../InFile/Swiss/PubMed201401"));
		PrintWriter Fout =   new PrintWriter("PubMedProteinRelationShip");
		String access = new String();
		int count = 0;
		int PubMed = 0;
		ArrayList<ArrayList<String>> accList = new ArrayList<ArrayList<String>>();
		ArrayList<Integer> PubMedList = new ArrayList<Integer>();
		HashMap<Integer,Integer> PubMedIndex = new HashMap<Integer,Integer>();
		while (In.hasNext())
		{
			access = In.next();
			count = In.nextInt();
			
			for (int i = 1;i<=count;i++)
			{
				PubMed = In.nextInt();
				if (PubMedIndex.containsKey(PubMed))
				{
					int index = PubMedIndex.get(PubMed);
					accList.get(index).add(access);
				}else
				{
					PubMedList.add(PubMed);
					accList.add(new ArrayList<String>());
					
					PubMedIndex.put(PubMed, PubMedList.size() - 1);
					
					int index = PubMedList.size() - 1;
					accList.get(index).add(access);
				}	
			}
		}
		
		for (int i = 0;i<PubMedList.size()-1;i++)
		{
			Fout.print(PubMedList.get(i) + "," + accList.get(i).size() + ",");
			for (int j = 0;j<accList.get(i).size();j++)
				Fout.print(" " + accList.get(i).get(j));
			Fout.println();
		}
		In.close();
		Fout.close();
	}

}

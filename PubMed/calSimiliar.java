package PubMed;

import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.Scanner;

import common.Pair;

public class calSimiliar {

	public static void loadFeature(Scanner In, ArrayList<String> Access,
			ArrayList<ArrayList<Pair<Integer,Double>>> Feature)
	{
		while(In.hasNext())
		{
			String access = In.next();
			//System.out.println(access);
			ArrayList<Pair<Integer,Double>> linshi  = new ArrayList<Pair<Integer,Double>>();
			int count = In.nextInt();
			for (int i=1;i<=count;i++)
			{
				String str = In.next();
				String[] strarr = str.split(",");
				int index = Integer.parseInt(strarr[0]);
				double dou = Double.parseDouble(strarr[1]);
				linshi.add(new Pair<Integer,Double>(index,dou));
			}
			Access.add(access);
			Feature.add(linshi);
		}
	}
	public static double calCosSimiliar(ArrayList<Pair<Integer,Double>> instance1,ArrayList<Pair<Integer,Double>> instance2)
	{
		double len1 = 0.0;
		double len2 = 0.0;
		double sim =  0.0;
		HashMap<Integer,Double> in1 = new HashMap<Integer,Double>();
		for (Pair<Integer,Double> pa:instance1)
		{
			len1 += (pa.getSecond() * pa.getSecond());
			in1.put(pa.getFirst(), pa.getSecond());
		}
		for (Pair<Integer,Double> pa:instance2)
			len2 += (pa.getSecond() * pa.getSecond());
		
		len1 = Math.sqrt(len1);
		len2 = Math.sqrt(len2);
		int index = 0;
		for (Pair<Integer,Double> pa:instance2)
		{
			if (in1.containsKey(pa.getFirst()))
				sim += pa.getSecond()*in1.get(pa.getFirst());
		}
		if (sim<0.00001) return 0.0;
			else
		return sim/(len1*len2);
	}
	public static void main(String[] args) throws FileNotFoundException 
	{
		// TODO Auto-generated method stub
		int k=10;
		PrintWriter Fout =   new PrintWriter("PubMedSimiliar");
		String measureDirectory = "../InFile/Measure/CAFA2/";
		String trainDirectory =   "../InFile/Train/201401/";
		Scanner InM = new Scanner(new FileInputStream(measureDirectory + "PubMedFeature"));
		Scanner InT = new Scanner(new FileInputStream(trainDirectory + "PubMedFeature"));
		ArrayList<String> measureAccess = new ArrayList<String>();
		ArrayList<ArrayList<Pair<Integer,Double>>> measureFeature = new ArrayList<ArrayList<Pair<Integer,Double>>>();
		ArrayList<String> trainAccess = new ArrayList<String>();
		ArrayList<ArrayList<Pair<Integer,Double>>> trainFeature = new ArrayList<ArrayList<Pair<Integer,Double>>>();
		
		loadFeature(InM,measureAccess,measureFeature);
		loadFeature(InT,trainAccess,trainFeature);
		
		double sim = calCosSimiliar(measureFeature.get(0),trainFeature.get(0));
		System.out.println(sim);
		
		for (int i = 0;i<measureAccess.size();i++)
		{
			ArrayList<Pair<String,Double>> simList = new ArrayList<Pair<String,Double>>();
			for (int j = 0;j<trainAccess.size();j++)
			{
				sim = calCosSimiliar(measureFeature.get(i),trainFeature.get(j));
				simList.add(new Pair<String,Double>(trainAccess.get(j),sim));
			}
			Comparator<Pair<String,Double>> compar = new Comparator<Pair<String,Double>>()
			{
				@Override
				public int compare(Pair<String, Double> o1, Pair<String, Double> o2) 
				{
					// TODO Auto-generated method stub
					return  o2.getSecond().compareTo(o1.getSecond());
				}
			};
			Collections.sort(simList, compar);
			for (int j=0;j<k;j++)
			{
				//Fout.println(measureAccess.get(i) +" "+ simList.get(j).getFirst() +" "+ simList.get(j).getSecond());
				Fout.format("%s %s %.3f\n", measureAccess.get(i),simList.get(j).getFirst(),simList.get(j).getSecond());
				System.out.format("%s %s %.3f\n", measureAccess.get(i),simList.get(j).getFirst(),simList.get(j).getSecond());
			}
		}
		InM.close();
		InT.close();
		Fout.close();

	}

}

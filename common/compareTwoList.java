package common;

import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Scanner;

public class compareTwoList 
{

	public static void main(String[] args) throws FileNotFoundException 
	{
		Scanner In = new Scanner(new FileInputStream("CompareKnnBlast/CC"));
		String line = In.nextLine();
		String[] strarr = line.split("\\s+",3);
		System.out.println(strarr[0]);
		System.out.println(strarr[1]);
		System.out.println(strarr[2]);
		double fre,auc1,auc2;
		ArrayList<Double> freList  = new ArrayList<Double>();
		ArrayList<Double> auc1List = new ArrayList<Double>();
		ArrayList<Double> auc2List = new ArrayList<Double>();
		while (In.hasNextDouble())
		{
			fre = In.nextDouble();    freList.add(fre);
			auc1 = In.nextDouble();   auc1List.add(auc1);
			auc2 = In.nextDouble();   auc2List.add(auc2);	
		}
		ArrayList<String> judgement = new ArrayList<String>(Arrays.asList
				(">0.1","0.1~0.01","0.01~0.001","<0.001"));		
		double[] average1 = new double[4];
		double[] average2 = new double[4];
		int[] bignum1 = new int[4];
		int[] bignum2 = new int[4];
		int[] number = new int[4];
		int range = 0;
		for (int i = 0; i<freList.size(); i++)
		{
			fre = freList.get(i);
			if (fre>0.1)
				range = 0;
			else if (fre>0.01)
				range = 1;
			else if (fre>0.001)
				range = 2;
			else 
				range = 3;
			
			number[range]++;
			average1[range]+=auc1List.get(i);
			average2[range]+=auc2List.get(i);
			if (auc1List.get(i)>auc2List.get(i)) bignum1[range]++;
			else bignum2[range]++;
		}
		System.out.println("range,Num,AverKNN,AverBlast,better");
		for (int i=0;i<=3;i++)
			System.out.format("%s,%d,%.3f,%.3f,%d | %d\n", judgement.get(i),number[i],average1[i]/number[i]
					,average2[i]/number[i],bignum1[i], bignum2[i]);

		In.close();
	}

}

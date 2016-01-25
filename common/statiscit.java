package common;

import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.PrintWriter;

public class statiscit 
{
	public static void statDistribution(int[] Arr,String OutFile) throws FileNotFoundException
	{
		PrintWriter Fout = new PrintWriter(new FileOutputStream(OutFile));
		int[] Dis = new int[100];
		for (int i = 0 ; i<Arr.length ; i++)
		{
			int num = Arr[i];
			if (num<10) Dis[num]++;
	        else if (num<100) Dis[num/10+9]++;
	        else if (num<1000) Dis[num/100+18]++;
	        else if (num<10000) Dis[num/1000+27]++;
	        else Dis[num/10000+36]++;
		}
	    for (int i=0;i<41;i++)
	    {
	        if (i<=9) Fout.println(i + "," + Dis[i]);
	        else
	        {
	        	int begin = (int) ((((i-1)%9)+1) * Math.pow(10,((i-1)/9)));
	        	int end = (int)   ((((i-1)%9)+2)*Math.pow(10,((i-1)/9))-1);
	        Fout.println(begin + "~" + end + "," + Dis[i]);
	    
	        }
	    }
	    Fout.close();
	}
}

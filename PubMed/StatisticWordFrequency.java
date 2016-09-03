package PubMed;

import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.PrintWriter;
import java.util.HashMap;
import java.util.Scanner;

public class StatisticWordFrequency {

	public static void main(String[] args) throws FileNotFoundException 
	{
		// TODO Auto-generated method stub
		Scanner In = new Scanner(new FileInputStream("../InFile/PubMedToken"));
		PrintWriter Fout =   new PrintWriter("wordFrequency");
		String publine = new String();
		String document = new String();
		HashMap<String,Integer> wordFrequency = new HashMap<String,Integer>();
		int count = 0;
		while(In.hasNext())
		{
			
			publine = In.nextLine();
			document = In.nextLine();
			count++;
			if (count%100 == 0) System.out.println(count);
			String[] strarr = document.split(" ");
			for (String word:strarr)
			{
				if (word.length()>2)
				{
					if ((Character.isLowerCase(word.charAt(1))) &&  (Character.isUpperCase(word.charAt(0))))
					{
						word = (Character.toLowerCase(word.charAt(0))) + word.substring(1);
					}
					if (!wordFrequency.containsKey(word))
					{
						wordFrequency.put(word, 1);
					}else
					{
						int num = wordFrequency.get(word);
						wordFrequency.put(word, num + 1);
					}
				}
			}
		}
		for (String key : wordFrequency.keySet())
		{
			if ((!key.contains(","))  && (Character.isLetter(key.charAt(0))))
				Fout.println(key + "," + wordFrequency.get(key));
		}
		In.close();
	}

}

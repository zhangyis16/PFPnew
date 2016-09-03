package PubMed;

import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Scanner;
import java.util.TreeMap;

public class proteinDocumentFeature 
{

	public static void main(String[] args) throws FileNotFoundException 
	{
		// TODO Auto-generated method stub
		Scanner In = new Scanner(new FileInputStream("../InFile/PubMed/Token"));
		Scanner InWord = new Scanner(new FileInputStream("../InFile/PubMed/selectWord"));
		String setDirectory = "../InFile/Train/201401/";
		Scanner InList = new Scanner(new FileInputStream(setDirectory + "PubMedList"));
		PrintWriter Fout =   new PrintWriter(setDirectory + "PubMedFeature");
		
		ArrayList<String> textList = new ArrayList<String>();
		HashMap<Integer,Integer> textIndex = new HashMap<Integer,Integer>();
		
		HashMap<Integer,Integer> wordFrequency = new HashMap<Integer,Integer>();
		HashMap<String,Integer> wordIndex = new HashMap<String,Integer>();
		while (In.hasNext())
		{
			String line = In.nextLine();
			int pubmed = Integer.parseInt(line.substring(9));
			String document = In.nextLine();
			textList.add(document);
			textIndex.put(pubmed, textList.size() - 1);
		}
		int count = 0;
		while (InWord.hasNext())
		{
			String word = InWord.next();
			int fre = InWord.nextInt();
			
			wordIndex.put(word, count);
			wordFrequency.put(count, fre);
			count++;
		}
		while (InList.hasNext())
		{
			TreeMap<Integer,Integer> wordFeature = new TreeMap<Integer,Integer>();
			String access = InList.next();
			count = InList.nextInt();
			for (int i = 1;i <= count;i++)
			{
				int pubmed = InList.nextInt();
				String document = "";
				if (textIndex.containsKey(pubmed))
					document = textList.get(textIndex.get(pubmed));
				else
					System.out.println(pubmed + " not exist");
				
				String[] strarr = document.split(" ");
				for (String word:strarr)
				{
					if (word.length()>2)
					if ((Character.isLowerCase(word.charAt(1))) &&  (Character.isUpperCase(word.charAt(0))))
					{
						word = (Character.toLowerCase(word.charAt(0))) + word.substring(1);
					}
					if (wordIndex.containsKey(word))
					{
						int index = wordIndex.get(word);
						if (wordFeature.containsKey(index))
						{
							int num = wordFeature.get(index);
							wordFeature.put(index, num + 1);
						}else
						{
							wordFeature.put(index, 1);
						}
						
					}
				}
			}
			
			Fout.println(access);
			Fout.print(wordFeature.size());
			for (Integer key : wordFeature.keySet())
			{
				double tfidf = wordFeature.get(key);
				int frequency = wordFrequency.get(key);
				
				tfidf = (1 + Math.log(tfidf)) * 200000/frequency;
				Fout.format(" %d,%.3f" , key,tfidf);
			}
			Fout.println();
			
		}
		In.close();
		InWord.close();
		InList.close();
		Fout.close();
	}

}

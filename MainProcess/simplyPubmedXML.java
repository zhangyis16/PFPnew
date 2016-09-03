package MainProcess;

import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.PrintWriter;
import java.util.Scanner;

public class simplyPubmedXML {

	public static void main(String[] args) throws FileNotFoundException 
	{
		// TODO Auto-generated method stub
		Scanner In = new Scanner(new FileInputStream("../InFile/pubmed.xml"));
		PrintWriter Fout = new PrintWriter(new FileOutputStream("pubmed"));
		String line = new String();
		while (In.hasNext())
		{
			line = In.nextLine();
			if (line.equals("<PubmedArticle>"))
			{
				In.nextLine();In.nextLine();In.nextLine();
				line = In.nextLine();
				int loc1 = line.indexOf('>');
				int loc2 = line.lastIndexOf('<');
				int PubMed = Integer.parseInt(line.substring(loc1 + 1, loc2));
				Fout.println("PubMed = " + PubMed);
			}
			if (line.contains("<ArticleTitle>"))
			{
				int loc1 = line.indexOf('>');
				int loc2 = line.lastIndexOf('<');
				String str = line.substring(loc1 + 1, loc2);
				Fout.println(str);
			}
			if (line.contains("<AbstractText>"))
			{
				int loc1 = line.indexOf('>');
				int loc2 = line.lastIndexOf('<');
				String str = line.substring(loc1 + 1, loc2);
				Fout.println(str);
			}	
		}
		In.close();
		Fout.close();
	}

}

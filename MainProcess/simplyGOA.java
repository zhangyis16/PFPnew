package MainProcess;

import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.PrintWriter;
import java.util.Scanner;
import protein.*;

public class simplyGOA {

	public static void main(String[] args) throws FileNotFoundException {
		// TODO Auto-generated method stub
		String data = args[1];
		String InFileName = "gene_association.goa_uniprot" + data;
		Scanner In = new Scanner(new FileInputStream(InFileName));
		PrintWriter Fout = new PrintWriter(new FileOutputStream("GoaSimply" + data));
		String line,evidence;
		while(In.hasNext())
		{
			line = In.nextLine();
			String[] strarr = line.split("\t",8);
			if (strarr.length > 6)
			{
				evidence = strarr[6];
				if (proteinCommon.EvidenceCode.contains(evidence))
				Fout.println(line);
			}
		}
		In.close();
		Fout.close();
	}

}

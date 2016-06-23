package Main;

import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.PrintWriter;
import java.util.Scanner;

import protein.proteinCommon;
import protein.proteinSet;

public class extractGOAIEA_Ann {

	public static void main(String[] args) throws FileNotFoundException {
		// TODO Auto-generated method stub
		proteinSet.LoadSwissMapAccess2UniAccess("../InFile/Swiss/ac2ac201401");
		

		String InFileName = args[1];
		Scanner In = new Scanner(new FileInputStream(InFileName));
		PrintWriter Fout = new PrintWriter(new FileOutputStream("GoaIEA_Ann"));
		String line,evidence;
		while(In.hasNext())
		{
			line = In.nextLine();
			String[] strarr = line.split("\t",8);
			if (strarr.length > 7)
			{
				String Access = strarr[1];
				String Ann = new String();
				Ann = strarr[4];
			    evidence = strarr[6];
				
				if (!proteinCommon.EvidenceCode.contains(evidence)&&
						(proteinSet.MapAccess2UniAccess.containsKey(Access)))
					Fout.println(Access + "\t" + Ann);
			}
		}
		In.close();
		Fout.close();
	}

}

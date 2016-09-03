package MainProcess;

import java.io.FileNotFoundException;

import protein.GoSet;
import protein.proteinSet;

public class getMeasure {
	
	public static void main(String[] args) throws FileNotFoundException
	{
		String Tstart = "201301";
		String Tend =   "201401";
		proteinSet.LoadAccess2NameMap("../InFile/Swiss/ac2Name" + Tend);
		proteinSet.LoadSwissMapAccess2UniAccess("../InFile/Swiss/ac2ac" + Tend);
		proteinSet measure = proteinSet.getNewProtein(Tstart, Tend);
		measure.UpdateIndex();
		
		
		


		measure.loadFastaSequence("../InFile/Swiss/Swiss" + Tend + ".fasta");
		measure.OutputAnnotation("../InFile/temporary/Ann");
		measure.OutputFastaAccessSequence("../InFile/temporary/Seq");
		System.out.println("measure.size = " + measure.size());
	}
}

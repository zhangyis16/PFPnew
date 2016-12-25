package MainProcess;

import java.io.IOException;

import Main.learningOfGO;
import protein.proteinSet;

public class filterSpecies {

	public static void main(String[] args) throws IOException 
	{
		// TODO Auto-generated method stub
		learningOfGO.ReadConfig();
		//learningOfGO.aGoSet.Load(learningOfGO.GoDirectory);
		proteinSet.LoadAccess2NameMap(learningOfGO.Access2NameDirectory);
		proteinSet measure =   new proteinSet();
		measure.AddAnnotation("../InFile/Measure/201509-201409/Ann");
		measure.loadFastaSequence("../InFile/Measure/201509-201409/Seq");
		System.out.println("Size = " + measure.size());
		
		measure.filterCAFA2Species();
		
		System.out.println("Size = " + measure.size());
		measure.filterProteinOnly5515();
		
		measure.OutputAnnotation("../InFile/Measure/201509-201409Species27/Ann");
		measure.OutputFastaAccessSequence("../InFile/Measure/201509-201409Species27/Seq");
	}
}

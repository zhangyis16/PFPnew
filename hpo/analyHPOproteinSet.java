package hpo;

import java.io.FileNotFoundException;
import java.util.ArrayList;
import java.util.Collection;
import java.util.Collections;

import common.Output;
import protein.proteinSet;

public class analyHPOproteinSet {

	public static void main(String[] args) throws FileNotFoundException {
		// TODO Auto-generated method stub
		
		proteinSet.LoadAccess2NameMap("../InFile/Swiss/ac2Name201409");
		learningHPO.aHPOSet.Load(learningHPO.HPODirectory);
		System.out.println(learningHPO.aHPOSet.sizeOfHPO());
		
		proteinSet train =   new proteinSet();
		train.AddHPOAnnotation("../InFile/HPO/20160125/Ann");
		
		
		
		//train.loadFastaSequence("../InFile/Swiss/Swiss201409.fasta");
		
		
		//train.OutputFastaAccessSequence("../InFile/HPO/HPO20160125/Seq");
		//train.filterNoAccess();
		//train.OutputHPOAnnotation("../InFile/HPO/201409/Ann");
		//train.OutputFastaAccessSequence("../InFile/HPO/201409/Seq");
		
		System.out.println(train.size());

		/*		proteinSet leave1 =   new proteinSet();
		proteinSet trainNew =   proteinSet.getNewProtein(leave1,train); 
		proteinSet trainCV =   new proteinSet();
		for (int i=1;i<=5;i++)
		{
			leave1 = trainNew.getRandomSubProteinSet(698);
			leave1.OutputHPOAnnotation("../InFile/HPO/Measure/HPO20160125CV" + i + "/Ann");
			leave1.OutputFastaAccessSequence("../InFile/HPO/Measure/HPO20160125CV" + i + "/Seq");
			
			trainNew = proteinSet.getNewProtein(leave1,trainNew); 
			System.out.println(trainNew.size());
			
			trainCV  = proteinSet.getNewProtein(leave1,train); 
			trainCV.OutputHPOAnnotation("../InFile/HPO/Train/HPO20160125CV" + i + "/Ann");
			trainCV.OutputFastaAccessSequence("../InFile/HPO/Train/HPO20160125CV" + i + "/Seq");
			System.out.println(trainCV.size());
		}
		*/
		
		
		System.out.println(train.sizeOfHPOAnn());
		train.removeHPONotIn(learningHPO.aHPOSet);
		train.addHPOFather(learningHPO.aHPOSet);
		System.out.println(train.sizeOfHPOAnn());
		
		
		int[] a = train.getEveryHPOFrequency(learningHPO.aHPOSet);
		//ArrayList<Integer> sizeList = train.getEveryProteinHPOAnnSizeList();
		ArrayList<Integer> sizeList = new ArrayList<Integer>();
		for (Integer e:a)
		{
			sizeList.add(e);
		}
		Collections.sort(sizeList);
		Output.OutputList("List", sizeList);
	}
}

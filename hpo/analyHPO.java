package hpo;

import java.io.FileNotFoundException;

public class analyHPO {

	public static void main(String[] args) throws FileNotFoundException {
		// TODO Auto-generated method stub
		HPOSet aHPOSet = new HPOSet();
		aHPOSet.Load("../InFile/HPO/Ontology/HPO2016-01-13.obo");
		System.out.println(aHPOSet.sizeOfHPO());
		System.out.println(aHPOSet.sizeOfRelationShip());
		
		aHPOSet.AddSon();
		aHPOSet.addDepth();
		
		//aHPOSet.OutputDepth("HPOdepth");
		
		aHPOSet.statMaxDepth(17);
		aHPOSet.statMinDepth(15);
	}

}

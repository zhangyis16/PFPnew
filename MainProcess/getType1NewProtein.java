package MainProcess;

import java.io.FileNotFoundException;

import Main.learningOfGO;
import protein.GoSet;
import protein.proteinSet;

public class getType1NewProtein {

	public static void main(String[] args) throws FileNotFoundException 
	{
		// TODO Auto-generated method stub
		String t1 = "201509";
		String t0 = "201409";
		learningOfGO.aGoSet = new GoSet("../InFile/GeneOntology/gene_ontology_edit.obo.2013-06-15");
		
		proteinSet AnnType1 = proteinSet.getCAFA2Protein(t0, t1, learningOfGO.aGoSet);
		System.out.println(AnnType1.size());
		AnnType1.OutputAnnotation("201509-201409Type1NewProtein");
	}

}

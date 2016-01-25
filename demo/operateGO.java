package demo;

import java.io.FileNotFoundException;

import protein.GoSet;

public class operateGO {

	public static void main(String[] args) throws FileNotFoundException 
	{
		// TODO Auto-generated method stub
		GoSet Go = new GoSet("../InFile/gene ontology/gene_ontology_edit.obo.2013-07-01");
		Go.AddAllFather();
		Go.OutputGO_Term("go201307Simply");
	}

}

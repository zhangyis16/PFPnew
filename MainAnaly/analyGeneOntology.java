package MainAnaly;

import java.io.FileNotFoundException;

import protein.GoSet;

public class analyGeneOntology 
{

	public static void main(String[] args) throws FileNotFoundException 
	{
		//GoSet Go20130615 = new GoSet("../InFile/GeneOntology/gene_ontology_edit.obo.2015-09-01");
		GoSet Go20130615 = new GoSet("../InFile/GeneOntology/gene_ontology_edit.obo.2016-06-01");
		
		Go20130615.AddSon();
		Go20130615.addDepth();
		Go20130615.statMaxDepth(4);
		Go20130615.statMinDepth(4);
		Go20130615.statSubGoNum();
	}
}

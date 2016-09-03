package MainAnaly;

import java.io.FileNotFoundException;
import java.lang.reflect.InvocationTargetException;

import Main.learning;
import protein.GoSet;
import protein.proteinSet;

public class analyIntegratedYear {

	public static void main(String[] args) throws FileNotFoundException, NoSuchMethodException, SecurityException, IllegalAccessException, IllegalArgumentException, InvocationTargetException 
	{
		proteinSet measure = new proteinSet();
		learning.aGoSet = new GoSet("../InFile/GeneOntology/gene_ontology_edit.obo.2013-06-15");
		//measure.loadSwissIntegratedYear("../InFile/SwissOri/uniprot_sprot201401.dat");
		//measure.invokeMethodEveryCell("outputIntegratedYear", "IntegratedYear201401");
		proteinSet train = new proteinSet();
		train.AddAnnotation("../InFile/Train/201401/Ann");
		train.loadIntegratedYear("../InFile/Swiss/IntegratedYear201401");
		train.statYear(1986, 2014, "statYear");
		
		//train.invokeMethodEveryCell("outputIntegratedYear", "train201401Year");

	}

}

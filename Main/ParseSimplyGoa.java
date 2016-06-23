package Main;
import protein.*;
import java.io.FileNotFoundException;

public class ParseSimplyGoa {

	public static void main(String[] args) throws FileNotFoundException 
	{
		//GoSet Go201307 = new GoSet("../InFile/gene ontology/gene_ontology_edit.obo.2013-07-01");
		//proteinSet.LoadSwissAccessMap("../InFile/Swiss/uniprot_sprot201307.dat");
		//proteinSet.OutputAccessMap("ac201307");
		
		/*proteinSet measure = new proteinSet();

		measure.addSwissProteinSequence("../InFile/Swiss/uniprot_sprot201401.dat");
		measure.OutputFastaSequence("Sequence201401.fasta");   */
		proteinSet measure = new proteinSet();
		measure.addGoaAnnotation("../InFile/Goa/GoaSimply201301");
		measure.OutputAnnotation("Annotation201301");
		System.out.println(measure.size());
	}

}

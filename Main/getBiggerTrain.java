package Main;

import java.io.FileNotFoundException;

import protein.proteinSet;

public class getBiggerTrain {

	public static void main(String[] args) throws FileNotFoundException {
		// TODO Auto-generated method stub
		proteinSet Train201401 = new proteinSet();
/*		Train201401.AddAnnotation("../InFile/Trembl/AnnGoa201401");
		Train201401.AddAnnotation("../InFile/Trembl/AnnGODB201401");
		Train201401.AddAnnotation("../InFile/Trembl/AnnSwiss201401");
		Train201401.AddAnnotation("../InFile/Trembl/AnnTrembl201401");
		System.out.println(Train201401.size());
		Train201401.loadFastaSequence("../InFile/Trembl/SeqSwiss201401.fasta");
		Train201401.loadTableSequence("../InFile/Trembl/SeqTrembl201401");
		Train201401.removeNoSeqProtein();
		System.out.println(Train201401.size());
		Train201401.OutputAnnotation("AnnAll");
		Train201401.OutputFastaSequence("SeqALL");    */
		Train201401.AddAnnotation("../InFile/Trembl/AnnAll");
		
	}

}

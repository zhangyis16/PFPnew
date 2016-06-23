package Main;

import java.io.FileNotFoundException;

import protein.proteinSet;

public class extractSeq {

	public static void main(String[] args) throws FileNotFoundException {
		// TODO Auto-generated method stub
		proteinSet Ann = new proteinSet();
		Ann.AddAnnotation(args[0]);
		Ann.loadFastaSequence(args[1]);
		Ann.OutputFastaSequence(args[2]);
	}

}

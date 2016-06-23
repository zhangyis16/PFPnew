package Main;
import java.io.FileNotFoundException;

import protein.proteinSet;
public class extractGOA_noIEAAnn {

	public static void main(String[] args) throws FileNotFoundException {
		// TODO Auto-generated method stub
		proteinSet GOA = new proteinSet();
		GOA.addGoaAnnotation("../InFile/Goa/GoaSimply201301");
		GOA.OutputAnnotation("../InFile/Goa/Ann201301");
	}

}

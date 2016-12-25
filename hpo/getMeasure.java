package hpo;

import java.io.FileNotFoundException;

import protein.proteinSet;

public class getMeasure {

	public static void main(String[] args) throws FileNotFoundException {
		// TODO Auto-generated method stub
		proteinSet T0set = new proteinSet();
		proteinSet T1set = new proteinSet();
		T0set.AddHPOAnnotation("../InFile/HPO/HPO-t0");
		T1set.AddHPOAnnotation("../InFile/HPO/HPO-t1");
		proteinSet Measure = proteinSet.getNewProtein(T0set, T1set);
		System.out.println(T0set.size());
		System.out.println(T0set.sizeOfHPOAnn());
		
		System.out.println(T1set.size());
		System.out.println(T1set.sizeOfHPOAnn());
		
		System.out.println(Measure.size());
		System.out.println(Measure.sizeOfHPOAnn());
		Measure.OutputHPOAnnotation("../InFile/HPO/CAFA2Ann");
	}

}

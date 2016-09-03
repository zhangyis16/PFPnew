package demo;
import java.util.ArrayList;

import java.util.Collections;
import java.util.Comparator;

import common.Pair;
public class sortDemo {

	public static void main(String[] args) {
		// TODO Auto-generated method stub
		ArrayList<Pair<Integer,Double>> arr = new ArrayList<Pair<Integer,Double>>();
		arr.add(new Pair<Integer,Double>(40, 1.0));
		arr.add(new Pair<Integer,Double>(50, 2.0));
		arr.add(new Pair<Integer,Double>(10, 3.0));
		arr.add(new Pair<Integer,Double>(20, 4.0));
		arr.add(new Pair<Integer,Double>(30, 5.0));
		Comparator<Pair<Integer,Double>> compar = new Comparator<Pair<Integer,Double>>()
				{
					@Override
					public int compare(Pair<Integer, Double> o1, Pair<Integer, Double> o2) {
						// TODO Auto-generated method stub
						return  o2.getSecond().compareTo(o1.getSecond());
					}
				};
				
		Collections.sort(arr, compar);
		for (Pair<Integer,Double> cell:arr)
		{
			System.out.println(cell.getFirst() +"   "+ cell.getSecond());
		}
	}
}

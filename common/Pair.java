package common;

import java.util.*;

public class Pair<T extends Comparable<T>,U> implements Comparable<Pair<T,U>>
{
	public static <T extends Comparable<T>,U> Map<T,U> TranPair2Map(ArrayList<Pair<T,U>> list)
	{
		Map<T,U> returnMap = new TreeMap<T,U>();
		for (Pair<T,U> e: list)
			returnMap.put(e.getFirst(), e.getSecond());
		return null;
	}
	
	private T first;
	private U second;
	public Pair() {first = null; second = null;}
	public Pair(T first,U second) {this.first = first; this.second = second;}
	public T getFirst() {return this.first;}
	public U getSecond() {return this.second;}
	public void setFirst(T first) {this.first = first;}
	public void setSecond(U second) {this.second = second;}
	@Override
	public int compareTo(Pair<T,U> other) 
	{
		// TODO Auto-generated method stub
		T i = this.getFirst();
		T j = other.getFirst();

		return  i.compareTo(j);
	}
}


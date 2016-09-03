package demo;

import java.util.HashSet;

public class SetHandle 
{
	public static double interRatioUnion(HashSet<Integer> Set1,HashSet<Integer> Set2)
	{
		HashSet<Integer> jiao = new HashSet<Integer>();
		HashSet<Integer> bing = new HashSet<Integer>();
		for (Integer ann:Set1)
		{
			if (Set2.contains(ann)) jiao.add(ann);
			bing.add(ann);
		}
		for (Integer ann:Set2)
			bing.add(ann);
		double s1 = jiao.size();
		double s2 = bing.size();
		return s1/s2;
	}
}

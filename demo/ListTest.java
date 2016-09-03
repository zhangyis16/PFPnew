package demo;

import java.util.*;

public class ListTest {

	public static void main(String[] args) {
		// TODO Auto-generated method stub
		List<String> staff = new LinkedList<String>();
		staff.add("Amd");
		staff.add("Bob");
		staff.add("Carl");
		Iterator<String> iter = staff.iterator();
		String first =  (String) iter.next();
		String second = (String) iter.next();
		iter.remove();
		for (String e: staff)
		{
			System.out.println(e);
		}
	}

}

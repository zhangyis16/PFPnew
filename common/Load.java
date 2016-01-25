package common;

import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.util.HashMap;
import java.util.Scanner;

public class Load {
	public static void LoadMap(HashMap<String,String> Map,String InFile) throws FileNotFoundException
	{
		Scanner In = new Scanner(new FileInputStream(InFile));
		String key = new String();
		String value = new String();
		Map.clear();
		while (In.hasNext())
		{
			key = In.next();
			value = In.next();
			Map.put(key, value);
		}
		In.close();
	}
}

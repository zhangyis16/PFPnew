package demo;

import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.Scanner;

import common.StringHandle;
import protein.proteinCommon;


public class demo {

	public static void main(String[] args) throws NoSuchMethodException, SecurityException, IOException, InterruptedException 
	{
		Scanner In = new Scanner(new FileInputStream("../InFile/Trembl/ac2ac201401"));
		PrintWriter Fout = new PrintWriter(new FileOutputStream("../InFile/Trembl/ac2ac201401new"));
		while(In.hasNext())
		{
			String line = In.nextLine();
			if (!line.contains(";"))
				Fout.println(line);
		}
		Fout.close();
		In.close();
	}

}

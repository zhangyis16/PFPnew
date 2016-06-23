package Main;

import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.PrintWriter;
import java.util.HashMap;
import java.util.Map;
import java.util.Scanner;

public class idmapGoDb {

	public static void main(String[] args) throws FileNotFoundException {
		// TODO Auto-generated method stub
		Scanner In = new Scanner(new FileInputStream("../InFile/GODB/go201409result.txt"));
		Scanner IDmap = new Scanner(new FileInputStream("../InFile/SimplyIDmap/SimplyIDmap"));
		PrintWriter Fout = new PrintWriter(new FileOutputStream("../InFile/GODB/Ann201409"));
		PrintWriter FoutNoFind = new PrintWriter(new FileOutputStream("../error/GOnofind"));
		String db,key,go;
		String uniprotAcc;
		Map<String,String> idmap = new HashMap<String,String>();
		while (IDmap.hasNext())
		{
			uniprotAcc = IDmap.next();
			db = IDmap.next();
			key = IDmap.next();
			key = db + key;
			idmap.put(key, uniprotAcc);
		}
		while (In.hasNext())
		{
			db = In.next();
			key = In.next();
			go = In.next();
			if (db.equals("UniProtKB"))
			{
				Fout.println(key + '\t' + go);
			}
			else
			{
				if (db.equals("FB")) db = "FlyBase";
				if (db.equals("WB")) db = "WormBase";	
				String dbkey = db + key;
				if (idmap.containsKey(dbkey))
				{
					Fout.println(idmap.get(dbkey) + '\t' + go);
				}
				else
					FoutNoFind.println(db + '\t' + key + '\t' + go);
			}
		}
		IDmap.close();
		In.close();
		Fout.close();
		FoutNoFind.close();
	}

}

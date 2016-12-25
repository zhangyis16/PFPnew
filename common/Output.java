package common;

import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.PrintWriter;
import java.util.ArrayList;

public class Output 
{
	public static <T extends Comparable<T>,U> void OutputPairList(String OutFile,ArrayList<Pair<T,U>> PairList) throws FileNotFoundException
	{
		PrintWriter Fout = new PrintWriter(new FileOutputStream(OutFile));
		for (Pair<T,U> node: PairList)
		{
			Fout.println(node.getFirst() + "\t" + node.getSecond());
		}
		Fout.close();
	}
	
	public static <T extends Comparable<T>,U> void OutputFmaxPair(String OutFile,ArrayList<Pair<T,U>> PairList) throws FileNotFoundException
	{
		PrintWriter Fout = new PrintWriter(new FileOutputStream(OutFile));
		for (Pair<T,U> node: PairList)
		{
			Fout.println(node.getFirst() + "\t" + node.getSecond());
		}
		Fout.close();
	}
	
	public static <T extends Comparable<T>,U> void OutputPair(String OutFile,ArrayList<Pair<T,U>> PairList,String message) throws FileNotFoundException
	{
		PrintWriter Fout = new PrintWriter(new FileOutputStream(OutFile));
		Fout.println(message);
		for (Pair<T,U> node: PairList)
		{
			Fout.println(node.getFirst() + "\t" + node.getSecond());
		}
		Fout.close();
	}
	public static <T extends Comparable<T>,U> void OutputPair(PrintWriter Fout,ArrayList<Pair<T,U>> PairList) throws FileNotFoundException
	{
		for (Pair<T,U> node: PairList)
		{
			Fout.println(node.getFirst() + "\t" + node.getSecond());
		}
	}
	public static <T> void OutputList(PrintWriter Fout,ArrayList<T> List) throws FileNotFoundException
	{
		for (T node: List)
		{
			Fout.print(node + "\t");
		}
	}
	
	public static <T> void OutputList(String OutFile ,ArrayList<T> List) throws FileNotFoundException
	{
		PrintWriter Fout = new PrintWriter(new FileOutputStream(OutFile));
		for (T node: List)
		{
			Fout.println(node + "\t");
		}
		Fout.close();
	}
}

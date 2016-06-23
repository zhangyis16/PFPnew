package common;

import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.PrintWriter;
import java.util.ArrayList;

public class FeatureTran 
{
	public static ArrayList<Double> VecMuitiMatrix(ArrayList<Double> Fea,ArrayList<ArrayList<Pair<Integer,Double>>> Matrix)
	{
		int length = Fea.size();
		ArrayList<Double> newFea = new ArrayList<Double>(length);
		for (int i = 0;i<length; i++)
			newFea.add(0.0);
		for (int i = 1;i<Math.min(Matrix.size() ,length)    ;i++)
		{
			if (Fea.get(i) > 0.001)
			for (int j = 0;j<Matrix.get(i).size();j++)
			{
				
				double sim = Matrix.get(i).get(j).getSecond();
				int Ob = Matrix.get(i).get(j).getFirst();
				newFea.set(Ob, newFea.get(Ob) + Fea.get(i) * sim);
			}
		}
		return newFea;
	}
	public static void OutputSparseMatrix(String OutFile,ArrayList<ArrayList<Pair<Integer,Double>>> Matrix) throws FileNotFoundException
	{
		PrintWriter Fout = new PrintWriter(new FileOutputStream(OutFile));
		for (int i=0;i<Matrix.size();i++)
		{
			for (int j=0;j<Matrix.get(i).size();j++)
			{
				Pair<Integer,Double> pair = Matrix.get(i).get(j);
				Fout.print(pair.getFirst() + ":" + pair.getSecond()+"    ");
			}
			Fout.println();
		}
		Fout.close();
	}
}

package atry;

import java.io.File;
import java.io.IOException;
import  java.lang.*;
import java.util.ArrayList;

import common.Pair;
import liblinear.Feature;
import liblinear.FeatureNode;
import liblinear.InvalidInputDataException;
import liblinear.Linear;
import liblinear.Model;
import liblinear.Parameter;
import liblinear.Problem;
import liblinear.SolverType;
import liblinear.Train;
public class atry 
{
	private static Parameter param            = null;
    private static Problem   prob             = null;
	public static void main(String[] args) throws IOException, InterruptedException, InvalidInputDataException 
	{
		// TODO Auto-generated method stub
		/*String cmd = "makeblastdb -in blastdb/train.fasta -parse_seqids -hash_index -dbtype prot";
		Process process1 = Runtime.getRuntime().exec(cmd);
		process1.waitFor();
		cmd = "blastp -task blastp -query blastdb/sample.fasta -db blastdb/train.fasta -out blastdb/blastResult.out -outfmt 6";
		Process process2 = Runtime.getRuntime().exec(cmd);
		process2.waitFor();    */
		
		//new Train().run(args);
		//measure.naiveBaseline(train);
		//measure.evalution();
		
/*		measure.addBlastResult("../InFile/blastResult/blastResult.out");   
		
		measure.GOtchaBaseline(train);
		measure.evalution();     */
		
		
/*		ArrayList<Pair<Integer,Double>> predList = train.getNaiveList();
		proteinCommon.NaivePred = Pair.TranPair2Map(predList);
		
		Set<Integer> goset = new TreeSet<Integer>();
		for (Pair<Integer,Double> e:predList)
		{
			if (e.getSecond()< (double)10.0/50000.0) goset.add(e.getFirst());
		}

		
		

		measure.removeAnnotation(goset);
		train.removeAnnotation(goset);         */

        param   = new Parameter(SolverType.L2R_LR, 1, 1000, 0.01);
        
        prob    = new Problem();

        prob.bias = -1;
        prob.l = 2;
        prob.n = 2;
        
        prob.x = new Feature[prob.l][];
        prob.y = new double[prob.l];
        
        prob.y[0] = 1;
        prob.y[1] = -1;
        for (int i = 0; i < prob.l; i++)
        {
        	Feature[] x = new Feature[2];
        	x[0] = new FeatureNode(1,i);
        	x[1] = new FeatureNode(2,i);
        	prob.x[i] = x;
        }
    	Model model = Linear.train(prob, param);
        Linear.saveModel(new File("Label" + 10 + ".model"), model);
	}

}

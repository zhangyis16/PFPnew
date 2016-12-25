package Main;

import java.io.IOException;

public class blastCommand {

	public static void main(String[] args) throws IOException, InterruptedException {
		// TODO Auto-generated method stub
		String cmd = "makeblastdb -in ../InFile/Sequence/Train201401.fasta -parse_seqids -hash_index -dbtype prot";
		Process process1 = Runtime.getRuntime().exec(cmd);
		process1.waitFor();
		cmd = "blastp -task blastp -query measure -db train.fasta -out blastResult -outfmt 6";
		cmd = "blastp -task blastp -query ../InFile/Sequence/CAFA2type1Seq.fasta -db ../InFile/Sequence/Train201401.fasta "
				+ "-out ../InFile/blastResult/blastResult -outfmt 6";
		String psicmd = "psiblast -query Q2TA -db Seq -out psiresult -outfmt 7 -num_iterations 2 -evalue 100";
		Process process2 = Runtime.getRuntime().exec(psicmd);
		process2.waitFor();
	}
}
package Main;

import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.PrintWriter;
import java.util.Scanner;

import javax.print.attribute.standard.PrinterLocation;

import protein.GoSet;

public class hehe {

	public static void main(String[] args) throws FileNotFoundException {
		// TODO Auto-generated method stub
		Scanner In = new Scanner(new FileInputStream("../bin/everyLabel"));
		PrintWriter Fout = new PrintWriter(new FileOutputStream("labelResult"));
		learning.aGoSet = new GoSet("../InFile/gene ontology/gene_ontology_edit.obo." + "2013-06-15");
		int label;
		double fre;
		double score1,score2;
		char sp;
		int knnF = 0,knnP =0,knnC=0;
		int llF=0,llP=0,llC=0;
		double AUCF1=0,AUCP1=0,AUCC1=0;
		double AUCF2=0,AUCP2=0,AUCC2=0;
		double SUMF=0,SUMP=0,SUMC=0;
		while (In.hasNext())
		{
			label = In.nextInt();
			fre = In.nextDouble();
			score1 = In.nextDouble();
			score2 = In.nextDouble();
			sp = learning.aGoSet.getSpace(label);
			if (fre>(double)10/60000) 
			{
				if ((score1<0.99) && (score2<0.99))
				{
					Fout.print(label + "," + (int)(fre*60000) + "," + sp);
					if (score1>score2) Fout.println(",Knn"); else Fout.println(",LL");
					if ((fre>(double)10/60000) &&(fre<=(double)50/60000))
					{
						if (sp == 'F')
						{
							AUCF1+=score1;
							AUCF2+=score2;
							SUMF+=1;
						}
						if (sp == 'P')
						{
							AUCP1+=score1;
							AUCP2+=score2;
							SUMP+=1;
						}
						if (sp == 'C')
						{
							AUCC1+=score1;
							AUCC2+=score2;
							SUMC+=1;
						}
							
						if (sp == 'F') if (score1>score2) knnF++; else llF++;
						if (sp == 'P') if (score1>score2) knnP++; else llP++;
						if (sp == 'C') if (score1>score2) knnC++; else llC++;
					}
			
				}
			}
		}
		System.out.println(knnF + ","+ llF);
		System.out.println(knnP + ","+ llP);
		System.out.println(knnC + ","+ llC);
		System.out.printf("%.4f /%.4f\n",AUCF1/SUMF,AUCF2/SUMF);
		System.out.printf("%.4f /%.4f\n",AUCP1/SUMP,AUCP2/SUMP);
		System.out.printf("%.4f /%.4f\n",AUCC1/SUMC,AUCC2/SUMC);

		In.close();
		Fout.close();
	}

}

package protein;
import java.util.*;
import java.io.*;

class GO_Term
{
    private int ID;
    private char Space;			//F表示MF   P表示BP  C表示CC
    private HashSet<Integer> father = new HashSet<Integer>();
    private HashSet<Integer> son = new HashSet<>();
    public GO_Term(int id){	ID = id; }
    private int max_Depth = 0;
    private int min_Depth = 100;
    
    public int getID(){	return ID;}
    public char getSpace(){	return Space;}
    private boolean addAllFather = false;
    
    public Set<Integer> getSonList(){return this.son;}
    public boolean addAllFather()
    {
    	return addAllFather;
    }
    public void addAllFather(boolean bo)
    {
    	this.addAllFather = bo;
    }
    public int getFatherSize() {return father.size();}
    public int getSonSize() {return son.size();}
    public void setSpace(char aspace){ Space = aspace;}
    
    public int getMaxDepth() {return this.max_Depth;}
    public int getMinDepth() {return this.min_Depth;}
    public void setMaxDepth(int depth) {this.max_Depth = depth;}
    public void setMinDepth(int depth) {this.min_Depth = depth;}
    public HashSet<Integer> getFatherList()
    {
    	//HashSet<Integer> clone = (HashSet<Integer>) this.father.clone();
    	return  this.father;
    }
    public void addFather(String GoStr)
    {
    	if (proteinCommon.isGO_ID(GoStr)) father.add(proteinCommon.GOStr2Int(GoStr));	
    }
    public void addFather(int Go)
    {
    	father.add(Go);
    }
    public void addSon(int Gonum)
    {
    	son.add(Gonum);
    }
    public void setFatherSet(HashSet<Integer> aSet)
    {
    	this.father = aSet;
    }
    public void OutputFatherList(PrintWriter Fout)
    {
    	for (int fatherNode:this.father)
    	{
    		Fout.print(fatherNode+"\t");
    	}
    	Fout.println();
    }
}

public class GoSet 
{
    private ArrayList<GO_Term> Cell = new ArrayList<GO_Term>();
    private Map<Integer,Integer> MapGO_Index = new HashMap<Integer,Integer>();
    private int GO_Num = 0;
    private boolean addAllFather = false;
    
    public boolean isAddAllFather()
    {
    	return this.addAllFather;
    }
    public GoSet(String InFile) throws FileNotFoundException 
    {
    	Load(InFile);
    }
    public GoSet() {
		// TODO Auto-generated constructor stub
	}
	public int getIndex(int gonum)
    {
    	return MapGO_Index.get(gonum);
    }
    public Set<Integer> getFatherList(int gonum)
    {
    	int index = MapGO_Index.get(gonum);
    	return Cell.get(index).getFatherList();
    }
    
    public char getSpace(int gonum)
    {
    	if (MapGO_Index.containsKey(gonum)) 
    	return 
    		Cell.get(MapGO_Index.get(gonum)).getSpace(); else return 'X';
    }
    public char getIndexSpace(int index)
    {
    	if ((index >= 0) && (index<Cell.size()))
    	return Cell.get(index).getSpace(); else return 'X';
    }
    public boolean containNode(int node)
    {
    	return MapGO_Index.containsKey(node);
    }
    public int size()
    {
    	return GO_Num;
    }

    public void OutputGO_Term(String OutFile) throws FileNotFoundException
    {
    	PrintWriter Fout = new PrintWriter(new FileOutputStream(OutFile));
    	Fout.println("GOindex\t"+"Space\t"+"FatherList");
    	for(int i = 0; i<Cell.size(); i++)
    	{
    		int index = Cell.get(i).getID();
    		Fout.print(proteinCommon.GOInt2Str(index)+"\t"+getSpace(index)+"\t");
    		this.Cell.get(i).OutputFatherList(Fout);
    	}
    	Fout.close();
    }
    public void AddSon()
    {
    	int indexfather;
    	for(int i = 0; i<GO_Num; i++)
    	{
    			Set<Integer> fatherlist = Cell.get(i).getFatherList();
    			for (int father:fatherlist)
    			{
    				indexfather = MapGO_Index.get(father);
    				Cell.get(indexfather).addSon(Cell.get(i).getID());
    			}
    	}
    }
    public int checkFather(int gonum1,int gonum2)
    {
    	int index1 = MapGO_Index.get(gonum1);
    	int index2 = MapGO_Index.get(gonum2);
    	if (Cell.get(index1).getFatherList().contains(gonum2)) return gonum2;
    	if (Cell.get(index2).getFatherList().contains(gonum1)) return gonum1;
    	return 0;
    }
    public void AddAllFather()
    {
    	for (int i=0;i<Cell.size();i++)
    		this.AddAllFather(i);
    	this.addAllFather = true;
    }
    public void AddAllFather(int index)
    {
    	if (!Cell.get(index).addAllFather())
    	{
    		HashSet<Integer> fatherSet = Cell.get(index).getFatherList();
    		@SuppressWarnings("unchecked")
			HashSet<Integer> result = (HashSet<Integer>) Cell.get(index).getFatherList().clone();
    		if (fatherSet.size()>0)
    		for (Integer father:fatherSet)
    		{
    			int j = this.MapGO_Index.get(father);
    			if (!Cell.get(j).addAllFather())
    				this.AddAllFather(j);
    			Set<Integer> ancestorSet = Cell.get(j).getFatherList();
    			result.addAll(ancestorSet);
    		}
    		Cell.get(index).setFatherSet(result);
    		Cell.get(index).addAllFather(true);
    	}
    }
    public void Load(String InFile) throws FileNotFoundException 
    {
    	Scanner In = new Scanner(new FileInputStream(InFile));
    	String line;
    	System.out.println("Read Gene Ontology Now");
    	MapGO_Index.clear();
    	GO_Num = 0;
    	while (In.hasNext())
    	{
    		line = In.nextLine();
    		String[] strarr = line.split(" ",4);
    		if ((strarr[0].equals("id:")) && (proteinCommon.isGO_ID(strarr[1])))
    		{
    			int i = proteinCommon.GOStr2Int(strarr[1]);
    			this.MapGO_Index.put(i,GO_Num);
    			this.Cell.add(new GO_Term(i));
    			this.GO_Num++;
    		}
    		if (strarr[0].equals("namespace:"))
    		{
    			if (strarr[1].equals("molecular_function")) Cell.get(GO_Num-1).setSpace('F');
    			if (strarr[1].equals("biological_process")) Cell.get(GO_Num-1).setSpace('P');
    			if (strarr[1].equals("cellular_component"))	Cell.get(GO_Num-1).setSpace('C');		
    		}
    		if (strarr[0].equals("alt_id:"))
    		{
    			int i = proteinCommon.GOStr2Int(strarr[1]);
    			MapGO_Index.put(i,GO_Num-1);
    		}
    		if (strarr[0].equals("is_a:"))
    		{
    			Cell.get(GO_Num-1).addFather(strarr[1]);
    		}
    		if ((strarr[0].equals("relationship:")) && (strarr[1].equals("part_of")))
    		{
    			Cell.get(GO_Num-1).addFather(strarr[2]);
    		}
    		if (strarr[0].equals("is_obsolete:"))
    		{
    			Cell.remove(--GO_Num);	
    		}
    	}
    	System.out.println("Read Gene Ontology Finish");
    	In.close();
    }
    
    public void stat()
    {
    	System.out.println(this.Cell.size());
    	int MFnum = 0;
    	int BPnum = 0;
    	int CCnum = 0;
    	for (int i=0;i<Cell.size();i++)
    	{
    		if (Cell.get(i).getSpace() == 'F') MFnum++;
    		if (Cell.get(i).getSpace() == 'P') BPnum++;
    		if (Cell.get(i).getSpace() == 'C') CCnum++;
    	}
    	System.out.println("MFnum =" + MFnum);
    	System.out.println("BPnum =" + BPnum);
    	System.out.println("CCnum =" + CCnum);
    	for (int i=0;i<Cell.size();i++)
    	{
    		int gonum = Cell.get(i).getID();
    		Set<Integer> sonList = this.getSonList(gonum);
    		for (Integer son:sonList)
    		{
    			if (this.getSpace(gonum) != this.getSpace(son))
    				System.out.println(gonum + " and " + son + "have different space");
    		}
    	}
    	
    	int MFmaxDepth = 0;
    	int BPmaxDepth = 0;
    	int CCmaxDepth = 0;
    	for (int i=0;i<Cell.size();i++)
    	{
    		//System.out.println(Cell.get(i).getID() + "  "+ Cell.get(i).getMaxDepth());
    		if (Cell.get(i).getSpace() == 'F')
    			MFmaxDepth = (Cell.get(i).getMaxDepth()> MFmaxDepth) ? Cell.get(i).getMaxDepth():MFmaxDepth;
        	if (Cell.get(i).getSpace() == 'P')
        		BPmaxDepth = (Cell.get(i).getMaxDepth()> BPmaxDepth) ? Cell.get(i).getMaxDepth():BPmaxDepth;
            if (Cell.get(i).getSpace() == 'C')
            	CCmaxDepth = (Cell.get(i).getMaxDepth()> CCmaxDepth) ? Cell.get(i).getMaxDepth():CCmaxDepth;
    	}
    	System.out.println("MFmaxDepth = " + MFmaxDepth);
    	System.out.println("BPmaxDepth = " + BPmaxDepth);
    	System.out.println("CCmaxDepth = " + CCmaxDepth);
    	
    	int MFminDepth = 0;
    	int BPminDepth = 0;
    	int CCminDepth = 0;
    	for (int i=0;i<Cell.size();i++)
    	{
    		//System.out.println(Cell.get(i).getID() + "  "+ Cell.get(i).getMaxDepth());
    		if (Cell.get(i).getSpace() == 'F')
    			MFminDepth = (Cell.get(i).getMinDepth()> MFminDepth) ? Cell.get(i).getMinDepth():MFminDepth;
        	if (Cell.get(i).getSpace() == 'P')
        		BPminDepth = (Cell.get(i).getMinDepth()> BPminDepth) ? Cell.get(i).getMinDepth():BPminDepth;
            if (Cell.get(i).getSpace() == 'C')
            	CCminDepth = (Cell.get(i).getMinDepth()> CCminDepth) ? Cell.get(i).getMinDepth():CCminDepth;
    	}
    	System.out.println("MFminDepth = " + MFminDepth);
    	System.out.println("BPminDepth = " + BPminDepth);
    	System.out.println("CCminDepth = " + CCminDepth);
    	
    }
    public void statMinDepth()
    {
    	int[] MFdepth = new int[14];
    	int[] BPdepth = new int[14];
    	int[] CCdepth = new int[14];
    	for (int i=0;i<Cell.size();i++)
    	{
    		char sp = Cell.get(i).getSpace();
    		if (sp == 'F') MFdepth[Cell.get(i).getMinDepth()]++;
    		if (sp == 'P') BPdepth[Cell.get(i).getMinDepth()]++;
    		if (sp == 'C') CCdepth[Cell.get(i).getMinDepth()]++;
    	}
    	System.out.println("depth,MF,BP,CC");
    	for (int i=1;i<=13;i++)
    	{
    		System.out.println(i + "," + MFdepth[i] + "," + BPdepth[i] + "," + CCdepth[i]);
    	}
    	
    }
    public void statMaxDepth()
    {
    	int[] MFdepth = new int[18];
    	int[] BPdepth = new int[18];
    	int[] CCdepth = new int[18];
    	for (int i=0;i<Cell.size();i++)
    	{
    		char sp = Cell.get(i).getSpace();
    		if (sp == 'F') MFdepth[Cell.get(i).getMaxDepth()]++;
    		if (sp == 'P') BPdepth[Cell.get(i).getMaxDepth()]++;
    		if (sp == 'C') CCdepth[Cell.get(i).getMaxDepth()]++;
    	}
    	System.out.println("depth,MF,BP,CC");
    	for (int i=1;i<=17;i++)
    	{
    		System.out.println(i + "," + MFdepth[i] + "," + BPdepth[i] + "," + CCdepth[i]);
    	}
    	
    }
    public void setMindepth(int gonum,int depth)
    {
    	int index = this.MapGO_Index.get(gonum);
    	Cell.get(index).setMinDepth(depth);
    }
    public void setMaxdepth(int gonum,int depth)
    {
    	int index = this.MapGO_Index.get(gonum);
    	Cell.get(index).setMaxDepth(depth);
    }
    public int getMindepth(int gonum)
    {
    	int index = this.MapGO_Index.get(gonum);
    	return Cell.get(index).getMinDepth();
    }
    public int getMaxdepth(int gonum)
    {
    	int index = this.MapGO_Index.get(gonum);
    	return Cell.get(index).getMaxDepth();
    }
    public Set<Integer> getSonList(int gonum)
    {
    	int index = this.MapGO_Index.get(gonum);
    	return Cell.get(index).getSonList();
    }
    
    public void outputDepth(String OutFile) throws FileNotFoundException
    {
    	PrintWriter Fout = new PrintWriter(new FileOutputStream(OutFile));
    	Fout.println("Label,maxDepth,minDepth,Spaec");
    	for (int i=0;i<Cell.size();i++)
    	{
    		int ID = Cell.get(i).getID();
    		int maxDepth = Cell.get(i).getMaxDepth();
    		int minDepth = Cell.get(i).getMinDepth();
    		Fout.println(ID + "," + maxDepth + "," + minDepth + "," + Cell.get(i).getSpace());	
    	}
    	Fout.close();
    }
    public void addDepth()
    {
    	this.setMaxdepth(8150, 1);
    	this.setMindepth(8150, 1);
    	this.setMaxdepth(3674, 1);
    	this.setMindepth(3674, 1);
    	this.setMaxdepth(5575, 1);
    	this.setMindepth(5575, 1);

    	Set<Integer> SetGonum = new HashSet<Integer>();
    	SetGonum.add(8150);  SetGonum.add(3674);  SetGonum.add(5575);
    	Queue<Integer> qc = new LinkedList<Integer>();
    	qc.offer(8150);  qc.offer(3674); qc.offer(5575);
    	
    	while (qc.peek() != null)
    	{
    		int gonum = (int) qc.remove();
    		SetGonum.remove(gonum);
    		int maxDepth = this.getMaxdepth(gonum) + 1;
    		int minDepth = this.getMindepth(gonum) + 1;
    		
    		Set<Integer> sonList = this.getSonList(gonum);
    		for (Integer son:sonList)
    		{
    			if (maxDepth>this.getMaxdepth(son))
    			{
    				this.setMaxdepth(son, maxDepth);
    				if (!SetGonum.contains(son))
    				{
    					SetGonum.add(son);
    					qc.offer(son);
    				}
    			}
    			if (minDepth<this.getMindepth(son))
    			{
    				this.setMindepth(son, minDepth);
    				if (!SetGonum.contains(son))
    				{
    					SetGonum.add(son);
    					qc.offer(son);
    				}
    			}
    		}
    		
    	}
    	
    }
}

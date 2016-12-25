package hpo;

import java.util.*;

import protein.proteinCommon;

import java.io.*;

class HPO_Term
{
    private int ID;
    private HashSet<Integer> father = new HashSet<Integer>();
    private HashSet<Integer> son = new HashSet<>();
    private String name = new String();
    private boolean addAllFather = false;
    private int max_Depth = 0;
    private int min_Depth = 100;
    
    
    public int getID(){	return ID;}
    public HPO_Term(int id){	ID = id; }
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
    
    public int getMaxDepth() {return this.max_Depth;}
    public int getMinDepth() {return this.min_Depth;}
    public void setMaxDepth(int depth) {this.max_Depth = depth;}
    public void setMinDepth(int depth) {this.min_Depth = depth;}
    public HashSet<Integer> getFatherList()
    {
    	//HashSet<Integer> clone = (HashSet<Integer>) this.father.clone();
    	return  this.father;
    }
    public void addFather(String HPOStr)
    {
    	if (proteinCommon.isHPO_ID(HPOStr)) father.add(proteinCommon.HPStr2Int(HPOStr));	
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
    	for (int fatherNode : this.father)
    	{
    		Fout.print("\t" + fatherNode);
    	}
    	Fout.println();
    }
}

public class HPOSet
{
    private ArrayList<HPO_Term> Cell = new ArrayList<HPO_Term>();
    private Map<Integer,Integer> MapHPO_Index = new HashMap<Integer,Integer>();
    
    private Map<Integer,Integer> MapHPOIndexToHPONum = new HashMap<Integer,Integer>();
    private int HPO_Num = 0;
    private boolean addAllFather = false;
    
    public boolean isAddAllFather()
    {
    	return this.addAllFather;
    }
    public HPOSet(String InFile) throws FileNotFoundException 
    {
    	Load(InFile);
    }
    public HPOSet() {
		// TODO Auto-generated constructor stub
	}
	public int getIndex(int HPOnum)
    {
    	return MapHPO_Index.get(HPOnum);
    }
	
	public String getGOFromIndex(int index)
    {
		int gonum = this.MapHPOIndexToHPONum.get(index);
    	return  proteinCommon.GOInt2Str(gonum);  
    }
	
    public Set<Integer> getFatherList(int HPOnum)
    {
    	int index = MapHPO_Index.get(HPOnum);
    	return Cell.get(index).getFatherList();
    }
    


    public boolean containNode(int node)
    {
    	return MapHPO_Index.containsKey(node);
    }
    public int sizeOfHPO()
    {
    	return HPO_Num;
    }
    
    public int sizeOfRelationShip()
    {
    	int sum = 0;
    	for(int i = 0; i<Cell.size(); i++)
    	{
    		sum += this.Cell.get(i).getFatherSize();
    	}
    	return sum;
    }

    public void OutputHPO_Term(String OutFile) throws FileNotFoundException
    {
    	PrintWriter Fout = new PrintWriter(new FileOutputStream(OutFile));
    	Fout.println("HPOindex\t"+"FatherList");
    	for(int i = 0; i<Cell.size(); i++)
    	{
    		int index = Cell.get(i).getID();
    		Fout.print(proteinCommon.HPInt2Str(index));
    		this.Cell.get(i).OutputFatherList(Fout);
    	}
    	Fout.close();
    }
    public void AddSon()
    {
    	int indexfather;
    	for(int i = 0; i<HPO_Num; i++)
    	{
    			Set<Integer> fatherlist = Cell.get(i).getFatherList();
    			for (int father:fatherlist)
    			{
    				indexfather = MapHPO_Index.get(father);
    				Cell.get(indexfather).addSon(Cell.get(i).getID());
    			}
    	}
    }
    public int checkFather(int gonum1,int gonum2)
    {
    	int index1 = MapHPO_Index.get(gonum1);
    	int index2 = MapHPO_Index.get(gonum2);
    	if (Cell.get(index1).getFatherList().contains(gonum2)) return gonum2;
    	if (Cell.get(index2).getFatherList().contains(gonum1)) return gonum1;
    	return 0;
    }
    public boolean checkHPO1_isHPO2father(int HPO1,int HPO2)
    {
    	int index2 = MapHPO_Index.get(HPO2);
    	if (Cell.get(index2).getFatherList().contains(HPO1)) return true;
    	return false;
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
    			int j = this.MapHPO_Index.get(father);
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
    	System.out.println("Read HPO Now");
    	MapHPO_Index.clear();
    	
    	this.MapHPOIndexToHPONum.clear();
    	HPO_Num = 0;
    	while (In.hasNext())
    	{
    		line = In.nextLine();
    		String[] strarr = line.split(" ",4);
    		if ((strarr[0].equals("id:")) && (proteinCommon.isHPO_ID(strarr[1])))
    		{
    			int i = proteinCommon.HPStr2Int(strarr[1]);
    			this.MapHPO_Index.put(i,HPO_Num);
    			this.MapHPOIndexToHPONum.put(HPO_Num,i);
    			
    			this.Cell.add(new HPO_Term(i));
    			this.HPO_Num++;
    		}
    		if (strarr[0].equals("alt_id:"))
    		{
    			int i = proteinCommon.HPStr2Int(strarr[1]);
    			MapHPO_Index.put(i,HPO_Num-1);
    		}
    		if (strarr[0].equals("is_a:"))
    		{
    			Cell.get(HPO_Num - 1).addFather(strarr[1]);
    		}
    		if ((strarr[0].equals("relationship:")) && (strarr[1].equals("part_of")))
    		{
    			Cell.get(HPO_Num - 1).addFather(strarr[2]);
    		}
    		if (strarr[0].equals("is_obsolete:"))
    		{
    			Cell.remove(--HPO_Num);	
    		}
    	}
    	System.out.println("Read HPO Finish");
    	In.close();
    }
    
    public void statMinDepth(int depth)
    {
    	int[] HPOdepth = new int[16];

    	for (int i=0;i<Cell.size();i++)
    	{
    		HPOdepth[Cell.get(i).getMinDepth()]++;
    	}
    	System.out.println("MinDepth,HPO");
    	for (int i=1;i<=depth;i++)
    	{
    		System.out.println(i + " & " + HPOdepth[i]);
    	}
    	
    }
    public void statMaxDepth(int depth)
    {
    	int[] HPOdepth = new int[21];

    	for (int i=0;i<Cell.size();i++)
    	{
    		HPOdepth[Cell.get(i).getMaxDepth()]++;
    	}
    	System.out.println("MaxDepth,HPO");
    	for (int i=1;i<=depth;i++)
    	{
    		System.out.println(i + " & " + HPOdepth[i]);
    	}
    }
    public void setMindepth(int HPOnum,int depth)
    {
    	int index = this.MapHPO_Index.get(HPOnum);
    	Cell.get(index).setMinDepth(depth);
    }
    public void setMaxdepth(int HPOnum,int depth)
    {
    	int index = this.MapHPO_Index.get(HPOnum);
    	Cell.get(index).setMaxDepth(depth);
    }
    public int getMindepth(int HPOnum)
    {
    	int index = this.MapHPO_Index.get(HPOnum);
    	return Cell.get(index).getMinDepth();
    }
    public int getMaxdepth(int HPOnum)
    {
    	int index = this.MapHPO_Index.get(HPOnum);
    	return Cell.get(index).getMaxDepth();
    }
    public Set<Integer> getSonList(int HPOnum)
    {
    	int index = this.MapHPO_Index.get(HPOnum);
    	return Cell.get(index).getSonList();
    }
    
    public void OutputDepth(String OutFile) throws FileNotFoundException
    {
    	PrintWriter Fout = new PrintWriter(new FileOutputStream(OutFile));
    	Fout.println("Label,maxDepth,minDepth");
    	for (int i=0;i<Cell.size();i++)
    	{
    		int ID = Cell.get(i).getID();
    		int maxDepth = Cell.get(i).getMaxDepth();
    		int minDepth = Cell.get(i).getMinDepth();
    		Fout.println(ID + "," + maxDepth + "," + minDepth);	
    	}
    	Fout.close();
    }
    public void addDepth()
    {
    	this.setMaxdepth(1, 1);
    	this.setMindepth(1, 1);


    	Set<Integer> SetHPOnum = new HashSet<Integer>();
    	SetHPOnum.add(1);  //1为HPO的root节点。
    	Queue<Integer> qc = new LinkedList<Integer>();
    	qc.offer(1);
    	
    	while (qc.peek() != null)
    	{
    		int HPOnum = (int) qc.remove();
    		SetHPOnum.remove(HPOnum);
    		int maxDepth = this.getMaxdepth(HPOnum) + 1;
    		int minDepth = this.getMindepth(HPOnum) + 1;
    		
    		Set<Integer> sonList = this.getSonList(HPOnum);
    		for (Integer son:sonList)
    		{
    			if (maxDepth>this.getMaxdepth(son))
    			{
    				this.setMaxdepth(son, maxDepth);
    				if (!SetHPOnum.contains(son))
    				{
    					SetHPOnum.add(son);
    					qc.offer(son);
    				}
    			}
    			if (minDepth<this.getMindepth(son))
    			{
    				this.setMindepth(son, minDepth);
    				if (!SetHPOnum.contains(son))
    				{
    					SetHPOnum.add(son);
    					qc.offer(son);
    				}
    			}
    		}
    		
    	}
    	
    }
}


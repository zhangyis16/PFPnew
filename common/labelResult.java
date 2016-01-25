package common;


public class labelResult implements Comparable<labelResult>
{
    private int Index;
    private boolean Result;
    private double Value;
    public labelResult(int aIndex ,boolean aResult,double aValue)
    {
    	this.Index = aIndex;
    	this.Result = aResult;
    	this.Value = aValue;
    }
    public int getIndex() {return Index;}
    public boolean getResult() {return Result;}
    public Double getValue() {return Value;}
	@Override
	public int compareTo(labelResult other) 
	{
		return other.getValue().compareTo(this.getValue());
	}
}

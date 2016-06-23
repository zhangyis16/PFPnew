package common;

public class Tri <F,S,T>
{
	private F first;
	private S second;
	private T thrid;
	public Tri() 
	{
		first = null; second = null; this.thrid = null;
	}
	public Tri(F first,S second,T thrid) 
	{
		this.first = first; this.second = second; this.thrid = thrid;
	}
	public F getFirst() {return this.first;}
	public S getSecond() {return this.second;}
	public T getThrid() {return this.thrid;}
	
	public void setFirst(F first) {this.first = first;}
	public void setSecond(S second) {this.second = second;}
	public void setThrid(T thrid) {this.thrid = thrid;}
}

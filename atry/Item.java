package atry;

import java.util.Objects;

public class Item implements Comparable<Item>
{
	private String description;
	private int partNumber;
	public Item(String aDescription,int apartNumber)
	{
		description = aDescription;
		partNumber = apartNumber;
	}
	
	public String getDescription()
	{
		return this.description;
	}
	
	public String toString()
	{
		return "[description=" + this.description +", partNumber=" + partNumber +"]"; 
	}
	
	public boolean equals(Object otherObject)
	{
		if (this == otherObject) return true;
		if (otherObject == null) return false;
		if (getClass() != otherObject.getClass()) return false;
		Item other = (Item) otherObject;
		return Objects.equals(description, other.description) && partNumber == other.partNumber;
	}
	public int compareTo(Item other) {
		// TODO Auto-generated method stub
		return Integer.compare(partNumber, other.partNumber);
	}
	
	public int hashCode()
	{
		return Objects.hash(description,partNumber);
	}
}

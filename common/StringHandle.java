package common;

public class StringHandle 
{
	public static String removeChar(String str,char ch)
	{
		String reStr = new String();
		for (int i=0;i<str.length();i++)
			if ( str.charAt(i) != ch)
			{
				reStr = reStr.concat(str.substring(i, i+1));
			}
		return reStr;
	}
}

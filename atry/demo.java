package atry;

import java.awt.*;
import java.awt.event.*;
import java.util.*;
import javax.swing.*;
import javax.swing.Timer;

public class demo {

	public static void main(String[] args) {
		// TODO Auto-generated method stub
		ActionListener listener = new TimePrinter();
		Timer t = new Timer(2000,listener);
		t.start();
		JOptionPane.showMessageDialog(null, "Quit program");
		System.exit(0);
	}

}

class TimePrinter implements ActionListener
{
	public void actionPerformed(ActionEvent event)
	{
		Date now = new Date();
		System.out.println("At the tone,the time is" + now);
		
	}
}
package atry;

import java.awt.Toolkit;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.util.Date;

import javax.swing.JOptionPane;
import javax.swing.Timer;

public class InnerClassTest
{
	public static void main (String[] args)
	{
		TalkingClock clock = new TalkingClock(1000,true);
		clock.start();
		JOptionPane.showMessageDialog(null, "Quit program");
		System.exit(0);
	}
}

class TalkingClock
{
	private int interval;
	private boolean beep;
	public TalkingClock(int Ainterval,boolean Abeep)
	{
		this.interval = Ainterval;
		this.beep = Abeep;
	}
	public void start()
	{
		ActionListener listener = new TimePrinter();
		Timer t = new Timer(interval, listener);
		t.start();	
	}
	public  class TimePrinter implements ActionListener
	{
		public void actionPerformed(ActionEvent event)
		{
			Date now = new Date();
			System.out.println(now);
			if (TalkingClock.this.beep) Toolkit.getDefaultToolkit().beep();
		}
	}
}

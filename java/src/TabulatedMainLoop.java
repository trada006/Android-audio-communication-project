
public class TabulatedMainLoop 
{

	/** Runs the main filter loop, this is where most time is spent */
	public static Complex mainFilterLoop(int mfrom, int mto, int mstep, int ifrom, int istep, Complex[] rmem, Complex[] pulsetable) 
	{
	    Complex sum=new Complex(0, 0);
	    //mainloop, this will greatly benefit from a mac instruction
	    for (int m = mfrom, i = ifrom; m <= mto; m += mstep, i += istep)
	    {
	    	sum=sum.add(rmem[m].multiply(pulsetable[i]));
	    }
	        
	    return sum;
	}

	/** 
	 * Runs the main filter loop. This assumes you have arranged memory so that the pulsetable 
	 * can be stepped through with unit increments.
	 */
	public static Complex mainBankedFilterLoop(int mfrom, int mto, int mstep, Complex[] rmem, Complex[] pulsetable) 
	{
	    Complex sum=new Complex(0, 0);
	    //mainloop, this will greatly benefit from a mac instruction
	    for (int m = mfrom, i = 0; m <= mto; m += mstep, i++)
	    {
	    	sum=sum.add(rmem[m].multiply(pulsetable[i]));
	    }
	        
	    return sum;
	}
	
	
}

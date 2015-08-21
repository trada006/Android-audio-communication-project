

public class Naive<Pulse extends FinitePulse> implements TimeOffsetEstimator 
{
	
	public Pulse g;
	public double T;
	public double Ts;
	public double taumin;
	public double taumax;
	public int c;
	public int L;
	public double Delta;
	public Complex[] rmem;
	public int rmemsize;
	public double DEFAULTBRENTTOL;
	public int[] P;
	public int[] D;
	public Complex[] pilots;
	double brenttol = DEFAULTBRENTTOL;
	
	

	
	public Naive(int[] P_, int[] D_, Complex[] pilots_, Pulse g_, double T_, double Ts_, double taumin_, double taumax_, int c_, double brenttol_) 
    {
  
    	
		L=D_.length + P_.length;
		Delta=T_ / c_;
		P=P_;
		D=D_;
		pilots=pilots_;
		g=g_;
		T=T_;
		Ts=Ts_;
		taumin=taumin_;
		taumax=taumax_;
		c=c_;
		brenttol=brenttol_*T;
    	//if(taumax <= taumin) throw "Maximum time offset taumax must be larger than minimum time offset taumin";
    }
 
    
    
    public double estimate(Complex[] r) 
    {
        //setupr(r);
        double tautilde = coarseMaximiseSS();
        double tauhat = refineCoarseEstimate(tautilde);
        return tauhat;
    }
    
    
    /** Obtain a coarse estimate of the time offset */
    public double coarseMaximiseSS() 
    {
        double maxSS = -1.0;
        double taubest = taumin;
        for (double tau = taumin; tau <= taumax; tau += Delta) {
            double thisSS = SS(tau);
            if (maxSS < thisSS) {
                maxSS = thisSS;
                taubest = tau;
            }
        }
        return taubest;
    }
    
    
    
    
    
    /** Refines the coarse estimate. */
    public double refineCoarseEstimate(double tautilde)
    {
        double a = tautilde - Delta;
        double c = tautilde + Delta;
             
        SingleVariateFunction f=new SingleVariateFunction()
        {
        	public double value(double tau)
        	{
        		return SS(tau);
        	}
        };
        
        Brent opt=new Brent(f, a, tautilde, c, brenttol);
        return opt.xmin();
    }
    
    
    
    
    /** Inner product between the received signal and pulse g */
    public final Complex m(double tau)
    {
        int A = (int) Math.ceil((g.tmin() + tau) / Ts); //ceil(x) returns the smallest integer greater than the value x, eg ceil(2.3) = 3
        int B = (int) Math.floor((g.tmax() + tau) / Ts); //floor(x) returns the largest integer not greater than x, eg floor(5.6) = 5
        Complex sum=new Complex(0, 0);
        for (int n = A; n <= B; n++) sum=sum.add( r(n).multiply((g.pulse(n * Ts - tau).conjugate())));
        return sum.times(Ts); //rescale discrete inner products by sample period.  This is not in the original paper, but it results in the same estimator
    }
    
    
    
    
    /** The Y function computing a correlation with the data and pilots */
    public final Complex Y(double tau)
    {
        Complex sum=new Complex(0, 0);
        
        for (int i = 0; i < P.length; i++) 
    	{
    		sum.add(m((P[i] + 1) * T + tau).multiply(pilots[i].conjugate()));
    	}
        return sum;
    }
    
    
    
    
    
    
    
    
    /** The amplitude accumulating Z function */
    public final double Z(double tau)
    {
        double sum = 0.0;
        for (int i : D) {
            sum += m((i + 1) * T + tau).abs();
        }
        return sum;
    }
    
    
    
    
    
    
    
    
    
    /** The objective function */
    public final double SS(double tau)
    {
        return Z(tau) + Y(tau).abs();
    }

    public final String name()
    {
        return "Naive";
    }
    
    
    
    /** The received sequence, zeros returned for unknown values */
    public final Complex r(int n)
    {
        if (n > 0 && n <= rmemsize) 
    	{
    		return rmem[n - 1];
    	}
        else return new Complex(0, 0);
    }
    
    
    
    public void setupr(Complex[] r) 
    {
        rmem = r;
        rmemsize = r.length;
    }
    
    
    
    
    
}

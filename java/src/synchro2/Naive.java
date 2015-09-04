package synchro2;

/** Naive implementation of the time offset estimator that computes the objective function directly */
 class Naive<Pulse extends FinitePulse> implements TimeOffsetEstimator {

	public final double DEFAULTBRENTTOL = 1e-7;
	
	public double brenttol_ = DEFAULTBRENTTOL; 	
	protected static int[] P=new int[0];
	protected static int[] D=new int[0];
	protected static Complex[] pilots=new Complex[0];    
	//public Pulse g;
	public static double T;
	public static double Ts;
	public static double taumin;
	public double taumax;
	public static int c;
	public double brenttol;
    /** Total number of transmitted symbols */
	public int L;
    /** Grid search width */
	public static double Delta;
    
	public TruncatedSincPulse<Pulse> g=new TruncatedSincPulse<Pulse>();


    /** Pointer to current received signal */
    protected static Complex[] rmem=new Complex[0];
    /** size of current received signal */
    protected static int rmemsize;
    
	
	public Naive()
	{}
	
    @SuppressWarnings("unchecked")
	public Naive(int[] P_,
            int[] D_,
            Complex[] pilots_,
            Pulse g_,
            double T_,
            double Ts_,
            double taumin_,
            double taumax_,
            int c_) throws Exception
            {
            	g=(TruncatedSincPulse<Pulse>) g_;
				L=D_.length + P_.length;
				Delta=T_ / c_;
				P=P_;
				D=D_;
				pilots=pilots_;
				//this.g=g_;
				T=T_;
				Ts=Ts_;
				taumin=taumin_;
				taumax=taumax_;
				c=c_;
				brenttol=brenttol_*T;
				if(taumax <= taumin) throw new Exception("Maximum time offset taumax must be larger than minimum time offset taumin");
			}
    
    

     public double estimate(final Complex[] r) {
        setupr(r);
        double tautilde = coarseMaximiseSS();
        double tauhat = refineCoarseEstimate(tautilde);
        return tauhat;
    }

    /** Obtain a coarse estimate of the time offset */
     double coarseMaximiseSS() {
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
     public double refineCoarseEstimate(double tautilde){
        double a = tautilde - Delta;
        double c = tautilde + Delta;
        
        SingleVariateFunction f=new SingleVariateFunction()
        {
        	public double value(double tau)
        	{
        		return -SS(tau);//function to minimise
        	}
        };
        Brent opt=new Brent(f, a, tautilde, c, brenttol);
        return opt.xmin();
    }

     
 /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
     
    /** Inner product between the received signal and pulse g */
     public Complex m(double tau){
    	int A = (int) Math.ceil((g.tmin() + tau) / Ts);    
        int B = (int) Math.floor((g.tmax() + tau) / Ts);
        Complex sum=new Complex(0,0);
        
        for (int n = A; n <= B; n++) 
        {
    		sum=sum.add(r(n).multiply((g.pulse(n * Ts - tau).conjugate())));	
    	}
        return sum.times(Ts); //rescale discrete inner products by sample period.  This is not in the original paper, but it results in the same estimator
    }
     
   ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////  
     

    /** The Y function computing a correlation with the data and pilots */
     public final Complex Y(double tau){
        Complex sum=new Complex(0, 0);
        for ( int i = 0; i < P.length; i++) sum=sum.add(m((P[i] + 1) * T + tau).multiply(pilots[i].conjugate()));
        return sum;
    }

    /** The amplitude accumulating Z function */
     public double Z(double tau){
        double sum = 0.0;
        for (int i : D)sum += (m((i + 1) * T + tau).abs());            
        return sum;
    }

    /** The objective function */
     public final double SS(double tau){
        return Z(tau) + Y(tau).abs();
    }

     public String name(){
        return "Naive";
    }

    

    /** The received sequence, zeros returned for unknown values */
     protected final Complex r(int n){ 
        if (n > 0 && n <= rmemsize) return new Complex(rmem[n - 1]);
        else return new Complex(0, 0);	
    }
    
     protected void setupr(final Complex[] r) {
        rmem = r;
        rmemsize = r.length;
    }

    
    
}

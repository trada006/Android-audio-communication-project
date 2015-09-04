package synchro2;

public class TabulatedFastAlgorithm 
{
	public final static int DEFAULTTABLESIZE = 100000;
	
	static class TabulatedDirect<Pulse extends FinitePulse> extends FastAlgorithmNoStorev.FastAlgorithm<Pulse> {

	    /** Table stores values of the transmit pulse g.pulse */
	    protected static Complex[] pulsetable;
	    /**  Width in t between elements in the table. */
	    protected static double stepwidth;
	    /** Multiplier for indexing the table */
	    protected static int tabls;
	    
	    public FastAlgorithmNoStorev.FastAlgorithm<Pulse> FA;

	    
	    public TabulatedDirect()
	    {}
	    
	    
	    TabulatedDirect( int[] P,
	             int[] D,
	             Complex[] pilots,
	             Pulse g,
	             double T,
	             double Ts,
	             double taumin,
	             double taumax,
	              int c,
	              int p,
	              int q) throws Exception{
	              int mintabs = DEFAULTTABLESIZE; 
	              FA = new FastAlgorithmNoStorev.FastAlgorithm<Pulse>(P, D, pilots, g, T, Ts, taumin, taumax, c, p, q); 
	        //int d = lcm(q,c);
	        int d = c * FA.b; //we require the number of table elements per symbol to be a multiple of d
	        tabls = (mintabs + 1) / d;
	        int tablesizepersymbol = tabls*d;
	        stepwidth = T / tablesizepersymbol;
	        for (double t = g.tmin(); t <= g.tmax(); t += stepwidth) pulsetable=Utili.Append(pulsetable, g.pulse(t).conjugate());
	    }

	    /** Inner product between the received signal and pulse g */
	      public final Complex m(double tau)  {

	        //step loop counters.
	        int A = (int) Math.ceil((FA.g.tmin() + tau) / FA.Ts);
	        int B = (int) Math.floor((FA.g.tmax() + tau) / FA.Ts);
	        int nfrom = Math.max(1, A);
	        int nto = Math.min((int) FA.rmemsize, B);
	        double startt = nfrom * FA.Ts - tau;
	        int ifrom = (int) Math.round((startt - FA.g.tmin()) / this.stepwidth);
	        int istep = FA.a*tabls;
	        int mfrom = nfrom - 1;
	        int mto = nto - 1;
	        int mstep = 1;

	        
	        return TabulatedMainLoop.mainFilterLoop(mfrom, mto, mstep, ifrom, istep, FA.rmem, pulsetable).times(FA.Ts); //rescale discrete inner products by sample period.  This is not in the original paper, but it results in the same estimator
	    }

	    /** The h_{\ell,k} sequence from the paper, computed by convolution*/
	      public final Complex h(int ell, int k)  {

	        //step loop counters.
	        double A = 1 - (FA.g.tmax() + FA.taumin) / FA.Delta + ((double) ell) / FA.b;
	        double B = 1 - (FA.g.tmin() + FA.taumin) / FA.Delta + ((double) ell) / FA.b;
	        int Bprime = (int) Math.ceil((k - B + ell * FA.n0) / FA.a);
	        int Aprime = (int) Math.floor((k - A + ell * FA.n0) / FA.a);
	        int Bpp = FA.a * Bprime - ell * FA.n0;
	        int App = FA.a * Aprime - ell * FA.n0;
	        int nfrom = Math.max((int) (ell / FA.b + 1), Bpp);
	        //int nto = min((int) (a * (rmemsize + 1) / b), App);
	        int nto = App;
	        double startt = -(k - nfrom - 1) * FA.Delta - FA.taumin + (ell * FA.Ts) / FA.a;
	        int ifrom = (int) Math.round((startt - FA.g.tmin()) / stepwidth);
	        int istep = FA.a*FA.b*tabls;
	        int mfrom = (FA.b * nfrom + ell) / FA.a - 1;
	        //int mto = (b * nto + ell) / a - 1;
	        int mto = Math.min((int)FA.rmemsize,(int)((FA.b*nto+ell)/FA.a))-1;
	        int mstep = FA.b;

	        return TabulatedMainLoop.mainFilterLoop(mfrom, mto, mstep, ifrom, istep, FA.rmem, pulsetable);
	    }

	    /** Return the values of Z computed on the grid taumin to taumax by Ts/c.  Direct computation. */
	     void fillZgrid() {
	        for ( int k = 1; k <= FA.K; k++) {
	            double sum = 0.0;
	            for (int i = 0; i < FA.D.length; i++) sum += FA.v(k + FA.c * (FA.D[i] + 1)).abs();
	            //for (int i : this.D) sum += Math.abs(bfunc(k + this.c * (i + 1)));
	            FA.Zgrid[k - 1] = sum;
	        }
	    }

	    /* Return the values of Y computed on the grid taumin to taumax by T/c. Direct convolution */
	     void fillYgrid() {
	        for ( int k = 1; k <= FA.K; k++) {
	            Complex sum = new Complex(0.0);
	            for ( int i = 0; i < FA.P.length; i++)
	                sum =sum.add(FA.v(k + FA.c * (FA.P[i] + 1)).multiply(FA.pilots[i].conjugate()));
	            FA.Ygridmag[k - 1] = sum.abs();
	        }
	    }

	     public String name()  {
	        return "TabulatedDirect";
	    }
	}
	     
	     
	     
	     
	     
	     
	     
	     
	     /**
	      * Recursive estimator of time offset.  The data symbol indices must now be 
	      * contiguous.  Exception is thrown if they are not.
	      */
	     class TabulatedRecursive<Pulse extends FinitePulse> extends TabulatedDirect<Pulse> {
	    	 
	    	 /** Minimum data symbol index */
	         protected int Dmin;
	         /** Maximum data symbol index */
	         protected int Dmax;
	    	 TabulatedDirect<Pulse> TD;
	    	 
	    	 public TabulatedRecursive()
	    	 {}
	    	 
	         public TabulatedRecursive( int[] P,
	                  int[] D,
	                  Complex[] pilots,
	                  Pulse g,
	                  double T,
	                  double Ts,
	                  double taumin,
	                  double taumax,
	                   int c,
	                   int p,
	                   int q) throws Exception
	         {
	        	 
	         TD = new TabulatedDirect<Pulse>(P, D, pilots, g, T, Ts, taumin, taumax, c, p, q);
	         Dmin=Utili.min(D);
	         Dmax=Utili.max(D); 
	             //check D is contiguous and sorted
	             for ( int i = 0; i < D.length - 1; i++) {
	                 if (D[i] != (D[i + 1] - 1)) {
	                     String str="Exception constructing Recursive. The data symbols must be contiguous.";
	                     throw new Exception(str);
	                 }
	             }
	         }

	         /* 
	          * Return the values of Z computed on the grid taumin to taumax by T/c. 
	          * Uses recursive algorithm that only applies when D is contiguous.
	          */
	          public void fillZgrid() {
	     		//cout << " yeyysudsudh   hello!!!" << endl;
	             for ( int k = 1; k <= this.c; k++) {
	                 double sum = 0.0;
	                 for (int i = 0; i < this.D.length; i++) sum += this.v(k + this.c * (this.D[i] + 1)).abs();
	                 this.Zgrid[k - 1] = sum;
	                 for (int m = 0; k + (m + 1) * this.c <= this.K; m++)
	                     this.Zgrid[k - 1 + (m + 1) * this.c] = this.Zgrid[k - 1 + m * this.c] - this.v(k + m * this.c + (Dmin + 1) * this.c).abs() + this.v(k + m * this.c + (Dmax + 1) * this.c + this.c).abs();
	             }
	         }

	          public final String TRname()  {
	             return "TabulatedRecursive";
	         }    
	     }
	     
	     
	     
	     
	     
	     
	     
	     
	     
	     
	     
	     /**
	      * Recursive estimator of time offset.  Does no refinement, i.e. no Brent's method.  Just here
	      * for benchmarking.
	      */
	
	     class TabulatedRecursiveNoRefine<Pulse extends FinitePulse> extends TabulatedRecursive<Pulse> {

	    	 public TabulatedRecursiveNoRefine()
	    	 {}
	    	 
	    	 private TabulatedRecursive<Pulse> TR;
	    	 
	         public TabulatedRecursiveNoRefine(int[] P,
	                  int[] D,
	                  Complex[] pilots,
	                  Pulse g,
	                  double T,
	                  double Ts,
	                  double taumin,
	                  double taumax,
	                   int c,
	                   int p,
	                   int q) throws Exception{
	         TR = new TabulatedRecursive<Pulse>(P, D, pilots, g, T, Ts, taumin, taumax, c, p, q); 
	         }

	   
	         public final double refineCoarseEstimate(double tautilde)  {
	             return tautilde;
	         }
	         
	         public final String name()  {
	             return "TabulatedRecursiveNoRefine";
	         }

	     }
	     

	}
	


package synchro2;

/* 
 * File:   FastAlgorithm.h
 * Author: Robby McKilliam
 *
 * Created on 4 February 2013, 8:13 PM
 */

/**
 * Override Naive with a faster way to compute the v vectors.  This does not store the v_k sequence.
 */
//P should be a FinitePulse
class FastAlgorithmNoStorev<Pulse extends FinitePulse> extends Naive<Pulse> {


	
	public int p;
	public int q;
	public static int a;
	public static int b;
	public static int n0;
	public int m0;
	/** Size of the search grid */
	public static int K;
	public static double[] Zgrid;//=new double[0];
	public static double[] Ygridmag;//=new double[0];
	public static double[] SSgrid;//=new double[0];
	public static Naive<?> naive;
	
	
	public FastAlgorithmNoStorev()
	{}
	
    public FastAlgorithmNoStorev( final int[] P,
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
          
    	
    naive=new Naive<Pulse>(P, D, pilots, g, T, Ts, taumin, taumax, c);
    
	//new Naive<Pulse>(P, D, pilots, g, T, Ts, taumin, taumax, c);
    this.p=p;
    this.q=q;
    K=(int)Math.floor((taumax - taumin) / Delta); //2147483647 from testdirect
    Zgrid=Utili.initializeDoubleArray(Zgrid, K, 0.0);
    this.Ygridmag=Utili.initializeDoubleArray(Ygridmag, K, 0.0);
    this.SSgrid=Utili.initializeDoubleArray(SSgrid, K, 0.0);
    
    
    
        //setup a, b, n0 and m0 for the polyphase filter
        int d = Utili.gcd(q, c * p);
        a = c * p / d;
        b = q / d;
        Utili.extended_gcd(b, a, d, n0, m0);
        n0=Utili.x_;
        m0=Utili.y_;
        d=Utili.gcd_;
    }

     public double coarseMaximiseSS() {
        fillSSgrid(); //compute objective function on the grid
        int khat = Utili.indexOfLargest(SSgrid); //get position of max element
        return taumin + khat*Delta;
    }

    /** fills the vector SSgrid with values of the objective function SS */
    public void fillSSgrid() {
        Direct.fillZgrid();
        Direct.fillYgrid();
        for ( int i = 0; i < K; i++) SSgrid[i] = Zgrid[i] + Ygridmag[i];
    }

    /** zero filled received signal z_n in the paper */
     public Complex zn(int n)  {
        if (n % a == 0) return r(n / a);
        else return new Complex(0, 0);
    }

    /** banked received signal z_{ell,n} from the paper*/
     public Complex z(int ell, int n)  { 
        return zn(b * n + ell);
    }

    /** The sequence g_{ell,n} from the paper */
     public Complex gfunc(int ell, int n)  {
        return naive.g.pulse(-(n - 1) * Delta - taumin + (ell * Ts) / a).conjugate();
    }

    /** The sequence bk from the paper */
     public Complex v(int k)  {
    	 
        Complex sum=new Complex(0, 0);
        for ( int ell = 0; ell < b; ell++) 
    	{
    		sum=sum.add(h(ell, k));
    	}
        return sum.times(Ts); //rescale discrete inner products by sample period.  This is not in the original paper, but it results in the same estimator
    }
     
     
     
     
     
     
     //ERROR generic g is null
     ///////////////////////////////////////////////////////////////////////////////////////////////////////
     

    /** The h_{\ell,k} sequence from the paper, computed by convolution*/
    public Complex h(int ell, int k)  {
    	
    	
        double A = 1 - (naive.g.tmax() + taumin) / Delta + ((double) ell) / b;
        double B = 1 - (naive.g.tmin() + taumin) / Delta + ((double) ell) / b;
        
        
        
        
        int Bprime = (int) Math.ceil((k - B + ell * n0) / a);
        
        
        int Aprime = (int) Math.floor((k - A + ell * n0) / a);
        
        
        Complex sum=new Complex(0, 0);
        int Bpp = a * Bprime - ell * n0;
        int App = ((int) a) * Aprime - ell * n0;
        for (int n = Bpp; n <= App; n += a) sum =sum.add( z(ell, n).multiply(gfunc(ell, k - n)));
        return sum;
    }

    
    
    
    
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////
    
    
    
    
    
    
    
    
    
    
    
    
    /** fill vector Zgrid with discretized Z function */
    public static void fillZgrid2(){};

    /** fill vector Ygridmag with discretized magnitude of the Y function */
    public static void fillYgrid2(){};

    //output for testing

    public double[] getZgrid() {
        return Zgrid;
    }

    public double[] getYgrid() {
        return Ygridmag;
    }

    public double[] getSSgrid() {
        return SSgrid;
    }

   
    
    
    /**
     * Fast polyphase computer for bk sequence.  Stores the bk sequence for fast access.
     */

//=================================================================================================================================================================

    static class FastAlgorithm<Pulse extends FinitePulse> extends FastAlgorithmNoStorev<Pulse> {

    	
    	static protected Complex[] vstore=new Complex[0];
        protected static int boffset;
    	
    	public FastAlgorithm()
    	{}
    	
        public FastAlgorithm( int[] P,
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
        	FastAlgorithmNoStorev<Pulse> fans = new FastAlgorithmNoStorev<Pulse>(P, D, pilots, g, T, Ts, taumin, taumax, c, p, q); 
    	    int Dmin = Utili.min(D);
    		int Dmax = Utili.max(D);
    		int Pmin = Utili.min(P);
    		int Pmax = Utili.max(P);
    		
    		vstore=Utili.resize(vstore, K + c*(Math.max(Dmax, Pmax) - Math.min(Dmin, Pmin)), 0.0);
    		boffset = 1 + c*(Math.min(Dmin, Pmin) + 1);
        }

         public double estimate(Complex[] r) {   	 
            setupr(r);
            for ( int k = 0; k < vstore.length; k++) 
            	{
            		vstore[k] = v(k + boffset);
            	}
            
            
            double tautilde = this.coarseMaximiseSS();
            double tauhat = this.refineCoarseEstimate(tautilde);   		
            return tauhat;
        }

        //b takes from array
         public static Complex vv(int k)  {
            //if(k - boffset < 0 || k - boffset >= bstore.length) std::cout << "Indexing out of bounds!" << std::endl; 
            return vstore[k - boffset];
            
        }

        

    }
    
    //=======================================================================================================================================================

    static class Direct<Pulse extends FinitePulse> extends FastAlgorithm<Pulse> {


    	public Direct()
    	{}
    	
        public Direct( int[] P,
                 int[] D,
                 Complex[] pilots,
                 Pulse g,
                 double T,
                 double Ts,
                 double taumin,
                 double taumax,
                  int c,
                  int p,
                  int q) throws Exception {
        
        	new FastAlgorithm<Pulse>(P, D, pilots, g, T, Ts, taumin, taumax, c, p, q); 

    	}	

        

         /**Return the values of Z computed on the grid taumin to taumax by Ts/c.  Direct computation.*/ 
         public static void fillZgrid() {

            for ( int k = 1; k <= K; k++) {
                double sum = 0.0;    
                for (int i : D)
                {
                	sum += vv(k + c * (i + 1)).abs();
                }          
                Zgrid[k - 1] = sum;
            }
        }

        /** Return the values of Y computed on the grid taumin to taumax by T/c. Direct convolution */
         public static void fillYgrid() {
            for ( int k = 1; k <= K; k++) {
                Complex sum = new Complex(0.0);
                for ( int i = 0; i < P.length; i++)
                    sum =sum.add(vv(k + c * (P[i] + 1)).multiply(pilots[i].conjugate()));
                Ygridmag[k - 1] = sum.abs();
            }
        }

         public String name()  {
            return "Direct";
        }
    }



    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    /**
     * Recursive estimator of time offset.  The data symbol indices must now be 
     * contiguous.  Exception is thrown if they are not.
     */

    class Recursive<Pulse extends FinitePulse> extends Direct<Pulse> {
   
    	
    	  /** Minimum data symbol index */
        final protected int Dmin;
        /** Maximum data symbol index */
        final protected int Dmax;
    	
    	
        public Recursive(int P,
                int[] D,
                double[] pilots,
                Pulse g,
                double T,
                double Ts,
                double taumin,
                double taumax,
                int c,
                int p,
                int q) throws Exception{
     
        //new Direct<Pulse>(P, D, pilots, g, T, Ts, taumin, taumax, c, p, q);
        Dmin=Utili.min(D);
        Dmax=Utili.max(D); 
            //check D is contiguous and sorted
            for (int i = 0; i < D.length - 1; i++) {
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
        public void DirectfillZgrid() {
            for (int k = 1; k <= this.c; k++) {
                double sum = 0.0;
                for (int i : this.D) sum += this.v(k + this.c * (i + 1)).abs();
                this.Zgrid[k - 1] = sum;
                for (int m = 0; k + (m + 1) * this.c <= this.K; m++)
                    this.Zgrid[k - 1 + (m + 1) * this.c] = this.Zgrid[k - 1 + m * this.c] - this.v(k + m * this.c + (Dmin + 1) * this.c).abs() + this.v(k + m * this.c + (Dmax + 1) * this.c + this.c).abs();
            }
        }

        public final String name(){
            return "Recursive";
        }
        
        
        
        
        
        
        
        
        
        
        /**
         * Estimator of time offset that only uses pilot symbols
         */
        class PilotsOnly<Pulse extends FinitePulse> extends Direct<Pulse> {
      
            public PilotsOnly(int[] P,
                    Complex[] pilots,
                    Pulse g,
                    double T,
                    double Ts,
                    double taumin,
                    double taumax,
                    int c,
                    int p,
                    int q) throws Exception{
            new Direct<Pulse>(P, P, pilots, g, T, Ts, taumin, taumax, c, p, q); 
            }

            /** Does nothing! There are only pilots */
        	public void PilotsOnlyfillZgrid(){}
            
            /** Return zero! There are only pilots */
            public final double Z(double tau){
                return 0;
            }
            
            /** 
             * This replaces the pilots, copy new ones.  This is really just a hack to make simulations
             * easier!  Ideally this would not exist and pilots would be const
             * @throws Exception 
             */
            public void setPilots(Complex[] pilots) throws Exception{
                if(this.pilots.length != pilots.length) throw new Exception("You shouldn't be changing the number of pilots");
                for(int i = 0; i < pilots.length; i++) this.pilots[i] = pilots[i];
            }
            
            public final String name(){
                return "PilotsOnly";
            }
           
        }
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
          
    }
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    

}




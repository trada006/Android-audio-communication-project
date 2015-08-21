public class TruncatedSincPulse implements FinitePulse 
{
	public double tmin_;
    public double tmax_;
    public double T_;
    public double sqrtT;
	

    public TruncatedSincPulse(double T, int numzeros) 
    {
		 tmin_=-T*numzeros;
		 tmax_=T*numzeros;
		 T_=T;
		 sqrtT=Math.sqrt(T); 
    }

    /** The truncated sinc pulse */
    public final Complex pulse(double t)
    {
        if (t > tmin_ && t < tmax_) 
    	{
    		return new Complex(Utili.sinc(t / T_) / sqrtT);
    	}
        else return new Complex(0.0);
    }

    public final double tmin()
    {
        return tmin_;
    }

    public final double tmax()  
    {
        return tmax_;
    }

    public final double T()  {
        return T_;
    }


    
    public static class NormalisedFinitePulse<P> implements FinitePulse 
    {
    	
    	P p;
        double tmin_;
        double tmax_;
        double T_;
        double normalisingconstant;
    	
    	
        public NormalisedFinitePulse(P p_)
        { 
        	p=p_;
        	tmin_=((FinitePulse) p_).tmin();
        	tmax_=((FinitePulse) p_).tmax();
        	T_=((FinitePulse) p_).T();
        	
   
            SingleVariateFunction f=new SingleVariateFunction()
            {
            	public double value(double x)
            	{
            		return ((FinitePulse) p).pulse(x).abs();
            	}
            };
            
            double energy = Utili.trapezoidal(f, tmin_, tmax_, 100000); //last number is integration steps
            normalisingconstant = Math.sqrt(energy);
        }
       

        public final double  tmin()
        {
            return tmin_;
        }

        public final double tmax()
        {
            return tmax_;
        }

        public final double T()
        {
            return T_;
        }

        //The normalised pulse 
        public final Complex pulse(double t)
        {
            return new Complex(((FinitePulse) p).pulse(t).divide(normalisingconstant));
        }
    }
    
    
    
    
    public static class TruncatedRootRaisedCosine implements FinitePulse 
    {
    	protected double tmin_;
        protected double tmax_;
        protected double T_;
        protected double beta;
    	

        public TruncatedRootRaisedCosine(double T, double beta, double duration, boolean dodgyflag)
        { 
        	T_=T;
        	this.beta=beta;  
            final double stepsize = 0.01; //step taken whilst looking for zeros (this will work only if zeros are atleast stepsize apart)
            double c = 0.0;
            int csign = 1;
            while(c < duration/T/2) 
            {
                while (Utili.signum(rootraisedcosine(c)) == csign) 
            	{
            		c += stepsize;
            	}
                csign = -csign;
            }
            
            SingleVariateFunction f = new SingleVariateFunction()
            {

    			@Override
    			public double value(double x) 
    			{
    				return rootraisedcosine(x);
    			}
            };
     
            double nt=new Utili.Bisection(f, c-stepsize, c, 1e-7).zero();
            
      
            tmin_ = -nt*T;
            tmax_ = nt*T;
            double dur = tmax_ - tmin_;
            //std::cout << dur << " but request duration is " << duration << std::endl;
            if( dur > 2*duration || dur < duration/2 ) 
            { //check that the duration obtained is reasonable
                //StringBuffer str;
                System.out.println( "Something when wrong with the bisection method, duration is ");
                System.out.println(" but request duration is " + duration);
                //throw str.toString();
            }
        }
        
        public TruncatedRootRaisedCosine(double T, double beta, int numzeros) 
        {
    	    this.T_=T;
    	    this.beta=beta; 
    	    final double stepsize = 0.01; //step taken whilst looking for zeros (this will work only if zeros are atleast stepsize apart)
            double c = 0.0;
            int csign = 1;
            for (int i = 1; i <= numzeros; i++) 
            {
                while (Utili.signum(rootraisedcosine(c)) == csign) c += stepsize;
                csign = -csign;
            }

            SingleVariateFunction f = new SingleVariateFunction()
            {
    			@Override
    			public double value(double x) 
    			{
    				return rootraisedcosine(x);
    			}
            };
            
            double nt = new Utili.Bisection(f, c - stepsize, c, 1e-7).zero();
            tmin_ = -nt*T;
            tmax_ = nt*T;
        }
        

        public final Complex pulse(double t)
        {
            if (t > tmin_ && t < tmax_) return new Complex(rootraisedcosine(t / T_));
            else return new Complex(0.0);
        }

        public final double tmin() {
            return tmin_;
        }

        public final double tmax(){
            return tmax_;
        }

        public final double T() {
            return T_;
        }

        /** A root raised cosine with period 1 and rolloff beta */
        public final double rootraisedcosine(double t)
        {
        	double pi=Utili.pi;
        	
            double abst = Math.abs(t); //pulse is symmetric about zero, so just use magnitude of t
            if (abst < 5e-3) { //second order expansion if t is near zero
                double term0 = 1 + beta * (4 / pi - 1);
                double term2 = (Utili.cub((beta - 1) * pi) + 96 * beta * beta * (4 * beta + pi - beta * pi) - 12 * beta * Utili.sqr(pi + beta * pi)) / 6 / pi;
                return term0 + term2 * abst*abst;
            }
            if (Math.abs(abst - 1.0 / 4 / beta) < 5e-4) { //first order expansion if t is near 1/(4beta)
                double a = (1 + beta) * pi / 4 / beta;
                double term0 = beta * (Math.sin(a) - 2 * Math.cos(a) / pi);
                double term1 = beta * ((12 * beta + pi * pi) * Math.cos(a) + 2 * (1 - 2 * beta) * pi * Math.sin(a)) / pi;
                return term0 + term1 * (abst - 1.0 / 4 / beta);
            } else { //otherwise use direct formula
                double a = pi*t;
                double b = beta*t;
                return (Math.sin(a * (1 - beta)) + 4 * b * Math.cos(a * (1 + beta))) / (a * (1 - 16 * b * b));
            }
        }


        
    }
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
}


import org.junit.After;
import org.junit.AfterClass;
import org.junit.Before;
import org.junit.BeforeClass;
import org.junit.Test;



import static org.junit.Assert.*;



public class FinitePulseTest
{
	static final double pi = 3.141592653589793238463;
	
	
	public FinitePulseTest() {
    }
    
    @BeforeClass
    public static void setUpClass() {
    }
    
    @AfterClass
    public static void tearDownClass() {
    }
    
    @Before
    public void setUp() {
    }
    
    @After
    public void tearDown() {
    }
	
	
    @Test
    public void finiteSinc4zeros() 
    {
        double tol = 1e-9;
        double nzs = 4;
        double T = 1;
        
        TruncatedSincPulse sincpulse=new TruncatedSincPulse(T, (int)nzs);

        assertTrue((sincpulse.tmin() - -T * nzs) < tol);
        assertTrue((sincpulse.tmax() - T * nzs) < tol);
        assertTrue((sincpulse.pulse(-4.4).subtract(new Complex(0.0, 0.0)).abs()) < tol); //output from scala
        assertTrue(((sincpulse.pulse(-3.9).subtract(new Complex(-0.025221324181627227, 0.0))).abs()) < tol);
        assertTrue((sincpulse.pulse(-3.4).subtract(new Complex(-0.08903843866360674, 0.0)).abs()) < tol);
        assertTrue((sincpulse.pulse(-2.9).subtract(new Complex(0.03391833252011936, 0.0)).abs()) < tol);
        assertTrue((sincpulse.pulse(-2.4).subtract(new Complex(0.12613778810677617, 0.0)).abs()) < tol);
        assertTrue((sincpulse.pulse(-1.9).subtract(new Complex(-0.05177008647807704, 0.0)).abs()) < tol);
        assertTrue((sincpulse.pulse(-1.9).subtract(new Complex(-0.05177008647807704, 0.0)).abs()) < tol);
        assertTrue((sincpulse.pulse(0.0).subtract(new Complex(1.0, 0.0)).abs())<tol);
        assertTrue((sincpulse.pulse(4e-3).subtract(new  Complex(0.9999973333354667, 0.0)).abs()) < tol);
        assertTrue((sincpulse.pulse(1.4).subtract(new Complex(-0.21623620818304484, 0.0)).abs()) < tol);       
    }
    
    
    
    @Test
    public void testRootRaisedCosine() 
    {
 
        double tol = 1e-6;
        double nzs = 4;
        double T = 1.0;
        double beta = 0.5;
        
        TruncatedSincPulse.TruncatedRootRaisedCosine p=new TruncatedSincPulse.TruncatedRootRaisedCosine(T, beta, (int)nzs);

        double[] expected = {2.0 / 15 / pi, -1.0 / 3 / Math.sqrt(2.0) / pi, -1.0 / 3 / pi, (pi + 2) / 2 / Math.sqrt(2.0) / pi, 0.5 + 2 / pi,
            (pi + 2) / 2 / Math.sqrt(2.0) / pi, -1.0 / 3 / pi, -1.0 / 3 / Math.sqrt(2.0) / pi, 2.0 / 15 / pi}; //from Mathematica
        
        double[] test=new double[9];
        int i=0;
  
    	for (double t = -2.0; t <= 2.0; t =t+ 0.5) 
    	{
    		test[i]=p.rootraisedcosine(t);
    		i++;
    	}

        boolean pass = false;
        for (int j = 0; j < test.length; j++)
        {
    		if (Math.abs(test[j] - expected[j]) < tol)
    		{
    			pass = true;
    			continue;
    		}
    		else pass = false;	
        }
        assertTrue(pass);
    }
    
    
    
    
    
    @Test
    public void testRootRaisedCosineFindZeros() 
    {
        final double tol = 1e-4;
        final double T = 1.0;
        final double beta = 0.5;
 
        assertTrue(Math.abs(0.873584 - new TruncatedSincPulse.TruncatedRootRaisedCosine(T, beta, 1).tmax()) < tol);
        assertTrue(Math.abs(1.69555 - new TruncatedSincPulse.TruncatedRootRaisedCosine(T, beta, 2).tmax()) < tol);
        assertTrue(Math.abs(-1.69555 - new TruncatedSincPulse.TruncatedRootRaisedCosine(T, beta, 2).tmin()) < tol);
        assertTrue(Math.abs(2.35734 - new TruncatedSincPulse.TruncatedRootRaisedCosine(T, beta, 3).tmax()) < tol);
        assertTrue(Math.abs(2.96409 - new TruncatedSincPulse.TruncatedRootRaisedCosine(T, beta, 4).tmax()) < tol);
        assertTrue(Math.abs(3.68054 - new TruncatedSincPulse.TruncatedRootRaisedCosine(T, beta, 5).tmax()) < tol);
        
    }
    
    
    @Test
    public void testRootRaisedCosineNearZero() 
    {
        double tol = 1e-7;
        double nzs = 4;
        double T = 1.0;
        double beta = 0.5;
        TruncatedSincPulse.TruncatedRootRaisedCosine p=new TruncatedSincPulse.TruncatedRootRaisedCosine(T, beta, (int)nzs);
       
        assertTrue(Math.abs(1.136619772 - p.rootraisedcosine(0.0)) < tol); //from mathematica
        assertTrue(Math.abs(1.136617045 - p.rootraisedcosine(1e-3)) < tol); //from mathematica
        assertTrue(Math.abs(1.136617045 - p.rootraisedcosine(-1e-3)) < tol); //from mathematica
        assertTrue(Math.abs(1.136576129 - p.rootraisedcosine(-4e-3)) < tol); //from mathematica
        assertTrue(Math.abs(1.136576129 - p.rootraisedcosine(4e-3)) < tol); //from mathematica
        
    }
    
    
    
    /*
    public void normalisedTruncatedSinc() 
    {
        double tol = 1e-4;
        int nzs = 4;
        double T = 1;
        
        
        SingleVariateFunction f = new SingleVariateFunction()
        {
        	@Override
			public double value(double x) 
			{
        		return Math.abs(p.pulse(t));
			}
        }
        
        NormalisedFinitePulse<TruncatedSincPulse> p=new NormalisedFinitePulse<TruncatedSincPulse>(new TruncatedSincPulse(T, nzs));
        
        
        
        auto f = [&p](double t){return std::norm(p.pulse(t));};
        
        
        
        
        double energy1 = trapezoidal(f, p.tmin() - 0.3, p.tmax() + 0.3, 100000); //compute pulse energy
        double energy2 = trapezoidal(f, p.tmin(), p.tmax(), 100000); //compute pulse energy
        bool pass = true;
        assertTrue(Math.abs(1.0 - energy1) < tol;
        assertTrue(Math.abs(1.0 - energy2) < tol;
    }
    */
    
    
    
    
    @Test
    public void finiteSinc2zeros() 
    {
        double tol = 1e-9;
        double nzs = 2;
        double T = 1;
        TruncatedSincPulse sincpulse=new TruncatedSincPulse(T, (int)nzs);
        assertTrue((sincpulse.tmin() - -T * nzs) < tol);
        assertTrue((sincpulse.tmax() - T * nzs) < tol);
        assertTrue((sincpulse.pulse(1.2).subtract(new Complex(-0.15591488063143982, 0.0))).abs() < tol); //output from scala
        assertTrue((sincpulse.pulse(1.0).subtract(new Complex(0.0, 0.0))).abs() < tol); //output from scala
        assertTrue((sincpulse.pulse(0.8).subtract(new Complex(0.23387232094715982, 0.0))).abs() < tol); //output from scala
        assertTrue((sincpulse.pulse(0.6).subtract(new Complex(0.5045511524271047, 0.0))).abs() < tol);
        assertTrue((sincpulse.pulse(0.4).subtract(new Complex(0.756826728640657, 0.0))).abs() < tol);
       
    }
    
    
    @Test
    public void testRootRaisedCosineNearBeta4() 
    {
        double tol = 1e-6;
        double nzs = 4;
        double T = 1.0;
        double beta = 0.5;
        double t = 1.0 / 4 / beta;
        TruncatedSincPulse.TruncatedRootRaisedCosine p=new TruncatedSincPulse.TruncatedRootRaisedCosine(T, beta, (int)nzs);
        
        assertTrue(Math.abs(0.5786324696 - p.rootraisedcosine(t)) < tol); //from Mathematica
        assertTrue(Math.abs(0.5786324696 - p.rootraisedcosine(-t)) < tol); //from Mathematica

        assertTrue(Math.abs(0.5784538720 - p.rootraisedcosine(t + 1e-4)) < tol); //from Mathematica
        assertTrue(Math.abs(0.5788110635 - p.rootraisedcosine(t - 1e-4)) < tol); //from Mathematica
        assertTrue(Math.abs(0.5784538720 - p.rootraisedcosine(-t - 1e-4)) < tol); //from Mathematica
        assertTrue(Math.abs(0.5788110635 - p.rootraisedcosine(-t + 1e-4)) < tol); //from Mathematica

        assertTrue(Math.abs(0.5793468225 - p.rootraisedcosine(t - 4e-4)) < tol); //from Mathematica
        assertTrue(Math.abs(0.5779180564 - p.rootraisedcosine(t + 4e-4)) < tol); //from Mathematica
        assertTrue(Math.abs(0.5793468225 - p.rootraisedcosine(-t + 4e-4)) < tol); //from Mathematica
        assertTrue(Math.abs(0.5779180564 - p.rootraisedcosine(-t - 4e-4)) < tol); //from Mathematica
    }
    
    
    
    @Test
    public void durationConstructedRRC() 
    {
   
        double tol = 1e-9;
        final double beta = 1.0/3.0;
        double duration = 10;
        int nzs = 7; //obtained by inspection.
        double T = 1;
        TruncatedSincPulse.TruncatedRootRaisedCosine  rrcpulse = new TruncatedSincPulse.TruncatedRootRaisedCosine(T, beta, duration, true);
        TruncatedSincPulse.TruncatedRootRaisedCosine rrcpulsetester= new TruncatedSincPulse.TruncatedRootRaisedCosine(T, beta, nzs);
        boolean pass = false;
        for(double t = -duration; t < duration; t+=0.01) 
        {
        	if ((rrcpulse.pulse(t).subtract(rrcpulsetester.pulse(t))).abs() < tol)
        	{
        		pass = true;
        		continue;
        	}
        	else pass = false;	
        }
        assertTrue(pass);
    }
    
    
    
    
    
    
    
    
    
}

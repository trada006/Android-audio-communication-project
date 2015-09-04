package synchro2;

public class FinitePulseTest 
{

	/* 
	 * File:   FinitePulseTest.cpp
	 * Author: Robby McKilliam
	 *
	 * Created on 23/01/2013, 2:19:26 PM
	 */
	public final double pi=Math.PI;


	void finiteSinc4zeros() {
	    System.out.println( "finiteSinc with 4 zeros ... ");
	    double tol = 1e-9;
	    double nzs = 4;
	    double T = 1;
	    TruncatedSincPulse sincpulse=new TruncatedSincPulse(T, (int)nzs);
	    boolean pass = true;
	    pass &= (sincpulse.tmin() - -T * nzs) < tol;
	    pass &= (sincpulse.tmax() - T * nzs) < tol;
	    pass &= (sincpulse.pulse(-4.4).subtract(new Complex(0.0, 0.0)).abs()) < tol; //output from scala
	    pass &= (sincpulse.pulse(-3.9).subtract(new Complex(-0.025221324181627227, 0.0)).abs()) < tol;
	    pass &= (sincpulse.pulse(-3.4).subtract(new Complex(-0.08903843866360674, 0.0)).abs()) < tol;
	    pass &= (sincpulse.pulse(-2.9).subtract(new Complex(0.03391833252011936, 0.0)).abs()) < tol;
	    pass &= (sincpulse.pulse(-2.4).subtract(new Complex(0.12613778810677617, 0.0)).abs()) < tol;
	    pass &= (sincpulse.pulse(-1.9).subtract(new Complex(-0.05177008647807704, 0.0)).abs()) < tol;
	    pass &= (sincpulse.pulse(-1.9).subtract(new Complex(-0.05177008647807704, 0.0)).abs()) < tol;
	    pass &= (sincpulse.pulse(0.0).subtract(new Complex(1.0, 0.0)).abs()) < tol;
	    pass &= (sincpulse.pulse(4e-3).subtract(new Complex(0.9999973333354667, 0.0)).abs()) < tol;
	    pass &= (sincpulse.pulse(1.4).subtract(new Complex(-0.21623620818304484, 0.0)).abs()) < tol;
	    if (pass) System.out.println( "PASS" );
	    else System.out.println( "FAIL" );
	}

	void normalisedTruncatedSinc() {
	    System.out.println( "normalised truncated sinc ... ");
	    double tol = 1e-4;
	     int nzs = 4;
	    double T = 1;
	    TruncatedSincPulse.NormalisedFinitePulse<TruncatedSincPulse> p=new TruncatedSincPulse.NormalisedFinitePulse<TruncatedSincPulse>(new TruncatedSincPulse(T, nzs));

	    
	    SingleVariateFunction f = new SingleVariateFunction()
        {
        	@Override
			public double value(double x) 
			{
        		return (p.pulse(x)).abs2();
			}
        };
	    
        
        
	    double energy1 = Utili.trapezoidal(f, p.tmin() - 0.3, p.tmax() + 0.3, 100000); //compute pulse energy
	    double energy2 = Utili.trapezoidal(f, p.tmin(), p.tmax(), 100000); //compute pulse energy
	    boolean pass = true;
	    pass &= (Math.abs(1.0 - energy1)) < tol;
	    pass &= (Math.abs(1.0 - energy2)) < tol;
	    if (pass) System.out.println( "PASS" ); 
	    else 
	    	{
	    		System.out.println( "FAIL" );
	    	}
	}

	void finiteSinc2zeros() {
	    System.out.println( "finiteSinc with 2 zeros ... ");
	    double tol = 1e-9;
	    double nzs = 2;
	    double T = 1;
	    TruncatedSincPulse sincpulse=new TruncatedSincPulse(T, (int)nzs);
	    boolean pass = true;
	    pass &= Math.abs(sincpulse.tmin() - -T * nzs) < tol;
	    pass &= Math.abs(sincpulse.tmax() - T * nzs) < tol;
	    pass &= Math.abs(sincpulse.pulse(1.2).subtract(new Complex(-0.15591488063143982, 0.0)).abs()) < tol; //output from scala
	    pass &= (sincpulse.pulse(1.0).subtract(new Complex(0.0, 0.0)).abs()) < tol; //output from scala
	    pass &= (sincpulse.pulse(0.8).subtract(new Complex(0.23387232094715982, 0.0)).abs()) < tol; //output from scala
	    pass &= (sincpulse.pulse(0.6).subtract(new Complex(0.5045511524271047, 0.0)).abs()) < tol;
	    pass &= (sincpulse.pulse(0.4).subtract(new Complex(0.756826728640657, 0.0)).abs()) < tol;
	    if (pass) System.out.println( "PASS" );
	    else System.out.println( "FAIL" );
	}

	void testRootRaisedCosineNearZero() {
	    System.out.println( "root raised cosine near zero... ");
	    double tol = 1e-7;
	    double nzs = 4;
	    double T = 1.0;
	    double beta = 0.5;
	    TruncatedSincPulse.TruncatedRootRaisedCosine p=new TruncatedSincPulse.TruncatedRootRaisedCosine(T, beta, (int)nzs);
	    boolean pass = true;
	    pass &= Math.abs(1.136619772 - p.rootraisedcosine(0.0)) < tol; //from mathematica
	    pass &= Math.abs(1.136617045 - p.rootraisedcosine(1e-3)) < tol; //from mathematica
	    pass &= Math.abs(1.136617045 - p.rootraisedcosine(-1e-3)) < tol; //from mathematica
	    pass &= Math.abs(1.136576129 - p.rootraisedcosine(-4e-3)) < tol; //from mathematica
	    pass &= Math.abs(1.136576129 - p.rootraisedcosine(4e-3)) < tol; //from mathematica
	    if (pass) System.out.println( "PASS" );
	    else System.out.println( "FAIL" );
	}

	void testRootRaisedCosineNearBeta4() {
	    System.out.println( "root raised cosine near beta/4... ");
	    double tol = 1e-6;
	    double nzs = 4;
	    double T = 1.0;
	    double beta = 0.5;
	    double t = 1.0 / 4 / beta;
	    TruncatedSincPulse.TruncatedRootRaisedCosine p=new TruncatedSincPulse.TruncatedRootRaisedCosine(T, beta, (int)nzs);
	    boolean pass = true;
	    pass &= Math.abs(0.5786324696 - p.rootraisedcosine(t)) < tol; //from Mathematica
	    pass &= Math.abs(0.5786324696 - p.rootraisedcosine(-t)) < tol; //from Mathematica

	    pass &= Math.abs(0.5784538720 - p.rootraisedcosine(t + 1e-4)) < tol; //from Mathematica
	    pass &= Math.abs(0.5788110635 - p.rootraisedcosine(t - 1e-4)) < tol; //from Mathematica
	    pass &= Math.abs(0.5784538720 - p.rootraisedcosine(-t - 1e-4)) < tol; //from Mathematica
	    pass &= Math.abs(0.5788110635 - p.rootraisedcosine(-t + 1e-4)) < tol; //from Mathematica

	    pass &= Math.abs(0.5793468225 - p.rootraisedcosine(t - 4e-4)) < tol; //from Mathematica
	    pass &= Math.abs(0.5779180564 - p.rootraisedcosine(t + 4e-4)) < tol; //from Mathematica
	    pass &= Math.abs(0.5793468225 - p.rootraisedcosine(-t + 4e-4)) < tol; //from Mathematica
	    pass &= Math.abs(0.5779180564 - p.rootraisedcosine(-t - 4e-4)) < tol; //from Mathematica
	    if (pass) System.out.println( "PASS" );
	    else System.out.println( "FAIL" );
	}

	void testRootRaisedCosine() {
	    System.out.println( "root raised cosine... ");
	    double tol = 1e-6;
	    double nzs = 4;
	    double T = 1.0;
	    double beta = 0.5;
	    TruncatedSincPulse.TruncatedRootRaisedCosine p=new TruncatedSincPulse.TruncatedRootRaisedCosine(T, beta, (int)nzs);
	    double[] expected = {2.0 / 15 / pi, -1.0 / 3 / Math.sqrt(2.0) / pi, -1.0 / 3 / pi, (pi + 2) / 2 / Math.sqrt(2.0) / pi, 0.5 + 2 / pi,
	        (pi + 2) / 2 / Math.sqrt(2.0) / pi, -1.0 / 3 / pi, -1.0 / 3 / Math.sqrt(2.0) / pi, 2.0 / 15 / pi}; //from Mathematica
	    double[] test = null;
	    for (double t = -2.0; t <= 2.0; t += 0.5) 
    	{
	    	test=Utili.Append(test, p.rootraisedcosine(t));
	    	//test.push_back(p.rootraisedcosine(t));
    	}
	    boolean pass = true;
	    for (int i = 0; i < test.length; i++)
	        pass &= Math.abs(test[i] - expected[i]) < tol; //from Mathematica
	    if (pass) System.out.println( "PASS" );
	    else System.out.println( "FAIL" );
	}

	void testRootRaisedCosineFindZeros() {
	    System.out.println( "root raised cosine find zeros ... ");
	    final double tol = 1e-4;
	    final double T = 1.0;
	    final double beta = 0.5;
	    //test against roots found by Mathematica
	    boolean pass = true;
	    pass &= Math.abs(0.873584 - new TruncatedSincPulse.TruncatedRootRaisedCosine(T, beta, 1).tmax()) < tol;
	    pass &= Math.abs(1.69555 - new TruncatedSincPulse.TruncatedRootRaisedCosine(T, beta, 2).tmax()) < tol;
	    pass &= Math.abs(-1.69555 - new TruncatedSincPulse.TruncatedRootRaisedCosine(T, beta, 2).tmin()) < tol;
	    pass &= Math.abs(2.35734 - new TruncatedSincPulse.TruncatedRootRaisedCosine(T, beta, 3).tmax()) < tol;
	    pass &= Math.abs(2.96409 - new TruncatedSincPulse.TruncatedRootRaisedCosine(T, beta, 4).tmax()) < tol;
	    pass &= Math.abs(3.68054 - new TruncatedSincPulse.TruncatedRootRaisedCosine(T, beta, 5).tmax()) < tol;
	    if (pass) System.out.println( "PASS" );
	    else System.out.println( "FAIL" );
	}

	void durationConstructedRRC() throws Exception {
	    System.out.println( "construct rrc pulse with duration ... ");
	    double tol = 1e-9;
	    final double beta = 1.0/3.0;
	    double duration = 10;
	    int nzs = 7; //obtained by inspection.
	    double T = 1;
	    
				
		TruncatedSincPulse.TruncatedRootRaisedCosine rrcpulse=new TruncatedSincPulse.TruncatedRootRaisedCosine(T, beta, duration, true);
	    TruncatedSincPulse.TruncatedRootRaisedCosine rrcpulsetester=new TruncatedSincPulse.TruncatedRootRaisedCosine(T, beta, nzs);
	    boolean pass = true;
	    for(double t = -duration; t < duration; t+=0.01) 
	        pass &= (rrcpulse.pulse(t).subtract(rrcpulsetester.pulse(t)).abs()) < tol;
	    if (pass) System.out.println( "PASS" );
	    else System.out.println( "FAIL" );
	}

	public static void main(String[] args) {
		
		FinitePulseTest x=new FinitePulseTest();

	    x.finiteSinc4zeros();
	    x.finiteSinc2zeros();
	    x.normalisedTruncatedSinc();

	    x.testRootRaisedCosineNearZero();
	    x.testRootRaisedCosineNearBeta4();
	    x.testRootRaisedCosine();
	    x.testRootRaisedCosineFindZeros();
	    try {
			x.durationConstructedRRC();
		} catch (Exception e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}

	  
	}


	
	
	
	
	
	
	
}

package synchro2;

public class DirectTest 
{

	/* 
	 * File:   DirectTest.cpp
	 * Author: Robby McKilliam
	 *
	 * Created on 25 January 2013, 7:31 PM
	 */


	static double tol = 1e-3;
	static double T = 1.0; //symbol period
	static  int p = 3;
	static  int q = 21;
	static  int c = 11; //new c, its will divide T
	static double Ts = p * T / q; //sample period
	static double Delta = T/c;
	static double taumin = 10.0;
	static double taumax = 30.0;
	static double tau0 = 20.111;
	static  int L = 100;
	static  int N = (int) Math.ceil((T * L + taumax) / Ts);
	static int numzeros = 2;
	static int[] P=new int[0];
	static int[] D=new int[0];
	static Complex[] pilots=new Complex[0];
	static Complex[] s=new Complex[0];
	static Complex[] r=new Complex[0];
	static  int M = 4; //QPSK
	static int numpilots = 30;
	static TruncatedSincPulse tx=new TruncatedSincPulse(T, numzeros);
	static Complex a0 = Complex.polar(1,0.1*Math.PI); //phase and amplitude

	void testDirect() 
	{

		//FastAlgorithmNoStorev.Direct<TruncatedSincPulse> dir=null;
		try 
		{
			FastAlgorithmNoStorev.Direct<TruncatedSincPulse> dir = new FastAlgorithmNoStorev.Direct<TruncatedSincPulse>(P, D, pilots, tx, T, Ts, taumin, taumax, c, p, q);
			Naive<TruncatedSincPulse> nai = new Naive<TruncatedSincPulse>(P, D, pilots, tx, T, Ts, taumin, taumax, c);
			double tauhatdir = dir.estimate(r);//wrong from direct
		    double tauhatnai = nai.estimate(r);  //correct
		    double[] zgrid = dir.getZgrid();

		  
		    
		    System.out.println( "Direct Zgrid test ... ");
		    boolean pass = true;


		    for(int i = 0; i < zgrid.length; i++ ) {
		        
		        pass &= Math.abs(zgrid[i] - nai.Z(taumin + i*Delta)) < tol;
		        pass &= Math.abs(zgrid[i] - dir.Z(taumin + i*Delta)) < tol;
		    }
		    if (pass) System.out.println( "PASS" );
		    else System.out.println( "FAIL" );
		    

		    System.out.println( "Direct Ygrid test ... ");
		    pass = true;
		    double[] ygrid = dir.getYgrid();

		    
		    for(int i = 0; i < ygrid.length; i++ ) {
		        pass &= Math.abs(ygrid[i] - nai.Y(taumin + i*Delta).abs()) < tol;
		        pass &= Math.abs(ygrid[i] - dir.Y(taumin + i*Delta).abs()) < tol;
		    }
		    if (pass) System.out.println( "PASS" );
		    else System.out.println( "FAIL" );
		    
		    System.out.println( "Direct SSgrid test ... ");
		    pass = true;
		    double[] ssgrid = dir.getSSgrid();

		    

		    for(int i = 0; i < ssgrid.length; i++ ) {
		        pass &= Math.abs(ssgrid[i] - Math.abs(nai.SS(taumin + i*Delta))) < tol;
		        pass &= Math.abs(ssgrid[i] - Math.abs(dir.SS(taumin + i*Delta))) < tol;
		    }
		    if (pass) System.out.println( "PASS" );
		    else System.out.println( "FAIL" );
		    
		    
		    
		    System.out.println( "Direct coarse estimate... ");
		    pass = Math.abs(dir.coarseMaximiseSS() - nai.coarseMaximiseSS()) < tol;
		    if (pass) System.out.println( "PASS" );
		    else 
	    	{
	    		System.out.println( "FAIL" );
	    	}
		    

		    System.out.println( "Direct fine estimate... ");
		    pass = Math.abs(tauhatdir - tauhatnai) < tol;
		    if (pass) System.out.println( "PASS" );
		    else System.out.println( "FAIL" );
		    

	    
		} 
		catch (Exception e) 
		{
			e.printStackTrace();
		}
	     
	}
	
	

		void testTabulatedDirect() {
			
	   try {
		   
		   System.out.println( "TabulatedDirect Zgrid test ... ");
		   
			TabulatedFastAlgorithm.TabulatedDirect<TruncatedSincPulse> dir=new TabulatedFastAlgorithm.TabulatedDirect<TruncatedSincPulse>(P, D, pilots, tx, T, Ts, taumin, taumax, c, p, q);
			Naive<TruncatedSincPulse> nai = new Naive<TruncatedSincPulse>(P, D, pilots, tx, T, Ts, taumin, taumax, c);
		/*
	    double tauhatdir = dir.estimate(r);
	    double tauhatnai = nai.estimate(r);
	    
	    System.out.println( "TabulatedDirect Zgrid test ... ");
	    boolean pass = true;
	    double[] zgrid = dir.getZgrid();
	    for(int i = 0; i < zgrid.length; i++ ) {
	        //System.out.println( zgrid[i] << ", " << nai.Z(taumin + i*Delta) );
	        pass &= Math.abs(zgrid[i] - nai.Z(taumin + i*Delta)) < tol;
	        pass &= Math.abs(zgrid[i] - dir.Z(taumin + i*Delta)) < tol;
	    }
	    if (pass) System.out.println( "PASS" );
	    else System.out.println( "FAIL" );
	    
	    System.out.println( "TabulatedDirect Ygrid test ... ");
	    pass = true;
	    double[] ygrid = dir.getYgrid();
	    for(int i = 0; i < ygrid.length; i++ ) {
	        pass &= Math.abs(ygrid[i] - nai.Y(taumin + i*Delta).abs()) < tol;
	        pass &= Math.abs(ygrid[i] - dir.Y(taumin + i*Delta).abs()) < tol;
	    }
	    if (pass) System.out.println( "PASS" );
	    else System.out.println( "FAIL" );
	    
	    System.out.println( "TabulatedDirect SSgrid test ... ");
	    pass = true;
	    double[] ssgrid = dir.getSSgrid();
	    for(int i = 0; i < ssgrid.length; i++ ) {
	        //System.out.println( taumin + i*Delta << ", " << ssgrid[i] << ", " << Math.abs(nai.SS(taumin + i*Delta)) );
	        pass &= Math.abs(ssgrid[i] - Math.abs(nai.SS(taumin + i*Delta))) < tol;
	        pass &= Math.abs(ssgrid[i] - Math.abs(dir.SS(taumin + i*Delta))) < tol;
	    }
	    if (pass) System.out.println( "PASS" );
	    else System.out.println( "FAIL" );
	    
	    System.out.println( "TabulatedDirect coarse estimate... ");
	    //System.out.println( dir.coarseMaximiseSS() << ", " << nai.coarseMaximiseSS() );
	    pass = Math.abs(dir.coarseMaximiseSS() - nai.coarseMaximiseSS()) < tol;
	    if (pass) System.out.println( "PASS" );
	    else System.out.println( "FAIL" );
	    
	    System.out.println( "TabulatedDirect fine estimate... ");
	    pass = Math.abs(tauhatdir - tauhatnai) < tol;
	    if (pass) System.out.println( "PASS" );
	    else System.out.println( "FAIL" );
*/
	    } catch (Exception e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		} 
	    
	}
	   
		 
	   
	   

	public static void main(String[] args) {
	    

		for (int m = 0; m < numpilots; m++) P=Utili.Append(P,m);
	    for (int m = numpilots; m < L; m++) D=Utili.Append(D,m);
	    for (int m = 1; m <= L; m++) s=Utili.Append(s,Complex.polar(1.0, (2*Math.PI*(Math.random()%M))/M));//s.push_back(polar<double>(1.0, (2*pi*(rand()%M))/M)); 
	    for ( int i = 0; i < P.length; i++ ) pilots=Utili.Append(pilots,s[P[i]]);

	    //function generates transmitted signal
	    ComplexVariateFunction f =new ComplexVariateFunction()
	    {
			@Override
			public Complex value(double x) 
			{						
				int mini = Math.max((int) 1, (int) Math.ceil((x - tx.tmax()) / T));
		        int maxi = Math.min((int) L, (int) Math.floor((x - tx.tmin()) / T));
		        Complex sum=new Complex(0, 0);

		        for (int i = mini; i <= maxi; i++) 
	        	{
	        		sum =sum.add(s[i - 1].multiply(tx.pulse(x - i * T))); //this here is already fairly time consuming!
	        	}
		        return sum;
			}
	    };
	    for(int n = 1; n <= N; n++) r=Utili.Append(r,a0.multiply(f.value(n*Ts-tau0))); //generate noiseless received signal

	    DirectTest test=new DirectTest();
	    
	    try {
	    	test.testDirect();
		} catch (Exception e) {
			e.printStackTrace();
		}
	    
	    
	    
	    
	    
	    for (int m = 0; m < numpilots; m++) P=Utili.Append(P,m);
	    for (int m = numpilots; m < L; m++) D=Utili.Append(D,m);
	    for (int m = 1; m <= L; m++) s=Utili.Append(s,Complex.polar(1.0, (2*Math.PI*(Math.random()%M))/M));//s.push_back(polar<double>(1.0, (2*pi*(rand()%M))/M)); 
	    for ( int i = 0; i < P.length; i++ ) pilots=Utili.Append(pilots,s[P[i]]);
	    for(int n = 1; n <= N; n++) r=Utili.Append(r,a0.multiply(f.value(n*Ts-tau0))); //generate noiseless received signal
	    
	    try {
			test.testTabulatedDirect();
		} catch (Exception e) {
			e.printStackTrace();
		}
	    
	    
	    
	    
	    
	    
	    
	    
	
	}
}

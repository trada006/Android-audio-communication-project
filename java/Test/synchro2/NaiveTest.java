package synchro2;

public class NaiveTest<Pulse extends FinitePulse> 
{

	/* 
	 * File:   TimeOffsetEstimatorTest.cpp
	 * Author: Robby McKilliam
	 *
	 * Created on 24/01/2013, 1:32:09 PM
	 */


	static double tol = 1e-6;
	static double T = 1.0; //symbol period
	static  int p = 1;
	static  int q = 5;
	static  int c = 5; //new c, its will divide T
	static double Ts = p * T / q; //sample period
	static double taumin = 10.0;
	static double taumax = 30.0;
	static double step = 5.0;
	static  int L = 10;
	static  int N = (int) Math.ceil((T * L + taumax) / Ts);
	static int numzeros = 2;
	static int[] P;
	static int[] D;
	static Complex[] pilots;
	static Complex[] r;
	static  int M = 4; //QPSK
	static int numpilots = 4;

	void NaiveDotr(Naive<Pulse> est) {
	    System.out.println( "NaiveDotr ... ");
	    boolean pass = true;
	    pass &= est.m(taumin).subtract(new Complex(1.3499373963348864, 4.563483625803526).times(Ts)).abs() < tol;
	    pass &= est.m(taumin+step).subtract(new Complex(1.6496241036561894, -4.463905699456121).times(Ts)).abs() < tol;
	    pass &= est.m(taumin+2*step).subtract(new Complex(-3.9931090337278663, 2.5889754772421516).times(Ts)).abs() < tol;
	    pass &= est.m(taumin+3*step).subtract(new Complex(4.748483513451538, 0.31562335065586683).times(Ts)).abs() < tol;
	    pass &= est.m(taumin+4*step).subtract(new Complex(-3.6153254669352815, -3.0946947418331097).times(Ts)).abs() < tol;
	    if (pass) System.out.println( "PASS" );
	    else System.out.println( "FAIL" );
	}

	void NaiveY(Naive<Pulse> est) {
	    System.out.println( "NaiveY ... ");
	    boolean pass = true;
	    pass &= est.Y(taumin).subtract(new Complex(14.657852409140666, -3.228790691768866).times(Ts)).abs() < tol;
	    pass &= est.Y(taumin+step).subtract(new Complex(-13.67538616337809, -6.185591310609769).times(Ts)).abs() < tol;
	    pass &= est.Y(taumin+2*step).subtract(new Complex(7.254044220717808, 13.139884665524072).times(Ts)).abs() < tol;
	    pass &= est.Y(taumin+3*step).subtract(new Complex(2.0523237347316847, -14.868278107005565).times(Ts)).abs() < tol;
	    if (pass) System.out.println( "PASS" );
	    else {
	        System.out.println( "FAIL" );
	        System.out.println( est.Y(taumin) );
	        System.out.println( est.Y(taumin+step) );
	        System.out.println( est.Y(taumin+2*step) );
	    }
	}

	void NaiveZ(Naive<Pulse> est) {
	    System.out.println( "NaiveZ ... ");
	    boolean pass = true;
	    pass &= Math.abs(est.Z(taumin) - 28.553768507361884*Ts) < tol;
	    pass &= Math.abs(est.Z(taumin+step) - 28.553768507361887*Ts) < tol;
	    pass &= Math.abs(est.Z(taumin+2*step) - 28.55376850736188*Ts) < tol;
	    pass &= Math.abs(est.Z(taumin+3*step) - 28.553768507361884*Ts) < tol;
	    pass &= Math.abs(est.Z(taumin+4*step) - 27.183844897660062*Ts) < tol;
	    if (pass) System.out.println( "PASS" );
	    else {
	        System.out.println( "FAIL" );
	        System.out.println( est.Z(taumin) );
	        System.out.println( est.Z(taumin+4*step) );
	    }
	}


	void NaiveCoarseEst(Naive<Pulse> est) {
	    System.out.println( "Naive coarse estimate ... ");
	    boolean pass = Math.abs(29.199999999999932 - est.coarseMaximiseSS()) < tol;
	    if (pass) System.out.println( "PASS" );
	    else {
	        System.out.println( "FAIL" );
	        System.out.println( est.coarseMaximiseSS() );
	    }
	}

	void NaiveFineEst(double tauhat) {
	    System.out.println( "Naive fine estimate ... ");
	    boolean pass = Math.abs(29.199999999999932 - tauhat) < tol;
	    if (pass) System.out.println( "PASS" );
	    else {
	        System.out.println( "FAIL" );
	        System.out.println(  tauhat );
	    }
	}

	public static void main(String[] args) {

	    for (int m = 0; m < numpilots; m++) P=Utili.Append(P,  m);
	    for (int m = numpilots; m < L; m++) D=Utili.Append(D,  m);
	    
	    for (int m = 1; m <= numpilots; m++) pilots=Utili.Append(pilots,Complex.polar(1.0, 0.1 * m));
	    
	    
	    
	    for (int n = 1; n <= N; n++) r=Utili.Append(r,Complex.polar(1.0, -0.1 * n));
	    
	    Naive<TruncatedSincPulse> est = null;
		try {
			 est = new Naive<TruncatedSincPulse>(P, D, pilots, new TruncatedSincPulse(T, numzeros), T, Ts, taumin, taumax, c);
			//est.testNaive();
			
		} catch (Exception e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
	    double tauhat = est.estimate(r);
	    
	    NaiveTest<TruncatedSincPulse> x=new NaiveTest<TruncatedSincPulse>();
	    
	    x.NaiveDotr(est);
	    x.NaiveY(est);
	    x.NaiveZ(est);
	    x.NaiveCoarseEst(est);
	    x.NaiveFineEst(tauhat);
	    
	}


	
	
	
	
	
	
}

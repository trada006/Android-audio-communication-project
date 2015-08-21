import org.junit.After;
import org.junit.AfterClass;
import org.junit.Before;
import org.junit.BeforeClass;
import org.junit.Test;



import static org.junit.Assert.*;
public class NaiveTest 
{
	public NaiveTest() {
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
	//static  int N = (int) ceil((T * L + taumax) / Ts);
	static int numzeros = 2;
	static int[] P;
	static int[] D;
	static Complex[] pilots;
	static Complex[] r;
	static  int M = 4; //QPSK
	static int numpilots = 4;
	
	
	
	/*
	
	public void NaiveDotr(Naive<Pulse> est) {
	    assertTrue(abs(est.m(taumin) - new Complex(1.3499373963348864, 4.563483625803526).times(Ts)) < tol);
	    assertTrue(abs(est.m(taumin+step) - new Complex(1.6496241036561894, -4.463905699456121).times(Ts)) < tol);
	    assertTrue(abs(est.m(taumin+2*step) - new Complex(-3.9931090337278663, 2.5889754772421516).times(Ts)) < tol);
	    assertTrue(abs(est.m(taumin+3*step) - new Complex(4.748483513451538, 0.31562335065586683).times(Ts)) < tol);
	    assertTrue(abs(est.m(taumin+4*step) - new Complex(-3.6153254669352815, -3.0946947418331097).times(Ts)) < tol);
	}
	
	
	
	*/
	
	
	
	
	
	
	
	
	
	void NaiveFineEst(double tauhat) 
	{
	    assertTrue(Math.abs(29.199999999999932 - tauhat) < tol);
	}
	
	
	
	
	
}

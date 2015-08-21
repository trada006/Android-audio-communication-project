import java.util.Arrays;

import org.junit.After;
import org.junit.AfterClass;
import org.junit.Before;
import org.junit.BeforeClass;
import org.junit.Test;





import static org.junit.Assert.*;

public class UtiliTest extends Utili 
{
	public UtiliTest() {
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
	
	
	public double sqrint(double a, double b)
	{
		return b*b*b/3 - a*a*a/3;
	}
	
	
	/*@Test
	public void trapezoidalSquareTest()
	{
		double tol = 1e-6;
		int N = 5000;
		SingleVariateFunction sqr = new SingleVariateFunction() 
		{
	         public double value(double x) 
	         {
	             return x;
	         }
	    };
	    System.out.println(trapezoidal(sqr, -3.0, 4.0, N) + "     " +sqrint(-3.0,4.0));
	    
	    
		assertTrue( Math.abs(trapezoidal(sqr, -3.0, 4.0, N) - sqrint(-3.0,4.0)) < tol);
		
		
	}*/
	
	@Test
	public void gcdTest()
	{
		assertTrue(gcd(10,2) == 2);
	}
	
	//int gcd, x, y;
	@Test
	public void extended_gcd_Test() 
	{
		//array elements
	    //1st element = gcd
		//2nd element = x
		//3rd element = y
		int gcd = 0, x = 0, y = 0;
	    int[] array = extended_gcd(10,2, gcd, x, y);
	    assertTrue (array[0]==2 && (10*array[1] + 2*array[2]) == array[0]);
	    array = extended_gcd(3458,4864,gcd,x, y);
	    assertTrue (array[0]==38 && (3458*array[1] + 4864*array[2]) == array[0]);
	}
	
	
	/*@Test
	public void msequenceTest() 
	{ 
		int[] s = msequence(4);
	    int[] exp4 = {1,1,1,0,0,0,1,0,0,1,1,0,1,0,1};
	    assertTrue(Arrays.equals(exp4, s));
	    s = msequence(5);
	    int[] exp5 = {1,1,1,1,0,0,0,1,1,0,1,1,1,0,1,0,1,0,0,0,0,1,0,0,1,0,1,1,0,0,1}; 	
	    assertTrue(Arrays.equals(exp5, s));
	}*/
	
	
	@Test
	public void fzeroLinearTest() 
	{
	    double tol = 1e-7;
		SingleVariateFunction f = new SingleVariateFunction() 
		{
	         public double value(double x) 
	         {
	             return x;
	         }
	    };
	    double x=new Bisection(f, -11.0,8.0, tol).zero();	
	    assertTrue(Math.abs(0.0 - f.value(x)) <  tol);  
	    assertTrue(Math.abs(0.0 - x) <  tol);
	  }
	
	
	
	@Test
	public void fzeroCubicTest() 
	{
	    double tol = 1e-7;
	    
	    SingleVariateFunction f = new SingleVariateFunction() 
		{
	         public double value(double x) 
	         {
	             return x*x*(x-1);
	         }
	    };
	 
	    double x=new Bisection(f, 0.5,1.7,tol).zero();	  
	    assertTrue(Math.abs(1.0 - x) <  tol);
	  }
	
	
	

	@Test
	public void mpilotsequenceTest() 
	{
	    double tol = 1e-8;
	    Complex[] pilots = pilotmsequence(3);
	    Complex[] exp3=new Complex[7];
	    
	    exp3[0]=new Complex(1,1);
	    exp3[1]=new Complex(1,1);
	    exp3[2]=new Complex(-1,-1);
	    exp3[3]=new Complex(-1,-1);
	    exp3[4]=new Complex(1,1);
	    exp3[5]=new Complex(-1,-1);
	    exp3[6]=new Complex(1,1);
	   
	    boolean pass=false;
	    for(int i = 0; i < pilots.length; i++) 
    	{
	    	if ((exp3[i].subtract(pilots[i])).abs()   < tol)
	    	{
	    		pass = true;
	    		continue;
	    	}
	    	else 
    		{
    			pass = false;
    		}
    	}
	    assertTrue(pass==true);
	}
	

	@Test
	public void innerproductTest() 
	{
	    Complex[] x = {new Complex(1),new Complex(2),new Complex(3),new Complex(4)};
	    Complex[] y = {new Complex(1),new Complex(0,1),new Complex(1),new Complex(0,1)};
	    
	    
	    
	    
	    Complex expected = new Complex(4,-6);
	    //System.out.println((expected.subtract(inner_product(x,y))).abs());
	    
	    assertTrue((expected.subtract(inner_product(x,y))).abs() < 1e-8);
	}
	
	
	
}

package synchro2;





public class UtiliTest extends Utili 
{
	
	public UtiliTest()
	{}
	
	
	SingleVariateFunction sqr=new SingleVariateFunction()
	{
		@Override
		public double value(double x) {
			return x*x;
		}
	};
	

	double sqrint(double a, double b) { return b*b*b/3 - a*a*a/3; }
	
	void trapezoidalSquareTest() {
		

	    double tol = 1e-6;
	    System.out.println( "Trapezoidal Square Test ... ");
	    int N = 5000;
	    boolean pass = Math.abs(trapezoidal(sqr, -3.0, 4.0, N) - sqrint(-3.0,4.0)) < tol;   
	    if(pass) System.out.println( "PASS" );
	    else {
	        System.out.println( "FAIL" );
	        System.out.println( "expected " + sqrint(-3.0,4.0) + ", but was " + trapezoidal(sqr, -3.0, 4.0, N) );
	    }
	}

	void gcdTest() {
	    System.out.println( "gcd ... ");
	    boolean pass = gcd(10,2)==2;
	    pass &= gcd(3458,4864)==38;
	    if(pass) System.out.println( "PASS" );
	    else System.out.println( "FAIL" );
	}

	void extended_gcd_Test() {
	    System.out.println( "extended gcd ... ");
	    int gcd = 0, x = 0, y = 0;
	    extended_gcd(10,2,gcd,x,y);
	    int numCorrect=0;
	    if ((gcd_==2)&&((10*x_ + 2*y_) == gcd_))
	    {
	    	numCorrect++;
	    }
	    extended_gcd(3458,4864,gcd,x,y); 
	    if ((gcd_==38)&&((3458*x_ + 4864*y_) == gcd_))
	    {
	    	numCorrect++;
	    }

	    if(numCorrect==2) System.out.println( "PASS" );
	    else 
    	{
    		System.out.println( "FAIL" );
    	}
	}

	void fzeroLinearTest() {
	    System.out.println( "Bisection linear ... ");
	    double tol = 1e-7;
	    SingleVariateFunction f = new SingleVariateFunction() 
		{
	         public double value(double x) 
	         {
	             return x;
	         }
	    };
	    double x = new Utili.Bisection(f, -11.0,8.0,tol).zero();
	    boolean pass = Math.abs(0.0 - f.value(x)) <  tol;
	    pass &= Math.abs(0.0 - x) <  tol;
	    if(pass) System.out.println( "PASS" );
	    else System.out.println( "FAIL" );
	  }

	void fzeroCubicTest() {
	    System.out.println( "Bisection cubic ... ");
	    double tol = 1e-7;
	    SingleVariateFunction f = new SingleVariateFunction() 
		{
	         public double value(double x) 
	         {
	             return x*x*(x-1);
	         }
	    };
	    double x = new Utili.Bisection(f, 0.5,1.7,tol).zero();
	    boolean pass = true;
	    pass &= Math.abs(1.0 - x) <  tol;
	    if(pass) System.out.println( "PASS" );
	    else System.out.println( "FAIL" );
	  }

	void msequenceTest() {
	    System.out.println( "m-sequence ... ");
	    boolean pass = true;
	    int[] s = msequence(3);
	    int exp3[] = {1,1,0,0,1,0,1};
	    for(int i = 0; i < s.length; i++) pass &= exp3[i] == s[i];
	    s = msequence(4);
	    int exp4[] = {1,1,1,0,0,0,1,0,0,1,1,0,1,0,1};
	    for(int i = 0; i < s.length; i++) pass &= exp4[i] == s[i];
	    s = msequence(5);
	    int exp5[] = {1,1,1,1,0,0,0,1,1,0,1,1,1,0,1,0,1,0,0,0,0,1,0,0,1,0,1,1,0,0,1};
	    for(int i = 0; i < s.length; i++) pass &= exp5[i] == s[i];
	    if(pass) System.out.println( "PASS" );
	    else System.out.println( "FAIL" );
	}

	void mpilotsequenceTest() {
	    System.out.println( "pilot m-sequence ... ");
	    double tol = 1e-8;
	    boolean pass = true;
	    Complex[] pilots = pilotmsequence(3);
	    Complex exp3[] = {new Complex(1),new Complex(1),new Complex(-1),new Complex(-1),new Complex(1),new Complex(-1),new Complex(1)};
	    for(int i = 0; i < pilots.length; i++) pass &= (exp3[i].subtract(pilots[i]).abs()) < tol;
	    if(pass) System.out.println( "PASS" );
	    else System.out.println( "FAIL" );
	}

	void innerproductTest() {
	    System.out.println( "inner product ... ");
	    Complex[] x = {new Complex(1),new Complex(2),new Complex(3),new Complex(4)};
	    Complex[] y = {new Complex(1), new Complex(0,1),new Complex(1) ,new Complex(0,1)};
	    Complex expected = new Complex(4,-6);
	    boolean pass=false;
		try {
			pass = expected.subtract(inner_product(x,y)).abs() < 1e-8;
		} catch (Exception e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
	    if(pass) System.out.println( "PASS" );
	    else System.out.println( "FAIL" );
	}

	public static void main(String[] args) {

		
		UtiliTest x=new UtiliTest();
		
		x.trapezoidalSquareTest();
	    x.extended_gcd_Test();
	    x.gcdTest();
	    x.fzeroLinearTest();
	    x.fzeroCubicTest();
	    x.msequenceTest();
	    x.mpilotsequenceTest();
	    x.innerproductTest(); 
	    
	    
	}


}

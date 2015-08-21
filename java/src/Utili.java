


public class Utili
{
	
	static final double pi = 3.141592653589793238463;
	
	public static double sqr(double x)
	{
		return x*x;
	}
	
	
	/** Cube of x */
	public static double cub(double x) 
	{
	    return x*x*x;
	}
	
	
	
	/** The sign of x, zero if x is zero. */
	public static double signum(double x)
    {
        if (x > 0)
        {
            return 1;
        }
        else if (x < 0)
        {
            return -1;
        }
        else return 0;
    }
	
	
	
	/** sinc function */
	public static double sinc(double t) 
	{
	    if (Math.abs(t) < 5e-3) 
    	{
    		return 1.0 - t * t * (1.0 / 6 - 1.0 / 120 * t * t);
    	}
	    else return Math.sin(pi * t) / (pi * t);
	}
	
	
	
	public static int gcd(int a, int b)
	{
		if (b ==0) return Math.abs(a);
		else return gcd(Math.abs(b), Math.abs(a) % Math.abs(b));
		
	}
	
	
	
	
	/** Trapezoidal integration of a function f mapping double to double */
	public static double trapezoidal(SingleVariateFunction f, double a, double b, int N) 
	{
	    double del = (b - a) / N;
	    double inner = 0.0;
	    for (int n = 1; n <= N - 1; n++) 
    	{
    		inner += 2 * f.value(a + n * del);
    	}
	    return del / 2 * (inner + f.value(a) + f.value(b));
	}
	
	
	
	public static int[] extended_gcd(int a, int b, int gcd, int x, int y) 
	{
	    x = 0;
	    y = 1;
	    int u = 1, v = 0, m, n, q, r;
	    gcd = b;
	    while (a != 0) 
	    {
	        q = gcd / a;
	        r = gcd % a;
	        m = x - u*q;
	        n = y - v*q;
	        gcd = a;
	        a = r;
	        x = u;
	        y = v;
	        u = m;
	        v = n;
	    }
	    
	    int[] array = new int[]{gcd, x, y};
	    return array;
	}
	
	
	
	public static int[] msequence(int N) throws IllegalArgumentException
	{
	if(N < 3 || N > 14) 
	   {
		throw new IllegalArgumentException("Input argument mismatch: N<3 || N>14");
	   }
	    
	    //initial taps
	    int[] taps = null;
	    if(N == 3) {
	        int[] array = {0,1,1};
	        taps=array;
	    }
	    
	    else if(N == 4)
	    {	
	        int[] array = {0,0,1,1};
	        taps = array; 
	    }
	    
	    
	    else if(N == 5){
	        int[] array = {0,0,1,0,1};
    		taps=array;        	
	    }

	    else if(N == 6){
	        int[] array = {0,0,0,0,1,1};
	        taps=array;
	    }
	    else if(N == 7){
	        int[] array = {0,0,0,1,0,0,1};
	        taps=array;
	    }
	    else if(N == 8){
	        int[] array = {0,0,0,1,1,1,0,1};
	        taps=array;
	    }
	    else if(N == 9){
	        int[] array = {0,0,0,0,1,0,0,0,1};
	        taps=array;
	    }
	    else if(N == 10){
	      int[] array = {0,0,0,0,0,0,1,0,0,1};
	      taps=array;
	    }
	    else if(N == 11){
	      int[] array = {0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 1};
	      taps=array;
	    }
	    else if(N == 12){
	      int[] array = {0, 0, 0, 0, 0, 1, 0, 1, 0, 0, 1, 1};
	      taps=array;
	    }
	    else if(N == 13){
	      int[] array = {0,0,0,0,0,0,0,0,1,1,0,1,1};
	      taps=array;
	    }
	    else if(N == 14){
	      int[] array = {0,0,0,1,0,0,0,1,0,0,0,0,1,1};
	      taps=array;
	    }
	    
	    int M = (1<<N) - 1; //M = 2^N-1
	    
	    
	    
	    int[] m = new int[N]; //vector of 1's of length N
	    for (int i=0; i<m.length; i++)
	    {
	    	m[i] = 1;
	    }
	    

	    
	    int[] regout = new int[M]; //output m sequence
	    for(int ind = 0; ind < M; ind++)
	    {
	        int buf = 0;
	        for(int i = 0; i < N; i++) 
        	{
        		buf += taps[i]*m[i];
        		
        	}
	        
	        buf = buf % 2; //xor
	        
	        //left shift
	        int i=0;
	        int a = m[0];
	        int NumOfShift = m.length - 1;
	        int j=0;
	        while (j < NumOfShift)
	        {
	        	i=0;
	        	a = m[0];
	        	while (i < m.length-1)
		        {
		        	m[i] = m[++i];
		        }
	        	m[i] = a;
	        	j++;
	        }
	        m[0] = buf;
	        regout[ind] = m[N-1];
	    }
    
	    return regout;
	    
	    
	    
	    
	}
	
	
	
	public static class Bisection 
	{
		public Bisection()
		{}
		
	    public Bisection(SingleVariateFunction f, double ax, double bx, double tol) 
	    {
	    	int ITRMAX = 100;
	        double a = ax;
	        double b = bx;
	        
	        
	        for (int i = 1; i <= ITRMAX; i++) 
	        {
	            double c = (a + b) / 2;
	            double fc = f.value(c);
	            if (fc == 0 || Math.abs(a - b) / 2 < tol) 
	            {
	                xzero = c;
	                return;
	            }
	            if (signum(fc) == signum(f.value(a)))
            	{
            		a = c;
            	}
	            else b = c;
	        }
	        System.out.println("Warning: Bisection failed. Reached the maximum number " + ITRMAX + " iterations");
	    }

	    /** The value of x where f(x) = 0 */
	    double zero() 
	    {
	        return xzero;
	    }
	    protected double xzero;
	}
	
	
	
	
	public double sign(double a, double b) 
	{
        return Math.abs(a) * signum(b);
    }
	
	
	
	
	
	//need concatenate
	
	/*
	static T[] concatenate(T[] A, T[] B) 
	{
	    T[] C;
	    C.reserve(A.size() + B.size());
	    for (T a : A) C.push_back(a);
	    for (T b : B) C.push_back(b);
	    return C;
	}
	
	*/
	
	
	
	
	/** The Hermitian inner product between two complex vectors. */

	public Complex inner_product(Complex[] x, Complex[] y)
	{
	    int N = x.length;
	    //if(N != y.length) throw "Vectors must be the same size for inner product";
	    Complex sum = new Complex();
	    for(int i = 0; i < N; i++)
	    {
	      sum=sum.add(x[i].multiply(y[i].conjugate()));
	    }
	    return sum;
	}
	
	
	
	
	
	public Complex[] pilotmsequence(int N)
	{
		int[] ms = msequence(N);
		Complex[] pilots = new Complex[ms.length];		
		for(int i = 0; i < ms.length; i++)
		{
	        if(ms[i]==0) 
        	{
        		pilots[i] = new Complex(-1,-1);
        	}
	        else if(ms[i]==1)
        	{
        		pilots[i] = new Complex(1,1);
        	}
	    }	
	    return pilots;
		
	}
	
	
	
}

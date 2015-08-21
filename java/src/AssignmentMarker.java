import org.junit.runner.JUnitCore;
import org.junit.runner.Result;
import org.junit.runner.notification.Failure;

public class AssignmentMarker { 
	
	
	//private static java.util.ArrayList<Failure> failures;
	
	
	
	

	
	// More complex test information
	private static void runATest(String name, Class c) {
		System.out.println("\n" + name);
		
		System.out.println();
	
	}
	
	public static void main(String[] args) {
	
	
		
		//Complex a=new Complex(1,2);
		//Complex b=new Complex(3,4);
		
		//Complex[] t = new Complex[]{a,b};
		
		//Complex[] t = new Complex[1];
	
		
		
		
		//t[0].add(a);
		
		//System.out.println(t[0].toString());
		
		runATest("BrentTest", BrentTest.class);
		
		runATest("UtilTest", UtiliTest.class);

		
		
		
}  
	
}

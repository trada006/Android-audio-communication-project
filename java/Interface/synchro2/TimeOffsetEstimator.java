package synchro2;

public interface TimeOffsetEstimator
{
	    /** Run the estimator. Return time offset estimate */
	   public double estimate(Complex[] r);
	    /** Return a string with the name of this estimator */
	   public String name();
}

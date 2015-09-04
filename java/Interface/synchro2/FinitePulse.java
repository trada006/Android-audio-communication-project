package synchro2;

public interface FinitePulse {

	 /** pulse is zero for t less than tmin */
    public double tmin();
    /** pulse is zero for t more than tmax */
    public double tmax();
    /** pulse period */
    public double T();
    /** The transmit pulse */
    public Complex pulse(double t); 

}

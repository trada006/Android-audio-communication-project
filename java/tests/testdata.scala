/**
* Generates some output to test the C implementation
*/

import numbers.finite.RectComplex
import numbers.finite.Complex
import numbers.finite.PolarComplex
import am.finitepulse.Recursive
import am.finitepulse.Naive
import am.finitepulse.PilotOnly
import am.OerderMeyerAndMassey
import am.TruncatedRootRaisedCosine
import am.GeneralFiniteDurationAmplitudeModulatedSignalFromFinitePulse
import am.FinitePulseWithLookupTable
import am.NormalisedFinitePulse
import am.PulseWithTappedFilter
import am.PulseWithFilter
import am.FinitePulse
import pubsim.distributions.GaussianNoise

val T = 1.0; //symbol period
  val p = 1
  val q = 5
  val c = 5 //new c, its will divide T
  val Ts = p*T/q; //sample period
  val taumin = 10.0 
  val taumax = 30.0
  val L = 10
  val N = scala.math.ceil((T*L+taumax)/Ts).toInt
  
  //transmitted symbols
  val M = 4 //QPSK
  val s = (1 to L).map(m => new PolarComplex(1, 0.1*m))
  val numpilots = 4
  val P = 0 until numpilots
  val D = numpilots until L

  //create the truncated sinc pulse we will use
  val numzeros = 2
  val tsincpulse = new NormalisedFinitePulse(new am.TruncatedSincPulse(T,numzeros))

//recieved signal 
  val r = (1 to N).map(n => new PolarComplex(1, -0.1*n))

  val est = new Naive(P,D,s,tsincpulse.pulse,tsincpulse.tmin,tsincpulse.tmax,T,Ts,taumin,taumax,c)
val tauhat =  est.estimate(r)

println("dotr")
for( t <- taumin to taumax by 5.0 ) println(est.dotr(t))
println("Y")
for( t <- taumin to taumax by 5.0 ) println(est.Y(t))
println("Z")
for( t <- taumin to taumax by 5.0 ) println(est.Z(t))
println("coarse est")
println(est.coarseMaximiseSS)
println("fine est")
println(tauhat)

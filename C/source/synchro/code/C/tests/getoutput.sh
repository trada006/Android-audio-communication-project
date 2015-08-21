#export JAVA_OPTS="-Xprof -d64 -server -Xms1g -Xmx1g"
export JAVA_OPTS="-d64 -server -Xms1g -Xmx1g"
scala -nocompdaemon -cp ../../lib/TimeOffsetEstimator.jar:../../lib/ScalaNumbers.jar:../../lib/bignums.jar:../../lib/PubSim.jar:../../lib/Jama-1.0.2.jar:../../lib/flanagan.jar:lib/colt.jar:../../lib/RngPack.jar:../../lib/jtransforms-2.4.jar testdata.scala

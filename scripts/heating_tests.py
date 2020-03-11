import heating_tests

gasT=100.0
gasDensity=100.0
hAbund=0.5
heAbund=0.1
electronAbund=1.0e-4
hxAbund=0.0
hexAbund=0.0
out=heating_tests.getatomiccooling(gasT,gasDensity,hAbund,heAbund,electronAbund,hxAbund,hexAbund)
print(out)
import heating_test


#gasT,gasDensity,hAbund,heAbund,electronAbund,hxAbund,hexAbund

for dens in [1.0e2,1.0e3,1.0e4]:
	for temp in [10,100,1000,10000,1.0e6]:
		a=heating_test.atomiccooling(temp,dens,0.5,0.1,1.0e-4,0.01,0.05)
		print(dens,temp,a)
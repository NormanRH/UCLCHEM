PROGRAM array_test
	IMPLICIT NONE
	DOUBLE PRECISION :: a(3),b(3),c(3)

	a=(/1,2,3/)
	b=2.0
	c=(/3.0,4.0,5.0/)
	a=a*b*c
	write(*,*) a
	write(*,*) sum(a)
END PROGRAM
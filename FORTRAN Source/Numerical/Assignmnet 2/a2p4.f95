Program Summation
IMPLICIT NONE
INTEGER :: I

DOUBLE PRECISION :: x,func,sumval
print*,"Enter the x value."
READ(*,*) x
print*,"Enter the I value:"
READ(*,*) I
sumval = -1.0
!Begin Do for the 
	do I = 1,I
	sumval = sumval + FUNC(x,I)
	end do
print*,"Sum:",sumval
End Program Summation


DOUBLE PRECISION FUNCTION FUNC(G,Q)
DOUBLE PRECISION :: G, Y
INTEGER :: Q

Y=( cos(G *dble(Q))/dble(Q) )

RETURN
END

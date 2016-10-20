PROGRAM Practice
IMPLICIT NONE 

DOUBLE PRECISION :: X=1.0,Y
INTEGER :: I
CHARACTER (len=60) :: scrpt_key = "set key inside left top vertical Right noreverse enhanced autotitle box lt black linewidth 1.0000"
OPEN (unit=17,file="results.gnu",action="write",status="replace")
write(17,*) "set term postscript eps font" ,' "Times-Roman,12"' 
write(17,*)"set output", "'ftran1.eps'"
write(17,*)scrpt_key
write(17,*)"set samples 50, 50"
!write(17,*)"dashtype solid"
write(17,*)"set title",' "Simple Plots" '
write(17,*)"set title  font ",'",20" norotate'
!write(17,*)"plot [-10:10] sin(x),atan(x),cos(atan(x))"
write(17,*) "plot [0:100] '1.dat' with lines "
CLOSE (17)

OPEN (unit=18,file= "1.dat", action="write", status="replace")

DO I=1,20
X=2.0 * X
write(18,*) X, Y(X)
END DO
CLOSE(18)



END PROGRAM Practice


DOUBLE PRECISION FUNCTION Y(X)
DOUBLE PRECISION :: OUTPUT,X


OUTPUT = X**2

RETURN
END FUNCTION

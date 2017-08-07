PROGRAM Practice
IMPLICIT NONE 

DOUBLE PRECISION :: X=1.0,Y,Xa=-2.0,Xb=2.0,SUMM=0,step=0
INTEGER :: I,j
!Generate Script.gnu-easier to make changes and generate each time
OPEN (unit=17,file="results.gnu",action="write",status="replace")
write(17,*) "set term postscript eps font" ,' "Times-Roman,12"' 
write(17,*)"set output", "'numerical1.eps'"
!write(17,*)"set samples 50, 50"
!write(17,*)"dashtype solid"
write(17,*)"set title",' "More Terms" '
write(17,*)"set title  font ",'",20" norotate'
write(17,*) "plot [",Xa,":",Xb,"] '1.dat' with lines"
CLOSE (17)

OPEN (unit=18,file= "1.dat", action="write", status="replace")
X=Xa
!Step in X from Xa to Xb
Do I=0,1000
step = step + 0.1d0
SUMM=0
	!Step in j however many terms are needed
	Do j=1,4
	SUMM = SUMM + Y(X+step,j)
	write(18,*) (X+step),SUMM
	End do
End do
CLOSE(18)



END PROGRAM Practice


DOUBLE PRECISION FUNCTION Y(X,j)
DOUBLE PRECISION :: OUTPUT,X,M
Integer:: j
M=j
OUTPUT = (cos(M*3.14*X)/M**2)

RETURN
END FUNCTION

PROGRAM Practice
IMPLICIT NONE 

DOUBLE PRECISION :: X=1.0,Y,Xa=-1.0,Xb=1.0,Ya=-1,Yb=1,SUMM=0,step=0,F
INTEGER :: I,j,S
!Generate Script.gnu-easier to make changes and generate each time
OPEN (unit=17,file="3dresults.gnu",action="write",status="replace")
write(17,*) "set term postscript eps font" ,' "Times-Roman,12"' 
write(17,*)"set output", "'numerical1.eps'"
write(17,*)"set mapping cartesian"
write(17,*)"set grid"
write(17,*)"set xrange [-1:1]"
write(17,*)"set yrange [-1:1]"
write(17,*)"set title",' "Surface Plot" '
write(17,*)"set title  font ",'",20" norotate'
write(17,*) 'splot "4.dat" u 1:2:3 with lines'
CLOSE (17)

OPEN (unit=18,file= "4.dat", action="write", status="replace")
X=Xa
Y=Ya
!Step in X from Xa to Xb
Do S=0,10
step = step + 0.1d0
SUMM=0
	!Step in j however many terms are needed
	Do j=1,4
		Do I=1,4
	SUMM = SUMM + F(X+step,Y+step,j,I)
		End Do
	write(18,*) (X+step),(Y+step),SUMM
	End do
End do
CLOSE(18)



END PROGRAM Practice


DOUBLE PRECISION FUNCTION F(X,Y,j,i)
DOUBLE PRECISION :: OUTPUT,X,Y,M,K
Integer:: j,i
M=j
K=i
OUTPUT = (cos(K*3.14*X)*sin(M*3.14*Y)/((K**2)*(M**3)))
RETURN
END FUNCTION

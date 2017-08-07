Program NewtonSummation
Implicit None
Double Precision :: FUNC,F,J,T,P,B,N,TOL,SUMM
Integer :: I,A
Real :: PI = 4.*atan(1.0)
TOL = 1.0E-8
N= 0.25
Do 10 A = 1,50 ! Max interval will be i=0,i=50 of sum ((N_-1) +1/5)
N = N + 1.0/5.0 !Step in P Guess
P = N
Print*,"PVAL:",P
SUMM=0

		Do 20 I=10000000,1,-1
B=P
SUMM = SUMM + F(P,I) 
	P = B - (F(B,I)*J(B,I))/((J(B,I)**2)-(F(B,I)*T(B,I))) !Numerical Method Based on Newton's Method
!Check For Cos(x) zeros
	If(FUNC(P)==0) Then !Test For Asymtotic Behavior
	
	P = P/4
Print*,"Approaching Discontinuity."
	
	End IF
If(abs(-1.0+SUMM).LT.TOL) Then
Print*,"Sum",SUMM
Print*,"Function Value:",F(P,I)
Print*,"P:",P
Print*,"I:",I
!Stop;
End IF


		20 END DO
!Print*,"Sum:",SUMM
10 END DO 
End Program NewtonSummation

!Original Function
Double Precision Function FUNC(X)
Double Precision X,Y
Y = cos(X)
Return
End

!MOD ORI FUNCT
Double Precision Function F(X,I)
Double Precision X,Y
Integer :: I
Y = (cos(dble(I)*x)/dble(I))
Return
End
!First Derivative
Double Precision Function J(X,I)
Double Precision X,Y
Integer :: I
Y = -(dble(I)*sin(dble(I)*x))
Return
End

!Second Derivative
Double Precision Function T(X,I)
Double Precision X,Y
Integer :: I
Y = dble(I)*cos((PI/2.0)+dble(I)*x)
Return 
End

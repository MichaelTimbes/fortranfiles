Program Secantp2
Implicit None
Double Precision :: A, B, P, Pold, F, TOL
Integer :: I
!I =0 !Initialize Iterator
A =1!Set to initial P to test tolerance
P=-1 !Left interval 
Pold = 1 !Outer interval expressed as Pold
Tol = 1.0E-8
DO! I = 1,20 
	B=P !B holds current P to store in Pold later
	P=P - (F(P) *((P-Pold))/(F(P)-F(Pold))) !Redefine P as outlined in the secant method
I = I + 1
	PRINT*,"Iteration:",I,"P:",P,"Pold:",Pold
	Pold=B !New Pold
	
If(abs(Pold-P) .LT. TOL) Then !If current P and Pold are very close- convergence is assumed
Print*,"Convergence is assumed:"
Print*,"Solution: ",Pold
Print*,"Iterations Needed: ",I
Stop
END If
END DO
End Program Secantp2

Double Precision Function F(G)
Double Precision :: G, Y

Y = ((-32.17)/(2.0 * G**2)) * (sinh(G) - sin(G))-1.7
Return 
End



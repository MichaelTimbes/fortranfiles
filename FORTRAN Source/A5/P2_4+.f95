Program P3_4
implicit none
integer:: n,filenum=7,itmax,i
double precision, allocatable,dimension (:,:) :: A
double precision, allocatable, dimension (:) :: b,x,l,u,d
double precision, dimension(:):: px(5)
print*,"n= "
read*,n
allocate(a(n,n),l(n),u(n),d(n), b(n),x(n))




!call Problem3_4v1(n,A)
call Problem3_4v2(n,L,U,D,B)
call Gaussv2(n,L,U,D,B,X,px)


open( unit=filenum, file='PROBLEM4.txt', status='replace')

write(filenum,*)"RESULT:"
write(filenum,*) "X:"
do i=1,5
write(filenum,'(10f10.6)')px(i)
enddo
write(filenum,*) 
write(filenum,*)
write(filenum,*)"ACTUAL RESULT:"
write(filenum,*) "X:"
do i=1,n
write(filenum,'(10f10.6)')x(i)
write(filenum,*)"AT:",i
enddo
write(filenum,*) 
write(filenum,*)

close(filenum)

deallocate(a,b,x)
End program P3_4






!2-D Matrix
!Subroutine Problem3_4v1(n,A)
!IMPLICIT NONE
!integer:: n,i,j,k,itmax
!DOUBLE PRECISION :: a(n,n),b(n),x(n)
!do i=1,n
!	do j=1,n
!	A(i,j)=0
!	
!	enddo
!b(i)=1
!enddo

!do i=1,n
!a(i,i)=5

!enddo
!do i=1,n-1
!	
!	a(i,i+1)=2
!	a(i+1,i)=2
!	
!enddo
!
!return 
!end

!Broken apart LUD matricies that make up A
Subroutine Problem3_4v2(n,L,U,D,B)
IMPLICIT NONE
integer:: n,i,j,k,itmax
DOUBLE PRECISION :: l(n),u(n),d(n),b(n)
do i=1,n
d(i)=5
l(i)=2
u(i)=2
b(i)=1
enddo
return 
end
!Based off a variation of Gauss Elimination for tridiagonal matricies (Thomas's Method)
Subroutine Gaussv2(n,L,U,D,B,X,px)
implicit none
integer :: n,itmax,i,j,k
double precision :: l(n),u(n),d(n),b(n),x(n),summ,tol=.0000000001,cp(n),dp(n),px(5)
double precision :: v1=0,v2=0



!Find the c' values for back substitution
do i=1,n
	IF(i.eq.1)THEN
	cp(i) = u(i)/d(i)
!print*,"C'",cp(i)
ELSE
cp(i) = u(i)/(d(i)-(l(i)*cp(i-1)))
!print*,"C'",cp(i)

ENDIF
enddo
!Find the d' values for back substitution
do i=1,n
IF(i.eq.1)THEN
	dp(i) = b(i)/d(i)
!print*,"D'",dp(i)
ELSE
dp(i) = (b(i)-(l(i)*dp(i-1)))/(d(i)-(cp(i-1)*l(i)))
!print*,"D'",dp(i)
ENDIF
enddo

!Now back solve for X values
do i=n,1,-1
IF(i.eq.n)THEN
x(i)=dp(i)
px(1)=x(i)
!print*,"Init X:",x(i)
ELSE
x(i) = dp(i)-(cp(i)*x(i+1))
!print*,"X:",x(i)
if(i.eq.250000) THEN
px(2) = x(i)
else if(i.eq.500000) THEN
px(3) = x(i)
else if(i.eq.750000) THEN
px(4) = x(i)
else if(i.eq.100000) THEN
px(5) = x(i)
ENDIF
ENDIF

enddo
return 
end


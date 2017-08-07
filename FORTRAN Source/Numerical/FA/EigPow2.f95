Program EigPow_2
implicit none
integer:: n,filenum=7,itmax,i
double precision, allocatable, dimension (:) :: b,x,l,u,d
double precision, dimension(:):: px(10)
double precision:: mul
!print*,"n= "
!read*,n
n=1000
allocate(l(n-1),u(n-1),d(n), b(n),x(n))



call Problem3(n,L,U,D,B)
!write(*,*) "L:"
!do i=1,n-1
!write(*,'(10f10.6)')L(i)
!enddo
!write(*,*) 
!write(*,*) "U:"
!do i=1,n-1
!write(*,'(10f10.6)')U(i)
!enddo
!write(*,*) 
!write(*,*) "D:"
!do i=1,n
!write(*,'(10f10.6)')D(i)
!enddo
!write(*,*) 

call InvPower2(n,L,U,D,B,X,px,mul)

open( unit=filenum, file='EPow2Res.txt', status='replace')

write(filenum,*)"RESULT:"
write(filenum,*) "X:"
do i=1,5
write(filenum,'(10f10.6)')px(i)
enddo
write(filenum,*) 
write(filenum,*)"Eig Val:"
write(filenum,'(10f10.6)') mul
close(filenum)

deallocate(l,u,d,b,x)
End program EigPow_2

!Broken apart LUD matricies that make up A
Subroutine Problem3(n,L,U,D,B)
IMPLICIT NONE
integer:: n,i,j,k,itmax
DOUBLE PRECISION :: l(n),u(n),d(n),b(n)
do i=1,n
d(i)=2+(dble(i**2)/dble(n**4))
b(i)=1
enddo
do i=1,n-1,1
l(i)=-1
u(i)=-1
enddo
return 
end
!!!
!!!
!!!
Subroutine InvPower2(n,L,U,D,B,X,px,mul)
implicit none
integer :: n,itmax=100,i,j,k,maxr
double precision :: l(n-1),u(n-1),d(n),b(n),x(n),tmp(n),summ,tol=.0000000001,px(10)
double precision :: maxv,tmpa=0.0,tmpb=0.0,q2=0.0,c,mul,mulold
x(:)=1
b(:)=1
!Create Q
tmp(1)=d(1)+u(1)
do i=2,n
IF(i.lt.n) THEN
tmp(i)=l(i-1)+d(i)+u(i)
ELSE
tmp(i)=d(i)+l(i-1)
endif
enddo
do i=1,n
tmpa=tmpa+(tmp(i)*x(i))
tmpb=tmpb+(x(i)*x(i))
enddo
q2=tmpa/tmpb

!write(*,'(10f10.6)')q2
do i=1,n
d(i)=d(i)-q2
!write(*,'(10f10.6)')d(i)
enddo

do j=1,20
call Jacv2(n,l,u,d,b,x,px)
call Get_Max(x,n,maxr,maxv)
c=x(maxr)
!print*,c
do i=1,n 
b(i) = 1/c*x(i)
enddo
mulold=mul
mul=q2+(1/c)
!DIAG
!write(*,*)"RESULT:"
!write(*,*) "X:"
!do i=1,n
!write(*,'(10f10.6)')x(i)
!enddo
!write(*,*) 
!write(*,*)
!write(*,*) "C:"
!write(*,'(10f10.6)')c
!write(*,*) 
!write(*,*)
IF(abs(mul-mulold)<.000001) THEN
write(*,'(10f10.6)')mul
if(n.gt.10) THEN
px=x(1:10)
endif
goto 2
endif
enddo


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
2 return 
end
!!!
!!!
!!!
!!!
!Modified Jacobi Method for tri diagonal matricies

Subroutine Jacv2(n,L,U,D,B,X,px)
implicit none
integer :: n,itmax=10,i,j,k
double precision :: l(n-1),u(n-1),d(n),b(n),x(n),xo(n),summ,tol=.00001,px(5)
double precision :: v1=0,v2=0,a
!Initial Guess that all are zero
do i=1,n
x(i)=0

enddo

do k=1,itmax

x(1) = (b(1) - (u(1)*x(2)))/d(1)
!print*,"First X:",x(1)
do i=2,n
summ=0

summ = (l(i-1)*x(i-1)) 
if(i.lt.n) THEN
summ = summ+(u(i)*x(i+1))
endif

!print*,"sum ",summ,i
!enddo

x(i)=((b(i)-summ)/d(i))
!print*,"X:",x(i)
xo(i)=x(i)

enddo !end-i

do i=1,n
v1= v1+abs(xo(i)**2)
v2= v2+abs(x(i)**2)
enddo
v1 =sqrt(v1)
v2=sqrt(v2)

if(abs(v2-v1) < tol) THEN
print*,abs(v2-v1)
goto 1
ELSE 
v1=0
v2=0
END IF

enddo !end-k

! print*,"Convergence Found at:",k
1 px=xo(1:5)
return 
end


Subroutine Get_Max(b,n,maxr,maxv)
implicit none
integer::k,i,n,maxr
double precision::b(n),maxv 

maxv=abs(b(1))
maxr=1
do i=2,n-1
IF(maxv.lt.abs(b(i))) THEN
maxv=abs(b(i))
maxr=i
ENDIF
enddo
return
end

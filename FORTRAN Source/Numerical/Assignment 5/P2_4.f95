Program P3_4
implicit none
integer:: n,filenum=7,itmax,i
double precision, allocatable, dimension (:) :: b,x,l,u,d
double precision, dimension(:):: px(5)
print*,"n= "
read*,n
allocate(l(n-1),u(n-1),d(n), b(n),x(n))



call Problem3_4v2(n,L,U,D,B)
call Jacv2(n,L,U,D,B,X,px)

open( unit=filenum, file='junk.txt', status='replace')

write(filenum,*)"RESULT:"
write(filenum,*) "X:"
do i=1,5
write(filenum,'(10f10.6)')px(i)
enddo
write(filenum,*) 
write(filenum,*)
close(filenum)

deallocate(l,u,d,b,x)
End program P3_4

!Broken apart LUD matricies that make up A
Subroutine Problem3_4v2(n,L,U,D,B)
IMPLICIT NONE
integer:: n,i,j,k,itmax
DOUBLE PRECISION :: l(n),u(n),d(n),b(n)
do i=1,n
d(i)=5
b(i)=1
enddo
do i=1,n-1,1
l(i)=2
u(i)=2
enddo
return 
end

!Modified Jacobi Method for tri diagonal matricies

Subroutine Jacv2(n,L,U,D,B,X,px)
implicit none
integer :: n,itmax=100,i,j,k
double precision :: l(n-1),u(n-1),d(n),b(n),x(n),xo(n),summ,tol=.0000000001,px(5)
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
!Grab Certain Values
if(i.eq.1) THEN
px(1) = x(i)
endif
if(i.eq.250000) THEN
px(2) = x(i)
endif
if(i.eq.500000) THEN
px(3) = x(i)
endif
if(i.eq.750000) THEN
px(4) = x(i)
endif
if(i.eq.100000) THEN
px(5) = x(i)
endif


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

1 print*,"Convergence Found at:",k

return 
end

Program A06
	IMPLICIT NONE
DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:) :: B,B2,D,X,X2,X3
DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:) :: A,Z
INTEGER :: N,n2,I,J,K,FILENUM,IOSTATUS,itmax






open(7,file='matrixd3.txt')

 do
  read(7,*,iostat=iostatus) 
	if(iostatus/=0) THEN !Logic Check-EOF
		exit
	else
		n=n+1
	end if
 end do	
print*,"Size of matrix:"
n=n-1
print*,n
allocate(a(n,n),z(n,n),b(n),b2(n),x(n),x2(n),d(n))

!Return to BOF
REWIND 7

do i=1,n
read(7,*)(a(i,j),j=1,n)
end do
print*,A

read(7,*)b
print*,B
close(7)
print*,"Please Enter number of max iterations: "
read*,itmax

call Jac(n,itmax,A,b,x)

call Gseidel(n,itmax,A,b,x2)

open( unit=filenum, file='matrix3results.txt', status='replace')

write(filenum,*)"RESULT:"
write(filenum,*) "A:"
do i=1,n
write(filenum,'(10f10.6)')A(i,1:n)
enddo
write(filenum,*) 
write(filenum,*)
write(filenum,*) "B:"
do i=1,n
write(filenum,'(10f10.6)')B(i)
enddo
write(filenum,*) 
write(filenum,*)
write(filenum,*) "X:"
do i=1,n
write(filenum,'(10f10.6)')x(i)
enddo
write(filenum,*) 
write(filenum,*)"X2:"
do i=1,n
write(filenum,'(10f10.6)')x2(i)
enddo

close(filenum)

deallocate(a,b,x,x2,d)	
End Program A06



Subroutine Jac(n,itmax,A,b,x)
implicit none
integer::n,i,j,k,itmax
double precision:: a(n,n),x(n),currentx(n),b(n),summ=0,break=0,tol=0.0000001,v1,v2
print*,"Jacobi:"

!Initial Guess that all are zero
do i=1,n
x(i)=0
enddo

!DO K TIMES (or until convergence)
do k=1,itmax
!Reset Break Condition
break=0
!Copy Current X Approximations
do j=1,n
	currentx(j) = x(j)
enddo
!As follows in Jacobi
do i=1,n
summ=0
	do j=1,n
	if(j/=i) THEN	
	summ = summ+(a(i,j)*x(j))
	ENDIF
	enddo

	x(i) = ((b(i)-summ)/a(i,i))
print*,x(i)

enddo

 do i=1,n
	v1= v1+abs(currentx(i)**2)
	v2= v2+abs(x(i)**2)
 enddo
	v1 =sqrt(v1)
	v2=sqrt(v2)
	if((v2-v1) < tol) THEN
		goto 1
	ELSE 
		v1=0
		v2=0
	END IF
enddo

1 print*,"Convergence Found at:",k
return 
end

Subroutine Gseidel(n,itmax,A,b,x2)
IMPLICIT NONE
integer::n,i,j,k,itmax
double precision:: a(n,n),x2(n),xo(n),b(n)
double precision:: summ=0,summb=0,break=0,tol=0.0000001,v1,v2
print*,"Gauss-Seidel:"

!Initial Guess
do i=1,n
x2(i)=1
enddo
do k=1,itmax

	do i=1,n
	summ=0
		do j=1,n
		if(j/=i) THEN
		summ= summ+(a(i,j)*xo(j))
		endif
		x2(i)=((b(i)-summ)/a(i,i))
		print*,x2(i)
		enddo!end-j
	xo(i)=x2(i)
	enddo !end-i

 do i=1,n
	v1= v1+abs(xo(i)**2)
	v2= v2+abs(x2(i)**2)
 enddo
	v1 =sqrt(v1)
	v2=sqrt(v2)
	if((v2-v1) < tol) THEN
		goto 1
	ELSE 
		v1=0
		v2=0
	END IF

enddo !end-k

1 print*,"Convergence Found at:",k

return
end


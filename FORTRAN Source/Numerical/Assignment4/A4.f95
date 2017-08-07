!Michael Timbes
!Numerical Analysis 
!Assignment-4 
PROGRAM A4
IMPLICIT NONE 
DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:) :: B,BCOP,X,XPLU,XFLU
DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:) :: A,ACOP,ATRI,APAR,AFUL
INTEGER, ALLOCATABLE, DIMENSION (:) :: P,Q
INTEGER :: N=0,R,I,J,FILENUM,IOSTATUS

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!This fortran script contains a few subroutines to find solutions of a ! 
!matrix through various methods.(1) Through standard Gauss elimination !
!through upper triangularization.(2) ParPLU Factorization using matrix !
!'P' that stores the row pivioting. (3) PLU Factorization that also    !
!uses a 'P' matrix for row pivioting but also includes code for column ! 
!pivoting."OPS' is used to measure the operations needed,result in out.!


open(7,file='mdata2.txt')

 do
  read(7,*,iostat=iostatus) 
	if(iostatus/=0) THEN !Logic Check-EOF
		exit
	else
		n=n+1
	end if
 end do							
!Allocate Storage Needed (based on number of lines in data file)
print*,"Size of matrix:"
!read*,n

n=n-1
print*,n

allocate(a(n,n),b(n),x(n),p(n),q(n),xflu(n))
allocate(xplu(n),acop(n,n),bcop(n),atri(n,n),apar(n,n),aful(n,n))

!Return to BOF
REWIND 7

do i=1,n
read(7,*)(a(i,j),j=1,n)
end do
print*,A

read(7,*)b

close(7)

print*,"Enter the B:"
!read*,b !Reads the last values of each row
print*,b
ACOP=A
BCOP=B
call TRIU(n,a,b)
ATRI=A
call BackSolve(n,a,b,x)
A=ACOP
B=BCOP
call ParPLU(n,a,p)
APAR=A
call ParPLUSol(n,a,b,p,xplu)
A=ACOP
B=BCOP
call FullPLU(n,a,p,q)
AFUL=A
call ParPLUSol(n,a,b,p,xflu)
filenum=9
open( unit=filenum, file='results.txt', status='replace')

write(filenum,*)"RESULT:"
write(filenum,*) "A after TRIU call:"
do i=1,n
write(filenum,'(10f10.6)')ATRI(i,1:n)
enddo
write(filenum,*) 
write(filenum,*)
write(filenum,*)"X Values:"
write(filenum,'(10f10.6)')x(1:n)
write(filenum,*)
write(filenum,*)"Matrix for the PLU:"
do i=1,n
write(filenum,'(10f10.6)')APAR(i,1:n)
enddo 
write(filenum,*)
write(filenum,*)
write(filenum,*)"X Values from PLU:"
write(filenum,'(10f10.6)')xplu(1:n)
write(filenum,*)
write(filenum,*)"Matrix for the FPLU:"
do i=1,n
write(filenum,'(10f10.6)')AFUL(i,1:n)
enddo
write(filenum,*)
write(filenum,*)
write(filenum,*)"X Values from FPLU:"
write(filenum,'(10f10.6)')xflu(1:n)
write(filenum,*)



close(filenum)

!Deallocate matricies
deallocate (a,b,x,p,q,xplu,xflu,acop,bcop)

END PROGRAM A4



!Upper triangular using Gauss elimination and multipliers are stored within the
!A matrix zero entries
Subroutine TRIU(N,A,B)
implicit none
integer:: n,i,j,k
double precision:: a(n,n), b(n), mult
!A is the coefficient matrix, B are the given solutions to the equations
!mult is the multiplier variable
do k=1,n-1 !Need to have the n-1 since there are no zeros in the last column

	do i=k+1,n !Go along the i th row values

		mult = a(i,k)/a(k,k)
		do j=1,n
			a(i,j)= a(i,j)-(mult*a(k,j)) !Subtract off multiplier  
		end do
	b(i) = b(i) - mult*b(k) !Also subtract off in the solutions
	end do
	
end do
return 
end

!Using the upper triangular matrix stored in the coefficient matrix,
!solve backwards for the X variables which are stored in X in main
Subroutine BackSolve(N,A,B,X)
implicit none 
integer:: n,i,j
double precision :: a(n,n), b(n), x(n),sub

!Initialize first X value-make sure divide by zero doesn't happen
if(a(n,n).NE.0) THEN
x(n) = b(n)/a(n,n)
else 
print*,"A(n,n) = 0 not allowed."
stop
end if
!Going from the next to last element in matrix A
do i=n-1,1,-1
	sub= b(i) !sub is the intermediate variable to help solve X(i)
		do j = i+1,n
			sub = sub-a(i,j)*x(j)
		end do
	x(i) = sub/a(i,i)
end do 
return 
end




!Subroutine for partial pivoting using the matrix A, the size (n), 
!and the permutation matrix P
Subroutine ParPLU(n,a,p)
implicit none
double precision :: a(n,n),mul,maxv,temp(n)
integer :: n,i,j,k,p(n), maxr

!Initialize the permutation matrix
do i=1,n
	p(i) = i
end do
!Begin to find the max value in the matrix
do k= 1,n-1
	maxv= abs(a(k,k))
	maxr= k
		do i = k+1,n
		if (abs(a(i,k)) >maxv) THEN
			maxr=i
			maxv= abs(a(i,k))
		end if
		end do
		!Condition for row swap is when the row with the maximum 
		!value is not the current row in the loop
	if(maxr.GT.k) THEN 
		p(k) = maxr
		temp(k:n) = a(maxr,k:n)
		a(maxr,k:n) = a(k,k:n)
		a(k,k:n) = temp(k:n)
	end if
	!Compute the multipliers 
	a(k+1:n,k) = a(k+1:n,k)/a(k,k)
		!Elimination process with multiplier from above
		do i = k+1,n
			a(i,k+1:n) = a(i,k+1:n) - a(i,k)*a(k,k+1:n)
		end do
end do
return 
end

!Finishes the back substitution for the ParPLU
Subroutine ParPLUSol(n,a,b,p,x)
implicit none 
double precision :: a(n,n),b(n),x(n),s,y(n)
integer :: n,i,j,k,p(n)

y=b !Mutator matrix

do k= 1,n-1
	if(p(k) > k) THEN
		s=y(k)
		y(k) = y(p(k))
		y(p(k)) = s
	end if 
	do i = k+1,n
		y(i) = y(i) - a(i,k)*y(k) !Gauss elimination done to the B
	end do 
end do 
!Solve for coefficients 
x(n) = y(n) / a(n,n)
	do i= n-1,1,-1
		s=y(i)
			do j = i+1,n
				s= s-a(i,j)*x(j)
			end do 
		x(i) = s/a(i,i)
	end do 
return 
end

!Subroutine for complete pivoting using the matrix A, the size (n), 
!and the permutation matrix P.
!Note that this algorithm is n^2+n since two max functions are being used
!as a result, the matrix values are being compared n times for the row, 
!then n times for the column, then the Gauss elimination is being done.

Subroutine FullPLU(n,a,p,q)
implicit none
double precision :: a(n,n),mul,maxv,temp(n)
integer :: n,i,j,k,p(n),q(n), maxr,maxc

!Initialize the permutation matrix
do i=1,n
	p(i) = i
end do
!Begin to find the max value in the matrix
do k= 1,n-1
	maxv= abs(a(k,k))
	maxr= k
	maxc= k
		!Find max row index
		do i = k+1,n
			if (abs(a(i,k)) >maxv) THEN
				maxr=i
				maxv= abs(a(i,k))
			end if
		end do

		
	!If the current row is not the max value
	if(maxr.GT.k) THEN 
		p(k) = maxr !Store change in the P matrix
		temp(k:n) = a(maxr,k:n)
		a(maxr,k:n) = a(k,k:n)
		a(k,k:n) = temp(k:n)
	
	end if
		!Find max column index
		maxv=abs(a(1,k)) !Reset Max 
		do i = k+1,n
			if (abs(a(1,i)) >maxv) THEN
				maxc=i
				maxv= abs(a(1,i))
			end if
		end do

	!If the current column is not the max value
	 if(maxc.GT.k) THEN
		q(k) = maxc !Store change in the Q matrix (not really needed)
		temp(1:n) = a(1:n,maxc)
		a(1:n,maxc) = a(1:n,k)
		a(1:n,k) = temp(1:n)
	end if

	!Compute the multipliers 
	a(k+1:n,k) = a(k+1:n,k)/a(k,k)
		!Elimination process with multiplier from above
		do i = k+1,n
			a(i,k+1:n) = a(i,k+1:n) - a(i,k)*a(k,k+1:n)
		end do
end do
return 
end


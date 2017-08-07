Program EigPow
      IMPLICIT NONE
      DOUBLE PRECISION, ALLOCATABLE, DIMENSION (:,:):: A
      DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:)::B,B2
      DOUBLE PRECISION::large=0,eig=0
      INTEGER::I,J,N=0,IOSTATUS,filenum=11

      open(7,file='A.txt')
     do
       read(7,*,iostat=iostatus) 
         if(iostatus/=0) THEN !Logic Check-EOF
           exit
         else
           n=n+1
         end if
     end do
     !print*,"Size of matrix:"
     n=n-1
     !print*,n
      allocate(A(n,n),b(n),b2(n))
      !Return to BOF
      REWIND 7

     do i=1,n
       read(7,*)(a(i,j),j=1,n)
     end do
      ! print*,A

       read(7,*)b
       !print*,B
     close(7)
!do i=1,n
!write(*,'(10f10.6)')A(i,1:n)
!enddo 
!write(*,*)
!write(*,*)"------------------------------"
b2=b
call EPow(A,B,n,large)
call InvPower(A,B2,N,eig)


open( unit=filenum, file='eigvals.txt', status='replace')

write(filenum,*)"RESULT:"
write(filenum,*) "A:"
do i=1,n
write(filenum,'(10f10.6)')A(i,1:n)
enddo
write(filenum,*)"------------------------" 
write(filenum,*)
write(filenum,*) "V:"
do i=1,n
write(filenum,'(10f10.6)')B(i)
enddo
write(filenum,*)"--------------------------"
write(filenum,*)"Largest Eigenvalue: "
write(filenum,'(10f10.6)')large
write(filenum,*)
write(filenum,*)
write(filenum,*) "Smallest Eig Vector:"
do i=1,n
write(filenum,'(10f10.6)')B2(i)
enddo
write(filenum,*)"--------------------------"
write(filenum,*)"Eigenvalue: "
write(filenum,'(10f10.6)')eig








close(filenum)
deallocate(A,B)

End Program EigPow

! (TAG:A) Soubroutine to take matrix A and multiply it by any vector(default is where all components are 1).
! (TAG:B) Find the largest value in the current iteration of the vector 'b' normally known as x in texts.
! (TAG:C) Normalize the B vector to prevent overflow. Also check for convergence.
Subroutine EPow(A,B,N,large)
IMPLICIT NONE
INTEGER:: N,I,J,it=1
DOUBLE PRECISION:: A(n,n),B(n),TEMP(N),old,large
! TAG:A
        1 do i=1,n
        temp(i)=0
                do j=1,n
                temp(i)=temp(i)+(A(i,j)*b(j))

                enddo
        enddo
        do i=1,n
                b(i)=temp(i)
        enddo
                
! TAG:B
        old=large
        large=abs(b(1))
        do i=2,n
        IF(ABS(b(i)) .GT. large) THEN
        large=ABS(b(i))
        ENDIF
        enddo
! TAG:C
        do i=1,n
        b(i)= b(i)/large
        enddo
        if(abs(large-old)>.000001)THEN
       it= it+1
         goto 1
        ENDIF

print*,"Iterations needed for problem 1:",it
RETURN
END

!Subroutine to solve problem two in the final HW assignment.
!(TAG:A) Create the AQ Matrix
!(TAG:B) Solve Linear System-Find largest value
!(TAG:C) Update the Xt vector-Reset AQ and note the iteration
Subroutine InvPower(A,B,N,eig)
IMPLICIT NONE
INTEGER:: N,I,J,itmax=1,P(n),ipos,it
DOUBLE PRECISION:: A(n,n),AQS(n,n),B(n),bsafe(n),xo(n),TMat(1,n),xoT(1,n),TMatA(n),AQ(n,n),E(n,n)
DOUBLE PRECISION::Q,tmp,tmpb,C=0,mul=0,mulold,eig
write(*,*)"Max iterations: "
read*,itmax
!TAG:A
	do i=1,n
		xo(i)=1
		
	enddo
xoT(1,:)=xo(:)
bsafe=b
tmp=0
tmpb=0
	do i=1,n
		TMatA(i)=0
		do j=1,n
			TMatA(i)=TMatA(i)+(A(i,j)*xo(j))
		enddo
	tmp=tmp+(TMatA(i)*xoT(1,i))
	tmpb=tmpb+(xoT(1,i)*xo(i))
	enddo
Q=tmp/tmpb
call EYE(E,N,Q) !Creates Identity matrix and handles multiplication of I*Q
!Creates the AQ matrix as in pg.380 in the Numerical Analysis book
do i=1,n
	do j=1,n
	IF(i.eq.j) THEN
	AQ(i,j)=A(i,j)-E(i,j)
	ELSE
	AQ(i,j)=A(i,j)
	ENDIF
	enddo
enddo
AQS=AQ
do j= 0,itmax
!Call Partial Pivoting Routine
call ParPLU(n,AQ,p,xo,b)

!Find Maximum Abs Value-Store into C
call Get_Max(b,n,ipos,c)
c=b(ipos)

!TAG:C
!B is the dummy matrix for the solution of the subroutine call, xo is being updated here.
do i=1,n 
xo(i) = 1/c*b(i)
enddo
mulold=mul
mul=q+(1/c)
IF(abs(mul-mulold)<.000001) THEN
goto 1
endif
it=j
AQ=AQS
enddo
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
1 eig=mulold
write(*,*)"Iterations needed for problem 2: ",it
RETURN
END


! Some utility subroutines to make things easier to read:
!Subroutine that creates the 2-d identiy matrix
Subroutine EYE(E,N,Q)
implicit none
integer::N,I,J
double precision::E(N,N),Q
	do i=1,n
	do j=1,n
	IF(i.eq.j) THEN 
		E(i,j)=1*Q
	ELSE
		E(i,j)=0
	END IF
	enddo
	enddo
return 
end

!PARTIAL PIVOTING SUBROUTINE FROM EARLIER ASSIGNMENT
Subroutine ParPLU(n,a,p,xo,b)
implicit none
double precision :: a(n,n),mul,maxv,temp(n),xo(n),b(n)
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

call ParPLUSol(n,a,xo,p,b) !Solves from within to maintain matrix information- also a bit easier to read.
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
IF(a(n,n)/=0) THEN 
x(n) = y(n) / a(n,n)
ELSE 
x(n)=0
endif


	do i= n-1,1,-1
		s=y(i)
			do j = i+1,n
				s= s-a(i,j)*x(j)
			end do 
		x(i) = s/a(i,i)

	end do 
!do i=1,n
!write(*,'(10f10.6)')x(i)
!enddo
!write(*,*)
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

!Michael Timbes
!A4-Prob:2-4
PROGRAM A4_CFACT
IMPLICIT NONE
DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:) :: B,BCOP,X
DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:) :: A,ACOP,L,GT,G
CHARACTER :: PD='F'
INTEGER :: N=0,I,J,FILENUM,IOSTATUS
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
allocate(a(n,n),l(n,n),gt(n,n),b(n),acop(n,n),bcop(n),x(n),g(n,n))
!Return to BOF
REWIND 7

do i=1,n
read(7,*)(a(i,j),j=1,n)
end do
do i=1,n
print'(10f10.6)',A(i,1:n)
enddo

ACOP=A
BCOP=B
close(7)
call CHO_F(n,a,g)
call TPOSE(n,g,gt,l,1)
A=ACOP
call PDEF(n,a,PD)
!WRITE OUT FILE
filenum=9
open( unit=filenum, file='CFACTRES.txt', status='replace')

write(filenum,*)"RESULT:"
write(filenum,*) "A:"
do i=1,n
write(filenum,'(15f15.6)')ACOP(i,1:n)
enddo
write(filenum,*)
write(filenum,*) "G:"
do i=1,n
write(filenum,'(15f15.10)')G(i,1:n)
enddo
write(filenum,*)
write(filenum,*) "G^T:"
do i=1,n
write(filenum,'(15f15.10)')GT(i,1:n)
enddo
write(filenum,*) 
write(filenum,*)
write(filenum,*) "G*G^T:"
do i=1,n
write(filenum,'(15f15.10)')l(i,1:n)
enddo
write(filenum,*)
write(filenum,*)"Positive Definate?"
write(filenum,*)PD 

 
CLOSE(FILENUM)

deallocate(a,b,acop,bcop,x,g,gt,l)
end program A4_CFACT

Subroutine CHO_F(n,a,G)
IMPLICIT NONE
INTEGER :: N,I,K,J 
DOUBLE PRECISION :: A(N,N),G(N,N)

!BEGIN THE FACTORIZATION

DO J=1,N
	IF(A(J,J).GT.0) THEN
	G(J,J) = SQRT(A(J,J))
	ELSE 
	PRINT*,"ERROR- SQRT OF A(J,J) NOT POSITIVE DEFINITE. END OF PROGRAM."
	STOP
	ENDIF
	!
	DO I=J+1,N
		G(I,J) = (A(I,J)/G(J,J))
		DO K= J+1,I
			A(I,K)=A(I,K) - (G(I,J)*G(K,J))
		ENDDO
	ENDDO
ENDDO


RETURN
END
Subroutine TPOSE(n,G,GT,l,comm)
IMPLICIT NONE
INTEGER :: N,I,K,J,comm
DOUBLE PRECISION :: g(n,n),gt(n,n),l(n,n),summ=0

do i=1,n
	gt(i,1:n)=g(1:n,i)
enddo
	!If the command var=1 then do 
	IF(COMM.eq.1) THEN
	do k=1,n 
		do i=1,n
		summ = 0 
			do j=1,n
			summ = summ +( g(i,j)*gt(j,k))
			enddo
		l(i,k)=summ
		enddo 
	enddo

	endif
return 
end
Subroutine PDEF(n,a,PD)
IMPLICIT NONE 
INTEGER :: N,I
DOUBLE PRECISION :: a(n,n)
CHARACTER:: PD
	do i=1,n
	if(a(i,i).LT.0) THEN
		PD='F'
		stop
	else
	PD='T'
	endif
	enddo
return 
end 
	

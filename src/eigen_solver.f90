!-----------------------------------------------------------------------
!Module: eigen_solver
!-----------------------------------------------------------------------
!! By: Benjamin Pieczynski
!!
!! This module solves for the eigenvalues(energies) and eigenvectors(wavefunctions).
!! This module contains more than just solve_eigenproblem, but those are
!! unused because the program uses lapack to solve the system of equations.
!!----------------------------------------------------------------------
!! Included subroutines:
!!
!! solve_eigenproblem
!!----------------------------------------------------------------------
module eigen_solver
use types
implicit none

private
public solve_eigenproblem

contains

! Here you can define a subroutine, like the one in the example we did 
! in class, that takes the diagonal and off-diagonal elements of a matrix,
! and returns the corresponding eigenvalues and eigenvectors. Whether 
! you use lapacks' dstev or the tqli subroutine provided below, be sure 
! to read the corresponding documentation.
! The arguments for each approach are slightly different.

! Again, notice how uninformative names and spaghetti code practices
! make it very hard to figure out what a subroutine or function does.
! If there was a bug in these pieces of code it would be very hard to 
! find.

!-----------------------------------------------------------------------
!! Subroutine: solve_eigenproblem
!-----------------------------------------------------------------------
!! by: Benjamin Pieczynski
!!
!! This subroutine determines the eigen values and eigen functions
!! for an input of diagaonals and offdiagonals using lapack. In this 
!! program, this subroutine is used to solve the schrodinger equation for 
!! differenr potential.
!!
!!----------------------------------------------------------------------
!! Input:
!! 
!!  d       real        array containing the diagonals
!!  e       real        array containing the off diagonals
!-----------------------------------------------------------------------
!! Output:
!!
!! lambdas      real        array containing the eigenvalues
!! vectors      real        array containing the eigenvectors
!-----------------------------------------------------------------------

subroutine solve_eigenproblem(d, e, lambdas, vectors)
    implicit none
    real(dp), intent(in) :: d(:), e(:)
    real(dp), intent(out) :: lambdas(:), vectors(:,:)

    character(len=1) :: jobz = 'V'
    integer :: n, info
    real(dp), allocatable :: e_work(:), work(:)

    n = size(d)

    lambdas = d
    allocate(e_work(1:size(e)))
    allocate(work(1:2*n-2))
    e_work = e
    call dstev(jobz, n, lambdas, e_work, vectors, n, work, info)

    if (info /= 0) then
        print*, 'Error while using dstev'
        print*, 'info = ', info
        stop
    endif

end subroutine solve_eigenproblem

!-----------------------------------------------------------------------
!! Subroutine: tqli
!-----------------------------------------------------------------------
!! by: Rodrigo Navarro Perez
!!
!! Modified from Numerical Recipes. Returns the eigenvalues and 
!! eigenvectors of a symmetric tridiagonal matrix.
!!
!! QL algorithm with implicit shifts, to determine the eigenvalues and 
!! eigenvectors of a real, symmetric, tridiagonal matrix.
!!
!! Unlike lapack's dstev, the eigenvalues are NOT ordered on output.
!! The subroutine eigsrt needs to be called after tqli to order
!! the eigenvalues and eigenvectors in ascending order.
!!
!!----------------------------------------------------------------------
!! Input:
!! d(:)     real(dp)    rank 1 array of size n. Contains the n diagonal 
!!                      elements of the tridiagonal matrix
!! e(:)     real(dp)    rank 1 array of size n. e(2:n) contains the
!!                      off-diagonal elements of the matrix.
!!                      e(1) is arbitrary.
!! z(:,:)   real(dp)    rank 2 array of size n by n. Contains the identity
!!                      matrix
!-----------------------------------------------------------------------
!! Output:
!! d(:)     real(dp)    Contains the eigenvalues of the tridiagonal matrix
!! e(:)     real(dp)    e is destroyed on output
!! z(:,:)   real(dp)    The k th column of z (i.e. z(:,k)) returns the 
!!                      normalized eigenvector corresponding to d(k)
!-----------------------------------------------------------------------
subroutine tqli(d, e, z)
    implicit none
    real(dp), intent(inout) :: d(:), e(:), z(:,:)
    integer :: n,i,iter,k,l,m
    real(dp) :: b,c,dd,f,g,p,r,s

    n = size(d)
    do i=2,n
        e(i-1)=e(i)
    enddo
    e(n)=0._dp

    do l=1,n
        iter=0
1       continue
        do m=l,n-1
            dd=abs(d(m))+abs(d(m+1))
            if (abs(e(m))+dd.eq.dd) goto 2
        enddo
        m=n
2       continue
        if (m.ne.l)then
            if(iter.eq.30) then
                print*, 'too many iterations in tqli'
                stop
            endif
            iter=iter+1
            g=(d(l+1)-d(l))/(2._dp*e(l))
            r=pythag(g,1._dp)
            g=d(m)-d(l)+e(l)/(g+sign(r,g))
            s=1._dp
            c=1._dp
            p=0._dp
            do i=m-1,l,-1
                f=s*e(i)
                b=c*e(i)
                r=pythag(f,g)
                e(i+1)=r
                if(r.eq.0._dp)then
                    d(i+1)=d(i+1)-p
                    e(m)=0._dp
                    goto 1
                endif
                s=f/r
                c=g/r
                g=d(i+1)-p
                r=(d(i)-g)*s+2._dp*c*b
                p=s*r
                d(i+1)=g+p
                g=c*r-b
                do k=1,n
                    f=z(k,i+1)
                    z(k,i+1)=s*z(k,i)+c*f
                    z(k,i)=c*z(k,i)-s*f
                enddo
            enddo
            d(l)=d(l)-p
            e(l)=g
            e(m)=0._dp
            goto 1
        endif
    enddo
    return
end subroutine tqli


!-----------------------------------------------------------------------
!! Function: pythag
!-----------------------------------------------------------------------
!! by: Rodrigo Navarro Perez
!!
!! Modified from Numerical Recipes. Computes (a^2 + b^2 )^(1/2) without
!! destructive underflow or overflow.
!!
!!----------------------------------------------------------------------
!! Input:
!! a        real(dp)    a real number
!! b        real(dp)    another real number
!-----------------------------------------------------------------------
!! Output:
!! r        real(dp)    (a^2 + b^2)^(1/2)
!-----------------------------------------------------------------------
real(dp) function pythag(a, b) result(r)
    implicit none
    real(dp), intent(in) :: a, b
    real(dp) absa,absb
    absa=abs(a)
    absb=abs(b)
    if(absa.gt.absb)then
        r=absa*sqrt(1._dp+(absb/absa)**2)
    else
        if(absb.eq.0._dp)then
            r=0._dp
        else
            r=absb*sqrt(1._dp+(absa/absb)**2)
        endif
    endif
    return
end function pythag

! 
!-----------------------------------------------------------------------
!! Subroutine: eigsrt
!-----------------------------------------------------------------------
!! by: Rodrigo Navarro Perez
!!
!! Modified from Numerical Recipes. Given the eigenvalues d and 
!! eigenvectors v as output from tqli, this routine sorts the eigenvalues
!! into ascending order, and rearranges the columns of v correspondingly.
!! The method is straight insertion.
!!
!!----------------------------------------------------------------------
!! Input:
!! d(:)     real(dp)    rank 1 array of size n. Contains the eigenvalues
!!                      of a matrix
!! z(:,:)   real(dp)    rank 2 array of size n by n. The k th column of
!!                      z (i.e. z(:,k)) contains the eigenvector
!!                      corresponding to d(k)
!-----------------------------------------------------------------------
!! Output:
!! d(:)     real(dp)    Ordered eigenvalues in ascending order
!! z(:,:)   real(dp)    Ordered eigenvectors according to the order of d
!-----------------------------------------------------------------------
subroutine eigsrt(d, v)
    implicit none
    real(dp), intent(inout) :: d(:), v(:,:)
    integer :: i,j,k,n
    real(dp) :: p
    n = size(d)
    do i=1,n-1
        k=i
        p=d(i)
        do j=i+1,n
            if(d(j).le.p)then
                k=j
                p=d(j)
            endif
        enddo
        if(k.ne.i)then
            d(k)=d(i)
            d(i)=p
            do j=1,n
                p=v(j,i)
                v(j,i)=v(j,k)
                v(j,k)=p
            enddo
        endif
    enddo
end subroutine eigsrt
    
end module eigen_solver
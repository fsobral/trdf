!     functions.f
!
!     This file contains the functions called by the TRDF algorithm
!     
!     - calfun: evaluates the objective function calobjf given by the
!               user
!     
!     and those that should be used by the nonlinear solver:
!
!     - evalf: evaluates the quadratic model
!     - evalg: evaluates the gradient of the quadratic model
!     - evalh: evaluates the Hessian of the quadratic model
!     - evalc: evaluates the user provided constraints (given by
!              user subroutine calcon)
!     - evaljac: evaluates the Jacobian of the constraints (given
!                by user subroutine caljac)
!     - evalhc: evaluates the Hessian of the constraints (given
!               by user subroutine calhc)


      subroutine calfun(n,x,f,flag)

!     SCALAR ARGUMENTS
      integer flag,n
      double precision f
      
!     ARRAY ARGUMENTS
      double precision x(n)

!     COMMON BLOCKS
      integer IC,MAXIC
      
      common /CONTA1/ IC, MAXIC

      flag = 0

C     Returns flag = 3 if reached the maximum of function evaluations
      if ( IC .eq. MAXIC ) then 
         flag = 3
         return
      end if

C     TODO: add 'flag' to calobjf
      call calobjf(n,x,f)

      IC = IC + 1

      end

      subroutine mevalf(n,x,f,flag)
      
      implicit none

#include "tr_params.par"

C     SCALAR ARGUMENTS
      integer flag,n
      double precision f

C     ARRAY ARGUMENTS
      double precision x(n)

!     COMMON BLOCKS
      double precision gopt_a(NMAX),hq_a(NMAX * NMAX)

      common /nometeste/ gopt_a,hq_a

!     LOCAL ARRAYS
      double precision hqd(n)

!     LOCAL SCALARS
      double precision gradd,dhqd

      flag = 0
      
!     avalia o modelo quadratico f, no ponto d
      
!     definindo o modelo quadratico f
      call mvv(GOPT_a,x, n, gradd)
      call mmv(HQ_a, x, n, hqd)
      call mvv(x, hqd, n, dhqd)

      f = gradd + dhqd / 2.0D0

      end

!     ******************************************************************
!     ******************************************************************

      subroutine mevalg(n,x,g,flag)
      
      implicit none

#include "tr_params.par"

C     SCALAR ARGUMENTS
      integer flag,n

C     ARRAY ARGUMENTS
      double precision g(n),x(n)

!     COMMON BLOCKS
      double precision gopt_a(NMAX),hq_a(NMAX * NMAX)

      common /nometeste/ gopt_a,hq_a

!     LOCAL ARRAYS
      double precision hqd(n)

!     LOCAL SCALARS
      integer i

      flag = 0
      
!     avalia o gradiente do modelo quadratico g, no ponto d
      call mmv(HQ_a, x, n, hqd)
      do i=1, n
         g(i)=  GOPT_a(i) + hqd(i)
      end do

      end

!     ******************************************************************
!     ******************************************************************

      subroutine mevalh(n,x,hrow,hcol,hval,hnnz,lim,lmem,flag)

      implicit none
      
#include "tr_params.par"
      
C     SCALAR ARGUMENTS
      logical lmem
      integer flag,n,hnnz,lim
      
C     ARRAY ARGUMENTS
      integer hcol(lim),hrow(lim)
      double precision hval(lim),x(n)

!     COMMON BLOCKS
      double precision gopt_a(NMAX),hq_a(NMAX * NMAX)

      common /nometeste/ gopt_a,hq_a

!     LOCAL SCALARS
      integer i,iii,j

      flag = 0
      lmem = .false.
      
      hnnz = (N + 1) * N / 2
      
!     cria hlin e hcol com as coordenadas da triangular inferior
      
      III=1
      DO WHILE (III .LE.  (N+1)*N/2 )
         do I=1, N
            DO J= 1, I
               hrow(III) = I
               hcol(III) = J
               III= III + 1
            end do
         end do
      END DO
  
      do i=1,(N+1)*N/2
         hval(i) = HQ_a(i)
      end do
  
      end

!     ******************************************************************
!     ******************************************************************

      subroutine mevalc(n,x,ind,c,flag)
      
      implicit none
      
#include "tr_params.par"
      
C     SCALAR ARGUMENTS
      integer ind,flag,n
      double precision c
      
C     ARRAY ARGUMENTS
      double precision x(n)

C     COMMON BLOCKS
      double precision XBASE_A(NMAX),XOPT_A(NMAX)

      COMMON /XBASEA/ XBASE_A
      COMMON /XOPTA/ XOPT_A

C     LOCAL ARRAYS
      double precision XA(NMAX)

C     LOCAL SCALARS
      integer i

      flag = 0
      
      DO I=1, N
         XA(I) = X(I) + XOPT_A(I) + XBASE_A(I)
      END DO
      
      call calcon(n,XA,ind,c)

      end

!     ******************************************************************
!     ******************************************************************

      subroutine mevaljac(n,x,ind,jcvar,jcval,jcnnz,lim,lmem,flag)

      implicit none

#include "tr_params.par"
      
C     SCALAR ARGUMENTS
      logical lmem
      integer flag,ind,jcnnz,lim,n

C     ARRAY ARGUMENTS
      integer jcvar(lim)
      double precision x(n),jcval(lim)

C     COMMON BLOCKS
      double precision XBASE_A(NMAX),XOPT_A(NMAX)

      COMMON /XBASEA/ XBASE_A
      COMMON /XOPTA/ XOPT_A

C     LOCAL ARRAYS
      double precision XA(NMAX)

C     LOCAL SCALARS
      integer i

      DO I=1, N
         XA(I) = X(I) + XOPT_A(I) + XBASE_A(I)
      END DO
      
      call caljac(n,XA,ind,jcvar,jcval,jcnnz,lim,lmem,flag)

      end

!     ******************************************************************
!     ******************************************************************

      subroutine mevalhc(n,x,ind,hcrow,hccol,hcval,hcnnz,lim,lmem,flag)

      implicit none

#include "tr_params.par"

C     SCALAR ARGUMENTS
      logical lmem
      integer flag,hcnnz,ind,lim,n

C     ARRAY ARGUMENTS
      integer hccol(lim),hcrow(lim)
      double precision hcval(lim),x(n)

C     COMMON BLOCKS
      double precision XBASE_A(NMAX),XOPT_A(NMAX)

      COMMON /XBASEA/ XBASE_A
      COMMON /XOPTA/ XOPT_A

C     LOCAL ARRAYS
      double precision XA(NMAX)

C     LOCAL SCALARS
      integer i

      DO i = 1,n
         XA(I) = X(I) + XOPT_A(I) + XBASE_A(I)
      END DO

      call calhc(n,XA,ind,hcrow,hccol,hcval,hcnnz,lim,lmem,flag)

      end

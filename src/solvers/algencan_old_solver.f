
      SUBROUTINE SOLVER(N, L, U, X, M, EQUATN, LINEAR, CCODED,
     +                  PHASE0, EPS, CNORM, FLAG)

!     This subroutine calls the user-supplied nonlinear programming
!     solver, in order to solve the trust-region subproblems. The user
!     who wants to use another solver, should write its own file
!     containing an implementation of this subroutine.
!
!     This subroutine calls the nonlinear solver ALGENCAN version ????.
!
!     The algorithm should be able to find only feasible points or to
!     solve the trust-region subproblems, if the value of
!     variable 'phase0' is .true. .

      implicit none

#include "tr_params.par"

!     SCALAR ARGUMENTS
      integer flag,n,m
      logical phase0
      double precision cnorm,eps

!     ARRAY ARGUMENTS
      logical ccoded(2),equatn(m),linear(m)
      double precision l(n),u(n),x(n)

!     LOCAL SCALARS
      logical checkder
      integer i,inform,iprint,ncomp
      double precision epsfeas,epsopt,f,nlpsupn,snorm

!     LOCAL ARRAYS
      logical coded(10)
      double precision lambda(MMAX)

C     EXTERNAL SUBROUTINES
      EXTERNAL ALGENCAN

      do i = 1,m
         lambda(i) = 0.0D0
      end do

!     Coded subroutines

      do i = 1,4
         coded(i) =    .true.
      end do
      coded(5)    = ccoded(1)
      coded(6)    = ccoded(2)
      do i = 7,10
         coded(i) =   .false.
      end do

      checkder = .false.

!     Detects if the algorithm should ignore the objective function in
!     order to find only a feasible point

      if (phase0) then
         coded(1) = .false.
         coded(2) = .false.
         coded(3) = .false.
      end if

!     Parameters setting

      EPSFEAS  = eps

      EPSOPT   = eps

      IPRINT   = 0

      NCOMP    = 6

!     Optimize

      flag = 0

      CALL ALGENCAN(EPSFEAS,EPSOPT,IPRINT,NCOMP,N,X,L,U,M,LAMBDA,EQUATN,
     +LINEAR,CODED,CHECKDER,F,CNORM,SNORM,NLPSUPN,INFORM)

      if ( inform .ne. 0 ) flag = 2

      end

C     ******************************************************************
C     ******************************************************************

      subroutine evalf(n,x,f,flag)

      implicit none

C     SCALAR ARGUMENTS
      integer flag,n, i
      double precision f

C     ARRAY ARGUMENTS
      double precision x(n)

      call mevalf(n,x,f,flag)

      end

C     ******************************************************************
C     ******************************************************************

      subroutine evalg(n,x,g,flag)

      implicit none

C     SCALAR ARGUMENTS
      integer flag, n, i

C     ARRAY ARGUMENTS
      double precision g(n),x(n)

      call mevalg(n,x,g,flag)

      end      

C     ******************************************************************
C     ******************************************************************

      subroutine evalh(n,x,hlin,hcol,hval,hnnz,flag)

      implicit none

#include "tr_params.par"

C     SCALAR ARGUMENTS
      integer flag,n,hnnz, i, j, III

C     ARRAY ARGUMENTS
      integer hcol(HCNNZMAX),hlin(HCNNZMAX)
      double precision hval(HCNNZMAX),x(n)

C     LOCAL SCALARS
      logical lmem

      call mevalh(n,x,hlin,hcol,hval,hnnz,HCNNZMAX,lmem,flag)

      end

C     ******************************************************************
C     ******************************************************************

      subroutine evalc(n,x,ind,c,flag)

      implicit none

C     SCALAR ARGUMENTS
      integer ind,flag,n
      double precision c

C     ARRAY ARGUMENTS
      double precision x(n)

      call mevalc(n,x,ind,c,flag)

      end

C     ******************************************************************
C     ******************************************************************

      subroutine evaljac(n,x,ind,jcvar,jcval,jcnnz,flag)

      implicit none

C     SCALAR ARGUMENTS
      integer flag,ind,jcnnz,n

C     ARRAY ARGUMENTS
      integer jcvar(n)
      double precision x(n),jcval(n)

C     LOCAL SCALARS
      logical lmem

      call mevaljac(n,x,ind,jcvar,jcval,jcnnz,n,lmem,flag)

      end

C     ******************************************************************
C     ******************************************************************

      subroutine evalhc(n,x,ind,hclin,hccol,hcval,hcnnz,flag)

      implicit none

#include "tr_params.par"

C     SCALAR ARGUMENTS
      integer flag,hcnnz,ind,n

C     ARRAY ARGUMENTS
      integer hccol(HCNNZMAX),hclin(HCNNZMAX)
      double precision hcval(HCNNZMAX),x(n)

C     LOCAL SCALARS
      logical lmem

      call mevalhc(n,x,ind,hclin,hccol,hcval,hcnnz,HCNNZMAX,lmem,flag)

      end

C     ******************************************************************
C     ******************************************************************

      subroutine evalfc(n,x,f,m,c,flag)

      implicit none

C     SCALAR ARGUMENTS
      integer flag,m,n
      double precision f

C     ARRAY ARGUMENTS
      double precision c(m),x(n)

      flag = - 1

      end

C     ******************************************************************
C     ******************************************************************

      subroutine evalgjac(n,x,g,m,jcfun,jcvar,jcval,jcnnz,lim,lmem,
     +flag)

      implicit none

C     SCALAR ARGUMENTS
      logical lmem
      integer flag,jcnnz,lim,m,n

C     ARRAY ARGUMENTS
      integer jcfun(lim),jcvar(lim)
      double precision g(n),jcval(lim),x(n)

      flag = - 1

      end

C     ******************************************************************
C     ******************************************************************

      subroutine evalgjacp(n,x,g,m,p,q,work,gotj,flag)

      implicit none

C     SCALAR ARGUMENTS
      logical gotj
      integer flag,m,n
      character work

C     ARRAY ARGUMENTS
      double precision g(n),p(m),q(n),x(n)

      flag = - 1

      end

C     ******************************************************************
C     ******************************************************************

      subroutine evalhl(n,x,m,lambda,sf,sc,hlrow,hlcol,hlval,hlnnz,
     +lim,lmem,flag)

      implicit none

C     SCALAR ARGUMENTS
      logical lmem
      integer flag,hlnnz,lim,m,n
      double precision sf

C     ARRAY ARGUMENTS
      integer hlcol(lim),hlrow(lim)
      double precision hlval(lim),lambda(m),sc(m),x(n)

      flag = - 1

      end

C     ******************************************************************
C     ******************************************************************

      subroutine evalhlp(n,x,m,lambda,sf,sc,p,hp,goth,flag)

      implicit none

C     SCALAR ARGUMENTS
      logical goth
      integer flag,m,n
      double precision sf

C     ARRAY ARGUMENTS
      double precision hp(n),lambda(m),p(n),sc(m),x(n)

      flag = - 1

      end

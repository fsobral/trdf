
      SUBROUTINE SOLVER(N, L, U, X, M, EQUATN, LINEAR, CCODED,
     +     MEVALF, MEVALG, MEVALH, MEVALC, MEVALJAC, MEVALHC,
     +     PHASE0, EPS, CNORM, FLAG)

!     This subroutine calls the user-supplied nonlinear programming
!     solver, in order to solve the trust-region subproblems. The user
!     who wants to use another solver, should write its own file
!     containing an implementation of this subroutine.
!
!     This subroutine calls the nonlinear solver ALGENCAN.
!
!     The algorithm should be able to find only feasible points or to
!     solve the trust-region subproblems, if the value of
!     variable 'phase0' is .true. .
!
!     This subroutine should return FLAG = 0 if it has terminated
!     correctly, or FLAG = 2 otherwise.

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
      integer i,inform,hnnzmax,nvparam
      double precision efacc,efstain,eoacc,eostain,epsfeas,epsopt,
     +     f,nlpsupn,snorm

!     LOCAL ARRAYS
      character * 80 specfnm,outputfnm,vparam(10)
      logical coded(11)
      double precision lambda(MMAX)

!     EXTERNAL SUBROUTINES
      external mevalf,mevalg,mevalh,mevalc,mevaljac,mevalhc,
     +     evalfc,evalgjac,evalgjacp,evalhl,evalhlp

      do i = 1,m
         lambda(i) = 0.0D0
      end do

!     Coded subroutines

      do i = 1,4
         coded(i) =    .true.
      end do
      coded(5)    = ccoded(1)
      coded(6)    = ccoded(2)
      do i = 7,11
         coded(i) =   .false.
      end do

      checkder = .false.

!     Upper bounds on the number of sparse-matrices non-null elements

      hnnzmax  = ((n + 1) * n / 2) * (1 + m)

!     Parameters setting

      epsfeas   = eps
      epsopt    = eps
      
      efstain   = sqrt( epsfeas )
      eostain   = epsopt ** 1.5d0

      efacc     = sqrt( epsfeas )
      eoacc     = sqrt( epsopt )

      outputfnm = ''
      specfnm   = ''
      
      nvparam = 0

      nvparam = nvparam + 1
      vparam(nvparam) = 'SAFEMODE'

!     Detects if the algorithm should ignore the objective function in
!     order to find only a feasible point
      if (phase0) then
         nvparam         =                 nvparam + 1
         vparam(nvparam) = 'IGNORE-OBJECTIVE-FUNCTION'
      end if

!     Optimize

      flag = 0

      call algencan(mevalf,mevalg,mevalh,mevalc,mevaljac,mevalhc,
     +     evalfc,evalgjac,evalgjacp,evalhl,evalhlp,JCNNZMAX,
     +     hnnzmax,epsfeas,epsopt,efstain,eostain,efacc,eoacc,outputfnm,
     +     specfnm,nvparam,vparam,n,x,l,u,m,lambda,equatn,linear,coded,
     +     checkder,f,cnorm,snorm,nlpsupn,inform)

      if ( inform .ne. 0 ) flag = 2

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

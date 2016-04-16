SUBROUTINE EASYCTRDF(N,X,XL,XU,M,EQUATN,LINEAR,CCODED,F,FEAS,FCNT) &
     BIND(C, name="easytrdf")
  
  use iso_c_binding, only : c_bool, c_double, c_int

  implicit none
  
  ! SCALAR ARGUMENTS
  integer(kind=c_int), value :: m,n
  integer(kind=c_int)        :: fcnt
  real(kind=c_double)        :: f,feas
  
  ! ARRAY ARGUMENTS
  real(kind=c_double) :: X(N),XL(N),XU(N)
  logical(kind=c_bool) :: ccoded(2),equatn(m),linear(m)
  
  ! LOCAL ARRAYS
  logical :: ccoded_(2),equatn_(m),linear_(m)
  
  ccoded_(1:2) = logical( ccoded(1:2) )
  equatn_(1:m) = logical( equatn(1:m) )
  linear_(1:m) = logical( linear(1:m) )

  CALL EASYTRDF(N,X,XL,XU,M,EQUATN_,LINEAR_,CCODED_,F,FEAS,FCNT)
  
END SUBROUTINE EASYCTRDF

SUBROUTINE CTRDF(N,NPT,X,XL,XU,M,EQUATN,LINEAR,CCODED,MAXFCNT,RBEG, &
     REND,XEPS,F,FEAS,FCNT) BIND(C, name="trdf")
  
  use iso_c_binding, only : c_bool, c_double, c_int

  implicit none
  
  ! SCALAR ARGUMENTS
  integer(kind=c_int), value :: m,maxfcnt,n,npt
  integer(kind=c_int)        :: fcnt
  real(kind=c_double), value :: rbeg,rend,xeps
  real(kind=c_double)        :: f,feas
  
  ! ARRAY ARGUMENTS
  real(kind=c_double) :: X(N),XL(N),XU(N)
  logical(kind=c_bool) :: ccoded(2),equatn(m),linear(m)
  
  ! LOCAL ARRAYS
  logical :: ccoded_(2),equatn_(m),linear_(m)
  
  ccoded_(1:2) = logical( ccoded(1:2) )
  equatn_(1:m) = logical( equatn(1:m) )
  linear_(1:m) = logical( linear(1:m) )

  CALL TRDF(N,NPT,X,XL,XU,M,EQUATN_,LINEAR_,CCODED_,MAXFCNT,RBEG, &
            REND,XEPS,F,FEAS,FCNT)
  
END SUBROUTINE CTRDF

! ******************************************************************
! ******************************************************************

subroutine calobjf(n,x,f)
  
  use iso_c_binding, only : c_double, c_int
  
  implicit none
  
  ! SCALAR ARGUMENTS
  integer(kind=c_int) :: n
  real(kind=c_double) :: f
  
  ! ARRAY ARGUMENTS
  real(kind=c_double) :: x(n)
  
  ! INTERFACES
  interface
     subroutine c_calobjf(n,x,f) bind(C)
       import :: c_double, c_int
       
       integer(kind=c_int), value :: n
       real(kind=c_double)        :: f,x(n)
     end subroutine c_calobjf
  end interface
  
  call c_calobjf(n,x,f)
  
end subroutine calobjf

! ******************************************************************
! ******************************************************************

subroutine calcon(n,x,ind,c)

  use iso_c_binding, only : c_double, c_int
  
  implicit none

  ! SCALAR ARGUMENTS
  integer(kind=c_int) :: ind,n
  real(kind=c_double) :: c

  ! ARRAY ARGUMENTS
  real(kind=c_double) :: x(n)

  ! INTERFACES
  interface
     subroutine c_calcon(n,x,ind,c) bind(C)
       import :: c_double, c_int
       
       integer(kind=c_int), value :: ind,n
       real(kind=c_double)        :: c,x(n)
     end subroutine c_calcon
  end interface
  
  call c_calcon(n,x,ind - 1,c)

end subroutine calcon

! ******************************************************************
! ******************************************************************

subroutine caljac(n,x,ind,jcvar,jcval,jcnnz,lim,lmem,flag)

  use iso_c_binding, only : c_bool, c_double, c_int
  
  implicit none

  ! SCALAR ARGUMENTS
  logical :: lmem
  integer(kind=c_int) :: flag,ind,jcnnz,lim,n

  ! ARRAY ARGUMENTS
  integer(kind=c_int) :: jcvar(lim)
  real(kind=c_double) :: x(n),jcval(lim)

  ! INTERFACES
  interface
     subroutine c_caljac(n,x,ind,jcvar,jcval,jcnnz,lim,lmem,flag) bind(C)
       import :: c_bool, c_double, c_int

       integer(kind=c_int),  value :: n,ind,lim
       integer(kind=c_int)         :: flag,jcnnz,jcvar(lim)
       logical(kind=c_bool)        :: lmem
       real(kind=c_double)         :: x(n),jcval(lim)
     end subroutine c_caljac
  end interface

  ! LOCAL SCALARS
  logical(kind=c_bool) :: lmem_
  
  call c_caljac(n,x,ind - 1,jcvar,jcval,jcnnz,lim,lmem_,flag)

  jcvar(1:jcnnz) = jcvar(1:jcnnz) + 1

  lmem = logical(lmem_)
  
end subroutine caljac


subroutine calhc(n,x,ind,hcrow,hccol,hcval,hcnnz,lim,lmem,flag)

  use iso_c_binding, only : c_bool, c_double, c_int
  
  implicit none

  ! SCALAR ARGUMENTS
  logical :: lmem
  integer(kind=c_int) :: flag,hcnnz,ind,lim,n

  ! ARRAY ARGUMENTS
  integer(kind=c_int) :: hccol(lim),hcrow(lim)
  real(kind=c_double) :: hcval(lim),x(n)

  ! INTERFACES
  interface
     subroutine c_calhc(n,x,ind,hcrow,hccol,hcval,hcnnz,lim,lmem,flag) bind(C)
       import :: c_bool, c_double, c_int

       integer(kind=c_int),  value :: n,ind,lim
       integer(kind=c_int)         :: flag,hcnnz,hccol(lim),hcrow(lim)
       logical(kind=c_bool)        :: lmem
       real(kind=c_double)         :: hcval(lim),x(n)
     end subroutine c_calhc
  end interface

  ! LOCAL SCALARS
  logical(kind=c_bool) :: lmem_
  
  call c_calhc(n,x,ind - 1,hcrow,hccol,hcval,hcnnz,lim,lmem_,flag)

  hcrow(1:hcnnz) = hcrow(1:hcnnz) + 1
  hccol(1:hcnnz) = hccol(1:hcnnz) + 1

  lmem = logical(lmem_)

end subroutine calhc
  

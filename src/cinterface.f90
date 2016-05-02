module cinterface

  ! Interface between user-defined C subroutines and the Fortran
  ! subroutines.
  !
  !
  ! TODO: A special attention should be given to the last 4 interface
  ! subroutines. They use the special types defined in the
  ! iso_c_binding module. Since they are actually Fortran subroutines,
  ! the correct way of doing this is to use tradicional Fortran sizes
  ! (or types) and converting them to the C types, as it is being done
  ! with 'logical' types.

  use iso_c_binding, only : c_bool, c_char, c_double, c_int, &
       c_funptr, c_f_procpointer

  implicit none

  ! INTERFACES

  abstract interface

     subroutine evalf(n,x,f,flag) bind(C)
       import :: c_double,c_int

       ! SCALAR ARGUMENTS
       integer(kind=c_int), value :: n
       integer(kind=c_int)        :: flag
       real(kind=c_double)        :: f
       ! ARRAY ARGUMENTS
       real(kind=c_double)        :: x(n)

       intent(in ) :: n,x
       intent(out) :: f,flag
     end subroutine evalf

     subroutine evalc(n,x,ind,c,flag) bind(C)
       import :: c_double,c_int

       ! SCALAR ARGUMENTS
       integer(kind=c_int), value :: ind,n
       integer(kind=c_int)        :: flag
       real(kind=c_double)        :: c
       ! ARRAY ARGUMENTS
       real(kind=c_double) :: x(n)

       intent(in ) :: ind,n,x
       intent(out) :: c,flag
     end subroutine evalc

     subroutine evaljac(n,x,ind,jcvar,jcval,jcnnz,lim,lmem,flag) &
          bind(C)
       import :: c_double,c_int,c_bool
       
       ! SCALAR ARGUMENTS
       integer(kind=c_int), value :: ind,lim,n
       integer(kind=c_int)        :: flag,jcnnz
       real(kind=c_double)        :: f
       logical(kind=c_bool)       :: lmem
       ! ARRAY ARGUMENTS
       integer(kind=c_int)        :: jcvar(lim)
       real(kind=c_double)        :: jcval(lim),x(n)

       intent(in ) :: ind,lim,n,x
       intent(out) :: flag,jcnnz,jcval,jcvar,lmem
     end subroutine evaljac

     subroutine evalhc(n,x,ind,hcrow,hccol,hcval,hcnnz,lim,lmem, &
          flag) bind(C)
       import :: c_double,c_int,c_bool
       
       ! SCALAR ARGUMENTS
       integer(kind=c_int), value :: ind,lim,n
       integer(kind=c_int)        :: flag,hcnnz
       real(kind=c_double)        :: f
       logical(kind=c_bool)       :: lmem
       ! ARRAY ARGUMENTS
       integer(kind=c_int)        :: hccol(lim),hcrow(lim)
       real(kind=c_double)        :: hcval(lim),x(n)
       
       intent(in ) :: ind,lim,n,x
       intent(out) :: flag,hccol,hcrow,hcval,hcnnz,lmem
     end subroutine evalhc
     
  end interface

  ! GLOBAL PROCEDURE POINTERS TO C USER PROVIDED ROUTINES
  procedure(evalf),     pointer :: c__evalf
  procedure(evalc),     pointer :: c__evalc
  procedure(evaljac),   pointer :: c__evaljac
  procedure(evalhc),    pointer :: c__evalhc

  private

  public ::  easyctrdf,fullctrdf

contains

  SUBROUTINE EASYCTRDF(N,X,XL,XU,M,EQUATN,LINEAR,CCODED,EVALF_,EVALC_, &
       EVALJAC_,EVALHC_,F,FEAS,FCNT) BIND(C, name="easytrdf")

    implicit none

    ! SCALAR ARGUMENTS
    integer(kind=c_int), value :: m,n
    integer(kind=c_int)        :: fcnt
    real(kind=c_double)        :: f,feas

    ! ARRAY ARGUMENTS
    real(kind=c_double) :: X(N),XL(N),XU(N)
    logical(kind=c_bool) :: ccoded(2),equatn(m),linear(m)

    ! EXTERNAL C SUBROUTINES
    type(c_funptr), value :: evalf_,evalc_,evaljac_,evalhc_

    ! LOCAL ARRAYS
    logical :: ccoded_(2),equatn_(m),linear_(m)

    call c_f_procpointer(evalf_,  c__evalf  )
    call c_f_procpointer(evalc_,  c__evalc  )
    call c_f_procpointer(evaljac_,c__evaljac)
    call c_f_procpointer(evalhc_, c__evalhc )

    ccoded_(1:2) = logical( ccoded(1:2) )
    equatn_(1:m) = logical( equatn(1:m) )
    linear_(1:m) = logical( linear(1:m) )

    CALL EASYTRDF(N,X,XL,XU,M,EQUATN_,LINEAR_,CCODED_, &
         CALOBJF,CALCON,CALJAC,CALHC,F,FEAS,FCNT)

  END SUBROUTINE EASYCTRDF

  SUBROUTINE FULLCTRDF(N,NPT,X,XL,XU,M,EQUATN,LINEAR,CCODED, &
       EVALF,EVALC_,EVALJAC_,EVALHC_,MAXFCNT,RBEG,REND,XEPS,&
       F,FEAS,FCNT) BIND(C, name="fulltrdf")

    implicit none

    ! SCALAR ARGUMENTS
    integer(kind=c_int), value :: m,maxfcnt,n,npt
    integer(kind=c_int)        :: fcnt
    real(kind=c_double), value :: rbeg,rend,xeps
    real(kind=c_double)        :: f,feas

    ! ARRAY ARGUMENTS
    real(kind=c_double) :: X(N),XL(N),XU(N)
    logical(kind=c_bool) :: ccoded(2),equatn(m),linear(m)

    ! EXTERNAL C SUBROUTINES
    type(c_funptr), value :: evalf,evalc_,evaljac_,evalhc_

    ! LOCAL ARRAYS
    logical :: ccoded_(2),equatn_(m),linear_(m)

    call c_f_procpointer(evalf,  c__evalf  )
    call c_f_procpointer(evalc_,  c__evalc  )
    call c_f_procpointer(evaljac_,c__evaljac)
    call c_f_procpointer(evalhc_, c__evalhc )

    ccoded_(1:2) = logical( ccoded(1:2) )
    equatn_(1:m) = logical( equatn(1:m) )
    linear_(1:m) = logical( linear(1:m) )

    CALL FULLTRDF(N,NPT,X,XL,XU,M,EQUATN_,LINEAR_,CCODED_,CALOBJF, &
         CALCON,CALJAC,CALHC,MAXFCNT,RBEG,REND,XEPS,F,FEAS,FCNT)

  END SUBROUTINE FULLCTRDF

  ! ******************************************************************
  ! ******************************************************************

  subroutine calobjf(n,x,f,flag)

    implicit none

    ! SCALAR ARGUMENTS
    integer(kind=c_int) :: flag,n
    real(kind=c_double) :: f

    ! ARRAY ARGUMENTS
    real(kind=c_double) :: x(n)
    
    intent(in ) :: n,x
    intent(out) :: f,flag

    call c__evalf(n,x,f,flag)

  end subroutine calobjf

  ! ******************************************************************
  ! ******************************************************************

  subroutine calcon(n,x,ind,c,flag)

    implicit none

    ! SCALAR ARGUMENTS
    integer(kind=c_int) :: ind,flag,n
    real(kind=c_double) :: c

    ! ARRAY ARGUMENTS
    real(kind=c_double) :: x(n)

    intent(in ) :: ind,n,x
    intent(out) :: c,flag

    call c__evalc(n,x,ind - 1,c,flag)

  end subroutine calcon

  ! ******************************************************************
  ! ******************************************************************

  subroutine caljac(n,x,ind,jcvar,jcval,jcnnz,lim,lmem,flag)

    implicit none

    ! SCALAR ARGUMENTS
    integer(kind=c_int) :: ind,flag,jcnnz,lim,n
    real(kind=c_double) :: f
    logical             :: lmem
    ! ARRAY ARGUMENTS
    integer(kind=c_int) :: jcvar(lim)
    real(kind=c_double) :: jcval(lim),x(n)

    intent(in ) :: ind,lim,n,x
    intent(out) :: flag,jcnnz,jcval,jcvar,lmem

    ! LOCAL SCALARS
    logical(kind=c_bool) :: lmem_

    call c__evaljac(n,x,ind - 1,jcvar,jcval,jcnnz,lim,lmem_,flag)

    jcvar(1:jcnnz) = jcvar(1:jcnnz) + 1

    lmem = logical(lmem_)

  end subroutine caljac

  ! ******************************************************************
  ! ******************************************************************

  subroutine calhc(n,x,ind,hcrow,hccol,hcval,hcnnz,lim,lmem,flag)

    implicit none

    ! SCALAR ARGUMENTS
    integer(kind=c_int) :: ind,lim,n
    integer(kind=c_int) :: flag,hcnnz
    real(kind=c_double) :: f
    logical             :: lmem
    ! ARRAY ARGUMENTS
    integer(kind=c_int) :: hccol(lim),hcrow(lim)
    real(kind=c_double) :: hcval(lim),x(n)

    intent(in ) :: ind,lim,n,x
    intent(out) :: flag,hccol,hcrow,hcval,hcnnz,lmem

    ! LOCAL SCALARS
    logical(kind=c_bool) :: lmem_

    call c__evalhc(n,x,ind - 1,hcrow,hccol,hcval,hcnnz,lim,lmem_,flag)

    hcrow(1:hcnnz) = hcrow(1:hcnnz) + 1
    hccol(1:hcnnz) = hccol(1:hcnnz) + 1

    lmem = logical(lmem_)

  end subroutine calhc

end module cinterface

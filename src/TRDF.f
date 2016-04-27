      SUBROUTINE EASYTRDF(N,X,XL,XU,M,EQUATN,LINEAR,CCODED,EVALF,EVALC,
     +     EVALJAC,EVALHC,F,FEAS,FCNT)

      use trdf

      IMPLICIT NONE

C     This subroutine calls the true TRDF subroutine using default
C     configuration parameters

C     SCALAR ARGUMENTS
      integer M,N,FCNT
      double precision F,FEAS

C     ARRAY ARGUMENTS
      double precision  X(N),XL(N),XU(N)
      logical CCODED(2),EQUATN(M),LINEAR(M)

C     LOCAL SCALARS
      integer MAXFCNT,NPT
      double precision RBEG,REND,XEPS

C     EXTERNAL SUBROUTINES
      external EVALF,EVALC,EVALJAC,EVALHC

C     NUMBER OF INTERPOLATION POITNS. 

      if ( N .le. 2 ) then
         NPT = 2 * N + 1
      else
         NPT = 2 * N + 3
      end if

C     SETS DEFAULT PARAMETERS
      
      MAXFCNT = 100000
      
      RBEG = 1.0D-1
      REND = 1.0D-4
      XEPS = 1.0D-8

C     CALLS THE TRUE SUBROUTINE

      CALL TRDFSUB(N,NPT,X,XL,XU,M,EQUATN,LINEAR,CCODED,EVALF,EVALC,
     +             EVALJAC,EVALHC,MAXFCNT,RBEG,REND,XEPS,F,FEAS,FCNT)

      END

C     ******************************************************************
C     ******************************************************************

      SUBROUTINE FULLTRDF(N,NPT,X,XL,XU,M,EQUATN,LINEAR,CCODED,EVALF,
     +                    EVALC,EVALJAC,EVALHC,MAXFCNT,RBEG,REND,XEPS,F,
     +                    FEAS,FCNT)     

C     This subroutine is the implementation of the Derivative-free
C     Trust-region algorithm for constrained optimization described in
C
C     P. D. Conejo, E. W. Karas, L. G. Pedroso, "A trust-region
C     derivative-free algorithm for constrained optimization".
C     Optimization Methods & Software, to appear, 2015.
C
C     For more information on configurations, recommended values of the
C     parameters and examples, please read the documentation provided
C     with the method.
C
C     The INPUT parameters are
C
C     N - dimension of the problem
C     
C     NPT - number of points used in the quadratic interpolation
C
C     X(N) - inital point
C
C     XL(N) - lower bounds on the variables
C
C     XU(N) - upper bounds on the variables
C
C     M - number of constraints
C
C     EQUATN(M) - a logical vector such that EQUATN(i) = .true. if and
C                 only if constraint 'i' is an equality constraint
C
C     LINEAR(M) - a logical vector such that LINEAR(i) = .true. if and
C                 only if constraint 'i' is a linear constraint
C
C     CCODED(2) - a logical vector indicating whether or not the
C                 Jacobian (CCODED(1)) and the Hessian (CCODED(2))
C                 of the constraints are provided by the user
C
C     MAXFCNT - maximum number of function evaluations
C
C     RBEG - initial value for the trust region radius and interpolation
C
C     REND - smallest trust-region radius for declaring convergence of
C            the method
C
C     XEPS - feasibility/optimality tolerance for solving the subproblems
C
C
C     The OUTPUT parameters are
C
C     X(N) - the final point
C
C     F - objective function value at the final point
C
C     FEAS - sup-norm of the infeasibility at the final point
C
C     FCNT - number of function evaluations

      use trdf

      IMPLICIT NONE

c$$$#include "tr_params.par"

C     SCALAR ARGUMENTS
      integer m,maxfcnt,N,NPT,FCNT
      double precision F,FEAS,RBEG,REND,XEPS

C     ARRAY ARGUMENTS
      DOUBLE PRECISION  X(N),XL(N),XU(N)
      logical ccoded(2),equatn(m),linear(m)

C     EXTERNAL SUBROUTINES
      external EVALF,EVALC,EVALJAC,EVALHC


      CALL TRDFSUB(N,NPT,X,XL,XU,M,EQUATN,LINEAR,CCODED,EVALF,EVALC,
     +             EVALJAC,EVALHC,MAXFCNT,RBEG,REND,XEPS,F,FEAS,FCNT)

      END


c$$$C     COMMON SCALARS
c$$$      integer IC,MAXIC
c$$$      double precision VQUAD,VQUAD_A
c$$$
c$$$C     COMMON ARRAYS
c$$$      DOUBLE PRECISION XBASE_A(INN), GOPT_A(INN), HQ_A(INN**2)
c$$$
c$$$      COMMON /VQUADA/ VQUAD_A, VQUAD
c$$$      COMMON /XBASEA/ XBASE_A 
c$$$      COMMON / NOMETESTE /  GOPT_A, HQ_A    
c$$$      COMMON /CONTA1/ IC, MAXIC 
c$$$
c$$$C     LOCAL ARRAYS
c$$$      DOUBLE PRECISION FF(NPT),D(INN),Y(NPT,N),Q(1+N+N*(N+1)/2),
c$$$     1     H(NPT+N+1,NPT+N+1), 
c$$$     1     XNOVO(INN), SL(INN),
c$$$     1     SU(INN), VETOR1(NPT+N+1) 
c$$$
c$$$!     LOCAL SCALARS
c$$$      integer i,it,j,k,kn,flag
c$$$      double precision alfa,beta,c,cnorm,delta,distsq,dsq,fopt,gama,
c$$$     +     mindelta,rho,rhobeg,rhoend,sigm,sum,tau
c$$$      double precision tempofinal,tempoinicial
c$$$
c$$$      IF ( OUTPUT ) WRITE(*,3000)
c$$$
c$$$      F = 1.0D300
c$$$      IC = 0
c$$$      MAXIC = maxfcnt ! MAXIMUM NUMBER OF FUNCTION EVALUATIONS
c$$$      RHOBEG = RBEG
c$$$      RHO = RHOBEG
c$$$      DELTA = RHO
c$$$      RHOEND = REND
c$$$      GAMA = 0.1D0
c$$$
c$$$!     ---------------------------
c$$$!     Feasibility phase - Phase 0
c$$$!     ---------------------------
c$$$
c$$$      CALL SOLVER(N, XL, XU, X, M, EQUATN, LINEAR, CCODED, .true.,
c$$$     +            XEPS, CNORM, FLAG)
c$$$
c$$$      IF ( FLAG .NE. 0 ) GOTO 31
c$$$      IF (OUTPUT) WRITE(*,1000) CNORM,X
c$$$
c$$$!     ------------------------
c$$$!     End of feasibility phase
c$$$!     ------------------------
c$$$
c$$$      IF (OUTPUT) WRITE(*,1001)
c$$$
c$$$      DO I=1,N
c$$$         XNOVO(I) = X(I)
c$$$         XBASE_A(I)= X(I)            
c$$$      END DO 
c$$$      
c$$$      GO TO 5
c$$$ 4    CONTINUE     
c$$$      DO I=1, N
c$$$         X(I) = XNOVO(I)                     
c$$$         XBASE_A(I)=  XNOVO(I)                  
c$$$      END DO 
c$$$
c$$$ 5    continue
c$$$      
c$$$      CALL  PRIMEIROMODELO1 (N,X,Q,H, NPT,RHO,Y,FF,FLAG) 
c$$$
c$$$      IF ( OUTPUT ) WRITE(*,1002) RHO,DELTA,FF(1),IC,MIN(N,MAXXEL),
c$$$     +              (X(I), I=1,MIN(N,MAXXEL))
c$$$      IF ( FLAG .NE. 0 ) GOTO 31
c$$$
c$$$      FOPT = FF(1)         
c$$$
c$$$ 11   CALL   SUBPROBLEMA(N,NPT,Q,DELTA,D, X, XL, XU, DSQ,
c$$$     +                   M, EQUATN, LINEAR, CCODED, XEPS, FLAG) 
c$$$
c$$$      IF ( OUTPUT ) WRITE(*,1003) RHO,DELTA,Q(1),FOPT,IC
c$$$      IF ( FLAG .NE. 0 ) GOTO 31
c$$$      
c$$$      DISTSQ=(10.D0*RHO)**2                 
c$$$      IF (SQRT(DSQ) .LT. 0.5D0*RHO) THEN 
c$$$         
c$$$         KN=0
c$$$         DO   K=1,NPT
c$$$            SUM=0D0
c$$$            DO   J=1,N
c$$$               SUM=SUM+(Y(K,J)-X(J))**2
c$$$            END DO 
c$$$            IF (SUM .GT. DISTSQ) THEN
c$$$               KN=K
c$$$               exit
c$$$            END IF
c$$$         END DO
c$$$                    
c$$$         IF (RHO .LE. RHOEND) GO TO 31
c$$$         IF (KN .EQ. 0)   RHO = GAMA * RHO              
c$$$         
c$$$         GO TO 4
c$$$      END IF                     
c$$$      
c$$$      CALL CALFUN(N,X,F,FLAG)
c$$$
c$$$      IF ( FLAG .NE. 0 ) GOTO 31
c$$$      
c$$$      IF ((F-FOPT) .GT. (0.1D0*VQUAD)) THEN
c$$$         DELTA= 0.5D0*DELTA                 
c$$$      ELSE IF ((F-FOPT) .GE. (0.7D0*VQUAD)) THEN
c$$$         DELTA=DELTA  
c$$$      ELSE
c$$$         DELTA =   DELTA + DELTA  
c$$$      END IF      
c$$$
c$$$C     CHOOSE WHO LEAVE Y CALCULATING THE VALUE OF SIGMA. THE VARIABLE
c$$$C     'IT' IS CHOOSEN FOR DEFINE WHO LEAVE.
c$$$
c$$$      CALL SIGMA(H,N,NPT,Y,X,VETOR1,SIGM,ALFA,BETA,TAU,IT,DELTA)         
c$$$                 
c$$$C     IF ANY REDUCTION IN F, PUT X IN INTERPOLATION SET.
c$$$      IF (F .LE. FOPT) THEN  
c$$$         IF ( OUTPUT ) WRITE(*,1005) IT
c$$$         DO I=1, N            
c$$$            Y(IT,I) = X(I) 
c$$$         END DO
c$$$      ELSE
c$$$         GO TO 23
c$$$      END IF 
c$$$      
c$$$C     UPDATE H              
c$$$      CALL INVERSAH(H, N, NPT,VETOR1, SIGM, IT, ALFA, BETA,TAU)
c$$$      
c$$$      CALL ATUALIZAQ(H, N, NPT, Q, DELTA, Y, X, F, IT) 
c$$$      
c$$$ 23   IF (F  .LE.  FOPT + 0.1D0*VQUAD) THEN                 
c$$$         FOPT = F
c$$$         DO I=1, N
c$$$            XNOVO(I) = X(I) 
c$$$         END DO                      
c$$$         GO TO 11  
c$$$      END IF
c$$$      IF (sigm .le. 0d0) go to 4
c$$$      IF ((RHO .LE. RHOEND) .OR. (IC == MAXIC)) GO TO 31
c$$$      KN=0
c$$$      DO   K=1,NPT
c$$$         SUM=0D0
c$$$         DO   J=1,N
c$$$            SUM=SUM+(Y(K,J)-X(J))**2
c$$$         END DO 
c$$$         IF (SUM .GT. DISTSQ) THEN
c$$$            KN=K
c$$$            exit
c$$$         END IF
c$$$      END DO    
c$$$      
c$$$      IF (KN .GT. 0)  GOTO 4       
c$$$      
c$$$      DELTA = RHO  
c$$$      RHO = GAMA * RHO      
c$$$      
c$$$      GO TO 11    
c$$$C******************************************************************************
c$$$    
c$$$C     OUTPUT DATA
c$$$ 31   continue           
c$$$
c$$$      if ( OUTPUT .and. flag .eq.  0 ) write(*,1020)
c$$$      if ( OUTPUT .and. flag .eq. -1 ) write(*,1021)
c$$$      if ( OUTPUT .and. flag .eq.  2 ) write(*,1022)
c$$$      if ( OUTPUT .and. flag .eq.  3 ) write(*,1023) MAXIC
c$$$
c$$$      F = FOPT
c$$$
c$$$      do i = 1,n
c$$$         x(i) = XNOVO(i)
c$$$      end do
c$$$
c$$$      FEAS = 0.0D0
c$$$      do I = 1,M
c$$$         CALL CALCON(N,X,I,C)
c$$$         IF ( EQUATN(I) ) THEN
c$$$            FEAS = MAX(FEAS,ABS(C))
c$$$         ELSE
c$$$            FEAS = MAX(FEAS,MAX(0.0D0,C))
c$$$         END IF
c$$$      end do
c$$$
c$$$      FCNT = IC
c$$$
c$$$      IF ( OUTPUT ) THEN
c$$$         call cpu_time(tempofinal)
c$$$         write(*,2000) F,FEAS,RHO,DELTA,IC,(tempofinal - tempoinicial),
c$$$     +                 MIN(N,MAXXEL),(X(I), I=1,MIN(N,MAXXEL))
c$$$c$$$         print '("Time = ",1PD23.8," seconds.")',tempofinal-tempoinicial
c$$$c$$$      
c$$$c$$$         PRINT*, "NUMBER OF CALL OBJECTIVE FUNCTION =", IC 
c$$$c$$$         DO I=1, N
c$$$c$$$            PRINT*, "X(",I,")=", XNOVO(I) 
c$$$c$$$         END DO
c$$$c$$$         PRINT* , " MIN OBJECTIVE FUNCTION ="  , DMIN1(F, FOPT)  
c$$$      END IF
c$$$
c$$$C     FORMATS
c$$$
c$$$ 1000 FORMAT(/,'PHASE 0',/,7('-'),/,/,'FEASIBILITY =',36X,D23.8,/,
c$$$     +       'NEW POINT',/,3(1X,D23.8))
c$$$ 1001 FORMAT(/,'PHASE 1',/,7('-'),/)
c$$$ 1002 FORMAT(/,'(RE)BUILDING MODEL from scratch.',/,
c$$$     +       5X,'RHO =',50X,D12.5,/,
c$$$     +       5X,'Delta =',48X,D12.5,/,
c$$$     +       5X,'Objective function =',24X,D23.8,/,
c$$$     +       5X,'Function evaluations =',35X,I10,/,
c$$$     +       5X,'Current model center (first ',I3,' elements)',/,6X,
c$$$     +       3(1X,D21.8))
c$$$ 1003 FORMAT(/,'SOLVED TR SUBPROBLEM.',/,
c$$$     +       5X,'RHO =',50X,D12.5,/,
c$$$     +       5X,'Delta =',48X,D12.5,/,
c$$$     +       5X,'Model value =',31X,D23.8,/,
c$$$     +       5X,'Objective function =',24X,D23.8,/,
c$$$     +       5X,'Function evaluations =',35X,I10)
c$$$ 1004 FORMAT(5X,'Objective function =',24X,D23.8)
c$$$ 1005 FORMAT(/,'REMOVING sampling point',1X,I4,'.')
c$$$
c$$$ 1020 FORMAT(/,'Solution was found!',/)
c$$$ 1021 FORMAT(/,'Flag -1: Error while evaluating functions.',/)
c$$$ 1022 FORMAT(/,'Flag 2: Error in the internal solver.',/)
c$$$ 1023 FORMAT(/,'Flag 3: Reached the maximum of',1X,I10,1X,
c$$$     +     'function evaluations.',/)
c$$$
c$$$ 2000 FORMAT(/,'Final Iteration',/,15('-'),2/,
c$$$     +       'Objective function =',29X,D23.8,/,
c$$$     +       'Feasibility =',36X,D23.8,/,
c$$$     +       'RHO =',55X,D12.5,/,
c$$$     +       'Delta =',53X,D12.5,/,
c$$$     +       'Function evaluations =',40X,I10,/,
c$$$     +       'CPU time =',30X,1PD23.8,1X,'seconds.',/,
c$$$     +       'Solution (first ',I3,' elements)',/,3(1X,D23.8))
c$$$
c$$$ 3000 FORMAT(/,'Welcome to TRDF Algorithm!',/,
c$$$     +     'This algorithm was based on paper',/,
c$$$     +     'P.D. Conejo, E.W. Karas, and L.G. Pedroso',/,
c$$$     +     '"A trust-region derivative-free algorithm for',/,
c$$$     +     'constrained problems", to appear in Optimization',/,
c$$$     +     'Methods & Software.',/)
c$$$      END
c$$$
c$$$C     ******************************************************************
c$$$C     ******************************************************************
c$$$
c$$$C********************************  FIRST MODEL  *******************************
c$$$        SUBROUTINE  PRIMEIROMODELO1 (N,X,Q,H,NPT,DELTA,Y,FF,FLAG)
c$$$        IMPLICIT REAL*8 (A-H,O-Z)
c$$$#include "tr_params.par"
c$$$        integer flag
c$$$        DIMENSION Q(*), FF(*), x(*)
c$$$        DIMENSION GOPT_A(INN),HQ_A(INN**2),XBASE_A(INN)
c$$$        dimension  E(N+1,NPT),OMEGA(NPT,NPT),Y(NPT,N),
c$$$     1  GAMA(N+1,N+1),Z(NPT,NPT-N-1), H(NPT+N+1,NPT+N+1),
c$$$     1  YY(N), HQ(n,n),  FFTEMP(npt)                        
c$$$        INTEGER IP(npt), IQ(npt)
c$$$         COMMON / NOMETESTE /  GOPT_A, HQ_A  
c$$$         COMMON /XBASEA/ XBASE_A
c$$$C         NPT IS THE NUMBER INTERPOLATION POINTS.
c$$$C         Y IS THE INTERPOLATION SET.
c$$$C         FF KEEP IN VALUES OF F IN Y. 
c$$$C         Q STORES THE HESSIAN AND GRADIENT OF MODEL.
c$$$C         H  IS THE INVERSE ASSOCIATED WITH SYSTEM.
c$$$C          YY STORES EACH VECTOR OF Y. 
c$$$C         HQ IS THE HESSIAN IN MATRIX FORMAT.  
c$$$            
c$$$          DO I=1, 1+N+N*(N+1)/2
c$$$           Q(I)=0.0D0
c$$$          END DO ! START Q
c$$$          DO I=1, NPT+N+1
c$$$           DO J=1, NPT+N+1
c$$$            H(I,J)=0.0D0
c$$$           END DO
c$$$         END DO
c$$$          
c$$$C         START HQ
c$$$         DO I=1, N
c$$$          DO J=1, N
c$$$           HQ(I,J)=0.0D0
c$$$          END DO
c$$$         END DO
c$$$           	            
c$$$C         NPT2N = 2*N+1
c$$$         DO I=1, N
c$$$         Y(1,I)=X(I) 
c$$$         END DO         
c$$$         DO I=1, N
c$$$           DO J=1, N 
c$$$            Y(I+1,J)= X(J)
c$$$            Y(I+1,I )= X(I )+DELTA 
c$$$            Y(I+N+1,J)= X(J)
c$$$            Y(I+N+1,I)= X(I )-DELTA                 
c$$$           END DO
c$$$         END DO
c$$$        DO I=1,  2*N+1 
c$$$         DO J=1, N
c$$$              YY(J) = Y(I,J) 
c$$$         END DO            
c$$$         CALL CALFUN(N,YY,FF(I),FLAG)
c$$$
c$$$         IF ( FLAG .NE. 0 ) RETURN
c$$$        
c$$$       END DO 
c$$$                      
c$$$C******************* MODEL ***************************
c$$$         Q(1)=FF(1)          
c$$$         ! DEFINE THE GRADIENT GOPT OF THE FIRST MODEL
c$$$         DO I=1, N 
c$$$         Q(I+1)=(1D0/(2*DELTA)) * (FF(I+1)-FF(I+1+N))      
c$$$         END DO     
c$$$C           DEFINE THE DIAGONAL OF THE HESSIAN MODEL   
c$$$         DO I=1, N 
c$$$           HQ(I,I)=(1D0/(DELTA**2))*(FF(I+1)+FF(I+1+N)-2*FF(1)) 
c$$$         END DO 
c$$$	         
c$$$C           NPT >= 2N+1       
c$$$ 
c$$$        IF (NPT .GT. 2*N+1) THEN             
c$$$C        SETTING THE POITS M-2N+1
c$$$          IF (NPT .GT. 2*N+1) THEN
c$$$           DO J= 2*N+2, NPT
c$$$              IF (J .LE. 3*N+1)  IP(J) = J-2*N-1
c$$$              IF (J .GE. 3*N+2)  IP(J) = IP(J-N) 
c$$$           END DO
c$$$          END IF 
c$$$       
c$$$              II =1 
c$$$              DO J= 2*N+2, NPT 
c$$$               IF (IP(J) + II .LE. N) THEN
c$$$                   IQ(J) = IP(J) + II
c$$$               ELSE 
c$$$                  IQ(J) = IP(J) + II - N 
c$$$               END IF
c$$$               IF (MOD(IP(J)  ,N) .EQ. 0) II = II+1 
c$$$               END DO
c$$$                                  
c$$$C        OBTAIN THE POINTS Y OF 2N+1 TO NPT.
c$$$             DO I=2*N+2, NPT
c$$$              DO J= 1, N 
c$$$               Y(I,J) = Y(IP(I)+1, J) + Y(IQ(I)+1, J) - Y(1,J)                
c$$$              END DO
c$$$             END DO
c$$$        DO I=2*N+2, NPT
c$$$         DO J=1, N
c$$$          YY(J) = Y(I,J) 
c$$$         END DO            
c$$$         CALL CALFUN(N,YY,FF(I))   
c$$$
c$$$         IF ( FLAG .NE. 0 ) RETURN
c$$$
c$$$        END DO     
c$$$                      
c$$$C        DEFINE OTHERS INPUTS OF HESSIAN FOR OVER 2N+1.     
c$$$         DO J=2*N+2, NPT
c$$$          HQ(IP(J),IQ(J))=(1.0D0/(DELTA**2))*(FF(J)-FF(IP(J)+1)
c$$$     1    -FF(IQ(J)+1)+ FF(1))  
c$$$          HQ(IQ(J),IP(J)) =  HQ(IP(J),IQ(J)) 
c$$$         END DO      
c$$$       END IF  
c$$$           K=1
c$$$           DO I=1 ,N
c$$$            DO J=1, I
c$$$               Q(K+N+1) = HQ(I,J)
c$$$               K = K+1
c$$$            END DO
c$$$           END DO    
c$$$                                    
c$$$C        UPDATE THE GRADIENT AND HESSIAN FOR THE FIRST MODEL. 
c$$$	          
c$$$               DO I=1, N
c$$$                GOPT_A(I) = Q(I+1)       
c$$$               END DO              
c$$$                DO J=1, (N+1)*N/2                  
c$$$                  HQ_A(J) =  Q(1+N+J)  
c$$$                  
c$$$               END DO    
c$$$ 	          
c$$$C******************* FIRST INVERSE***************************
c$$$                             
c$$$             DO I=1, NPT
c$$$              DO J=1, NPT
c$$$               OMEGA(I,J)=0D0
c$$$              END DO
c$$$             END DO
c$$$             DO I=1, N+1
c$$$              DO J=1, N+1
c$$$              GAMA(I,J)=0D0
c$$$              END DO
c$$$             END DO
c$$$              DO I=1, NPT  
c$$$              DO J=1, NPT-N -1
c$$$              Z(I,J)=0D0
c$$$              END DO
c$$$             END DO
c$$$C      MATRIX E( N+1 X NPT)
c$$$            DO I=1, N+1
c$$$              DO J=1,NPT  
c$$$               E(I,J)=0D0
c$$$               END DO
c$$$            END DO
c$$$                E(1,1)=1D0     
c$$$            DO I=2, N+1
c$$$               E(I,I)= 1.0D0/(2.0D0*DELTA)
c$$$               E(I,I+N)= - 1.0D0/(2.0D0*DELTA)
c$$$            END DO
c$$$C      MATRIX Z(NPT X NPT-N-1)            
c$$$            DO I=1, N 
c$$$              Z(1,I)= -SQRT(2.0D0)/(DELTA**2)
c$$$              Z(I+1,I)=  SQRT(2.0D0)/(2.0D0*DELTA**2)
c$$$              Z(N+I+1,I)= SQRT(2.0D0)/(2.0D0*DELTA**2)
c$$$            END DO  
c$$$           
c$$$C       THE NEW INVERSE FOR MORE OF 2N+1 POINTS IN Y
c$$$            IF (NPT .GT.  2*N+1) THEN             
c$$$               DO I=N+1, NPT-N-1 
c$$$                Z(1,I)= 1.0D0/(DELTA**2)
c$$$                Z(N+I+1, I) = Z(1,I) 
c$$$                Z(IP(N+I+1)+1,I) = -1.0D0/(DELTA**2)
c$$$                Z(IQ(N+I+1)+1,I) = -1.0D0/(DELTA**2)        
c$$$               END DO
c$$$            END IF 
c$$$         
c$$$C         MULTIPLYING ZZ^T FOR DETERMINE OMEGA          
c$$$         ACUMULADOR=0D0
c$$$         DO I=1, NPT   
c$$$            DO K=1, NPT              
c$$$               DO J=1, NPT-N-1
c$$$                ACUMULADOR =  ACUMULADOR + Z(I,J)*Z(K,J)
c$$$               END DO
c$$$               OMEGA(I,K) = ACUMULADOR               
c$$$               ACUMULADOR = 0D0              
c$$$             END DO
c$$$         END DO
c$$$             
c$$$C        THE MATRIX INVERSE H     
c$$$        DO I=1, NPT
c$$$          DO J=1, NPT
c$$$           H(I,J)=OMEGA(I,J)         
c$$$          END DO
c$$$        END DO
c$$$                            
c$$$C       THE N+1 LINES OF H                 
c$$$        DO I=NPT+1, NPT+N+1
c$$$         DO J= 1, NPT 
c$$$          H(I,J) = E(I-NPT,J)         
c$$$         END DO
c$$$        END DO
c$$$C       THE N+1 COLUMNS OF H 
c$$$                                       
c$$$        DO I=1, NPT
c$$$         DO J= NPT+1, NPT+N+1
c$$$          H(I,J) = H(J,I)        
c$$$         END DO
c$$$        END DO  
c$$$              
c$$$        RETURN
c$$$       END SUBROUTINE PRIMEIROMODELO1 
c$$$
c$$$C     ******************************************************************
c$$$C     ******************************************************************
c$$$
c$$$       SUBROUTINE SIGMA(H,N,NPT,Y,X,VETOR1,SIGM,ALFA,BETA,TAU,IT,DELTA)
c$$$         IMPLICIT REAL*8 (A-H,O-Z)
c$$$#include "tr_params.par"
c$$$         DIMENSION X(*), VETOR1(*) 
c$$$          DIMENSION   WW(NPT+N+1,1),AUXILIAR(4)           
c$$$          DIMENSION H(NPT+N+1,NPT+N+1), Y(NPT,N), XBASE_A(INN)
c$$$          COMMON /XBASEA/ XBASE_A
c$$$C        WW STORAGE THE VETOR1 IN W TO PRODUCE ALFA BETA TAU HOW IN DEFINITION.   
c$$$C        SIGMA = ALFA BETA + TAU**2. ALFA = ET^T H ET, NAMELY, HTT (T = IT)
c$$$C        IT  INDICATE THAT THE VETOR1 IN POSITION TWO (SECOND LINE OF Y) LEAVE OF Y
c$$$C        CHOOSE THE LARGEST SIGMA (AND THEREFORE IT) UNLESS DE CURRENT POINT (IT CURRENT)
c$$$          
c$$$         ITT=IT    ! STORAGE THE POSITION OF BEST ITERATING YET.          
c$$$         IT =1
c$$$         SIGMI = -1.0D100
c$$$         WW(NPT+1,1) = 1.0D0
c$$$         DO I=1, N
c$$$         WW(I+NPT+1, 1) =  X(I)- XBASE_A(I)
c$$$         END DO 
c$$$         DO I=1, NPT
c$$$          CONT = 0.0D0
c$$$         DO K=1, N         
c$$$         CONT = CONT +  (Y(I,K) - XBASE_A(K)) * (X(K) - XBASE_A(K))          
c$$$         END DO
c$$$         WW(I, 1) = 0.5D0*CONT**2
c$$$         END DO  
c$$$           
c$$$C       CALCULATING ALL SIGMA AND CHOOSE IT FOR THE GREATER SIGMA. MULTIPLY T-th LINE OF H FOR W FOR OBTAIN TAU. 
c$$$          DO KKK = 1, NPT-1            
c$$$             CONT=0.0D0
c$$$             DO I=1, NPT+ N +1
c$$$              CONT = CONT +  H(IT, I) * WW(I,1)
c$$$             END DO
c$$$             TAU = CONT
c$$$C      CALCULUS OF BETA = 0.5||X^+-XBASE||-WHW                  
c$$$            DO I=1, NPT+N+1
c$$$             CONT = 0.0D0 
c$$$             DO K=1, NPT+N+1
c$$$              CONT = CONT + H(I,K) * WW(K,1)
c$$$              VETOR1(I) = CONT                  
c$$$             END DO             
c$$$            END DO
c$$$               
c$$$            CONT = 0.0D0
c$$$            DO I=1,  NPT+N+1
c$$$             CONT = CONT+WW(I,1)*VETOR1(I)
c$$$             AGUARD = CONT
c$$$            END DO
c$$$                       
c$$$C       CALCULUS  OF X-XB^4
c$$$          CONT = 0.0D0
c$$$         DO I=1, N
c$$$          CONT = CONT + (X(I)-XBASE_A(I))**2
c$$$         END DO         
c$$$            BETA =   0.5D0 * CONT**2 - AGUARD
c$$$            ALFA = H(IT,IT) 
c$$$            SIGM = ALFA * BETA + TAU**2  
c$$$            IF (SIGM .GE. SIGMI .AND. IT .NE. ITT) THEN                      
c$$$                 SIGMI = SIGM
c$$$                 IAUXILIAR  =  IT                
c$$$		 AUXILIAR(1) = ALFA
c$$$ 		 AUXILIAR(2) = BETA
c$$$		 AUXILIAR(3) = TAU
c$$$                 AUXILIAR(4) = SIGM                           
c$$$             END IF 
c$$$           IT = IT + 1
c$$$       END DO
c$$$          
c$$$CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
c$$$                    IT   = IAUXILIAR  
c$$$		    ALFA = AUXILIAR(1) 
c$$$ 		    BETA = AUXILIAR(2)
c$$$		    TAU  = AUXILIAR(3) 
c$$$                    SIGM = AUXILIAR(4)
c$$$         
c$$$           CONT=0.0D0
c$$$           DO I=1, NPT+ N +1
c$$$            CONT = CONT +  H(IT, I) * WW(I,1)
c$$$           END DO
c$$$           TAU = CONT
c$$$          
c$$$           RETURN
c$$$         END SUBROUTINE SIGMA 
c$$$
c$$$C     ******************************************************************
c$$$C     ******************************************************************
c$$$
c$$$C **************** UPDAT THE INVERSE H ******************************
c$$$          SUBROUTINE INVERSAH(H, N, NPT,VETOR1,SIGM,IT,ALFA,BETA,TAU)
c$$$          IMPLICIT REAL*8 (A-H,O-Z)
c$$$#include "tr_params.par"
c$$$           DIMENSION VETOR1(*)
c$$$           DIMENSION P1(NPT+N+1,NPT+N+1) ,H(NPT+N+1,NPT+N+1)
c$$$           DIMENSION P2(NPT+N+1,NPT+N+1),P3(NPT+N+1,NPT+N+1)
c$$$C         ALFA*(E-MM) * (E-MM)'- BETA* H * E * E'*H+TAU*H*E*(E-MM)'+ (E-MM)* E'* H
c$$$C         MM = H*WW THAT IS STORED IN VETOR1 
c$$$           VETOR1(IT) = VETOR1(IT)-1.0D0              
c$$$          DO I=1, N+NPT+1
c$$$            DO J=1, N+NPT+1
c$$$              P1(I,J)= VETOR1(I) * VETOR1(J)
c$$$              P2(I,J)= H(I, IT) * H(J, IT)
c$$$              P3(I,J)= (H(IT, I)*(-VETOR1(J)))+(H(IT,J)*(-VETOR1(I)))
c$$$            END DO
c$$$          END DO
c$$$            DO I=1, N+NPT+1
c$$$             DO J=1, N+NPT+1
c$$$            if (sigm .eq. 0d0) return
c$$$         H(I,J)=H(I,J)+(1/SIGM)*(ALFA*P1(I,J)-BETA*P2(I,J)+TAU*P3(I,J)) 
c$$$             END DO
c$$$            END DO
c$$$          RETURN
c$$$          END SUBROUTINE INVERSAH 
c$$$ 
c$$$C     ******************************************************************
c$$$C     ******************************************************************
c$$$
c$$$       SUBROUTINE  SUBPROBLEMA(N, NPT, Q, DELTA, D, X, XL, XU, DSQ,
c$$$     +                         M, EQUATN, LINEAR, CCODED, XEPS, FLAG)
c$$$       IMPLICIT REAL*8 (A-H,O-Z)
c$$$#include "tr_params.par"
c$$$       DIMENSION Q(*), X(*), XL(*),XU(*)
c$$$       DOUBLE PRECISION XBASE_A(INN) 
c$$$       DOUBLE PRECISION VQUAD_A, GRADD , VQUAD , XANTIGO(INN),L(N)
c$$$       DOUBLE PRECISION H(NPT+N+1,NPT+N+1),XOPT_A(INN),D(INN),U(N) 
c$$$       INTEGER N
c$$$
c$$$!     SCALAR ARGUMENTS
c$$$       integer flag,m
c$$$
c$$$!     ARRAY ARGUMENTS
c$$$       logical ccoded(2),equatn(m), linear(m)
c$$$
c$$$       COMMON /XBASEA/ XBASE_A           
c$$$       COMMON /VQUADA/VQUAD_A, VQUAD   
c$$$
c$$$!     LOCAL SCALARS
c$$$       double precision cnorm
c$$$         
c$$$         DO I = 1,N
c$$$          L(I)=  DMAX1(XL(I) - XBASE_A(I),X(I) - XBASE_A(I)-DELTA) 
c$$$          U(I)=  DMIN1(XU(I) - XBASE_A(I),X(I) - XBASE_A(I)+DELTA)           
c$$$          XANTIGO(I) = X(I)       
c$$$         END DO
c$$$            
c$$$         DO I=1, N        
c$$$          D(I)=X(I)-XBASE_A(I)         
c$$$      
c$$$         END DO   
c$$$
c$$$         CALL MEVALF(N,D,F,FLAG)
c$$$
c$$$         IF ( FLAG .NE. 0 ) RETURN
c$$$
c$$$         VQUAD=  F + Q(1)                                 
c$$$         
c$$$         CALL SOLVER(N, L, U, D, M, EQUATN, LINEAR, CCODED, .false.,
c$$$     +               XEPS, CNORM, FLAG)
c$$$         
c$$$         IF ( FLAG .NE. 0 ) RETURN
c$$$                   
c$$$C CALCULUS THE STEP LENGTH 
c$$$           SUM = 0.0D0
c$$$           DO I=1, N
c$$$            SUM = SUM + (X(I)-(D(I)+XBASE_A(I)))**2
c$$$            END DO
c$$$            DSQ = SUM
c$$$          DO I=1, N
c$$$           X(I)= D(I) +  XBASE_A(I)   
c$$$          END DO   
c$$$
c$$$           CALL MEVALF(N,D,F,IFLAG)
c$$$
c$$$           IF ( FLAG .NE. 0 ) RETURN
c$$$
c$$$           VQUAD_A=  F + Q(1) !MODEL IN XNOVO 
c$$$           VQUAD=  -VQUAD + VQUAD_A  
c$$$         IF (VQUAD .GE. 0.D0 .or. cnorm .gt. xeps)  THEN
c$$$          DO I=1, N
c$$$           X(I) = XANTIGO(I) 
c$$$          END DO
c$$$          DSQ = 0D0
c$$$         END IF
c$$$
c$$$         RETURN
c$$$        END
c$$$
c$$$C     ******************************************************************
c$$$C     ****************************************************************** 
c$$$
c$$$       SUBROUTINE  ATUALIZAQ(H, N, NPT, Q, DELTA, Y, X, F, IT)
c$$$       IMPLICIT REAL*8 (A-H,O-Z)
c$$$#include "tr_params.par"
c$$$       DOUBLE PRECISION Q(*), X(*) 
c$$$     
c$$$       DIMENSION H(NPT+N+1,NPT+N+1),Y(NPT,N),VETORAUX(INN),DD(N,N)
c$$$       DOUBLE PRECISION  GOPT_A(INN), HQ_A(INN**2), XBASE_A(INN)  
c$$$       DIMENSION QQ((N+1)*(N+2)/2),  TEMP(1+N+NPT)     
c$$$             COMMON /VQUADA/ VQUAD_A, VQUAD 
c$$$             COMMON /XBASEA/ XBASE_A  
c$$$             COMMON / NOMETESTE /  GOPT_A, HQ_A    
c$$$               
c$$$              DO I=1, 1+N+NPT
c$$$               TEMP(I) =  (F - VQUAD_A)* H(I, IT) ! IS LAMBDA 
c$$$              END DO                     
c$$$               DO I=1, N 
c$$$                DO J=1, N
c$$$                 DD(I,J)=0.0D0
c$$$                END DO
c$$$               END DO
c$$$C             M=M+     LAMBCG(J)*((Y(:,J)-XB') * ( Y(: ,J)-XB')')  
c$$$             DO I=1, NPT
c$$$               DO JJ=1, N
c$$$                VETORAUX(JJ) =  Y(I,JJ)-XBASE_A(JJ)  
c$$$               END DO                      
c$$$                  DO K=1, N                    
c$$$                    DO J=1, N
c$$$                  DD(K,J) = DD(K,J)+ TEMP(I)* VETORAUX(K) * VETORAUX(J)  
c$$$C      DDEH IS THE GRADIENT OF QUADRATIC D                               
c$$$                    END DO                   
c$$$                  END DO                  
c$$$             END DO
c$$$C  N+1 FIRST ELEMENTS OF THE QQ (PARAMETER OF THE NEW MODEL)
c$$$            QQ(1) = TEMP(NPT+1)
c$$$            DO I=2, N+1 
c$$$             QQ(I) = TEMP(NPT+I)
c$$$            END DO   
c$$$C PUT IN THE QQ FOR SYMMETRY OF DD
c$$$              II=1
c$$$               J=1
c$$$C  DO WHILE (II .LE.  (N+1)*N/2 ) 
c$$$              DO WHILE (II .LE. N) 
c$$$              DO I=1, II 
c$$$               QQ(J+1+N) = DD(II,I) 
c$$$                J=J+1
c$$$               END DO
c$$$                II = II + 1
c$$$               END  DO   
c$$$C ADD Q TO QQ.  Q IS THE OLD MODEL, QQ IS THE MODEL D, AND Q + D = Q+    
c$$$              DO I=1, 1+N+ (N+1)*N/2
c$$$               Q(I) = Q(I) + QQ(I)
c$$$              END DO
c$$$C  UPDAT THE MODEL IN THE FILE INIP FOR COMMON FUNCTION 
c$$$               DO I=1, N      
c$$$               GOPT_A(I) = Q(I+1)
c$$$               END DO   
c$$$               DO J=1, (N+1)*N/2 
c$$$                  HQ_A(J) =  Q(1+N+J) 
c$$$               END DO  
c$$$           RETURN
c$$$         END SUBROUTINE  ATUALIZAQ  
c$$$
c$$$C     ******************************************************************
c$$$C     ******************************************************************
c$$$
c$$$      subroutine mvv(v1,v2,n, gradd)
c$$$      
c$$$      implicit none
c$$$
c$$$#include "tr_params.par"
c$$$
c$$$!     multiplica vetor por vetor
c$$$      double precision v1(INN), v2(INN), soma, gradd
c$$$      integer j, n
c$$$      
c$$$      soma=0
c$$$      do j=1 , n
c$$$         soma=soma +  v1(j)*v2(j)
c$$$         gradd=soma
c$$$      end do
c$$$      return
c$$$      end subroutine mvv
c$$$
c$$$C     ******************************************************************
c$$$C     ******************************************************************
c$$$
c$$$      subroutine  mmv(HQ, S,n, v)
c$$$
c$$$      implicit none
c$$$
c$$$#include "tr_params.par"
c$$$
c$$$!     multiplica matriz simetrica (dada como vetor) por vetor
c$$$!     estava com hs, e troquei para hss para nao atualizar hs desnec
c$$$      
c$$$      double precision S(INN), HQ(INN ** 2), HSS(INN), v(INN)
c$$$      
c$$$      integer n, i, j , IH
c$$$      IH=0
c$$$      DO J=1,N
c$$$         HSS(J)= 0
c$$$         DO I=1,J
c$$$            IH=IH+1
c$$$            IF (I .LT. J) HSS(J)=HSS(J)+HQ(IH)*S(I)
c$$$            HSS(I)=HSS(I)+HQ(IH)*S(J)
c$$$            v(I)=HSS(I)
c$$$         end DO
c$$$      end DO
c$$$      return
c$$$      end subroutine mmv
c$$$

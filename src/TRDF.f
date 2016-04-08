      SUBROUTINE EASYTRDF(N,X,XL,XU,M,EQUATN,LINEAR,CCODED,F,FEAS,
     +                    FCNT)

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

      CALL TRDF(N,NPT,X,XL,XU,M,EQUATN,LINEAR,CCODED,MAXFCNT,RBEG,REND,
     +          XEPS,F,FEAS,FCNT)

      END

C     ******************************************************************
C     ******************************************************************

      SUBROUTINE TRDF(N,NPT,X,XL,XU,M,EQUATN,LINEAR,CCODED,MAXFCNT,RBEG,
     +                REND,XEPS,F,FEAS,FCNT)     

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

      IMPLICIT NONE

#include "tr_params.par"

C     SCALAR ARGUMENTS
      integer m,maxfcnt,N,NPT,FCNT
      double precision F,FEAS,RBEG,REND,XEPS

C     ARRAY ARGUMENTS
      DOUBLE PRECISION  X(N),XL(N),XU(N)
      logical ccoded(2),equatn(m),linear(m)

C     COMMON SCALARS
      integer IC,MAXIC
      double precision VQUAD,VQUAD_A

C     COMMON ARRAYS
      DOUBLE PRECISION XBASE_A(INN), GOPT_A(INN), HQ_A(INN**2)

      COMMON /VQUADA/ VQUAD_A, VQUAD
      COMMON /XBASEA/ XBASE_A 
      COMMON / NOMETESTE /  GOPT_A, HQ_A    
      COMMON /CONTA1/ IC, MAXIC 

C     LOCAL ARRAYS
      DOUBLE PRECISION FF(NPT),D(INN),Y(NPT,N),Q(1+N+N*(N+1)/2),
     1     H(NPT+N+1,NPT+N+1), 
     1     XNOVO(INN), SL(INN),
     1     SU(INN), VETOR1(NPT+N+1) 

!     LOCAL SCALARS
      integer i,it,j,k,kn,flag
      double precision alfa,beta,c,cnorm,delta,distsq,dsq,fopt,gama,
     +     mindelta,rho,rhobeg,rhoend,sigm,sum,tau
      double precision tempofinal,tempoinicial

      IF ( OUTPUT ) WRITE(*,3000)

      F = 1.0D300
      IC = 0
      MAXIC = maxfcnt ! MAXIMUM NUMBER OF FUNCTION EVALUATIONS
      RHOBEG = RBEG
      RHO = RHOBEG
      DELTA = RHO
      RHOEND = REND
      GAMA = 0.1D0

!     ---------------------------
!     Feasibility phase - Phase 0
!     ---------------------------

      CALL SOLVER(N, XL, XU, X, M, EQUATN, LINEAR, CCODED, .true.,
     +            XEPS, CNORM, FLAG)

      IF ( FLAG .NE. 0 ) GOTO 31
      IF (OUTPUT) WRITE(*,1000) CNORM,X

!     ------------------------
!     End of feasibility phase
!     ------------------------

      IF (OUTPUT) WRITE(*,1001)

      DO I=1,N
         XNOVO(I) = X(I)
         XBASE_A(I)= X(I)            
      END DO 
      
      GO TO 5
 4    CONTINUE     
      DO I=1, N
         X(I) = XNOVO(I)                     
         XBASE_A(I)=  XNOVO(I)                  
      END DO 

 5    continue
      
      CALL  PRIMEIROMODELO1 (N,X,Q,H, NPT,RHO,Y,FF,FLAG) 

      IF ( OUTPUT ) WRITE(*,1002) RHO,DELTA,FF(1),IC,MIN(N,MAXXEL),
     +              (X(I), I=1,MIN(N,MAXXEL))
      IF ( FLAG .NE. 0 ) GOTO 31

      FOPT = FF(1)         

 11   CALL   SUBPROBLEMA(N,NPT,Q,DELTA,D, X, XL, XU, DSQ,
     +                   M, EQUATN, LINEAR, CCODED, XEPS, FLAG) 

      IF ( OUTPUT ) WRITE(*,1003) RHO,DELTA,Q(1),FOPT,IC
      IF ( FLAG .NE. 0 ) GOTO 31
      
      DISTSQ=(10.D0*RHO)**2                 
      IF (SQRT(DSQ) .LT. 0.5D0*RHO) THEN 
         
         KN=0
         DO   K=1,NPT
            SUM=0D0
            DO   J=1,N
               SUM=SUM+(Y(K,J)-X(J))**2
            END DO 
            IF (SUM .GT. DISTSQ) THEN
               KN=K
               exit
            END IF
         END DO
                    
         IF (RHO .LE. RHOEND) GO TO 31
         IF (KN .EQ. 0)   RHO = GAMA * RHO              
         
         GO TO 4
      END IF                     
      
      CALL CALFUN(N,X,F,FLAG)

      IF ( FLAG .NE. 0 ) GOTO 31
      
      IF ((F-FOPT) .GT. (0.1D0*VQUAD)) THEN
         DELTA= 0.5D0*DELTA                 
      ELSE IF ((F-FOPT) .GE. (0.7D0*VQUAD)) THEN
         DELTA=DELTA  
      ELSE
         DELTA =   DELTA + DELTA  
      END IF      

C     CHOOSE WHO LEAVE Y CALCULATING THE VALUE OF SIGMA. THE VARIABLE
C     'IT' IS CHOOSEN FOR DEFINE WHO LEAVE.

      CALL SIGMA(H,N,NPT,Y,X,VETOR1,SIGM,ALFA,BETA,TAU,IT,DELTA)         
                 
C     IF ANY REDUCTION IN F, PUT X IN INTERPOLATION SET.
      IF (F .LE. FOPT) THEN  
         IF ( OUTPUT ) WRITE(*,1005) IT
         DO I=1, N            
            Y(IT,I) = X(I) 
         END DO
      ELSE
         GO TO 23
      END IF 
      
C     UPDATE H              
      CALL INVERSAH(H, N, NPT,VETOR1, SIGM, IT, ALFA, BETA,TAU)
      
      CALL ATUALIZAQ(H, N, NPT, Q, DELTA, Y, X, F, IT) 
      
 23   IF (F  .LE.  FOPT + 0.1D0*VQUAD) THEN                 
         FOPT = F
         DO I=1, N
            XNOVO(I) = X(I) 
         END DO                      
         GO TO 11  
      END IF
      IF (sigm .le. 0d0) go to 4
      IF ((RHO .LE. RHOEND) .OR. (IC == MAXIC)) GO TO 31
      KN=0
      DO   K=1,NPT
         SUM=0D0
         DO   J=1,N
            SUM=SUM+(Y(K,J)-X(J))**2
         END DO 
         IF (SUM .GT. DISTSQ) THEN
            KN=K
            exit
         END IF
      END DO    
      
      IF (KN .GT. 0)  GOTO 4       
      
      DELTA = RHO  
      RHO = GAMA * RHO      
      
      GO TO 11    
C******************************************************************************
    
C     OUTPUT DATA
 31   continue           

      if ( OUTPUT .and. flag .eq.  0 ) write(*,1020)
      if ( OUTPUT .and. flag .eq. -1 ) write(*,1021)
      if ( OUTPUT .and. flag .eq.  2 ) write(*,1022)
      if ( OUTPUT .and. flag .eq.  3 ) write(*,1023) MAXIC

      F = FOPT

      do i = 1,n
         x(i) = XNOVO(i)
      end do

      FEAS = 0.0D0
      do I = 1,M
         CALL CALCON(N,X,I,C)
         IF ( EQUATN(I) ) THEN
            FEAS = MAX(FEAS,ABS(C))
         ELSE
            FEAS = MAX(FEAS,MAX(0.0D0,C))
         END IF
      end do

      FCNT = IC

      IF ( OUTPUT ) THEN
         call cpu_time(tempofinal)
         write(*,2000) F,FEAS,RHO,DELTA,IC,(tempofinal - tempoinicial),
     +                 MIN(N,MAXXEL),(X(I), I=1,MIN(N,MAXXEL))
c$$$         print '("Time = ",1PD23.8," seconds.")',tempofinal-tempoinicial
c$$$      
c$$$         PRINT*, "NUMBER OF CALL OBJECTIVE FUNCTION =", IC 
c$$$         DO I=1, N
c$$$            PRINT*, "X(",I,")=", XNOVO(I) 
c$$$         END DO
c$$$         PRINT* , " MIN OBJECTIVE FUNCTION ="  , DMIN1(F, FOPT)  
      END IF

C     FORMATS

 1000 FORMAT(/,'PHASE 0',/,7('-'),/,/,'FEASIBILITY =',36X,D23.8,/,
     +       'NEW POINT',/,3(1X,D23.8))
 1001 FORMAT(/,'PHASE 1',/,7('-'),/)
 1002 FORMAT(/,'(RE)BUILDING MODEL from scratch.',/,
     +       5X,'RHO =',50X,D12.5,/,
     +       5X,'Delta =',48X,D12.5,/,
     +       5X,'Objective function =',24X,D23.8,/,
     +       5X,'Function evaluations =',35X,I10,/,
     +       5X,'Current model center (first ',I3,' elements)',/,6X,
     +       3(1X,D21.8))
 1003 FORMAT(/,'SOLVED TR SUBPROBLEM.',/,
     +       5X,'RHO =',50X,D12.5,/,
     +       5X,'Delta =',48X,D12.5,/,
     +       5X,'Model value =',31X,D23.8,/,
     +       5X,'Objective function =',24X,D23.8,/,
     +       5X,'Function evaluations =',35X,I10)
 1004 FORMAT(5X,'Objective function =',24X,D23.8)
 1005 FORMAT(/,'REMOVING sampling point',1X,I4,'.')

 1020 FORMAT(/,'Solution was found!',/)
 1021 FORMAT(/,'Flag -1: Error while evaluating functions.',/)
 1022 FORMAT(/,'Flag 2: Error in the internal solver.',/)
 1023 FORMAT(/,'Flag 3: Reached the maximum of',1X,I10,1X,
     +     'function evaluations.',/)

 2000 FORMAT(/,'Final Iteration',/,15('-'),2/,
     +       'Objective function =',29X,D23.8,/,
     +       'Feasibility =',36X,D23.8,/,
     +       'RHO =',55X,D12.5,/,
     +       'Delta =',53X,D12.5,/,
     +       'Function evaluations =',40X,I10,/,
     +       'CPU time =',30X,1PD23.8,1X,'seconds.',/,
     +       'Solution (first ',I3,' elements)',/,3(1X,D23.8))

 3000 FORMAT(/,'Welcome to TRDF Algorithm!',/,
     +     'This algorithm was based on paper',/,
     +     'P.D. Conejo, E.W. Karas, and L.G. Pedroso',/,
     +     '"A trust-region derivative-free algorithm for',/,
     +     'constrained problems", to appear in Optimization',/,
     +     'Methods & Software.',/)
      END

C     ******************************************************************
C     ******************************************************************

C********************************  FIRST MODEL  *******************************
        SUBROUTINE  PRIMEIROMODELO1 (N,X,Q,H,NPT,DELTA,Y,FF,FLAG)
        IMPLICIT REAL*8 (A-H,O-Z)
#include "tr_params.par"
        integer flag
        DIMENSION Q(*), FF(*), x(*)
        DIMENSION GOPT_A(INN),HQ_A(INN**2),XBASE_A(INN)
        dimension  E(N+1,NPT),OMEGA(NPT,NPT),Y(NPT,N),
     1  GAMA(N+1,N+1),Z(NPT,NPT-N-1), H(NPT+N+1,NPT+N+1),
     1  YY(N), HQ(n,n),  FFTEMP(npt)                        
        INTEGER IP(npt), IQ(npt)
         COMMON / NOMETESTE /  GOPT_A, HQ_A  
         COMMON /XBASEA/ XBASE_A
C         NPT IS THE NUMBER INTERPOLATION POINTS.
C         Y IS THE INTERPOLATION SET.
C         FF KEEP IN VALUES OF F IN Y. 
C         Q STORES THE HESSIAN AND GRADIENT OF MODEL.
C         H  IS THE INVERSE ASSOCIATED WITH SYSTEM.
C          YY STORES EACH VECTOR OF Y. 
C         HQ IS THE HESSIAN IN MATRIX FORMAT.  
            
          DO I=1, 1+N+N*(N+1)/2
           Q(I)=0.0D0
          END DO ! START Q
          DO I=1, NPT+N+1
           DO J=1, NPT+N+1
            H(I,J)=0.0D0
           END DO
         END DO
          
C         START HQ
         DO I=1, N
          DO J=1, N
           HQ(I,J)=0.0D0
          END DO
         END DO
           	            
C         NPT2N = 2*N+1
         DO I=1, N
         Y(1,I)=X(I) 
         END DO         
         DO I=1, N
           DO J=1, N 
            Y(I+1,J)= X(J)
            Y(I+1,I )= X(I )+DELTA 
            Y(I+N+1,J)= X(J)
            Y(I+N+1,I)= X(I )-DELTA                 
           END DO
         END DO
        DO I=1,  2*N+1 
         DO J=1, N
              YY(J) = Y(I,J) 
         END DO            
         CALL CALFUN(N,YY,FF(I),FLAG)

         IF ( FLAG .NE. 0 ) RETURN
        
       END DO 
                      
C******************* MODEL ***************************
         Q(1)=FF(1)          
         ! DEFINE THE GRADIENT GOPT OF THE FIRST MODEL
         DO I=1, N 
         Q(I+1)=(1D0/(2*DELTA)) * (FF(I+1)-FF(I+1+N))      
         END DO     
C           DEFINE THE DIAGONAL OF THE HESSIAN MODEL   
         DO I=1, N 
           HQ(I,I)=(1D0/(DELTA**2))*(FF(I+1)+FF(I+1+N)-2*FF(1)) 
         END DO 
	         
C           NPT >= 2N+1       
 
        IF (NPT .GT. 2*N+1) THEN             
C        SETTING THE POITS M-2N+1
          IF (NPT .GT. 2*N+1) THEN
           DO J= 2*N+2, NPT
              IF (J .LE. 3*N+1)  IP(J) = J-2*N-1
              IF (J .GE. 3*N+2)  IP(J) = IP(J-N) 
           END DO
          END IF 
       
              II =1 
              DO J= 2*N+2, NPT 
               IF (IP(J) + II .LE. N) THEN
                   IQ(J) = IP(J) + II
               ELSE 
                  IQ(J) = IP(J) + II - N 
               END IF
               IF (MOD(IP(J)  ,N) .EQ. 0) II = II+1 
               END DO
                                  
C        OBTAIN THE POINTS Y OF 2N+1 TO NPT.
             DO I=2*N+2, NPT
              DO J= 1, N 
               Y(I,J) = Y(IP(I)+1, J) + Y(IQ(I)+1, J) - Y(1,J)                
              END DO
             END DO
        DO I=2*N+2, NPT
         DO J=1, N
          YY(J) = Y(I,J) 
         END DO            
         CALL CALFUN(N,YY,FF(I))   

         IF ( FLAG .NE. 0 ) RETURN

        END DO     
                      
C        DEFINE OTHERS INPUTS OF HESSIAN FOR OVER 2N+1.     
         DO J=2*N+2, NPT
          HQ(IP(J),IQ(J))=(1.0D0/(DELTA**2))*(FF(J)-FF(IP(J)+1)
     1    -FF(IQ(J)+1)+ FF(1))  
          HQ(IQ(J),IP(J)) =  HQ(IP(J),IQ(J)) 
         END DO      
       END IF  
           K=1
           DO I=1 ,N
            DO J=1, I
               Q(K+N+1) = HQ(I,J)
               K = K+1
            END DO
           END DO    
                                    
C        UPDATE THE GRADIENT AND HESSIAN FOR THE FIRST MODEL. 
	          
               DO I=1, N
                GOPT_A(I) = Q(I+1)       
               END DO              
                DO J=1, (N+1)*N/2                  
                  HQ_A(J) =  Q(1+N+J)  
                  
               END DO    
 	          
C******************* FIRST INVERSE***************************
                             
             DO I=1, NPT
              DO J=1, NPT
               OMEGA(I,J)=0D0
              END DO
             END DO
             DO I=1, N+1
              DO J=1, N+1
              GAMA(I,J)=0D0
              END DO
             END DO
              DO I=1, NPT  
              DO J=1, NPT-N -1
              Z(I,J)=0D0
              END DO
             END DO
C      MATRIX E( N+1 X NPT)
            DO I=1, N+1
              DO J=1,NPT  
               E(I,J)=0D0
               END DO
            END DO
                E(1,1)=1D0     
            DO I=2, N+1
               E(I,I)= 1.0D0/(2.0D0*DELTA)
               E(I,I+N)= - 1.0D0/(2.0D0*DELTA)
            END DO
C      MATRIX Z(NPT X NPT-N-1)            
            DO I=1, N 
              Z(1,I)= -SQRT(2.0D0)/(DELTA**2)
              Z(I+1,I)=  SQRT(2.0D0)/(2.0D0*DELTA**2)
              Z(N+I+1,I)= SQRT(2.0D0)/(2.0D0*DELTA**2)
            END DO  
           
C       THE NEW INVERSE FOR MORE OF 2N+1 POINTS IN Y
            IF (NPT .GT.  2*N+1) THEN             
               DO I=N+1, NPT-N-1 
                Z(1,I)= 1.0D0/(DELTA**2)
                Z(N+I+1, I) = Z(1,I) 
                Z(IP(N+I+1)+1,I) = -1.0D0/(DELTA**2)
                Z(IQ(N+I+1)+1,I) = -1.0D0/(DELTA**2)        
               END DO
            END IF 
         
C         MULTIPLYING ZZ^T FOR DETERMINE OMEGA          
         ACUMULADOR=0D0
         DO I=1, NPT   
            DO K=1, NPT              
               DO J=1, NPT-N-1
                ACUMULADOR =  ACUMULADOR + Z(I,J)*Z(K,J)
               END DO
               OMEGA(I,K) = ACUMULADOR               
               ACUMULADOR = 0D0              
             END DO
         END DO
             
C        THE MATRIX INVERSE H     
        DO I=1, NPT
          DO J=1, NPT
           H(I,J)=OMEGA(I,J)         
          END DO
        END DO
                            
C       THE N+1 LINES OF H                 
        DO I=NPT+1, NPT+N+1
         DO J= 1, NPT 
          H(I,J) = E(I-NPT,J)         
         END DO
        END DO
C       THE N+1 COLUMNS OF H 
                                       
        DO I=1, NPT
         DO J= NPT+1, NPT+N+1
          H(I,J) = H(J,I)        
         END DO
        END DO  
              
        RETURN
       END SUBROUTINE PRIMEIROMODELO1 

C     ******************************************************************
C     ******************************************************************

       SUBROUTINE SIGMA(H,N,NPT,Y,X,VETOR1,SIGM,ALFA,BETA,TAU,IT,DELTA)
         IMPLICIT REAL*8 (A-H,O-Z)
#include "tr_params.par"
         DIMENSION X(*), VETOR1(*) 
          DIMENSION   WW(NPT+N+1,1),AUXILIAR(4)           
          DIMENSION H(NPT+N+1,NPT+N+1), Y(NPT,N), XBASE_A(INN)
          COMMON /XBASEA/ XBASE_A
C        WW STORAGE THE VETOR1 IN W TO PRODUCE ALFA BETA TAU HOW IN DEFINITION.   
C        SIGMA = ALFA BETA + TAU**2. ALFA = ET^T H ET, NAMELY, HTT (T = IT)
C        IT  INDICATE THAT THE VETOR1 IN POSITION TWO (SECOND LINE OF Y) LEAVE OF Y
C        CHOOSE THE LARGEST SIGMA (AND THEREFORE IT) UNLESS DE CURRENT POINT (IT CURRENT)
          
         ITT=IT    ! STORAGE THE POSITION OF BEST ITERATING YET.          
         IT =1
         SIGMI = -1.0D100
         WW(NPT+1,1) = 1.0D0
         DO I=1, N
         WW(I+NPT+1, 1) =  X(I)- XBASE_A(I)
         END DO 
         DO I=1, NPT
          CONT = 0.0D0
         DO K=1, N         
         CONT = CONT +  (Y(I,K) - XBASE_A(K)) * (X(K) - XBASE_A(K))          
         END DO
         WW(I, 1) = 0.5D0*CONT**2
         END DO  
           
C       CALCULATING ALL SIGMA AND CHOOSE IT FOR THE GREATER SIGMA. MULTIPLY T-th LINE OF H FOR W FOR OBTAIN TAU. 
          DO KKK = 1, NPT-1            
             CONT=0.0D0
             DO I=1, NPT+ N +1
              CONT = CONT +  H(IT, I) * WW(I,1)
             END DO
             TAU = CONT
C      CALCULUS OF BETA = 0.5||X^+-XBASE||-WHW                  
            DO I=1, NPT+N+1
             CONT = 0.0D0 
             DO K=1, NPT+N+1
              CONT = CONT + H(I,K) * WW(K,1)
              VETOR1(I) = CONT                  
             END DO             
            END DO
               
            CONT = 0.0D0
            DO I=1,  NPT+N+1
             CONT = CONT+WW(I,1)*VETOR1(I)
             AGUARD = CONT
            END DO
                       
C       CALCULUS  OF X-XB^4
          CONT = 0.0D0
         DO I=1, N
          CONT = CONT + (X(I)-XBASE_A(I))**2
         END DO         
            BETA =   0.5D0 * CONT**2 - AGUARD
            ALFA = H(IT,IT) 
            SIGM = ALFA * BETA + TAU**2  
            IF (SIGM .GE. SIGMI .AND. IT .NE. ITT) THEN                      
                 SIGMI = SIGM
                 IAUXILIAR  =  IT                
		 AUXILIAR(1) = ALFA
 		 AUXILIAR(2) = BETA
		 AUXILIAR(3) = TAU
                 AUXILIAR(4) = SIGM                           
             END IF 
           IT = IT + 1
       END DO
          
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
                    IT   = IAUXILIAR  
		    ALFA = AUXILIAR(1) 
 		    BETA = AUXILIAR(2)
		    TAU  = AUXILIAR(3) 
                    SIGM = AUXILIAR(4)
         
           CONT=0.0D0
           DO I=1, NPT+ N +1
            CONT = CONT +  H(IT, I) * WW(I,1)
           END DO
           TAU = CONT
          
           RETURN
         END SUBROUTINE SIGMA 

C     ******************************************************************
C     ******************************************************************

C **************** UPDAT THE INVERSE H ******************************
          SUBROUTINE INVERSAH(H, N, NPT,VETOR1,SIGM,IT,ALFA,BETA,TAU)
          IMPLICIT REAL*8 (A-H,O-Z)
#include "tr_params.par"
           DIMENSION VETOR1(*)
           DIMENSION P1(NPT+N+1,NPT+N+1) ,H(NPT+N+1,NPT+N+1)
           DIMENSION P2(NPT+N+1,NPT+N+1),P3(NPT+N+1,NPT+N+1)
C         ALFA*(E-MM) * (E-MM)'- BETA* H * E * E'*H+TAU*H*E*(E-MM)'+ (E-MM)* E'* H
C         MM = H*WW THAT IS STORED IN VETOR1 
           VETOR1(IT) = VETOR1(IT)-1.0D0              
          DO I=1, N+NPT+1
            DO J=1, N+NPT+1
              P1(I,J)= VETOR1(I) * VETOR1(J)
              P2(I,J)= H(I, IT) * H(J, IT)
              P3(I,J)= (H(IT, I)*(-VETOR1(J)))+(H(IT,J)*(-VETOR1(I)))
            END DO
          END DO
            DO I=1, N+NPT+1
             DO J=1, N+NPT+1
            if (sigm .eq. 0d0) return
         H(I,J)=H(I,J)+(1/SIGM)*(ALFA*P1(I,J)-BETA*P2(I,J)+TAU*P3(I,J)) 
             END DO
            END DO
          RETURN
          END SUBROUTINE INVERSAH 
 
C     ******************************************************************
C     ******************************************************************

       SUBROUTINE  SUBPROBLEMA(N, NPT, Q, DELTA, D, X, XL, XU, DSQ,
     +                         M, EQUATN, LINEAR, CCODED, XEPS, FLAG)
       IMPLICIT REAL*8 (A-H,O-Z)
#include "tr_params.par"
       DIMENSION Q(*), X(*), XL(*),XU(*)
       DOUBLE PRECISION XBASE_A(INN) 
       DOUBLE PRECISION VQUAD_A, GRADD , VQUAD , XANTIGO(INN),L(N)
       DOUBLE PRECISION H(NPT+N+1,NPT+N+1),XOPT_A(INN),D(INN),U(N) 
       INTEGER N

!     SCALAR ARGUMENTS
       integer flag,m

!     ARRAY ARGUMENTS
       logical ccoded(2),equatn(m), linear(m)

       COMMON /XBASEA/ XBASE_A           
       COMMON /VQUADA/VQUAD_A, VQUAD   

!     LOCAL SCALARS
       double precision cnorm
         
         DO I = 1,N
          L(I)=  DMAX1(XL(I) - XBASE_A(I),X(I) - XBASE_A(I)-DELTA) 
          U(I)=  DMIN1(XU(I) - XBASE_A(I),X(I) - XBASE_A(I)+DELTA)           
          XANTIGO(I) = X(I)       
         END DO
            
         DO I=1, N        
          D(I)=X(I)-XBASE_A(I)         
      
         END DO   

         CALL MEVALF(N,D,F,FLAG)

         IF ( FLAG .NE. 0 ) RETURN

         VQUAD=  F + Q(1)                                 
         
         CALL SOLVER(N, L, U, D, M, EQUATN, LINEAR, CCODED, .false.,
     +               XEPS, CNORM, FLAG)
         
         IF ( FLAG .NE. 0 ) RETURN
                   
C CALCULUS THE STEP LENGTH 
           SUM = 0.0D0
           DO I=1, N
            SUM = SUM + (X(I)-(D(I)+XBASE_A(I)))**2
            END DO
            DSQ = SUM
          DO I=1, N
           X(I)= D(I) +  XBASE_A(I)   
          END DO   

           CALL MEVALF(N,D,F,IFLAG)

           IF ( FLAG .NE. 0 ) RETURN

           VQUAD_A=  F + Q(1) !MODEL IN XNOVO 
           VQUAD=  -VQUAD + VQUAD_A  
         IF (VQUAD .GE. 0.D0 .or. cnorm .gt. xeps)  THEN
          DO I=1, N
           X(I) = XANTIGO(I) 
          END DO
          DSQ = 0D0
         END IF

         RETURN
        END

C     ******************************************************************
C     ****************************************************************** 

       SUBROUTINE  ATUALIZAQ(H, N, NPT, Q, DELTA, Y, X, F, IT)
       IMPLICIT REAL*8 (A-H,O-Z)
#include "tr_params.par"
       DOUBLE PRECISION Q(*), X(*) 
     
       DIMENSION H(NPT+N+1,NPT+N+1),Y(NPT,N),VETORAUX(INN),DD(N,N)
       DOUBLE PRECISION  GOPT_A(INN), HQ_A(INN**2), XBASE_A(INN)  
       DIMENSION QQ((N+1)*(N+2)/2),  TEMP(1+N+NPT)     
             COMMON /VQUADA/ VQUAD_A, VQUAD 
             COMMON /XBASEA/ XBASE_A  
             COMMON / NOMETESTE /  GOPT_A, HQ_A    
               
              DO I=1, 1+N+NPT
               TEMP(I) =  (F - VQUAD_A)* H(I, IT) ! IS LAMBDA 
              END DO                     
               DO I=1, N 
                DO J=1, N
                 DD(I,J)=0.0D0
                END DO
               END DO
C             M=M+     LAMBCG(J)*((Y(:,J)-XB') * ( Y(: ,J)-XB')')  
             DO I=1, NPT
               DO JJ=1, N
                VETORAUX(JJ) =  Y(I,JJ)-XBASE_A(JJ)  
               END DO                      
                  DO K=1, N                    
                    DO J=1, N
                  DD(K,J) = DD(K,J)+ TEMP(I)* VETORAUX(K) * VETORAUX(J)  
C      DDEH IS THE GRADIENT OF QUADRATIC D                               
                    END DO                   
                  END DO                  
             END DO
C  N+1 FIRST ELEMENTS OF THE QQ (PARAMETER OF THE NEW MODEL)
            QQ(1) = TEMP(NPT+1)
            DO I=2, N+1 
             QQ(I) = TEMP(NPT+I)
            END DO   
C PUT IN THE QQ FOR SYMMETRY OF DD
              II=1
               J=1
C  DO WHILE (II .LE.  (N+1)*N/2 ) 
              DO WHILE (II .LE. N) 
              DO I=1, II 
               QQ(J+1+N) = DD(II,I) 
                J=J+1
               END DO
                II = II + 1
               END  DO   
C ADD Q TO QQ.  Q IS THE OLD MODEL, QQ IS THE MODEL D, AND Q + D = Q+    
              DO I=1, 1+N+ (N+1)*N/2
               Q(I) = Q(I) + QQ(I)
              END DO
C  UPDAT THE MODEL IN THE FILE INIP FOR COMMON FUNCTION 
               DO I=1, N      
               GOPT_A(I) = Q(I+1)
               END DO   
               DO J=1, (N+1)*N/2 
                  HQ_A(J) =  Q(1+N+J) 
               END DO  
           RETURN
         END SUBROUTINE  ATUALIZAQ  

C     ******************************************************************
C     ******************************************************************

      subroutine mvv(v1,v2,n, gradd)
      
      implicit none

#include "tr_params.par"

!     multiplica vetor por vetor
      double precision v1(INN), v2(INN), soma, gradd
      integer j, n
      
      soma=0
      do j=1 , n
         soma=soma +  v1(j)*v2(j)
         gradd=soma
      end do
      return
      end subroutine mvv

C     ******************************************************************
C     ******************************************************************

      subroutine  mmv(HQ, S,n, v)

      implicit none

#include "tr_params.par"

!     multiplica matriz simetrica (dada como vetor) por vetor
!     estava com hs, e troquei para hss para nao atualizar hs desnec
      
      double precision S(INN), HQ(INN ** 2), HSS(INN), v(INN)
      
      integer n, i, j , IH
      IH=0
      DO J=1,N
         HSS(J)= 0
         DO I=1,J
            IH=IH+1
            IF (I .LT. J) HSS(J)=HSS(J)+HQ(IH)*S(I)
            HSS(I)=HSS(I)+HQ(IH)*S(J)
            v(I)=HSS(I)
         end DO
      end DO
      return
      end subroutine mmv


$SIZES      LNP4=5000
$PROBLEM    Final model of glucagon, ACTH and cortisol
$INPUT      STD   ; Study: 1 - Lundqvist (hypo), 2 - Abrahamsson, 3 - Almby, 4 - Lundqvist (Hyper)
            ID    ; Identification of individual
            TIME  ; Time in hours
            DV    ; Logarithm of DV
            ODV   ; Original DV: glucagon [pmol/L], ACTH [pmol/L], cortisol [nmol/L]
            FLAG  ; Indicator of variable: 1 - glucagon, 2 - cortisol, 3 - ACTH
            EVID  ; Event identification: 0 - observation, 1 - dosing, 2 - missing
            L2    ; Indictor for observations collected at the same time-point
            GLU   ; Observed glucose [mmol/L]
            SLG   ; Slope between adjacent glucose observations
            INS   ; Observed insulin [mU/L]
            SLI   ; Slope between adjacent insulin observations
            SEX   ; Sex: 0 - man, 1 - woman
            AGE   ; Age [years]
            BMI   ; Body mass index [kg/m2]
            MV    ; M-value []
            FPG   ; Fasting plasma glucose [mmol/L]
            HBA1C ; Glycated hemoglobin [%]
            WHR   ; Waist-Hip ratio [-]
            BFP   ; Body fat percentage [%]
            HOMAIR; HOMA-IR [-]
$DATA       data.csv IGNORE=@
$SUBROUTINE ADVAN13 TOL=5
$MODEL      COMP=(INS)  ; linear interpolation of insulin 
            COMP=(GLUC) ; linear interpolation of glucose
            COMP=(GN)   ; glucagon cmt
            COMP=(ACTH) ; ACTH cmt
            COMP=(G2)   ; effect cmt of glucose
            COMP=(CORTISOL) ; cortisol cmt
$PK

OCC = 1
IF(STD.EQ.4) OCC = 2

;-------Setting variables for linear interpolation of insulin and glucose---------
IF(TIME.EQ.0) THEN
  PSLG  = 0
  PSLI  = 0
  PGLU  = GLU
  PINS  = INS
  PPINS = PINS
ENDIF

NSLG = PSLG
NSLI = PSLI
IF(SLG.NE.-99) NSLG = SLG
IF(SLI.NE.-99) NSLI = SLI

;------------PARAMETERS-----------------
ADDRES = THETA(18)

;-----------Covariate parameters---------------
IC50GICOV = EXP(THETA(20)*(MV - 8.13))
IC50ACOV  = EXP(THETA(19)*(MV - 8.13))

;------------Glucagon parameters----------------
TVKoutG  = THETA(1)
TVGn0    = THETA(2)
TVFrac   = THETA(3)
TVIC50GG = THETA(4)
TVHillG  = THETA(5)
TVIC50GI = THETA(6)*IC50GICOV

;----------ACTH parameters------------------------
TVKoutA  = THETA(7)
TVACTH0  = THETA(8)
TVt12    = THETA(9)
TVImaxA  = THETA(10)
TVIC50A  = THETA(11)*IC50ACOV
TVHillA  = THETA(12)

;--------Cortisol parameters----------------------
TVKoutC  = THETA(13)
TVCort0  = THETA(14)
TVEmaxC  = THETA(15)
TVIC50C  = THETA(16)
TVHillC  = THETA(17)

;-------INDIVIDUAL PARAMETERS---------------------

;-------Glucagon parameters--------------------  
IF(STD.NE.4) IOVG0 = ETA(8)
IF(STD.EQ.4) IOVG0 = ETA(9)

KoutG   = TVKoutG                      
Gn0     = TVGn0    * EXP(ETA(1) + IOVG0)
Frac    = TVFrac
IC50GG  = TVIC50GG * EXP(ETA(2))
HillG   = TVHillG 
IC50GI  = TVIC50GI * EXP(ETA(3))

;-------ACTH parameters----------------------
IF(STD.NE.4) IOVA0 = ETA(10)
IF(STD.EQ.4) IOVA0 = ETA(11)

KoutA  = TVKoutA
ACTH0  = TVACTH0 * EXP(ETA(4) + IOVA0)
t12    = TVt12    
ImaxA  = TVImaxA  
IC50A  = TVIC50A * EXP(ETA(5))
HillA  = TVHillA

;------Cortisol parameters------------------
KoutC  = TVKoutC  
Cort0  = TVCort0  * EXP(ETA(6))
EmaxC  = TVEmaxC  
IC50C  = TVIC50C  * EXP(ETA(7))
HillC  = TVHillC

;-----Baseline Glucose and Insulin ---------
IF(TIME.EQ.0) THEN
  INS0 = INS
  GLU0 = GLU
ENDIF

;------Derived parameters------------------------------------
KE0    = LOG(2)/t12

GnG0   = Frac*GLU0**HillG/(IC50GG**HillG + GLU0**HillG)
GnI0   = (1-Frac)*INS0/(IC50GI + INS0)

GnGI0  =  1 - (GnG0 + GnI0)
KinG   = Gn0*KoutG/GnGI0

AG0    = 1 - ImaxA * GLU0**HillA/(IC50A**HillA + GLU0**HillA)
KinA   = ACTH0*KoutA/AG0

CA0    = EmaxC * ACTH0**HillC/(IC50C**HillC + ACTH0**HillC)
KinC   = Cort0*KoutC/CA0

;-------Initialization---------------------------------------

A_0(1) = INS0
A_0(2) = GLU0
A_0(3) = Gn0
A_0(4) = ACTH0
A_0(5) = GLU0
A_0(6) = Cort0

;------------ODEs---------------------------------------------
$DES

INT = A(1)
IF(INT.LE.0) INT = 1E-6

IF(A(2).LE.0) THEN
   GnG = 0
ELSE
   GnG = Frac*A(2)**HillG/(IC50GG**HillG + A(2)**HillG)
ENDIF

IF(A(1).LE.0) THEN
   GnI = 0
ELSE
   GnI = (1-Frac)*A(1)/(IC50GI + A(1))
ENDIF

GnGI = 1 - (GnG + GnI)

IF(A(5).LE.0) THEN
   AG = 1
ELSE
   AG = 1 - ImaxA * A(5)**HillA / (IC50A**HillA + A(5)**HillA)
ENDIF

IF(A(4).LE.0) THEN
   CA = 0
ELSE
   CA  = EmaxC * A(4)**HillC / (IC50C**HillC + A(4)**HillC)
ENDIF

;---------------------------Differential equations----------------------------------------------
DADT(1) = NSLI                             ;Insulin
DADT(2) = NSLG                             ;Glucose
DADT(5) = KE0*(A(2) - A(5))                ;Effect compartment Glucose on ACTH

DADT(3) = KinG*GnGI - KoutG*A(3)           ;Glucagon
DADT(4) = KinA*AG   - KoutA*A(4)           ;ACTH
DADT(6) = KinC*CA   - KoutC*A(6)           ;Cortisol
;----------------Residual Unexplained Variability Model----------------------------------------
$ERROR

IF(FLAG.EQ.1) THEN
   LOQ = 1.7
   IPRED = A(3)                               ;Glucagon Proportional RUV
   RUV = EPS(1)
   SD = SQRT(SIGMA(1,1))
ENDIF

IF(FLAG.EQ.2) THEN                            ;Cortisol Proportional RUV
   LOQ = 70
   IPRED = A(6)
   RUV = EPS(3)
   SD = SQRT(SIGMA(3,3))
ENDIF

IF(FLAG.EQ.3) THEN                            ;ACTH Proportional RUV
  LOQ = 1.1
  IPRED = A(4)
  RUV = EPS(2)
  SD = SQRT(SIGMA(2,2))
ENDIF

IRES = EXP(DV)-IPRED
IF(SD.EQ.0) SD=1
IWRES = IRES/SD

IF(IPRED.LE.0) THEN 
   LGIPRD = -4.605
ELSE
   LGIPRD = LOG(IPRED)
END IF
IRES = EXP(DV)-IPRED

Y = LGIPRD + RUV
IF(IPRED.LE.LOQ) Y = LGIPRD + RUV*ADDRES

PSLG  = NSLG
PSLI  = NSLI
PGLU  = GLU 
PINS  = INS
PPINS = PINS

;---------------Initial estimates-----------------------------------------
$THETA  (0.01,0.0612001,1)  ; 1. Kout glucagon
$THETA  (6,8.34766,9)       ; 2. Baseline glucagon
$THETA  (0,0.750223,1)      ; 3. Frac
$THETA  (1,2.90692,5)       ; 4. IC50_GG
$THETA  (1,5.26063,10)      ; 5. HillG
$THETA  (4,5.51287,12)      ; 6. IC50_GI
$THETA  (0.01,0.0708762,1)  ; 7. Kout ACTH
$THETA  (1,2.94248,4)       ; 8. Baseline ACTH
$THETA  (2,8.67535,10)      ; 9. t1/2 glucose delay
$THETA  (0,0.941753,1)      ;10. ImaxA
$THETA  (1,2.85041,5)       ;11. IC50_A
$THETA  (1,20.8777,25)      ;12. HillA
$THETA  (0.01,0.0608134,1)  ;13. Kout cortisol
$THETA  (100,206.928,300)   ;14. Baseline cortisol
$THETA  (1,41.7789,60)      ;15. EmaxC
$THETA  (1,6.90405,10)      ;16. IC50_C
$THETA  (1,1.28996,25)      ;17. HillC
$THETA  (1,2.00765,3)       ;18. Add. uncertainty for BLQ
$THETA  (-0.5,-0.01551,0.5) ;19. MV on IC50_A
$THETA  (-0.2,-0.15697,0.2) ;20. MV on IC50_GI

$OMEGA  BLOCK(1) 0.156437   ; 1. IIV Baseline glucagon
$OMEGA  BLOCK(1) 0.0546883  ; 2. IIV IC50_GG
$OMEGA  BLOCK(1) 0.868819   ; 3. IIV IC50_GI
$OMEGA  BLOCK(1) 0.0862211  ; 4. IIV Baseline ACTH
$OMEGA  BLOCK(1) 0.00602187 ; 5. IIV IC50_A
$OMEGA  BLOCK(1) 0.0631253  ; 6. IIV Baseline cortisol
$OMEGA  BLOCK(1) 0 FIX      ; 6. IIV IC50_C
$OMEGA  BLOCK(1) 0.0401028  ; 7. IOV Baseline glucagon STD 1-3
$OMEGA  BLOCK(1) SAME       ; 8. IOV -"- STD 4
$OMEGA  BLOCK(1) 0.00656913 ; 9. IOV Baseline ACTH STD 1-3
$OMEGA  BLOCK(1) SAME       ;10. IOV -"- STD 4

$SIGMA  BLOCK(1) 0.114088   ; 1. RUV Glucagon
$SIGMA  BLOCK(1) 0.144903   ; 2. RUV ACTH
$SIGMA  BLOCK(1) 0.0856894  ; 3. RUV Cortisol

$ESTIMATION MAXEVAL=9990 PRINT=1 METHOD=1 INTER SADDLE_RESET=1
            SADDLE_HESS=1 NOABORT
$COVARIANCE UNCOND

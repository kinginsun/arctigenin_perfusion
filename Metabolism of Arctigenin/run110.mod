;; 1. Based on: run109
;; 2. Description: Arctigenin to AA + AG, IG, mean Cl 2, 2.4 mg
;; x1. Author: user

$PROBLEM  METABOLISM OF ARCTIGENIN IN RATS AFTER IV - ADVAN6
$DATA   Arctigenin_IG_for_nonmem_H.csv IGNORE=@
$INPUT  ID TIME DV DOSE=AMT MDV CMT DL
$SUBROUTINE ADVAN6 TOL=6
$MODEL
COMP=(LUMEN DEFDOSE)			; 1. intestine lumen
COMP=(AG_CENTR)				; 2. AG in central
COMP=(AA_CENTR)				; 3. AA in Central
COMP=(INTESTIN)				; 4. Intestine compartment
COMP=(AG_INTES)				; 5. AG in intestine
COMP=(AG_PER)				; 6. AG in perp
COMP=(AA_PER)				; 6. AA in perp

$PK
Ka = THETA(1)*EXP(ETA(1))		; absorption rate of AR
K20= THETA(2)				; elimination rate of AG from Central
K30= THETA(3)				; elimination rate of AA from central comp
;K52= THETA(4)				; AG from intestine to central
;K63= THETA(5)				; AR from intestine to central
V2=  THETA(4)				; volume of central
;V4=  THETA(7)				; volume of intestine
S2=  V2
S3=  V2
Ki= THETA(5)				; scale factor for blood hydrolysis
Kb= THETA(6)				; scale factor for liver glucuronidation
;Fi3= THETA(10)				; scale factor for intestine glucuronidation
;Km1= THETA(11)				; AA blood hydrolysis
;Vm1= THETA(12)
;Km2= THETA(13)				; AG liver glucuronidation
;Vm2= THETA(14)
;Km3= THETA(15)				; AG intestine glucuronidation
;Vm3= THETA(16)
K52= THETA(7)
K25= THETA(8)
K63= THETA(9)
K36= THETA(10)

$DES
DADT(1)=-Ka*A(1)
DADT(2)= Ki*A(4)-K20*A(2)-K25*A(2)+K52*A(5)
DADT(3)= Kb*A(4)-K30*A(3)-K36*A(3)+K63*A(6)
DADT(4)= Ka*A(1)-Ki*A(4)-Kb*A(4)
DADT(5)= K25*A(2)-K52*A(5)
DADT(6)= K36*A(3)-K63*A(6)

$ERROR
Y = F+ERR(1)

$THETA
(0, 17.5) ; Ka
(0, 0.3) ; K20
(0, 0.1) ; K30
;(0, 0.132) ; K52
;(0, 0.132) ; K63
(0, 70) ; V2
(0, 6.4) ; Ki
(0, 11) ; Kb
;(0, 0.0157) FIX ; Fi2
;(0, 1) FIX ; Fi3
;(13.5) FIX ; Km1
;(0.292) FIX ; Vm1
;(2.302) FIX ; Km2
;(1.231) FIX ; Vm2
;(20.38) FIX ; Km3
;(12.449) FIX ; Vm3
(0, 0.1) ; K52
(0, 0.1) ; K25
(0, 0.1) ; K63
(0, 0.1) ; K36

$OMEGA
0 FIX

$SIGMA
0.326 ; Residual error

$EST METHOD=0 MAXEVAL=999 NOABORT SIG=5 PRINT=1
$COV

; Xpose
$TABLE ID TIME DV MDV CMT DL ONEHEADER NOPRINT FILE=sdtab110.tab
$TABLE ID Ka K20 K30 V2 ONEHEADER NOPRINT FILE=patab110.tab


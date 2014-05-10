;; 1. Based on: run105
;; 2. Description: Arctigenin to AA + AG, IG, mean Cl
;; x1. Author: user

$PROBLEM  METABOLISM OF ARCTIGENIN IN RATS AFTER IV - ADVAN6
$DATA   Arctigenin_IG_for_nonmem.csv IGNORE=@
$INPUT  ID TIME DV DOSE=AMT MDV CMT DL
$SUBROUTINE ADVAN6 TOL=6
$MODEL
COMP=(LUMEN DEFDOSE)			; 1. intestine lumen
COMP=(AG_CENTR)				; 2. AG in central
COMP=(AA_CENTR)				; 3. AA in Central
COMP=(INTESTIN)				; 4. Intestine compartment
COMP=(AG_INTES)				; 5. AG in intestine
;COMP=(AA_INTES)				; 6. AA in intestine
COMP=(AR_CENTR)				; 6. AR in central

$PK
Ka = THETA(1)*EXP(ETA(1))		; absorption rate of AR
K20= THETA(2)				; elimination rate of AG from Liver, not considering the transfer from Liver to Central
K30= THETA(3)				; elimination rate of AA from central comp
K52= THETA(4)				; AG from intestine to central
K46= THETA(5)				; AR from intestine to central
V6=  THETA(6)				; volume of central
V4=  THETA(7)				; volume of intestine
S2=  V6
S3=  V6
Fi1= THETA(8)				; scale factor for blood hydrolysis
Fi2= THETA(9)				; scale factor for liver glucuronidation
Fi3= THETA(10)				; scale factor for intestine glucuronidation
Km1= THETA(11)				; AA blood hydrolysis
Vm1= THETA(12)
Km2= THETA(13)				; AG liver glucuronidation
Vm2= THETA(14)
Km3= THETA(15)				; AG intestine glucuronidation
Vm3= THETA(16)

$DES
DADT(1)=-Ka*A(1)
DADT(2)= K52*A(5)-K20*A(2)+Fi2*Vm2*A(6)/(Km2*V6+A(6))
DADT(3)=-K30*A(3)+Fi1*Vm1*A(6)/(Km1*V6+A(6))
DADT(4)= Ka*A(1)-Fi3*Vm3*A(4)/(Km3*V4+A(4))-K46*A(4)
DADT(5)= Fi3*Vm3*A(4)/(Km3*V4+A(4))-K52*A(5)
DADT(6)= K46*A(4)-Fi1*Vm1*A(6)/(Km1*V6+A(6))-Fi2*Vm2*A(6)/(Km2*V6+A(6))
$ERROR
Y = F+ERR(1)

$THETA
(0, 1) ; Ka
(0, 0.206) ; K20
(0, 0.132) ; K30
(0, 0.132) ; K52
(0, 0.132) ; K63
(0, 0.302) FIX ; V6
(0, 0.311) ; V4
(0, 1.8) FIX ; Fi1
(0, 0.0157) FIX ; Fi2
(0, 1) FIX ; Fi3
(13.5) FIX ; Km1
(0.292) FIX ; Vm1
(2.302) FIX ; Km2
(1.231) FIX ; Vm2
(20.38) FIX ; Km3
(12.449) FIX ; Vm3

$OMEGA
0 FIX

$SIGMA
0.326 ; Residual error

$EST METHOD=0 MAXEVAL=999 NOABORT SIG=5 PRINT=1
$COV

; Xpose
$TABLE ID TIME DV MDV CMT DL ONEHEADER NOPRINT FILE=sdtab108.tab
$TABLE ID Ka K20 K30 K52 K46 V4 Fi3 ONEHEADER NOPRINT FILE=patab108.tab


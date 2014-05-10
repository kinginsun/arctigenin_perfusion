;; 1. Based on: run103
;; 2. Description: prior method
;; x1. Author: user
;; 3. Label:
$PROB BAYES TEST 
$DATA data.csv IGNORE=# 
$INPUT ID TIME OBS=DV DVID
$ESTIMATE MAXEVALS=9990
METHOD=COND SLOW

$THETA 1 ; K
$OMEGA 0 FIX ; PPVK
$SIGMA 0.04 FIX ; RUV

$SIGMA 4 FIX ; Kprior_uncertainty
$PRED 
K = THETA(1) + ETA(1)
C = 10*EXP(-K*TIME)
IF (DVID.EQ.1) THEN ; prior
Y=THETA(1)+ERR(2)
ENDIF
IF (DVID.EQ.2) THEN ; obs
Y = C + ERR(1)
ENDIF
$TABLE ID TIME K ETA1 NOPRINT NOAPPEND ONEHEADER FILE=patab104.tab


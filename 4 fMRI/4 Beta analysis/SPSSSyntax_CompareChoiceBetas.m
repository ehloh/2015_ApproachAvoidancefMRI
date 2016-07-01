% * SPSS-SYNTAX COMMANDS FOR BETA ANALYSIS
% *	All statistical tests for comparing choice betas (not cluster models)
% *	Change roi name (default: HPC_L) preferably in Matlab
% *
% * ----------------------------------------------------------------------.

%% BETA ANALYSIS

* HPC_L.
COMPUTE mean_HPC_L=MEAN(HPC_LcF_Accept,HPC_LcF_Reject, HPC_LcF_Explore,HPC_Lct_NoBomb,HPC_Lct_Bomb,HPC_Lct_Explore).
COMPUTE HPC_LcF_Accept  =  HPC_LcF_Accept    -mean_HPC_L.
COMPUTE HPC_LcF_Reject   =  HPC_LcF_Reject    -mean_HPC_L.
COMPUTE HPC_LcF_Explore  =  HPC_LcF_Explore    -mean_HPC_L.
COMPUTE HPC_Lct_NoBomb    =  HPC_Lct_NoBomb    -mean_HPC_L.
COMPUTE HPC_Lct_Bomb    =  HPC_Lct_Bomb    -mean_HPC_L.
COMPUTE HPC_Lct_Explore  =  HPC_Lct_Explore    -mean_HPC_L.
EXECUTE.

GGRAPH
  /GRAPHDATASET NAME="graphdataset" VARIABLES=MEANSE(HPC_LcF_Accept, 1) MEANSE(HPC_LcF_Reject,1) 
    MEANSE(HPC_LcF_Explore, 1) MEANSE(HPC_Lct_NoBomb, 1) MEANSE(HPC_Lct_Bomb, 1) 
    MEANSE(HPC_Lct_Explore, 1) MISSING=LISTWISE REPORTMISSING=NO
    TRANSFORM=VARSTOCASES(SUMMARY="#SUMMARY" INDEX="#INDEX" LOW="#LOW" HIGH="#HIGH")
  /GRAPHSPEC SOURCE=INLINE.
BEGIN GPL
  SOURCE: s=userSource(id("graphdataset"))
  DATA: SUMMARY=col(source(s), name("#SUMMARY"))
  DATA: INDEX=col(source(s), name("#INDEX"), unit.category())
  DATA: LOW=col(source(s), name("#LOW"))
  DATA: HIGH=col(source(s), name("#HIGH"))
  GUIDE: axis(dim(2), label("Mean"))
  GUIDE: text.title(label("HPC_L"))
  GUIDE: text.footnote(label("Error Bars: +/- 1 SE"))
  SCALE: cat(dim(1), include("0", "1", "2", "2.5","3", "4", "5"))
  SCALE: linear(dim(2), include(0))
  ELEMENT: interval(position(INDEX*SUMMARY), shape.interior(shape.square))
  ELEMENT: interval(position(region.spread.range(INDEX*(LOW+HIGH))), shape.interior(shape.ibeam))
END GPL.

GLM HPC_LcF_Accept HPC_LcF_Reject HPC_LcF_Explore HPC_Lct_NoBomb HPC_Lct_Bomb HPC_Lct_Explore
  /WSFACTOR=TASK 2 Polynomial CHOICE 3 Polynomial 
  /METHOD=SSTYPE(3)
  /EMMEANS=TABLES(TASK*CHOICE)  COMPARE(CHOICE)
  /CRITERIA=ALPHA(.05)
  /WSDESIGN=TASK CHOICE TASK*CHOICE.

T-TEST
  /TESTVAL=0
  /MISSING=ANALYSIS
  /VARIABLES=HPC_LcF_Accept HPC_LcF_Reject HPC_LcF_Explore HPC_Lct_NoBomb HPC_Lct_Bomb 
    HPC_Lct_Explore
  /CRITERIA=CI(.95).

T-TEST PAIRS=HPC_LcF_Accept HPC_LcF_Reject HPC_LcF_Explore WITH HPC_Lct_NoBomb HPC_Lct_Bomb 
    HPC_Lct_Explore (PAIRED)
  /CRITERIA=CI(.9500)
  /MISSING=ANALYSIS.


  
  
  
  
  
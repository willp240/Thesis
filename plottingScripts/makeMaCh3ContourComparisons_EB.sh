#!/bin/bash
OUTFOLDER=compMaCh3Contours_EB_asimov
OUTNAMEBASE="compMaCh3Contours_EB"
rm -rf $OUTFOLDER
mkdir $OUTFOLDER
OUTNAME="$OUTFOLDER/$OUTNAMEBASE"

CONTOURFOLDER2020="/vols/build/t2k/wparker2/MaCh3_2020NewSplines/MaCh3/contours_AsimovA_poly"
CONTOURFOLDER2018="/vols/build/t2k/wparker2/MaCh3_2020NewSplines/MaCh3/contours_AsimovA_poly_EB"

WRCCONTOURFOLDER2020="/home/kwood/t2k/software/mach3-software/OA2020/jointFit_25Mar2020/paraWork/MaCh3/contours_wRC_AsimovA_poly"
WRCCONTOURFOLDER2018="/home/kwood/t2k/software/mach3-software/OA2020/jointFit_25Mar2020/paraWork/MaCh3/contours_15Apr2020_RCrw"


FILEBASE="contours_th23dm23_both.root"
NHFILEBASE="contours_th23dm23_NH.root"
IHFILEBASE="contours_th23dm23_IH.root"

APPFILEBASE="contours_th13dcp_both.root"
APPNHFILEBASE="contours_th13dcp_NH.root"
APPIHFILEBASE="contours_th13dcp_IH.root"

DCPFILEBASE="contours_app1D_both.root"
DCPNHFILEBASE="contours_app1D_NH.root"
DCPIHFILEBASE="contours_app1D_IH.root"

WRCFILEBASE="contours_th23dm23_both_RCrw.root"
WRCNHFILEBASE="contours_th23dm23_NH_RCrw.root"
WRCIHFILEBASE="contours_th23dm23_IH_RCrw.root"

WRCAPPFILEBASE="contours_th13dcp_both_RCrw.root"
WRCAPPNHFILEBASE="contours_th13dcp_NH_RCrw.root"
WRCAPPIHFILEBASE="contours_th13dcp_IH_RCrw.root"

WRCDCPFILEBASE="contours_app1D_both_RCrw.root"
WRCDCPNHFILEBASE="contours_app1D_NH_RCrw.root"
WRCDCPIHFILEBASE="contours_app1D_IH_RCrw.root"

CONTOURFILE2018="$CONTOURFOLDER2018/$FILEBASE"
CONTOURFILE2018NH="$CONTOURFOLDER2018/$NHFILEBASE"
CONTOURFILE2018IH="$CONTOURFOLDER2018/$IHFILEBASE"

CONTOURFILE2018APP="$CONTOURFOLDER2018/$APPFILEBASE"
CONTOURFILE2018APPNH="$CONTOURFOLDER2018/$APPNHFILEBASE"
CONTOURFILE2018APPIH="$CONTOURFOLDER2018/$APPIHFILEBASE"

CONTOURFILE2018DCP="$CONTOURFOLDER2018/$DCPFILEBASE"
CONTOURFILE2018DCPNH="$CONTOURFOLDER2018/$DCPNHFILEBASE"
CONTOURFILE2018DCPIH="$CONTOURFOLDER2018/$DCPIHFILEBASE"

CONTOURFILE2020="$CONTOURFOLDER2020/$FILEBASE"
CONTOURFILE2020NH="$CONTOURFOLDER2020/$NHFILEBASE"
CONTOURFILE2020IH="$CONTOURFOLDER2020/$IHFILEBASE"

CONTOURFILE2020APP="$CONTOURFOLDER2020/$APPFILEBASE"
CONTOURFILE2020APPNH="$CONTOURFOLDER2020/$APPNHFILEBASE"
CONTOURFILE2020APPIH="$CONTOURFOLDER2020/$APPIHFILEBASE"

CONTOURFILE2020DCP="$CONTOURFOLDER2020/$DCPFILEBASE"
CONTOURFILE2020DCPNH="$CONTOURFOLDER2020/$DCPNHFILEBASE"
CONTOURFILE2020DCPIH="$CONTOURFOLDER2020/$DCPIHFILEBASE"

# WRC
WRCCONTOURFILE2018APP="$WRCCONTOURFOLDER2018/$WRCAPPFILEBASE"
WRCCONTOURFILE2018APPNH="$WRCCONTOURFOLDER2018/$WRCAPPNHFILEBASE"
WRCCONTOURFILE2018APPIH="$WRCCONTOURFOLDER2018/$WRCAPPIHFILEBASE"

WRCCONTOURFILE2018="$WRCCONTOURFOLDER2018/$WRCFILEBASE"
WRCCONTOURFILE2018NH="$WRCCONTOURFOLDER2018/$WRCNHFILEBASE"
WRCCONTOURFILE2018IH="$WRCCONTOURFOLDER2018/$WRCIHFILEBASE"

WRCCONTOURFILE2018DCP="$WRCCONTOURFOLDER2018/$WRCDCPFILEBASE"
WRCCONTOURFILE2018DCPNH="$WRCCONTOURFOLDER2018/$WRCDCPNHFILEBASE"
WRCCONTOURFILE2018DCPIH="$WRCCONTOURFOLDER2018/$WRCDCPIHFILEBASE"

WRCCONTOURFILE2020APP="$WRCCONTOURFOLDER2020/$WRCAPPFILEBASE"
WRCCONTOURFILE2020APPNH="$WRCCONTOURFOLDER2020/$WRCAPPNHFILEBASE"
WRCCONTOURFILE2020APPIH="$WRCCONTOURFOLDER2020/$WRCAPPIHFILEBASE"

WRCCONTOURFILE2020="$WRCCONTOURFOLDER2020/$WRCFILEBASE"
WRCCONTOURFILE2020NH="$WRCCONTOURFOLDER2020/$WRCNHFILEBASE"
WRCCONTOURFILE2020IH="$WRCCONTOURFOLDER2020/$WRCIHFILEBASE"

WRCCONTOURFILE2020DCP="$WRCCONTOURFOLDER2020/$WRCDCPFILEBASE"
WRCCONTOURFILE2020DCPNH="$WRCCONTOURFOLDER2020/$WRCDCPNHFILEBASE"
WRCCONTOURFILE2020DCPIH="$WRCCONTOURFOLDER2020/$WRCDCPIHFILEBASE"



# 2018 vs 2020 woRC asimov
root -l -b -q 'utils/CompareContours.C("'$CONTOURFILE2018NH'","'$CONTOURFILE2020NH'","E_{b} Fixed", "E_{b} Free",false,"'$OUTNAME'_disapp_asimovA_NH",false,false,false)'
root -l -b -q 'utils/CompareContours.C("'$CONTOURFILE2018IH'","'$CONTOURFILE2020IH'","E_{b} Fixed", "E_{b} Free",false,"'$OUTNAME'_disapp_asimovA_IH",false,true,false)'

#root -l -b -q 'utils/CompareContours.C("'$CONTOURFILE2018'","'$CONTOURFILE2020'","TH2D","THPoly Preliminary",false,"'$OUTNAME'_disapp_asimovA_both_NHpost",true)'
#root -l -b -q 'utils/CompareContours.C("'$CONTOURFILE2018'","'$CONTOURFILE2020'","TH2D","THPoly Preliminary",false,"'$OUTNAME'_disapp_asimovA_both_IHpost",false)'

root -l -b -q 'utils/CompareContours.C("'$CONTOURFILE2018APPNH'","'$CONTOURFILE2020APPNH'","E_{b} Fixed", "E_{b} Free",true,"'$OUTNAME'_app_asimovA_NH",false,false,false)'
root -l -b -q 'utils/CompareContours.C("'$CONTOURFILE2018APPIH'","'$CONTOURFILE2020APPIH'","E_{b} Fixed", "E_{b} Free",true,"'$OUTNAME'_app_asimovA_IH",false,true,false)'
#root -l -b -q 'utils/CompareContours.C("'$CONTOURFILE2018APP'","'$CONTOURFILE2020APP'","TH2D","THPoly Preliminary",true,"'$OUTNAME'_app_asimovA_both",false)'

root -l -b -q 'utils/CompareContours1Ddcp.C("'$CONTOURFILE2018DCPNH'","'$CONTOURFILE2020DCPNH'","E_{b} Fixed", "E_{b} Free","'$OUTNAME'_dcp_asimovA_NH")'
root -l -b -q 'utils/CompareContours1Ddcp.C("'$CONTOURFILE2018DCPIH'","'$CONTOURFILE2020DCPIH'","E_{b} Fixed", "E_{b} Free","'$OUTNAME'_dcp_asimovA_IH")'
#root -l -b -q 'utils/CompareContours1Ddcp.C("'$CONTOURFILE2018DCP'","'$CONTOURFILE2020DCP'","TH2D","THPoly Preliminary","'$OUTNAME'_dcp_asimovA_both")'


# 2018 vs 2020 wRC asimov
# root -l -b -q 'utils/CompareContours.C("'$WRCCONTOURFILE2018NH'","'$WRCCONTOURFILE2020NH'","TH2D","THPoly",false,"'$OUTNAME'_wRC_disapp_asimovA_NH",true)'
# root -l -b -q 'utils/CompareContours.C("'$WRCCONTOURFILE2018IH'","'$WRCCONTOURFILE2020IH'","TH2D","THPoly",false,"'$OUTNAME'_wRC_disapp_asimovA_IH",true,true)'

# #root -l -b -q 'utils/CompareContours.C("'$CONTOURFILE2018'","'$CONTOURFILE2020'","TH2D","THPoly Preliminary",false,"'$OUTNAME'_disapp_asimovA_both_NHpost",true)'
# #root -l -b -q 'utils/CompareContours.C("'$CONTOURFILE2018'","'$CONTOURFILE2020'","TH2D","THPoly Preliminary",false,"'$OUTNAME'_disapp_asimovA_both_IHpost",false)'

# root -l -b -q 'utils/CompareContours.C("'$WRCCONTOURFILE2018APPNH'","'$WRCCONTOURFILE2020APPNH'","TH2D","THPoly",true,"'$OUTNAME'_wRC_app_asimovA_NH",true)'
# root -l -b -q 'utils/CompareContours.C("'$WRCCONTOURFILE2018APPIH'","'$WRCCONTOURFILE2020APPIH'","TH2D","THPoly",true,"'$OUTNAME'_wrC_app_asimovA_IH",true,true)'

# #root -l -b -q 'utils/CompareContours.C("'$CONTOURFILE2018APP'","'$CONTOURFILE2020APP'","TH2D","THPoly Preliminary",true,"'$OUTNAME'_app_asimovA_both")'

# root -l -b -q 'utils/CompareContours1Ddcp.C("'$WRCCONTOURFILE2018DCPNH'","'$WRCCONTOURFILE2020DCPNH'","TH2D","THPoly","'$OUTNAME'_wRC_dcp_asimovA_NH")'
# root -l -b -q 'utils/CompareContours1Ddcp.C("'$WRCCONTOURFILE2018DCPIH'","'$WRCCONTOURFILE2020DCPIH'","TH2D","THPoly","'$OUTNAME'_wRC_dcp_asimovA_IH")'
# #root -l -b -q 'utils/CompareContours1Ddcp.C("'$CONTOURFILE2018DCP'","'$CONTOURFILE2020DCP'","TH2D","THPoly Preliminary","'$OUTNAME'_dcp_asimovA_both")'

# # plot 1D dcp wRC for (i) 2018 (ii) 2020 (iii) 2020 with 2018 RC
# root -l -b -q 'utils/CompareContours1Ddcp.C("'$WRCCONTOURFILE2018DCPNH'","'$WRCCONTOURFILE2020DCPNH'","TH2D","THPoly","'$OUTNAME'_wRC2018_dcp_asimovA_NH","/home/kwood/t2k/software/mach3-software/OA2020/jointFit_25Mar2020/paraWork/MaCh3/contours_wRC2018_AsimovA_poly/'$WRCDCPNHFILEBASE'","THPoly, 2018 RC")'
# root -l -b -q 'utils/CompareContours1Ddcp.C("'$WRCCONTOURFILE2018DCPIH'","'$WRCCONTOURFILE2020DCPIH'","TH2D","THPoly","'$OUTNAME'_wRC2018_dcp_asimovA_IH","/home/kwood/t2k/software/mach3-software/OA2020/jointFit_25Mar2020/paraWork/MaCh3/contours_wRC2018_AsimovA_poly/'$WRCDCPIHFILEBASE'","THPoly, 2018 RC")'


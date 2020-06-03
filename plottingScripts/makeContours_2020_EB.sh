#!/bin/bash

for wrc in "false" # "true" -- note from Kevin - don't use this for wrc...separate script
do
for doasimov in "true" # "false"
do


if [[ "$doasimov" == "false" ]]
then
    trueth23=0
    truedm23=0
    trueth13=0
    truedcp=0
    if [[ "$wrc" == "true" ]]
    then
	    FILENAME="reduceddatasetsforTN320/wRC_data_joint_allto230717.root"
	    CONTOURDIR=contours_wRC_datafit_270717
	    #WRC Run  1-7C 1D dcp flat prior in sindcp limits
	    leftarr68=-0.66
	    rightarr68=-2.48
	    leftarr90=-0.28
	    rightarr90=-2.86
	    leftarr95=-0.1
	    rightarr95=-3.04
    else #worc
	    FILENAME="/storage/shared/kwood/MaCh3Data/OA2020/chains/jointDataFit_04May2020/MaCh3-jointDataFit2020_woRC_red.root"
	    CONTOURDIR=contours_woRC_datafit_18May2020
	    #WORC Run 1-7C 1D dcp flat prior in sindcp limits obtain from contours_app1D with (file,0,true,true) as arguments
	    burnin=120000
	    leftarr68=-0.27
	    rightarr68=-2.87
	    leftarr90=0.28
	    rightarr90=2.86
	    leftarr95=0.57
	    rightarr95=2.57
    fi
else #asimov
    trueth23=0.528
    truedm23=0.002509
    trueth13=0.0218
    truedcp=-1.601
    if [[ "$wrc" == "true" ]]
    then
	    FILENAME="reduced_25Mar2020_noBurnin_reweighted.root"
	    CONTOURDIR=contours_wRC_25Mar2020
	    burnin=150000
	    leftarr68=-0.426
	    rightarr68=-2.72
	    leftarr90=0.04
	    rightarr90=-3.10
	    leftarr95=0.28
	    rightarr95=2.86
    else #worc
	    #FILENAME="reduceddatasetsforTN320/woRC_asimovA_050717_100000burn.root"
	    FILENAME="/vols/build/t2k/wparker2/MaCh3_2020NewSplines/MaCh3/MaCh3_AsimovA_Poly_woRC_280520_EbSlice.root"
	    CONTOURDIR=contours_AsimovA_poly_EB
	    burnin=120000
	    #burnin=0
	    #WORC Run 1-7C 1D dcp flat prior in sindcp limits obtain from contours_app1D with (file,0,true,true) as arguments
	    leftarr68=-0.15
	    rightarr68=-2.99
	    leftarr90=0.52
	    rightarr90=2.63
	    leftarr95=0.86
	    rightarr95=2.28
    fi
fi


mkdir $CONTOURDIR

#Get 2D disappearance contours
echo "GETTING 2D DISAPPEARANCE CONTOURS:"
echo "  BOTH HIERARCHIES:"
root -l -q -b 'utils/contours_EB.C+("'$FILENAME'",'$trueth23','$truedm23',-999,0,'${burnin}','$wrc')'
mv contours_th23dm23.root $CONTOURDIR/contours_th23dm23_both.root
mv contours_th23dm23.pdf $CONTOURDIR/contours_th23dm23_both.pdf
mv contours_th23dm23_sigmas.pdf $CONTOURDIR/contours_th23dm23_sigmas_both.pdf
mv contours_th23dm23forroot.root $CONTOURDIR/contours_th23dm23forroot_both.root
#rm contours_th23dm23_conf.pdf
#rm contours_th23dm23forroot.root
#root -l -q -b 'utils/macros/makeOfficialPlotStyle.C+("'$CONTOURDIR'/contours_th23dm23forroot_both.root","'$CONTOURDIR'/contours_th23dm23_both_official.root","T2K Run 1-10b preliminary","c")'

echo "  NORMAL HIERARCHY:"
root -l -q -b 'utils/contours_EB.C+("'$FILENAME'",0,0,-999,1,'${burnin}','$wrc')'
mv contours_th23dm23.root $CONTOURDIR/contours_th23dm23_NH.root
mv contours_th23dm23.pdf $CONTOURDIR/contours_th23dm23_NH.pdf
mv contours_th23dm23_sigmas.pdf $CONTOURDIR/contours_th23dm23_sigmas_NH.pdf
mv contours_th23dm23forroot.root $CONTOURDIR/contours_th23dm23forroot_NH.root
#rm contours_th23dm23_conf.pdf
#rm contours_th23dm23forroot.root
#root -l -q -b 'utils/macros/makeOfficialPlotStyle.C+("'$CONTOURDIR'/contours_th23dm23forroot_NH.root","'$CONTOURDIR'/contours_th23dm23_NH_official.root","T2K Run 1-10b preliminary","c")'

echo "  INVERTED HIERARCHY:"
root -l -q -b 'utils/contours_EB.C+("'$FILENAME'",0,0,-999,-1,'${burnin}','$wrc')'
mv contours_th23dm23.root $CONTOURDIR/contours_th23dm23_IH.root
mv contours_th23dm23.pdf $CONTOURDIR/contours_th23dm23_IH.pdf
mv contours_th23dm23_sigmas.pdf $CONTOURDIR/contours_th23dm23_sigmas_IH.pdf
mv contours_th23dm23forroot.root $CONTOURDIR/contours_th23dm23forroot_IH.root
#rm contours_th23dm23_conf.pdf
#rm contours_th23dm23forroot.root
#root -l -q -b 'utils/macros/makeOfficialPlotStyle.C+("'$CONTOURDIR'/contours_th23dm23forroot_IH.root","'$CONTOURDIR'/contours_th23dm23_IH_official.root","T2K Run 1-10b preliminary","c")'



#get 2D appearance contours
echo "GETTING 2D APPEARANCE CONTOURS:"
echo "  BOTH HIERARCHIES:"
root -l -q -b 'utils/contours_app_EB.C+("'$FILENAME'",'$trueth13','$truedcp',-999,0,'${burnin}','$wrc')'
mv contours_th13dcp.root $CONTOURDIR/contours_th13dcp_both.root
mv contours_th13dcp.pdf $CONTOURDIR/contours_th13dcp_both.pdf
mv contours_th13dcp_sigmas.pdf $CONTOURDIR/contours_th13dcp_sigmas_both.pdf
mv contours_th13dcp_sigmas.root $CONTOURDIR/contours_th13dcp_sigmas_both.root
mv contours_th13dcp_conf.pdf $CONTOURDIR/contours_th13dcp_conf_both.pdf
mv contours_th13dcp_conf.root $CONTOURDIR/contours_th13dcp_conf_both.root
mv contours_th13dcpforroot.root $CONTOURDIR/contours_th13dcpforroot_both.root
#rm contours_th13dcp_conf.pdf
#rm contours_th13dcpforroot.root
#root -l -q -b 'utils/macros/makeOfficialPlotStyle.C+("'$CONTOURDIR'/contours_th13dcpforroot_both.root","'$CONTOURDIR'/contours_th13dcp_both_official.root","T2K Run 1-10b preliminary","c")'

echo "  NORMAL HIERARCHY:"
root -l -q -b 'utils/contours_app_EB.C+("'$FILENAME'",'$trueth13','$truedcp',-999,1,'${burnin}','$wrc')'
mv contours_th13dcp.root $CONTOURDIR/contours_th13dcp_NH.root
mv contours_th13dcp.pdf $CONTOURDIR/contours_th13dcp_NH.pdf
mv contours_th13dcp_sigmas.pdf $CONTOURDIR/contours_th13dcp_sigmas_NH.pdf
mv contours_th13dcp_sigmas.root $CONTOURDIR/contours_th13dcp_sigmas_NH.root
mv contours_th13dcp_conf.pdf $CONTOURDIR/contours_th13dcp_conf_NH.pdf
mv contours_th13dcp_conf.root $CONTOURDIR/contours_th13dcp_conf_NH.root
mv contours_th13dcpforroot.root $CONTOURDIR/contours_th13dcpforroot_NH.root
#rm contours_th13dcp_conf.pdf
#rm contours_th13dcpforroot.root
#root -l -q -b 'utils/macros/makeOfficialPlotStyle.C+("'$CONTOURDIR'/contours_th13dcpforroot_NH.root","'$CONTOURDIR'/contours_th13dcp_NH_official.root","T2K Run 1-10b preliminary","c")'
echo "  INVERTED HIERARCHY:"
root -l -q -b 'utils/contours_app_EB.C+("'$FILENAME'",'$trueth13','$truedcp',-999,-1,'${burnin}','$wrc')'
mv contours_th13dcp.root $CONTOURDIR/contours_th13dcp_IH.root
mv contours_th13dcp.pdf $CONTOURDIR/contours_th13dcp_IH.pdf
mv contours_th13dcp_sigmas.pdf $CONTOURDIR/contours_th13dcp_sigmas_IH.pdf
mv contours_th13dcp_sigmas.root $CONTOURDIR/contours_th13dcp_sigmas_IH.root
mv contours_th13dcp_conf.pdf $CONTOURDIR/contours_th13dcp_conf_IH.pdf
mv contours_th13dcp_conf.root $CONTOURDIR/contours_th13dcp_conf_IH.root
mv contours_th13dcpforroot.root $CONTOURDIR/contours_th13dcpforroot_IH.root
#rm contours_th13dcp_conf.pdf
#rm contours_th13dcpforroot.root
#root -l -q -b 'utils/macros/makeOfficialPlotStyle.C+("'$CONTOURDIR'/contours_th13dcpforroot_IH.root","'$CONTOURDIR'/contours_th13dcp_IH_official.root","T2K Run 1-10b preliminary","c")'

echo "GETTING 1D DCP CONTOURS:"
echo "  BOTH HIERARCHIES:"
#get 1D appearance contours
root -l -q -b 'utils/contours_app1D_EB.C+("'$FILENAME'",0,'${burnin}',false,false)'
mv contours_app1D.root $CONTOURDIR/contours_app1D_both.root
mv contours_1D_th13.pdf $CONTOURDIR/contours_1D_th13_both.pdf
mv contours_1D_dcp.pdf $CONTOURDIR/contours_1D_dcp_both.pdf
mv contours_1D_dcp.root $CONTOURDIR/contours_1D_dcp_both.root
mv contours_1D_dcp_sigmas.pdf $CONTOURDIR/contours_1D_dcp_sigmas_both.pdf
mv contours_1D_dcp_sigmas.root $CONTOURDIR/contours_1D_dcp_sigmas_both.root


echo "  NORMAL HIERARCHY:"
root -l -q -b 'utils/contours_app1D_EB.C+("'$FILENAME'",1,'${burnin}',false,false)'
mv contours_app1D.root $CONTOURDIR/contours_app1D_NH.root
mv contours_1D_th13.pdf $CONTOURDIR/contours_1D_th13_NH.pdf
mv contours_1D_dcp.pdf $CONTOURDIR/contours_1D_dcp_NH.pdf
mv contours_1D_dcp.root $CONTOURDIR/contours_1D_dcp_NH.root
mv contours_1D_dcp_sigmas.pdf $CONTOURDIR/contours_1D_dcp_sigmas_NH.pdf
mv contours_1D_dcp_sigmas.root $CONTOURDIR/contours_1D_dcp_sigmas_NH.root

echo "  INVERTED HIERARCHY:"
root -l -q -b 'utils/contours_app1D_EB.C+("'$FILENAME'",-1,'${burnin}',false,false)'
mv contours_app1D.root $CONTOURDIR/contours_app1D_IH.root
mv contours_1D_th13.pdf $CONTOURDIR/contours_1D_th13_IH.pdf
mv contours_1D_dcp.pdf $CONTOURDIR/contours_1D_dcp_IH.pdf
mv contours_1D_dcp.root $CONTOURDIR/contours_1D_dcp_IH.root
mv contours_1D_dcp_sigmas.pdf $CONTOURDIR/contours_1D_dcp_sigmas_IH.pdf
mv contours_1D_dcp_sigmas.root $CONTOURDIR/contours_1D_dcp_sigmas_IH.root

# LogY axis version:
root -l -q -b 'utils/contours_app1D_EB.C+("'$FILENAME'",0,'${burnin}',false,true)'
mv contours_app1D.root $CONTOURDIR/contours_app1D_LogY_both.root
mv contours_1D_th13.pdf $CONTOURDIR/contours_1D_th13_LogY_both.pdf
mv contours_1D_dcp.pdf $CONTOURDIR/contours_1D_dcp_LogY_both.pdf
mv contours_1D_dcp.root $CONTOURDIR/contours_1D_dcp_LogY_both.root
mv contours_1D_dcp_sigmas.pdf $CONTOURDIR/contours_1D_dcp_LogY_sigmas_both.pdf
mv contours_1D_dcp_sigmas.root $CONTOURDIR/contours_1D_dcp_LogY_sigmas_both.root

echo "  NORMAL HIERARCHY:"
root -l -q -b 'utils/contours_app1D_EB.C+("'$FILENAME'",1,'${burnin}',false,true)'
mv contours_app1D.root $CONTOURDIR/contours_app1D_LogY_NH.root
mv contours_1D_th13.pdf $CONTOURDIR/contours_1D_th13_LogY_NH.pdf
mv contours_1D_dcp.pdf $CONTOURDIR/contours_1D_dcp_LogY_NH.pdf
mv contours_1D_dcp.root $CONTOURDIR/contours_1D_dcp_LogY_NH.root
mv contours_1D_dcp_sigmas.pdf $CONTOURDIR/contours_1D_dcp_LogY_sigmas_NH.pdf
mv contours_1D_dcp_sigmas.root $CONTOURDIR/contours_1D_dcp_LogY_sigmas_NH.root

echo "  INVERTED HIERARCHY:"
root -l -q -b 'utils/contours_app1D_EB.C+("'$FILENAME'",-1,'${burnin}',false,true)'
mv contours_app1D.root $CONTOURDIR/contours_app1D_LogY_IH.root
mv contours_1D_th13.pdf $CONTOURDIR/contours_1D_th13_LogY_IH.pdf
mv contours_1D_dcp.pdf $CONTOURDIR/contours_1D_dcp_LogY_IH.pdf
mv contours_1D_dcp.root $CONTOURDIR/contours_1D_dcp_LogY_IH.root
mv contours_1D_dcp_sigmas.pdf $CONTOURDIR/contours_1D_dcp_LogY_sigmas_IH.pdf
mv contours_1D_dcp_sigmas.root $CONTOURDIR/contours_1D_dcp_LogY_sigmas_IH.root

echo "GETTING 1D CREDIBLE INTERVALS:"
echo "  BOTH HIERARCHIES:"
#get 1D appearance contours
root -l -q -b 'utils/contours_1D_EB.C+("'$FILENAME'",0,'${burnin}')'
mv contours_1D_th23.pdf $CONTOURDIR/contours_1D_th23_both.pdf
mv contours_1D_dm23.pdf $CONTOURDIR/contours_1D_dm23_both.pdf
mv contours_1D_th23.root $CONTOURDIR/contours_1D_th23_both.root
mv contours_1D_dm23.root $CONTOURDIR/contours_1D_dm23_both.root

echo "  NORMAL HIERARCHY:"
#get 1D appearance contours
root -l -q -b 'utils/contours_1D_EB.C+("'$FILENAME'",1,'${burnin}')'
mv contours_1D_th23.pdf $CONTOURDIR/contours_1D_th23_NH.pdf
mv contours_1D_dm23.pdf $CONTOURDIR/contours_1D_dm23_NH.pdf
mv contours_1D_th23.root $CONTOURDIR/contours_1D_th23_NH.root
mv contours_1D_dm23.root $CONTOURDIR/contours_1D_dm23_NH.root

echo "  INVERTED HIERARCHY:"
#get 1D appearance contours
root -l -q -b 'utils/contours_1D_EB.C+("'$FILENAME'",-1,'${burnin}')'
mv contours_1D_th23.pdf $CONTOURDIR/contours_1D_th23_IH.pdf
mv contours_1D_dm23.pdf $CONTOURDIR/contours_1D_dm23_IH.pdf
mv contours_1D_th23.root $CONTOURDIR/contours_1D_th23_IH.root
mv contours_1D_dm23.root $CONTOURDIR/contours_1D_dm23_IH.root




#get 2D th23 dcp contours
#echo "GETTING 2D th23 dcp CONTOURS:"
#echo "  BOTH HIERARCHIES:"
#root -l -q -b 'utils/contours_th23dcp.C+("'$FILENAME'",'$trueth23','$truedcp',-999,0,'${burnin}')'
#mv contours_th23dcp.root $CONTOURDIR/contours_th23dcp_both.root
#mv contours_th23dcp.pdf $CONTOURDIR/contours_th23dcp_both.pdf
#mv contours_th23dcp_conf.pdf $CONTOURDIR/contours_th23dcp_conf_both.pdf
#mv contours_th23dcpforroot.root $CONTOURDIR/contours_th23dcpforroot_both.root
#root -l -q -b 'utils/macros/makeOfficialPlotStyle.C+("'$CONTOURDIR'/contours_th23dcpforroot_both.root","'$CONTOURDIR'/contours_th23dcp_both_official.root","T2K Run 1-10b preliminary","c")'
#echo "  NORMAL HIERARCHY:"
#root -l -q -b 'utils/contours_th23dcp.C+("'$FILENAME'",0,0,-999,1,'${burnin}')'
#mv contours_th23dcp.root $CONTOURDIR/contours_th23dcp_NH.root
#mv contours_th23dcp.pdf $CONTOURDIR/contours_th23dcp_NH.pdf
#mv contours_th23dcp_conf.pdf $CONTOURDIR/contours_th23dcp_conf_NH.pdf
#mv contours_th23dcpforroot.root $CONTOURDIR/contours_th23dcpforroot_NH.root
#root -l -q -b 'utils/macros/makeOfficialPlotStyle.C+("'$CONTOURDIR'/contours_th23dcpforroot_NH.root","'$CONTOURDIR'/contours_th23dcp_NH_official.root","T2K Run 1-10b preliminary","c")'
#echo "  INVERTED HIERARCHY:"
#root -l -q -b 'utils/contours_th23dcp.C+("'$FILENAME'",0,0,-999,-1,'${burnin}')'
#mv contours_th23dcp.root $CONTOURDIR/contours_th23dcp_IH.root
#mv contours_th23dcp.pdf $CONTOURDIR/contours_th23dcp_IH.pdf
#mv contours_th23dcp_conf.pdf $CONTOURDIR/contours_th23dcp_conf_IH.pdf
#mv contours_th23dcpforroot.root $CONTOURDIR/contours_th23dcpforroot_IH.root
#root -l -q -b 'utils/macros/makeOfficialPlotStyle.C+("'$CONTOURDIR'/contours_th23dcpforroot_IH.root","'$CONTOURDIR'/contours_th23dcp_IH_official.root","T2K Run 1-10b preliminary","c")'




#root -l -q -b 'utils/macros/makeOfficialPlotStyle.C+("'$CONTOURDIR'/contours_1D_dcp_NH.root","'$CONTOURDIR'/contours_1D_dcp_NH_official.root","T2K Run 1-10b preliminary","c1")'
#root -l -q -b 'utils/macros/makeOfficialPlotStyle.C+("'$CONTOURDIR'/contours_1D_th13_NH.root","'$CONTOURDIR'/contours_1D_th13_NH_official.root","T2K Run 1-10b preliminary","c")'

#root -l -q -b 'utils/macros/makeOfficialPlotStyle.C+("'$CONTOURDIR'/contours_1D_dcp_both.root","'$CONTOURDIR'/contours_1D_dcp_both_official.root","T2K Run 1-10b preliminary","c1")'
#root -l -q -b 'utils/macros/makeOfficialPlotStyle.C+("'$CONTOURDIR'/contours_1D_th13_both.root","'$CONTOURDIR'/contours_1D_th13_both_official.root","T2K Run 1-10b preliminary","c")'

#root -l -q -b 'utils/macros/makeOfficialPlotStyle.C+("'$CONTOURDIR'/contours_1D_dcp_IH.root","'$CONTOURDIR'/contours_1D_dcp_IH_official.root","T2K Run 1-10b preliminary","c1")'
#root -l -q -b 'utils/macros/makeOfficialPlotStyle.C+("'$CONTOURDIR'/contours_1D_th13_IH.root","'$CONTOURDIR'/contours_1D_th13_IH_official.root","T2K Run 1-10b preliminary","c")'
#
##get 1D disappearance contours
#root -l -q -b 'utils/contours_1D.C+("'$FILENAME'",0)'
#mv contours_1D.root $CONTOURDIR/contours_1D_both.root
#mv contours_1D_th23.pdf $CONTOURDIR/contours_1D_th23_both.pdf
#mv contours_1D_dm23.pdf $CONTOURDIR/contours_1D_dm23_both.pdf
#mv contours_1D_th23.root $CONTOURDIR/contours_1D_th23_both.root
#mv contours_1D_dm23.root $CONTOURDIR/contours_1D_dm23_both.root
#root -l -q -b 'utils/macros/makeOfficialPlotStyle.C+("'$CONTOURDIR'/contours_1D_dm23_both.root","'$CONTOURDIR'/contours_1D_dm23_both_official.root","T2K Run 1-10b preliminary","c1")'
#root -l -q -b 'utils/macros/makeOfficialPlotStyle.C+("'$CONTOURDIR'/contours_1D_th23_both.root","'$CONTOURDIR'/contours_1D_th23_both_official.root","T2K Run 1-10b preliminary","c")'
#
##!!
#root -l -q -b 'utils/contours_1D.C+("'$FILENAME'",1)'
#mv contours_1D.root $CONTOURDIR/contours_1D_NH.root
#mv contours_1D_th23.pdf $CONTOURDIR/contours_1D_th23_NH.pdf
#mv contours_1D_dm23.pdf $CONTOURDIR/contours_1D_dm23_NH.pdf
#mv contours_1D_th23.root $CONTOURDIR/contours_1D_th23_NH.root
#mv contours_1D_dm23.root $CONTOURDIR/contours_1D_dm23_NH.root
#root -l -q -b 'utils/macros/makeOfficialPlotStyle.C+("'$CONTOURDIR'/contours_1D_dm23_NH.root","'$CONTOURDIR'/contours_1D_dm23_NH_official.root","T2K Run 1-10b preliminary","c1")'
#root -l -q -b 'utils/macros/makeOfficialPlotStyle.C+("'$CONTOURDIR'/contours_1D_th23_NH.root","'$CONTOURDIR'/contours_1D_th23_NH_official.root","T2K Run 1-10b preliminary","c")'
#
#root -l -q -b 'utils/contours_1D.C+("'$FILENAME'",-1)'
#mv contours_1D.root $CONTOURDIR/contours_1D_IH.root
#mv contours_1D_th23.pdf $CONTOURDIR/contours_1D_th23_IH.pdf
#mv contours_1D_dm23.pdf $CONTOURDIR/contours_1D_dm23_IH.pdf
#mv contours_1D_th23.root $CONTOURDIR/contours_1D_th23_IH.root
#mv contours_1D_dm23.root $CONTOURDIR/contours_1D_dm23_IH.root
#root -l -q -b 'utils/macros/makeOfficialPlotStyle.C+("'$CONTOURDIR'/contours_1D_dm23_IH.root","'$CONTOURDIR'/contours_1D_dm23_IH_official.root","T2K Run 1-10b preliminary","c1")'
#root -l -q -b 'utils/macros/makeOfficialPlotStyle.C+("'$CONTOURDIR'/contours_1D_th23_IH.root","'$CONTOURDIR'/contours_1D_th23_IH_official.root","T2K Run 1-10b preliminary","c")'
#
#
#
##plot with arrows for sindcp prior overlaid
#root -l -q -b 'utils/Overlay_CI_1Ddcp.C+("'$CONTOURDIR'/contours_1D_dcp_both.root",'$rightarr68','$leftarr68','$rightarr90','$leftarr90','$rightarr95','$leftarr95')'
#mv contours_1D_dcp_otherprior.pdf $CONTOURDIR/contours_1D_dcp_otherprior_both.pdf
#mv contours_1D_dcp_otherprior.png $CONTOURDIR/contours_1D_dcp_otherprior_both.png
#mv contours_1D_dcp_otherprior.root $CONTOURDIR/contours_1D_dcp_otherprior_both.root
#
#root -l -q -b 'utils/macros/makeOfficialPlotStyle.C+("'$CONTOURDIR'/contours_1D_dcp_otherprior_both.root","'$CONTOURDIR'/contours_1D_dcp_otherprior_official_both.root","T2K Run 1-10b preliminary","c1")'
#
##sindcp axis
#root -l -q -b 'utils/contours_app1D.C+("'$FILENAME'",0,true)'
#mv contours_app1D.root $CONTOURDIR/contours_app1D_sindcpaxis_both.root
#mv contours_1D_dcp.pdf $CONTOURDIR/contours_1D_dcp_sindcpaxis_both.pdf
#mv contours_1D_dcp.root $CONTOURDIR/contours_1D_dcp_sindcpaxis_both.root
#root -l -q -b 'utils/macros/makeOfficialPlotStyle.C+("'$CONTOURDIR'/contours_1D_dcp_sindcpaxis_both.root","'$CONTOURDIR'/contours_1D_dcp_sindcpaxis_both_official.root","T2K Run 1-10b preliminary","c1")'
#
#root -l -q -b 'utils/contours_app1D.C+("'$FILENAME'",1,true)'
#mv contours_app1D.root $CONTOURDIR/contours_app1D_sindcpaxis_NH.root
#mv contours_1D_dcp.pdf $CONTOURDIR/contours_1D_dcp_sindcpaxis_NH.pdf
#mv contours_1D_dcp.root $CONTOURDIR/contours_1D_dcp_sindcpaxis_NH.root
#root -l -q -b 'utils/macros/makeOfficialPlotStyle.C+("'$CONTOURDIR'/contours_1D_dcp_sindcpaxis_NH.root","'$CONTOURDIR'/contours_1D_dcp_sindcpaxis_NH_official.root","T2K Run 1-10b preliminary","c1")'
#root -l -q -b 'utils/contours_app1D.C+("'$FILENAME'",-1,true)'
#mv contours_app1D.root $CONTOURDIR/contours_app1D_sindcpaxis_IH.root
#mv contours_1D_dcp.pdf $CONTOURDIR/contours_1D_dcp_sindcpaxis_IH.pdf
#mv contours_1D_dcp.root $CONTOURDIR/contours_1D_dcp_sindcpaxis_IH.root
#root -l -q -b 'utils/macros/makeOfficialPlotStyle.C+("'$CONTOURDIR'/contours_1D_dcp_sindcpaxis_IH.root","'$CONTOURDIR'/contours_1D_dcp_sindcpaxis_IH_official.root","T2K Run 1-10b preliminary","c1")'
#
##sindcp axis and prior
#root -l -q -b 'utils/contours_app1D.C+("'$FILENAME'",0,true,true)'
#mv contours_app1D.root $CONTOURDIR/contours_app1D_sindcpaxisandprior_both.root
#mv contours_1D_dcp.pdf $CONTOURDIR/contours_1D_dcp_sindcpaxisandprior_both.pdf
#mv contours_1D_dcp.root $CONTOURDIR/contours_1D_dcp_sindcpaxisandprior_both.root
#root -l -q -b 'utils/macros/makeOfficialPlotStyle.C+("'$CONTOURDIR'/contours_1D_dcp_sindcpaxisandprior_both.root","'$CONTOURDIR'/contours_1D_dcp_sindcpaxisandprior_both_official.root","T2K Run 1-10b preliminary","c1")'
#
#root -l -q -b 'utils/contours_app1D.C+("'$FILENAME'",1,true,true)'
#mv contours_app1D.root $CONTOURDIR/contours_app1D_sindcpaxisandprior_NH.root
#mv contours_1D_dcp.pdf $CONTOURDIR/contours_1D_dcp_sindcpaxisandprior_NH.pdf
#mv contours_1D_dcp.root $CONTOURDIR/contours_1D_dcp_sindcpaxisandprior_NH.root
#root -l -q -b 'utils/macros/makeOfficialPlotStyle.C+("'$CONTOURDIR'/contours_1D_dcp_sindcpaxisandprior_NH.root","'$CONTOURDIR'/contours_1D_dcp_sindcpaxisandprior_NH_official.root","T2K Run 1-10b preliminary","c1")'
#root -l -q -b 'utils/contours_app1D.C+("'$FILENAME'",-1,true,true)'
#mv contours_app1D.root $CONTOURDIR/contours_app1D_sindcpaxisandprior_IH.root
#mv contours_1D_dcp.pdf $CONTOURDIR/contours_1D_dcp_sindcpaxisandprior_IH.pdf
#mv contours_1D_dcp.root $CONTOURDIR/contours_1D_dcp_sindcpaxisandprior_IH.root
#root -l -q -b 'utils/macros/makeOfficialPlotStyle.C+("'$CONTOURDIR'/contours_1D_dcp_sindcpaxisandprior_IH.root","'$CONTOURDIR'/contours_1D_dcp_sindcpaxisandprior_IH_official.root","T2K Run 1-10b preliminary","c1")'
#
##sindcp axis and prior
#root -l -q -b 'utils/contours_app1D.C+("'$FILENAME'",0,false,true)'
#mv contours_app1D.root $CONTOURDIR/contours_app1D_sindcpprior_both.root
#mv contours_1D_dcp.pdf $CONTOURDIR/contours_1D_dcp_sindcpprior_both.pdf
#mv contours_1D_dcp.root $CONTOURDIR/contours_1D_dcp_sindcpprior_both.root
#root -l -q -b 'utils/macros/makeOfficialPlotStyle.C+("'$CONTOURDIR'/contours_1D_dcp_sindcpprior_both.root","'$CONTOURDIR'/contours_1D_dcp_sindcpprior_both_official.root","T2K Run 1-10b preliminary","c1")'
#
#root -l -q -b 'utils/contours_app1D.C+("'$FILENAME'",1,false,true)'
#mv contours_app1D.root $CONTOURDIR/contours_app1D_sindcpprior_NH.root
#mv contours_1D_dcp.pdf $CONTOURDIR/contours_1D_dcp_sindcpprior_NH.pdf
#mv contours_1D_dcp.root $CONTOURDIR/contours_1D_dcp_sindcpprior_NH.root
#root -l -q -b 'utils/macros/makeOfficialPlotStyle.C+("'$CONTOURDIR'/contours_1D_dcp_sindcpprior_NH.root","'$CONTOURDIR'/contours_1D_dcp_sindcpprior_NH_official.root","T2K Run 1-10b preliminary","c1")'
#root -l -q -b 'utils/contours_app1D.C+("'$FILENAME'",-1,false,true)'
#mv contours_app1D.root $CONTOURDIR/contours_app1D_sindcpprior_IH.root
#mv contours_1D_dcp.pdf $CONTOURDIR/contours_1D_dcp_sindcpprior_IH.pdf
#mv contours_1D_dcp.root $CONTOURDIR/contours_1D_dcp_sindcpprior_IH.root
#root -l -q -b 'utils/macros/makeOfficialPlotStyle.C+("'$CONTOURDIR'/contours_1D_dcp_sindcpprior_IH.root","'$CONTOURDIR'/contours_1D_dcp_sindcpprior_IH_official.root","T2K Run 1-10b preliminary","c1")'
done
done

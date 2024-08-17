export DAT=/media/mayukh/data1/cov-test-with-swiftj-uttley/0311590901/ODF
#export DAT=/media/mayukh/data1/cov-test-with-ngc4593/data/2002/0059830101/ODF

cd $DAT

rm -rf *.cif *.SAS

export SAS_ODF=$DAT

gunzip *.gz

cifbuild

export SAS_CCF=$DAT/ccf.cif

odfingest

cp *.SAS copyodf.SAS

export SAS_ODF=$DAT/copyodf.SAS

rm -rf lightcurve*set/

read -p "What is the binsize of the lightcurve in sec ? " n

evselect table=pn_filt_src.fits withfilteredset=yes expression='(PATTERN <= 4)&&(PI in [300:500])&&#XMMEA_EP&&(RAWX in [27:47])' filteredset=pn_filt_src_pt3_pt5.fits filtertype=expression keepfilteroutput=yes updateexposure=yes filterexposure=yes

evselect table=pn_filt_src_pt3_pt5.fits withrateset=yes rateset=en1.fits maketimecolumn=yes timecolumn=TIME timebinsize="$n" makeratecolumn=yes

evselect table=pn_filt_src.fits withfilteredset=yes expression='(PATTERN <= 4)&&(PI in [500:700])&&#XMMEA_EP&&(RAWX in [27:47])' filteredset=pn_filt_src_pt5_pt7.fits filtertype=expression keepfilteroutput=yes updateexposure=yes filterexposure=yes

evselect table=pn_filt_src_pt5_pt7.fits withrateset=yes rateset=en2.fits maketimecolumn=yes timecolumn=TIME timebinsize="$n" makeratecolumn=yes

evselect table=pn_filt_src.fits withfilteredset=yes expression='(PATTERN <= 4)&&(PI in [700:900])&&#XMMEA_EP&&(RAWX in [27:47])' filteredset=pn_filt_src_pt7_pt9.fits filtertype=expression keepfilteroutput=yes updateexposure=yes filterexposure=yes

evselect table=pn_filt_src_pt7_pt9.fits withrateset=yes rateset=en3.fits maketimecolumn=yes timecolumn=TIME timebinsize="$n" makeratecolumn=yes

evselect table=pn_filt_src.fits withfilteredset=yes expression='(PATTERN <= 4)&&(PI in [900:1100])&&#XMMEA_EP&&(RAWX in [27:47])' filteredset=pn_filt_src_pt9_1pt1.fits filtertype=expression keepfilteroutput=yes updateexposure=yes filterexposure=yes

evselect table=pn_filt_src_pt9_1pt1.fits withrateset=yes rateset=en4.fits maketimecolumn=yes timecolumn=TIME timebinsize="$n" makeratecolumn=yes

evselect table=pn_filt_src.fits withfilteredset=yes expression='(PATTERN <= 4)&&(PI in [1100:1400])&&#XMMEA_EP&&(RAWX in [27:47])' filteredset=pn_filt_src_1pt1_1pt4.fits filtertype=expression keepfilteroutput=yes updateexposure=yes filterexposure=yes

evselect table=pn_filt_src_1pt1_1pt4.fits withrateset=yes rateset=en5.fits maketimecolumn=yes timecolumn=TIME timebinsize="$n" makeratecolumn=yes

evselect table=pn_filt_src.fits withfilteredset=yes expression='(PATTERN <= 4)&&(PI in [1400:1700])&&#XMMEA_EP&&(RAWX in [27:47])' filteredset=pn_filt_src_1pt4_1pt7.fits filtertype=expression keepfilteroutput=yes updateexposure=yes filterexposure=yes

evselect table=pn_filt_src_1pt4_1pt7.fits withrateset=yes rateset=en6.fits maketimecolumn=yes timecolumn=TIME timebinsize="$n" makeratecolumn=yes

evselect table=pn_filt_src.fits withfilteredset=yes expression='(PATTERN <= 4)&&(PI in [1700:2000])&&#XMMEA_EP&&(RAWX in [27:47])' filteredset=pn_filt_src_1pt7_2.fits filtertype=expression keepfilteroutput=yes updateexposure=yes filterexposure=yes

evselect table=pn_filt_src_1pt7_2.fits withrateset=yes rateset=en7.fits maketimecolumn=yes timecolumn=TIME timebinsize="$n" makeratecolumn=yes

evselect table=pn_filt_src.fits withfilteredset=yes expression='(PATTERN <= 4)&&(PI in [2000:2500])&&#XMMEA_EP&&(RAWX in [27:47])' filteredset=pn_filt_src_2_2pt5.fits filtertype=expression keepfilteroutput=yes updateexposure=yes filterexposure=yes

evselect table=pn_filt_src_2_2pt5.fits withrateset=yes rateset=en8.fits maketimecolumn=yes timecolumn=TIME timebinsize="$n" makeratecolumn=yes

evselect table=pn_filt_src.fits withfilteredset=yes expression='(PATTERN <= 4)&&(PI in [2500:3000])&&#XMMEA_EP&&(RAWX in [27:47])' filteredset=pn_filt_src_2pt5_3.fits filtertype=expression keepfilteroutput=yes updateexposure=yes filterexposure=yes

evselect table=pn_filt_src_2pt5_3.fits withrateset=yes rateset=en9.fits maketimecolumn=yes timecolumn=TIME timebinsize="$n" makeratecolumn=yes

evselect table=pn_filt_src.fits withfilteredset=yes expression='(PATTERN <= 4)&&(PI in [3000:4000])&&#XMMEA_EP&&(RAWX in [27:47])' filteredset=pn_filt_src_3_4.fits filtertype=expression keepfilteroutput=yes updateexposure=yes filterexposure=yes

evselect table=pn_filt_src_3_4.fits withrateset=yes rateset=en10.fits maketimecolumn=yes timecolumn=TIME timebinsize="$n" makeratecolumn=yes

evselect table=pn_filt_src.fits withfilteredset=yes expression='(PATTERN <= 4)&&(PI in [4000:5000])&&#XMMEA_EP&&(RAWX in [27:47])' filteredset=pn_filt_src_4_5.fits filtertype=expression keepfilteroutput=yes updateexposure=yes filterexposure=yes

evselect table=pn_filt_src_4_5.fits withrateset=yes rateset=en11.fits maketimecolumn=yes timecolumn=TIME timebinsize="$n" makeratecolumn=yes

evselect table=pn_filt_src.fits withfilteredset=yes expression='(PATTERN <= 4)&&(PI in [5000:7000])&&#XMMEA_EP&&(RAWX in [27:47])' filteredset=pn_filt_src_5_7.fits filtertype=expression keepfilteroutput=yes updateexposure=yes filterexposure=yes

evselect table=pn_filt_src_5_7.fits withrateset=yes rateset=en12.fits maketimecolumn=yes timecolumn=TIME timebinsize="$n" makeratecolumn=yes

evselect table=pn_filt_src.fits withfilteredset=yes expression='(PATTERN <= 4)&&(PI in [7000:9000])&&#XMMEA_EP&&(RAWX in [27:47])' filteredset=pn_filt_src_7_9.fits filtertype=expression keepfilteroutput=yes updateexposure=yes filterexposure=yes

evselect table=pn_filt_src_7_9.fits withrateset=yes rateset=en13.fits maketimecolumn=yes timecolumn=TIME timebinsize="$n" makeratecolumn=yes

evselect table=pn_filt_src.fits withfilteredset=yes expression='(PATTERN <= 4)&&(PI in [1000:4000])&&#XMMEA_EP&&(RAWX in [27:47])' filteredset=pn_filt_src_1_4.fits filtertype=expression keepfilteroutput=yes updateexposure=yes filterexposure=yes

evselect table=pn_filt_src_1_4.fits withrateset=yes rateset=en111.fits maketimecolumn=yes timecolumn=TIME timebinsize="$n" makeratecolumn=yes


### BACKGROUND LIGHTCURVE GENERATION AND SUBTRACTION ########

evselect table=pn_filt_bkg.fits withfilteredset=yes expression='(PATTERN <= 4)&&(PI in [300:500])&&#XMMEA_EP&&(RAWX in [3:5])' filteredset=pn_filt_bkg_pt3_pt5.fits filtertype=expression keepfilteroutput=yes updateexposure=yes filterexposure=yes

evselect table=pn_filt_bkg_pt3_pt5.fits withrateset=yes rateset=back_en1.fits maketimecolumn=yes timecolumn=TIME timebinsize="$n" makeratecolumn=yes

lcmath infile=en1.fits bgfile=back_en1.fits outfile=en-sub1.fits multi=1. multb=1. addsubr=no

evselect table=pn_filt_bkg.fits withfilteredset=yes expression='(PATTERN <= 4)&&(PI in [500:700])&&#XMMEA_EP&&(RAWX in [3:5])' filteredset=pn_filt_bkg_pt5_pt7.fits filtertype=expression keepfilteroutput=yes updateexposure=yes filterexposure=yes

evselect table=pn_filt_bkg_pt5_pt7.fits withrateset=yes rateset=back_en2.fits maketimecolumn=yes timecolumn=TIME timebinsize="$n" makeratecolumn=yes

lcmath infile=en2.fits bgfile=back_en2.fits outfile=en-sub2.fits multi=1. multb=1. addsubr=no

evselect table=pn_filt_bkg.fits withfilteredset=yes expression='(PATTERN <= 4)&&(PI in [700:900])&&#XMMEA_EP&&(RAWX in [3:5])' filteredset=pn_filt_bkg_pt7_pt9.fits filtertype=expression keepfilteroutput=yes updateexposure=yes filterexposure=yes

evselect table=pn_filt_bkg_pt7_pt9.fits withrateset=yes rateset=back_en3.fits maketimecolumn=yes timecolumn=TIME timebinsize="$n" makeratecolumn=yes

lcmath infile=en3.fits bgfile=back_en3.fits outfile=en-sub3.fits multi=1. multb=1. addsubr=no

evselect table=pn_filt_bkg.fits withfilteredset=yes expression='(PATTERN <= 4)&&(PI in [900:1100])&&#XMMEA_EP&&(RAWX in [3:5])' filteredset=pn_filt_bkg_pt9_1pt1.fits filtertype=expression keepfilteroutput=yes updateexposure=yes filterexposure=yes

evselect table=pn_filt_bkg_pt9_1pt1.fits withrateset=yes rateset=back_en4.fits maketimecolumn=yes timecolumn=TIME timebinsize="$n" makeratecolumn=yes

lcmath infile=en4.fits bgfile=back_en4.fits outfile=en-sub4.fits multi=1. multb=1. addsubr=no

evselect table=pn_filt_bkg.fits withfilteredset=yes expression='(PATTERN <= 4)&&(PI in [1100:1400])&&#XMMEA_EP&&(RAWX in [3:5])' filteredset=pn_filt_bkg_1pt1_1pt4.fits filtertype=expression keepfilteroutput=yes updateexposure=yes filterexposure=yes

evselect table=pn_filt_bkg_1pt1_1pt4.fits withrateset=yes rateset=back_en5.fits maketimecolumn=yes timecolumn=TIME timebinsize="$n" makeratecolumn=yes

lcmath infile=en5.fits bgfile=back_en5.fits outfile=en-sub5.fits multi=1. multb=1. addsubr=no

evselect table=pn_filt_bkg.fits withfilteredset=yes expression='(PATTERN <= 4)&&(PI in [1400:1700])&&#XMMEA_EP&&(RAWX in [3:5])' filteredset=pn_filt_bkg_1pt4_1pt7.fits filtertype=expression keepfilteroutput=yes updateexposure=yes filterexposure=yes

evselect table=pn_filt_bkg_1pt4_1pt7.fits withrateset=yes rateset=back_en6.fits maketimecolumn=yes timecolumn=TIME timebinsize="$n" makeratecolumn=yes

lcmath infile=en6.fits bgfile=back_en6.fits outfile=en-sub6.fits multi=1. multb=1. addsubr=no

evselect table=pn_filt_bkg.fits withfilteredset=yes expression='(PATTERN <= 4)&&(PI in [1700:2000])&&#XMMEA_EP&&(RAWX in [3:5])' filteredset=pn_filt_bkg_1pt7_2.fits filtertype=expression keepfilteroutput=yes updateexposure=yes filterexposure=yes

evselect table=pn_filt_bkg_1pt7_2.fits withrateset=yes rateset=back_en7.fits maketimecolumn=yes timecolumn=TIME timebinsize="$n" makeratecolumn=yes

lcmath infile=en7.fits bgfile=back_en7.fits outfile=en-sub7.fits multi=1. multb=1. addsubr=no

evselect table=pn_filt_bkg.fits withfilteredset=yes expression='(PATTERN <= 4)&&(PI in [2000:2500])&&#XMMEA_EP&&(RAWX in [3:5])' filteredset=pn_filt_bkg_2_2pt5.fits filtertype=expression keepfilteroutput=yes updateexposure=yes filterexposure=yes

evselect table=pn_filt_bkg_2_2pt5.fits withrateset=yes rateset=back_en8.fits maketimecolumn=yes timecolumn=TIME timebinsize="$n" makeratecolumn=yes

lcmath infile=en8.fits bgfile=back_en8.fits outfile=en-sub8.fits multi=1. multb=1. addsubr=no

evselect table=pn_filt_bkg.fits withfilteredset=yes expression='(PATTERN <= 4)&&(PI in [2500:3000])&&#XMMEA_EP&&(RAWX in [3:5])' filteredset=pn_filt_bkg_2pt5_3.fits filtertype=expression keepfilteroutput=yes updateexposure=yes filterexposure=yes

evselect table=pn_filt_bkg_2pt5_3.fits withrateset=yes rateset=back_en9.fits maketimecolumn=yes timecolumn=TIME timebinsize="$n" makeratecolumn=yes

lcmath infile=en9.fits bgfile=back_en9.fits outfile=en-sub9.fits multi=1. multb=1. addsubr=no

evselect table=pn_filt_bkg.fits withfilteredset=yes expression='(PATTERN <= 4)&&(PI in [3000:4000])&&#XMMEA_EP&&(RAWX in [3:5])' filteredset=pn_filt_bkg_3_4.fits filtertype=expression keepfilteroutput=yes updateexposure=yes filterexposure=yes

evselect table=pn_filt_bkg_3_4.fits withrateset=yes rateset=back_en10.fits maketimecolumn=yes timecolumn=TIME timebinsize="$n" makeratecolumn=yes

lcmath infile=en10.fits bgfile=back_en10.fits outfile=en-sub10.fits multi=1. multb=1. addsubr=no

evselect table=pn_filt_bkg.fits withfilteredset=yes expression='(PATTERN <= 4)&&(PI in [4000:5000])&&#XMMEA_EP&&(RAWX in [3:5])' filteredset=pn_filt_bkg_4_5.fits filtertype=expression keepfilteroutput=yes updateexposure=yes filterexposure=yes

evselect table=pn_filt_bkg_4_5.fits withrateset=yes rateset=back_en11.fits maketimecolumn=yes timecolumn=TIME timebinsize="$n" makeratecolumn=yes

lcmath infile=en11.fits bgfile=back_en11.fits outfile=en-sub11.fits multi=1. multb=1. addsubr=no

evselect table=pn_filt_bkg.fits withfilteredset=yes expression='(PATTERN <= 4)&&(PI in [5000:7000])&&#XMMEA_EP&&(RAWX in [3:5])' filteredset=pn_filt_bkg_5_7.fits filtertype=expression keepfilteroutput=yes updateexposure=yes filterexposure=yes

evselect table=pn_filt_bkg_5_7.fits withrateset=yes rateset=back_en12.fits maketimecolumn=yes timecolumn=TIME timebinsize="$n" makeratecolumn=yes

lcmath infile=en12.fits bgfile=back_en12.fits outfile=en-sub12.fits multi=1. multb=1. addsubr=no

evselect table=pn_filt_bkg.fits withfilteredset=yes expression='(PATTERN <= 4)&&(PI in [7000:9000])&&#XMMEA_EP&&(RAWX in [3:5])' filteredset=pn_filt_bkg_7_9.fits filtertype=expression keepfilteroutput=yes updateexposure=yes filterexposure=yes

evselect table=pn_filt_bkg_7_9.fits withrateset=yes rateset=back_en13.fits maketimecolumn=yes timecolumn=TIME timebinsize="$n" makeratecolumn=yes

lcmath infile=en13.fits bgfile=back_en13.fits outfile=en-sub13.fits multi=1. multb=1. addsubr=no

#########   REFERENCE BAND BACKGROUND LIGHTCURVE EXTRACTION AND SUBTRACTION   ######################

evselect table=pn_filt_bkg.fits withfilteredset=yes expression='(PATTERN <= 4)&&(PI in [1000:4000])&&#XMMEA_EP&&(RAWX in [3:5])' filteredset=pn_filt_bkg_1_4.fits filtertype=expression keepfilteroutput=yes updateexposure=yes filterexposure=yes

evselect table=pn_filt_bkg_1_4.fits withrateset=yes rateset=back_ref.fits maketimecolumn=yes timecolumn=TIME timebinsize="$n" makeratecolumn=yes

lcmath infile=en111.fits bgfile=back_ref.fits outfile=en-sub-ref.fits multi=1. multb=1. addsubr=no


mkdir lightcurve-"$n"sec-set
cp en-sub*.fits lightcurve-"$n"sec-set/
cp /home/mayukh/covar_specV1.py lightcurve-"$n"sec-set/
cd lightcurve-"$n"sec-set/
pwd
python3 covar_specV3.py

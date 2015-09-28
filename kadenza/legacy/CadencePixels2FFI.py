#!/usr/bin/env python

import PyMSG, PyIO
import getopt, sys, os, re
import datetime
import time as ti
import numpy, pyfits
from pyfits import *
from numpy import *

def main():

# input argument

    cpf = '/Volumes/data/Kepler/SecondLight/CPF/kplr2013304171549_lcs-targ.fits'
    pmrf = '/Volumes/data/Kepler/PMRF/oct2013_lcm.fits'
    try:
	opts, args = getopt.getopt(sys.argv[1:],"h:cpa",["help","cpf=","pmrf=","cal"])
    except getopt.GetoptError:
	usage()
    cal = False
    for o, a in opts:
	if o in ("-h", "--help"):
	    usage()
	if o in ("-c", "--cpf"):
            cpf = str(a)
	if o in ("-p", "--pmrf"):
            pmrf = str(a)
	if o in ("-a", "--cal"):
            cal = True

# output FFI filename

    if cal:
        outname = re.sub('_lcs-targ','_ffi-cal',cpf)
    else:
        outname = re.sub('_lcs-targ','_ffi-raw',cpf)

# creation date
            
    createdate = str(datetime.date.today())

# calibrated data?

    if cal:
        imagtype = 'SocCal'
        bunit = 'electrons/s'
    else:
        imagtype = 'raw'
        bunit = 'counts'

# create primary FITS extensions

    hdu0 = PrimaryHDU()
    hdu0.header.update('EXTEND',True,'file may contain standard extensions')
    hdu0.header.update('EXTNAME','PRIMARY','name of extension')
    hdu0.header.update('EXTVER',1,'extension version number')
    hdu0.header.update('ORIGIN','NASA/Ames','organization that generated this file')
    hdu0.header.update('DATE',createdate,'file creation date')
    hdu0.header.update('CREATOR','CadencePixels2FFI.py','SW used to create this file')
    hdu0.header.update('FILEVER','2.0','file format version')
    hdu0.header.update('PROCVER',1.0,'processing script version')
    hdu0.header.update('TIMVERSN','OGIP/93-003','OGIP memo number for file format')
    hdu0.header.update('TELESCOP','Kepler','telescope')
    hdu0.header.update('INSTRUME','Kepler photometer','detector type')
    hdu0.header.update('DATA_REL',-1,'version of data release notes describing data')
    hdu0.header.update('OBSMODE','full frame image','observing mode')
    hdu0.header.update('DATSETNM',os.path.basename(re.sub('_lcs-targ.fits','',cpf)),'data set name')
    hdu0.header.update('DCT_TIME',os.path.basename(re.sub('kplr','',re.sub('_lcs-targ.fits','',cpf))),
                       'data collection time: yyyydddhhss')
    hdu0.header.update('DCT_TYPE','FFI','data type')
    hdu0.header.update('DCT_PURP','cadence image','purpose of data')
    hdu0.header.update('IMAGTYPE',imagtype,'FFI image type: raw, SocCal, SocUnc')
    hdu0.header.update('QUARTER',-1,'mission quarter during which data was collected')
    hdu0.header.update('SEASON',-1,'mission season during which data was collected')
    hdu0.header.update('PIXELTYP','target','pixel type: target, bkgd, coll, all')
    hdu0.header.update('VSMRSROW',1046,'collateral virtual smear region start row')
    hdu0.header.update('VSMREROW',1057,'collateral virtual smear region row end')
    hdu0.header.update('NROWVSMR',12,'number of rows binned in virtual smear')
    hdu0.header.update('VSMRSCOL',12,'collateral virtual smear region start column')
    hdu0.header.update('VSMRECOL',1111,'collateral virtual smear region end column')
    hdu0.header.update('NCOLVSMR',1100,'number of columns in virtual smear region')
    hdu0.header.update('MASKSROW',6,'science collateral masked region start row')
    hdu0.header.update('MASKEROW',17,'science collateral masked region end row')
    hdu0.header.update('NROWMASK',12,'number of rows binned in masked region')
    hdu0.header.update('MASKSCOL',12,'science collateral masked region start')
    hdu0.header.update('MASKECOL',1111,'science collateral masked region end')
    hdu0.header.update('NCOLMASK',1100,'number of columns in masked region')
    hdu0.header.update('BLCKSROW',0,'science collateral black region start row')
    hdu0.header.update('BLCKEROW',1069,'science collateral black region end column')
    hdu0.header.update('NROWBLCK',1070,'number of rows in black region')
    hdu0.header.update('BLCKSCOL',1118,'science collateral black region start')
    hdu0.header.update('BLCKECOL',1131,'science collateral black region end column')
    hdu0.header.update('NCOLBLK',1070,'number of columns binned in black region')
    hdu0.header.update('RADESYS','ICRS','reference frame of celestial coordinates')
    hdu0.header.update('EQUINOX',2000.0,'equinox of celestial coordinate system')
    hdu0.header.update('RA_NOM',-1.0,'[deg] RA of spacecraft boresight')
    hdu0.header.update('DEC_NOM',-1.0,'[deg] declination of spacecraft boresight')
    hdu0.header.update('ROLL_NOM',-1.0,'[deg] roll angle of spacecraft')
    hdulist = HDUList(hdu0)
             
# open the Pixel Mapping Reference File

    try:
        inpmrf = open(pmrf,mode='readonly',memmap=True)
    except:
        message =  'ERROR -- Cannot read the pixel reference file'
        status = PyMSG.err(None,message)

# open the Cadence Pixel File

    try:
        incpf = open(cpf,mode='readonly',memmap=True)
    except:
        message =  'ERROR -- Cannot read the cadence pixel file'
        status = PyMSG.err(None,message)

# start and stop times of observations from the CPF

    card = incpf[0].header.ascardlist()
    startkey = card['STARTIME'].value
    endkey = card['END_TIME'].value
    dateobs = card['DATE-OBS'].value
    timeobs = card['TIME-OBS'].value
    tstart = startkey + 2400000.5
    tstop = endkey + 2400000.5
    deadc = 0.92063492
    telapse = tstop - tstart
    exposure = telapse * deadc

# iterate through the FITS extensions of the PMRF and read the data 

    for ext in range(1,85):
        start = ti.time()
        card = incpf[ext].header.ascardlist()
        gain = card['GAIN'].value
        readnois = card['READNOIS'].value
        timslice = card['TIMSLICE'].value
        meanblck = card['MEANBLCK'].value
        card = inpmrf[ext].header.ascardlist()
        chankey = card['CHANNEL'].value
        modkey = card['MODULE'].value
        outkey = card['OUTPUT'].value
        r = inpmrf[ext].data.field('row')[:]
        c = inpmrf[ext].data.field('column')[:]
        k = inpmrf[ext].data.field('target_id')[:]
        if cal:
            f = incpf[ext].data.field('cal_value')[:]
        else:
            f = incpf[ext].data.field('orig_value')[:]

# initialize FFI image

        ffimage = numpy.zeros((1070,1132),dtype=float32)
        ffimage[:,:] = numpy.nan

# populate FFI image

        for i in range(len(k)):
            ffimage[r[i],c[i]] = f[i]

# create image extension headers
   
        hdu1 = ImageHDU(ffimage)
        hdu1.header.update('INHERIT',True,'inherit primary keywords')
        hdu1.header.update('EXTNAME','MOD.OUT %d.%d' % (modkey,outkey),'extension name')
        hdu1.header.update('EXTVER',1,'extension version number')
        hdu1.header.update('TELESCOP','Kepler','telescope')
        hdu1.header.update('INSTRUME','Kepler photometer','detector type')
        hdu1.header.update('CHANNEL',chankey,'CCD channel')
        hdu1.header.update('SKYGROUP',-1,'roll-independent location of channel')
        hdu1.header.update('MODULE',modkey,'CCD module')
        hdu1.header.update('OUTPUT',outkey,'CCD output')
        hdu1.header.update('TIMEREF','SOLARSYSTEM','barycentric correction applied to times')
        hdu1.header.update('TASSIGN','SPACECRAFT','where time is assigned')
        hdu1.header.update('TIMESYS','TDB','time system is barycentric JD')
        hdu1.header.update('MJDSTART',startkey,'observation start time in MJD')
        hdu1.header.update('MJDEND',endkey,'observation stop time in MJD')
        hdu1.header.update('BJDREFI',2454833,'integer part of BJD reference date')
        hdu1.header.update('BJDREFF',0.0,'fraction of day in BJD reference date')
        hdu1.header.update('TIMEUNIT','d','time unit for TIME, TSTART and TSTOP')
        hdu1.header.update('TSTART',tstart,'observation start time in BJD-BJDREF')
        hdu1.header.update('TSTOP',tstop,'observation stop time in BJD-BJDREF')
        hdu1.header.update('TELAPSE',telapse,'TSTOP - TSTART')
        hdu1.header.update('EXPOSURE',exposure,'time on source')
        hdu1.header.update('LIVETIME',exposure,'TELAPSE multiplied by DEADC')
        hdu1.header.update('TIMEPIXR',0.5,'bin time beginning=0 middle=0.5 end=1')
        hdu1.header.update('TIERRELA',5.78E-07,'relative time error')
        hdu1.header.update('INT_TIME',6.019802903270,'[s] photon accumulation time per frame')
        hdu1.header.update('READTIME',0.518948526144,'readout time per frame')
        hdu1.header.update('FRAMETIM',6.5387514294144005,'[s] frame time (INT_TIME + READTIME)')
        hdu1.header.update('NUM_FRM',270,'number of frames per time stamp')
        hdu1.header.update('FGSFRPER',103.7897052288,'[ms] FGS frame period')
        hdu1.header.update('NUMFGSFP',58,'number of FGS frame periods per exposure')
        hdu1.header.update('TIMEDEL',telapse,'[d] time resolution of data')
        hdu1.header.update('DATE-OBS','','TSTART as UTC calendar date')
        hdu1.header.update('DATE-END','%sT%sZ' % (dateobs,timeobs),'TSTART as UTC calendar date')
        hdu1.header.update('BTC_PIX1',536,'reference col for barycentric time correction')
        hdu1.header.update('BTC_PIX2',567,'reference row for barycentric time correction')
        hdu1.header.update('BUNIT',bunit,'physical units of image data')
        hdu1.header.update('BARYCORR',-1.0,'[d] barycentric time correction')
        hdu1.header.update('BACKAPP',False,'background is subtracted')
        hdu1.header.update('DEADAPP',True,'deadtime applied')
        hdu1.header.update('VIGNAPP',False,'vignetting or collimator correction applied')
        hdu1.header.update('GAIN',gain,'[electrons/ADU] gain')
        hdu1.header.update('READNOIS',readnois,'[electrons] read noise')
        hdu1.header.update('NREADOUT',270,'number of read per cadence')
        hdu1.header.update('TIMSLICE',timslice,'time-slice readout sequence section')
        hdu1.header.update('MEANBLCK',meanblck,'[count] FSW mean black level')
        hdu1.header.update('RADESYS','ICRS','reference frame of celestial coordinates')
        hdu1.header.update('EQUINOX',2000.0,'equinox of celestial coordinate system')
        hdu1.header.update('WCSNAMEP','PHYSICAL','name of world coordinate system alternate P')
        hdu1.header.update('WCSAXESP',2,'number of WCS physical axes')
        hdu1.header.update('CTYPE1P','RAWX','physical WCS axis 1 type CCD col')
        hdu1.header.update('CUNIT1P','PIXEL','physical WCS axis 1 unit')
        hdu1.header.update('CRPIX1P',1,'reference CCD column')
        hdu1.header.update('CRVAL1P',0,'value at reference CCD column')
        hdu1.header.update('CDELT1P',1.0,'physical WCS axis 1 step')
        hdu1.header.update('CTYPE2P','RAWY','physical WCS axis 2 type CCD row')
        hdu1.header.update('CUNIT2P','PIXEL','physical WCS axis 2 units')
        hdu1.header.update('CRPIX2P',1,'reference CCD row')
        hdu1.header.update('CRVAL2P',0,'value at reference CCD row')
        hdu1.header.update('CDELT2P',1.0,'physical WCS axis 2 step')
        hdulist.append(hdu1)

# close the FFI file

    try:
        os.remove(outname)
    except:
        pass
    hdulist.writeto(outname,checksum=True)
        
# close the cadence pixel file

    try:
        incpf.close()
    except:
        message =  'ERROR -- Cannot close the cadence pixel file'
        status = PyMSG.err(None,message)
        
# close the pixel mapping reference file

    try:
        inpmrf.close()
    except:
        message =  'ERROR -- Cannot close the pixel mapping reference file'
        status = PyMSG.err(None,message)
        

#-------------------------------
def usage():

    print ' -----------------------------------------------------------------------'
    print ' Martin Still NASA ARC Oct 25, 2013'
    print ' '
    print ' Extract FFI from Cadence Pixel Files'
    print ' '
    print ' Typical usage:'
    print ' CadencePixels2FFI.py --cpf=kplr2010269175646_lcs-targ.fits --pmrf=/Volumes/data/Kepler/PMRF/q5-lc.fits'
    print ' '
    print '   --cpf        Cadence Pixel File name (FITS)'
    print '   --pmrf       Pixel Mapping Reference File (i.e. which pixels are collected)'
    print ' -----------------------------------------------------------------------'
    PyMSG.exit(None)

#-------------------------------
if __name__ == "__main__":
    main()

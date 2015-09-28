#!/usr/bin/env python

import PyMSG, PyIO
import getopt, sys, os, re
import time as ti
import numpy, pyfits
from pyfits import *
from numpy import *

def main():

# input argument

    cpflist = '/Volumes/data/Kepler/2W/1311/CPF/cpf.lis'
    pmrf = '/Volumes/data/Kepler/PMRF/nov2013_lcm.fits'
    try:
	opts, args = getopt.getopt(sys.argv[1:],"h:cpa",["help","cpflist=","pmrf=","cal"])
    except getopt.GetoptError:
	usage()
    calib = False
    for o, a in opts:
	if o in ("-h", "--help"):
	    usage()
	if o in ("-c", "--cpflist"):
            cpflist = str(a)
	if o in ("-p", "--pmrf"):
            pmrf = str(a)
	if o in ("-a", "--cal"):
            calib = True

# read the list of Cadence Pixel Files

    cpfs = []
    lines, status = PyIO.openascii(cpflist,'r',None)
    if status == 0:
        for line in lines:
            cpfs.append(line.strip())
        ncad = len(cpfs)
        lastcad = re.sub('_lcs-targ.fits','',os.path.basename(cpfs[-1])).strip('kplr')

# open the Pixel Mapping Reference File

    try:
        instr = open(pmrf,mode='readonly',memmap=True)
    except:
        message =  'ERROR -- Cannot read the pixel reference file'
        status = PyMSG.err(None,message)

# iterate through the FITS extensions of the PMRF and read the data 

    for ext in range(1,2):
        start = ti.time()
        card = instr[ext].header.ascardlist()
        chankey = card['CHANNEL'].value
        modkey = card['MODULE'].value
        outkey = card['OUTPUT'].value
        r = instr[ext].data.field('row')[:]
        c = instr[ext].data.field('column')[:]
        k = instr[ext].data.field('target_id')[:]

# populate dictionaries for channel and MORC defintions (module output row column)

        channel = []
        module = []
        output = []
        row = {}
        column = {}
        for line in range(len(r)):
            channel.append((k[line],chankey))
            module.append((k[line],modkey))
            output.append((k[line],outkey))
            try: 
                row[k[line]].append(r[line])
                column[k[line]].append(c[line])
            except:
                row[k[line]] = [r[line]]
                column[k[line]] = [c[line]]
        channel = dict(channel)
        module = dict(module)
        output = dict(output)

# physical wcs information for each target

        dcol = {}
        drow = {}
        eformat = {}
        jformat = {}
        kformat = {}
        coldim = {}
        crval1p = {}
        crval2p = {}
        raw_cnts = {}
        flux = {}
        flux_err = {}
        flux_bkg = {}
        flux_bkg_err = {}
        cosmic_rays = {}
        mask = {}
        for kepid in channel.keys():
            dcol[kepid] = max(column[kepid]) - min(column[kepid]) + 1
            drow[kepid] = max(row[kepid]) - min(row[kepid]) + 1
            coldim[kepid] = '(' + str(dcol[kepid]) + ',' + str(drow[kepid]) + ')'
            eformat[kepid] = str(drow[kepid]*dcol[kepid]) + 'E'
            jformat[kepid] = str(drow[kepid]*dcol[kepid]) + 'J'
            kformat[kepid] = str(drow[kepid]*dcol[kepid]) + 'K'
            crval1p[kepid] = min(column[kepid])
            crval2p[kepid] = min(row[kepid])
            
# declare pixel data array sizes

            raw_cnts[kepid] = zeros((ncad,drow[kepid],dcol[kepid]),dtype='int')
            flux[kepid] = zeros((ncad,drow[kepid],dcol[kepid]),dtype='float32')
            flux_err[kepid] = zeros((ncad,drow[kepid],dcol[kepid]),dtype='float32')
            flux_bkg[kepid] = zeros((ncad,drow[kepid],dcol[kepid]),dtype='float32')
            flux_bkg_err[kepid] = zeros((ncad,drow[kepid],dcol[kepid]),dtype='float32')
            cosmic_rays[kepid] = zeros((ncad,drow[kepid],dcol[kepid]),dtype='float32')
            mask[kepid] = ones((drow[kepid],dcol[kepid]),dtype='int')
            raw_cnts[kepid][:,:,:] = -1
            flux[kepid][:,:,:] = nan
            flux_err[kepid][:,:,:] = nan
            flux_bkg[kepid][:,:,:] = nan
            flux_bkg_err[kepid][:,:,:] = nan
            cosmic_rays[kepid][:,:,:] = nan

# declare  array sizes
            
        mjd = zeros((ncad),dtype='float64')
        time = zeros((ncad),dtype='float64')
        timecorr = zeros((ncad),dtype='float32')
        cadenceno = zeros((ncad),dtype='int')
        quality = zeros((ncad),dtype='int')
        pos_corr1 = zeros((ncad),dtype='float32')
        pos_corr2 = zeros((ncad),dtype='float32')

# open cadence pixel files
                    
        cadence = 0
        for cpf in cpfs:
            try:
                incpf = open(cpf,mode='readonly',memmap=True)
                card0 = incpf[0].header.ascardlist()
                carde = incpf[ext].header.ascardlist()
            except:
                message =  'ERROR -- Cannot read the cadence pixel file ' + cpf
                status = PyMSG.err(None,message)

# read cadence pixel file keywords

            dateobs = card0['DATE-OBS'].value
            timeobs = card0['TIME-OBS'].value
            lcfxdoff = card0['LCFXDOFF'].value
            scfxdoff = card0['SCFXDOFF'].value
            crpix1 = carde['CRPIX1'].value
            crpix2 = carde['CRPIX2'].value
            crval1 = carde['CRVAL1'].value
            crval2 = carde['CRVAL2'].value
            gain = carde['GAIN'].value
            readnois = carde['READNOIS'].value
            timslice = carde['TIMSLICE'].value           
            meanblck = carde['MEANBLCK'].value

# read cadence pixel file tables

            mjd[cadence] = card0['MID_TIME'].value
            time[cadence] = card0['MID_TIME'].value + 2400000.5 - 2454833.0
            timecorr[cadence] = 0.0
            cadenceno[cadence] = card0['LC_COUNT'].value
            raw = incpf[ext].data.field('orig_value')[:]
            if calib:
                cal = incpf[ext].data.field('cal_value')[:]
                unc = incpf[ext].data.field('cal_uncert')[:]
            else:
                cal = incpf[ext].data.field('orig_value')[:]
                unc = numpy.sqrt(cal)

# populate the pixel data arrays
            
            for line in range(len(r)):
                i,j = ccd2mask(1.0,1.0,crval1p[k[line]],crval2p[k[line]],1.0,1.0,c[line],r[line])
                raw_cnts[k[line]][cadence,j,i] = raw[line] 
                flux[k[line]][cadence,j,i] = cal[line] 
                flux_err[k[line]][cadence,j,i] = unc[line] 

# close the cadence pixel file

            try:
                incpf.close()
            except:
                message =  'ERROR -- Cannot close the cadence pixel file'
                status = PyMSG.err(None,message)
            cadence += 1

# iterate through each target in order to create a Target Pixel File

        for kepid in channel.keys():
#            if drow[kepid] * dcol[kepid] < 1700:
            if drow[kepid] * dcol[kepid] < 1e9:

# reshape the data arrays for the target pixel files

                raw_cnts[kepid] = numpy.reshape(raw_cnts[kepid],(ncad,dcol[kepid]*drow[kepid]))
                flux[kepid] = numpy.reshape(flux[kepid],(ncad,dcol[kepid]*drow[kepid]))
                flux_err[kepid] = numpy.reshape(flux_err[kepid],(ncad,dcol[kepid]*drow[kepid]))
                flux_bkg[kepid] = numpy.reshape(flux_bkg[kepid],(ncad,dcol[kepid]*drow[kepid]))
                flux_bkg_err[kepid] = numpy.reshape(flux_bkg_err[kepid],(ncad,dcol[kepid]*drow[kepid]))
                cosmic_rays[kepid] = numpy.reshape(cosmic_rays[kepid],(ncad,dcol[kepid]*drow[kepid]))

# create primary FITS extensions

                hdu0 = PrimaryHDU()
                hdu0.header.update('EXTEND',True,'file may contain standard extensions')
                hdu0.header.update('EXTNAME','PRIMARY','name of extension')
                hdu0.header.update('EXTVER',1,'extension version number')
                hdu0.header.update('ORIGIN','NASA/Ames','organization that generated this file')
                hdu0.header.update('DATE','2010-04-08','file creation date')
                hdu0.header.update('CREATOR','CadencePixels.py','SW version used to create this file')
                hdu0.header.update('PROCVER',1.0,'processing script version')
                hdu0.header.update('FILEVER','2.0','file format version')
                hdu0.header.update('TIMVERSN','OGIP/93-003','OGIP memo number for file format')
                hdu0.header.update('TELESCOP','Kepler','telescope')
                hdu0.header.update('INSTRUME','Kepler photometer','detector type')
                hdu0.header.update('OBJECT','KIC %s' % kepid,'string version of kepID')
                hdu0.header.update('KEPLERID',str(kepid),'unique Kepler target identifier')
                hdu0.header.update('CHANNEL',channel[kepid],'CCD channel')
                hdu0.header.update('SKYGROUP',0,'roll-independent location of channel')
                hdu0.header.update('MODULE',module[kepid],'CCD module')
                hdu0.header.update('OUTPUT',output[kepid],'CCD output')
                hdu0.header.update('QUARTER',0,'mission quarter during which data was collected')
                hdu0.header.update('SEASON',0,'mission season during which data was collected')
                hdu0.header.update('DATA_REL',0,'version of data release notes describing data')
                hdu0.header.update('OBSMODE','long cadence','observing mode')
                hdu0.header.update('RADESYS','ICRS','reference frame of celestial coordinates')
                hdu0.header.update('RA_OBJ',0.0,'[deg] right ascension from KIC')
                hdu0.header.update('DEC_OBJ',0.0,'[deg] declination from KIC')
                hdu0.header.update('EQUINOX',2000.0,'equinox of celestial coordinate system')
                hdu0.header.update('PMRA',0.0,'[arcsec/yr] RA proper motion')
                hdu0.header.update('PMDEC',0.0,'[arcsec/yr] Dec proper motion')
                hdu0.header.update('PMTOTAL',0.0,'[arcsec/yr] total proper motion')
                hdu0.header.update('PARALLAX',0.0,'[arcsec] parallax')
                hdu0.header.update('GLON',0.0,'[deg] galactic longitude')
                hdu0.header.update('GLAT',0.0,'[deg] galactic latitude')
                hdu0.header.update('GMAG',0.0,'[mag] SDSS g band magnitude from KIC')
                hdu0.header.update('RMAG',0.0,'[mag] SDSS r band magnitude from KIC')
                hdu0.header.update('IMAG',0.0,'[mag] SDSS i band magnitude from KIC')
                hdu0.header.update('ZMAG',0.0,'[mag] SDSS z band magnitude from KIC')
                hdu0.header.update('D51MAG',0.0,'[mag] D51 magnitude, from KIC')
                hdu0.header.update('JMAG',0.0,'[mag] J band magnitude from 2MASS')
                hdu0.header.update('HMAG',0.0,'[mag] H band magnitude from 2MASS')
                hdu0.header.update('KMAG',0.0,'[mag] K band magnitude from 2MASS')
                hdu0.header.update('KEPMAG',0.0,'[mag] Kepler magnitude (Kp) from KIC')
                hdu0.header.update('GRCOLOR',0.0,'[mag] (g-r) color, SDSS bands')
                hdu0.header.update('JKCOLOR',0.0,'[mag] (J-K) color, 2MASS bands')
                hdu0.header.update('GKCOLOR',0.0,'[mag] (g-K) color, SDSS g - 2MASS K')
                hdu0.header.update('TEFF',0.0,'[K] effective temperature from KIC')
                hdu0.header.update('LOGG',0.0,'[cm/s2] log10 surface gravity from KIC')
                hdu0.header.update('FEH',0.0,'[log10([Fe/H])] metallicity from KIC')
                hdu0.header.update('EBMINUSV',0.0,'[mag] E(B-V) redenning from KIC')
                hdu0.header.update('AV',0.0,'[mag] A_v extinction from KIC')
                hdu0.header.update('RADIUS',0.0,'[solar radii] stellar radius from KIC')
                hdu0.header.update('TMINDEX',0,'unique 2MASS catalog ID from KIC')
                hdu0.header.update('SCPID',0,'unique SCP processing ID from KIC') 
                hdulist = HDUList(hdu0)
                
# create the outfile table extension

                col1 = Column(name='TIME',format='D',unit='BJD - 2454833',array=time)
                col2 = Column(name='TIMECORR',format='E',unit='D',array=timecorr)
                col3 = Column(name='CADENCENO',format='J',array=cadenceno)
                col4 = Column(name='RAW_CNTS',format=jformat[kepid],unit='counts',dim=coldim[kepid],array=raw_cnts[kepid])
                col5 = Column(name='FLUX',format=eformat[kepid],unit='e-/s',dim=coldim[kepid],array=flux[kepid])
                col6 = Column(name='FLUX_ERR',format=eformat[kepid],unit='e-/s',dim=coldim[kepid],array=flux_err[kepid])
                col7 = Column(name='FLUX_BKG',format=eformat[kepid],unit='e-/s',dim=coldim[kepid],array=flux_bkg[kepid])
                col8 = Column(name='FLUX_BKG_ERR',format=eformat[kepid],unit='e-/s',dim=coldim[kepid],array=flux_bkg_err[kepid])
                col9 = Column(name='COSMIC_RAYS',format=eformat[kepid],unit='e-/s',dim=coldim[kepid],array=cosmic_rays[kepid])
                col10 = Column(name='QUALITY',format='J',array=quality)
                col11 = Column(name='POS_CORR1',format='E',unit='pixels',array=pos_corr1)
                col12 = Column(name='POS_CORR2',format='E',unit='pixels',array=pos_corr2)
                cols = ColDefs([col1,col2,col3,col4,col5,col6,col7,col8,col9,col10,col11,col12])
                hdu1 = new_table(cols)
                
                hdu1.header.update('TTYPE1','TIME','column title: data time stamps')
                hdu1.header.update('TFORM1','D','column format: 64-bit floating point')
                hdu1.header.update('TUNIT1','BJD - 2454833','column units: barycenter corrected JD')
                hdu1.header.update('TDISP1','D14.7','column display format')
                
                hdu1.header.update('TTYPE2','TIMECORR','column title: barycenter - timeslice correction')
                hdu1.header.update('TFORM2','E','column format: 32-bit floating point')
                hdu1.header.update('TUNIT2','d','column units: day')
                hdu1.header.update('TDISP2','E14.7','column display format')
                
                hdu1.header.update('TTYPE3','CADENCENO','column title: unique cadence number')
                hdu1.header.update('TFORM3','J','column format: signed 32-bit integer')
                
                hdu1.header.update('TTYPE4','RAW_CNTS','column title: raw pixel counts')
                hdu1.header.update('TFORM4',jformat[kepid],'column format: signed integer32')
                hdu1.header.update('TUNIT4','count','column units: count')
                hdu1.header.update('TDISP4','I8','column display format')
                hdu1.header.update('TDIM4',coldim[kepid],'column dimensions: pixel apeture array')
                hdu1.header.update('TNULL4',-1,'column null value indicator')
                hdu1.header.update('WCSN4P','PHYSICAL','table column WCS name')
                hdu1.header.update('WCAX4P',2,'table column physical WCS dimensions')
                hdu1.header.update('1CTY4P','RAWX','table column physical WCS axis 1 type, CCD col')
                hdu1.header.update('2CTY4P','RAWY','table column physical WCS axis 2 type, CCD row')
                hdu1.header.update('1CUN4P','PIXEL','table column physical WCS axis 1 unit')
                hdu1.header.update('2CUN4P','PIXEL','table column physical WCS axis 2 unit')
                hdu1.header.update('1CRV4P',crval1p[kepid],'table column physical WCS axis 1 ref value')
                hdu1.header.update('2CRV4P',crval2p[kepid],'table column physical WCS axis 2 ref value')
                hdu1.header.update('1CDL4P',1.0,'table column physical WCS axis 1 step')
                hdu1.header.update('2CDL4P',1.0,'table column physical WCS axis 2 step')
                hdu1.header.update('1CRP4P',1,'table column physical WCS axis 1 reference')
                hdu1.header.update('2CRP4P',1,'table column physical WCS axis 2 reference')
                hdu1.header.update('WCAX4',2,'number of WCS axes')
                hdu1.header.update('1CTYP4','RA---TAN','right ascension coordinate type')
                hdu1.header.update('2CTYP4','DEC--TAN','declination coordinate type')
                hdu1.header.update('1CRPX4',crpix1,'[pixel] reference pixel along image axis 1')
                hdu1.header.update('2CRPX4',crpix2,'[pixel] reference pixel along image axis 2')
                hdu1.header.update('1CRVL4',crval1,'[deg] right ascension at reference pixel')
                hdu1.header.update('2CRVL4',crval2,'[deg] delination at reference pixel')
                hdu1.header.update('1CUNI4','deg','physical unit in column dimension')
                hdu1.header.update('2CUNI4','deg','physical unit in row dimension')
                hdu1.header.update('1CDLT4',-0.001107921039197,'[deg] pixel scale in RA dimension')
                hdu1.header.update('2CDLT4',0.001107921039197,'[deg] pixel scale in Dec dimension')
                hdu1.header.update('11PC4',-0.48060068512850346,'linear transformation matrix element cos(th)')
                hdu1.header.update('12PC4',0.8815429040054924,'linear transformation matrix element -sin(th)')
                hdu1.header.update('21PC4',-0.8748562715817224,'linear transformation matrix element sin(th)')
                hdu1.header.update('22PC4',-0.4760223379649577,'linear transformation matrix element cos(th)')          
                hdu1.header.update('TTYPE5','FLUX','column title: calibrated pixel flux')
                hdu1.header.update('TFORM5',eformat[kepid],'column format: image of 32-bit floating point')
                hdu1.header.update('TUNIT5','e-/s','column units: electrons per second')
                hdu1.header.update('TDISP5','E14.7','column display format')
                hdu1.header.update('TDIM5',coldim[kepid],'column dimensions: pixel apeture array')
                hdu1.header.update('WCSN5P','PHYSICAL','table column WCS name')
                hdu1.header.update('WCAX5P',2,'table column physical WCS dimensions')
                hdu1.header.update('1CTY5P','RAWX','table column physical WCS axis 1 type, CCD col')
                hdu1.header.update('2CTY5P','RAWY','table column physical WCS axis 2 type, CCD row')
                hdu1.header.update('1CUN5P','PIXEL','table column physical WCS axis 1 unit')
                hdu1.header.update('2CUN5P','PIXEL','table column physical WCS axis 2 unit')
                hdu1.header.update('1CRV5P',crval1p[kepid],'table column physical WCS axis 1 ref value')
                hdu1.header.update('2CRV5P',crval2p[kepid],'table column physical WCS axis 2 ref value')
                hdu1.header.update('1CDL5P',1.0,'table column physical WCS axis 1 step')
                hdu1.header.update('2CDL5P',1.0,'table column physical WCS axis 2 step')
                hdu1.header.update('1CRP5P',1,'table column physical WCS axis 1 reference')
                hdu1.header.update('2CRP5P',1,'table column physical WCS axis 2 reference')
                hdu1.header.update('WCAX5',2,'number of WCS axes')
                hdu1.header.update('1CTYP5','RA---TAN','right ascension coordinate type')
                hdu1.header.update('2CTYP5','DEC--TAN','declination coordinate type')
                hdu1.header.update('1CRPX5',crpix1,'[pixel] reference pixel along image axis 1')
                hdu1.header.update('2CRPX5',crpix2,'[pixel] reference pixel along image axis 2')
                hdu1.header.update('1CRVL5',crval1,'[deg] right ascension at reference pixel')
                hdu1.header.update('2CRVL5',crval2,'[deg] delination at reference pixel')
                hdu1.header.update('1CUNI5','deg','physical unit in column dimension')
                hdu1.header.update('2CUNI5','deg','physical unit in row dimension')
                hdu1.header.update('1CDLT5',-0.001107921039197,'[deg] pixel scale in RA dimension')
                hdu1.header.update('2CDLT5',0.001107921039197,'[deg] pixel scale in Dec dimension')
                hdu1.header.update('11PC5',-0.48060068512850346,'linear transformation matrix element cos(th)')
                hdu1.header.update('12PC5',0.8815429040054924,'linear transformation matrix element -sin(th)')
                hdu1.header.update('21PC5',-0.8748562715817224,'linear transformation matrix element sin(th)')
                hdu1.header.update('22PC5',-0.4760223379649577,'linear transformation matrix element cos(th)')
                
                hdu1.header.update('TTYPE6','FLUX_ERR','column title: 1-sigma calibrated uncertainty')
                hdu1.header.update('TFORM6',eformat[kepid],'column format: image of 32-bit floating point')
                hdu1.header.update('TUNIT6','e-/s','column units: electrons per second')
                hdu1.header.update('TDISP6','E14.7','column display format')
                hdu1.header.update('TDIM6',coldim[kepid],'column dimensions: pixel apeture array')
                hdu1.header.update('WCSN6P','PHYSICAL','table column WCS name')
                hdu1.header.update('WCAX6P',2,'table column physical WCS dimensions')
                hdu1.header.update('1CTY6P','RAWX','table column physical WCS axis 1 type, CCD col')
                hdu1.header.update('2CTY6P','RAWY','table column physical WCS axis 2 type, CCD row')
                hdu1.header.update('1CUN6P','PIXEL','table column physical WCS axis 1 unit')
                hdu1.header.update('2CUN6P','PIXEL','table column physical WCS axis 2 unit')
                hdu1.header.update('1CRV6P',crval1p[kepid],'table column physical WCS axis 1 ref value')
                hdu1.header.update('2CRV6P',crval2p[kepid],'table column physical WCS axis 2 ref value')
                hdu1.header.update('1CDL6P',1.0,'table column physical WCS axis 1 step')
                hdu1.header.update('2CDL6P',1.0,'table column physical WCS axis 2 step')
                hdu1.header.update('1CRP6P',1,'table column physical WCS axis 1 reference')
                hdu1.header.update('2CRP6P',1,'table column physical WCS axis 2 reference')
                hdu1.header.update('WCAX6',2,'number of WCS axes')
                hdu1.header.update('1CTYP6','RA---TAN','right ascension coordinate type')
                hdu1.header.update('2CTYP6','DEC--TAN','declination coordinate type')
                hdu1.header.update('1CRPX6',crpix1,'[pixel] reference pixel along image axis 1')
                hdu1.header.update('2CRPX6',crpix2,'[pixel] reference pixel along image axis 2')
                hdu1.header.update('1CRVL6',crval1,'[deg] right ascension at reference pixel')
                hdu1.header.update('2CRVL6',crval2,'[deg] delination at reference pixel')
                hdu1.header.update('1CUNI6','deg','physical unit in column dimension')
                hdu1.header.update('2CUNI6','deg','physical unit in row dimension')
                hdu1.header.update('1CDLT6',-0.001107921039197,'[deg] pixel scale in RA dimension')
                hdu1.header.update('2CDLT6',0.001107921039197,'[deg] pixel scale in Dec dimension')
                hdu1.header.update('11PC6',-0.48060068512850346,'linear transformation matrix element cos(th)')
                hdu1.header.update('12PC6',0.8815429040054924,'linear transformation matrix element -sin(th)')
                hdu1.header.update('21PC6',-0.8748562715817224,'linear transformation matrix element sin(th)')
                hdu1.header.update('22PC6',-0.4760223379649577,'linear transformation matrix element cos(th)')
                
                hdu1.header.update('TTYPE7','FLUX_BKG','column title: calibrated background flux')
                hdu1.header.update('TFORM7',eformat[kepid],'column format: image of 32-bit floating point')
                hdu1.header.update('TUNIT7','e-/s','column units: electrons per second')
                hdu1.header.update('TDISP7','E14.7','column display format')
                hdu1.header.update('TDIM7',coldim[kepid],'column dimensions: pixel apeture array')
                hdu1.header.update('WCSN7P','PHYSICAL','table column WCS name')
                hdu1.header.update('WCAX7P',2,'table column physical WCS dimensions')
                hdu1.header.update('1CTY7P','RAWX','table column physical WCS axis 1 type, CCD col')
                hdu1.header.update('2CTY7P','RAWY','table column physical WCS axis 2 type, CCD row')
                hdu1.header.update('1CUN7P','PIXEL','table column physical WCS axis 1 unit')
                hdu1.header.update('2CUN7P','PIXEL','table column physical WCS axis 2 unit')
                hdu1.header.update('1CRV7P',crval1p[kepid],'table column physical WCS axis 1 ref value')
                hdu1.header.update('2CRV7P',crval2p[kepid],'table column physical WCS axis 2 ref value')
                hdu1.header.update('1CDL7P',1.0,'table column physical WCS axis 1 step')
                hdu1.header.update('2CDL7P',1.0,'table column physical WCS axis 2 step')
                hdu1.header.update('1CRP7P',1,'table column physical WCS axis 1 reference')
                hdu1.header.update('2CRP7P',1,'table column physical WCS axis 2 reference')
                hdu1.header.update('WCAX7',2,'number of WCS axes')
                hdu1.header.update('1CTYP7','RA---TAN','right ascension coordinate type')
                hdu1.header.update('2CTYP7','DEC--TAN','declination coordinate type')
                hdu1.header.update('1CRPX7',crpix1,'[pixel] reference pixel along image axis 1')
                hdu1.header.update('2CRPX7',crpix2,'[pixel] reference pixel along image axis 2')
                hdu1.header.update('1CRVL7',crval1,'[deg] right ascension at reference pixel')
                hdu1.header.update('2CRVL7',crval2,'[deg] delination at reference pixel')
                hdu1.header.update('1CUNI7','deg','physical unit in column dimension')
                hdu1.header.update('2CUNI7','deg','physical unit in row dimension')
                hdu1.header.update('1CDLT7',-0.001107921039197,'[deg] pixel scale in RA dimension')
                hdu1.header.update('2CDLT7',0.001107921039197,'[deg] pixel scale in Dec dimension')
                hdu1.header.update('11PC7',-0.48060068512850346,'linear transformation matrix element cos(th)')
                hdu1.header.update('12PC7',0.8815429040054924,'linear transformation matrix element -sin(th)')
                hdu1.header.update('21PC7',-0.8748562715817224,'linear transformation matrix element sin(th)')
                hdu1.header.update('22PC7',-0.4760223379649577,'linear transformation matrix element cos(th)')
                
                hdu1.header.update('TTYPE8','FLUX_BKG_ERR','column title: 1-sigma cal. background uncertain')
                hdu1.header.update('TFORM8',eformat[kepid],'column format: image of 32-bit floating point')
                hdu1.header.update('TUNIT8','e-/s','column units: electrons per second')
                hdu1.header.update('TDISP8','E14.7','column display format')
                hdu1.header.update('TDIM8',coldim[kepid],'column dimensions: pixel apeture array')
                hdu1.header.update('WCSN8P','PHYSICAL','table column WCS name')
                hdu1.header.update('WCAX8P',2,'table column physical WCS dimensions')
                hdu1.header.update('1CTY8P','RAWX','table column physical WCS axis 1 type, CCD col')
                hdu1.header.update('2CTY8P','RAWY','table column physical WCS axis 2 type, CCD row')
                hdu1.header.update('1CUN8P','PIXEL','table column physical WCS axis 1 unit')
                hdu1.header.update('2CUN8P','PIXEL','table column physical WCS axis 2 unit')
                hdu1.header.update('1CRV8P',crval1p[kepid],'table column physical WCS axis 1 ref value')
                hdu1.header.update('2CRV8P',crval2p[kepid],'table column physical WCS axis 2 ref value')
                hdu1.header.update('1CDL8P',1.0,'table column physical WCS axis 1 step')
                hdu1.header.update('2CDL8P',1.0,'table column physical WCS axis 2 step')
                hdu1.header.update('1CRP8P',1,'table column physical WCS axis 1 reference')
                hdu1.header.update('2CRP8P',1,'table column physical WCS axis 2 reference')
                hdu1.header.update('WCAX8',2,'number of WCS axes')
                hdu1.header.update('1CTYP8','RA---TAN','right ascension coordinate type')
                hdu1.header.update('2CTYP8','DEC--TAN','declination coordinate type')
                hdu1.header.update('1CRPX8',crpix1,'[pixel] reference pixel along image axis 1')
                hdu1.header.update('2CRPX8',crpix2,'[pixel] reference pixel along image axis 2')
                hdu1.header.update('1CRVL8',crval1,'[deg] right ascension at reference pixel')
                hdu1.header.update('2CRVL8',crval2,'[deg] delination at reference pixel')
                hdu1.header.update('1CUNI8','deg','physical unit in column dimension')
                hdu1.header.update('2CUNI8','deg','physical unit in row dimension')
                hdu1.header.update('1CDLT8',-0.001107921039197,'[deg] pixel scale in RA dimension')
                hdu1.header.update('2CDLT8',0.001107921039197,'[deg] pixel scale in Dec dimension')
                hdu1.header.update('11PC8',-0.48060068512850346,'linear transformation matrix element cos(th)')
                hdu1.header.update('12PC8',0.8815429040054924,'linear transformation matrix element -sin(th)')
                hdu1.header.update('21PC8',-0.8748562715817224,'linear transformation matrix element sin(th)')
                hdu1.header.update('22PC8',-0.4760223379649577,'linear transformation matrix element cos(th)')
                
                hdu1.header.update('TTYPE9','COSMIC_RAYS','column title: cosmic ray detections')
                hdu1.header.update('TFORM9',eformat[kepid],'column format: image of 32-bit floating point')
                hdu1.header.update('TUNIT9','e-/s','column units: electrons per second')
                hdu1.header.update('TDISP9','E14.7','column display format')
                hdu1.header.update('TDIM9',coldim[kepid],'column dimensions: pixel apeture array')
                hdu1.header.update('WCSN9P','PHYSICAL','table column WCS name')
                hdu1.header.update('WCAX9P',2,'table column physical WCS dimensions')
                hdu1.header.update('1CTY9P','RAWX','table column physical WCS axis 1 type, CCD col')
                hdu1.header.update('2CTY9P','RAWY','table column physical WCS axis 2 type, CCD row')
                hdu1.header.update('1CUN9P','PIXEL','table column physical WCS axis 1 unit')
                hdu1.header.update('2CUN9P','PIXEL','table column physical WCS axis 2 unit')
                hdu1.header.update('1CRV9P',crval1p[kepid],'table column physical WCS axis 1 ref value')
                hdu1.header.update('2CRV9P',crval2p[kepid],'table column physical WCS axis 2 ref value')
                hdu1.header.update('1CDL9P',1.0,'table column physical WCS axis 1 step')
                hdu1.header.update('2CDL9P',1.0,'table column physical WCS axis 2 step')
                hdu1.header.update('1CRP9P',1,'table column physical WCS axis 1 reference')
                hdu1.header.update('2CRP9P',1,'table column physical WCS axis 2 reference')
                hdu1.header.update('WCAX9',2,'number of WCS axes')
                hdu1.header.update('1CTYP9','RA---TAN','right ascension coordinate type')
                hdu1.header.update('2CTYP9','DEC--TAN','declination coordinate type')
                hdu1.header.update('1CRPX9',crpix1,'[pixel] reference pixel along image axis 1')
                hdu1.header.update('2CRPX9',crpix2,'[pixel] reference pixel along image axis 2')
                hdu1.header.update('1CRVL9',crval1,'[deg] right ascension at reference pixel')
                hdu1.header.update('2CRVL9',crval2,'[deg] delination at reference pixel')
                hdu1.header.update('1CUNI9','deg','physical unit in column dimension')
                hdu1.header.update('2CUNI9','deg','physical unit in row dimension')
                hdu1.header.update('1CDLT9',-0.001107921039197,'[deg] pixel scale in RA dimension')
                hdu1.header.update('2CDLT9',0.001107921039197,'[deg] pixel scale in Dec dimension')
                hdu1.header.update('11PC9',-0.48060068512850346,'linear transformation matrix element cos(th)')
                hdu1.header.update('12PC9',0.8815429040054924,'linear transformation matrix element -sin(th)')
                hdu1.header.update('21PC9',-0.8748562715817224,'linear transformation matrix element sin(th)')
                hdu1.header.update('22PC9',-0.4760223379649577,'linear transformation matrix element cos(th)')
                
                hdu1.header.update('TTYPE10','QUALITY','column title: pixel quality flags')
                hdu1.header.update('TFORM10','J','column format: signed 32-bit integer')
                hdu1.header.update('TDISP10','B16.16','column display format')
                
                hdu1.header.update('TTYPE11','POS_CORR1','column title: row position correction')
                hdu1.header.update('TFORM11','E','column format: 32-bit floating point')
                hdu1.header.update('TUNIT11','pixels','column units: pixels')
                hdu1.header.update('TDISP11','D14.7','column display format')
                
                hdu1.header.update('TTYPE12','POS_CORR2','column title: column position correction')
                hdu1.header.update('TFORM12','E','column format: 32-bit floating point')
                hdu1.header.update('TUNIT12','pixels','column units: pixels')
                hdu1.header.update('TDISP12','D14.7','column display format')
                
                hdu1.header.update('INHERIT',True,'inherit the primary header')
                hdu1.header.update('EXTNAME','TARGETTABLES','name of extension')
                hdu1.header.update('EXTVER',1,'extension version number')
                hdu1.header.update('TELESCOP','Kepler','telescope')
                hdu1.header.update('INSTRUME','Kepler photometer','detector type')
                hdu0.header.update('OBJECT','KIC %s' % kepid,'string version of kepID')
                hdu0.header.update('KEPLERID',str(kepid),'unique Kepler target identifier')
                hdu1.header.update('RADESYS','ICRS','reference frame of celestial coordinates')
                hdu1.header.update('RA_OBJ',0.0,'[deg] right ascension from KIC')
                hdu1.header.update('DEC_OBJ',0.0,'[deg] declination from KIC')
                hdu1.header.update('EQUINOX',2000.0,'equinox of celestial coordinate system')
                hdu1.header.update('EXPOSURE',(time[-1]-time[0])*0.92063492,'[d] time on source')
                hdu1.header.update('TIMEREF','SOLARSYSTEM','barycentric correction applied to times')
                hdu1.header.update('TASSIGN','SPACECRAFT','where time is assigned')
                hdu1.header.update('TIMESYS','TDB','time system is barycentric JD')
                hdu1.header.update('BJDREFI',2454833,'integer part of BJD reference date')
                hdu1.header.update('BJDREFF',0.0,'fraction of day in BJD reference date')
                hdu1.header.update('TIMEUNIT','d','time unit for TIME, TSTART and TSTOP')
                hdu1.header.update('TELAPSE',time[-1]-time[0],'[d] TSTOP - TSTART')
                hdu1.header.update('LIVETIME',(time[-1]-time[0])*0.92063492,'[d] TELAPSE multiplied by DEADC')
                hdu1.header.update('TSTART',time[0],'observation start time in BJD - BJDREF')
                hdu1.header.update('TSTOP',time[-1],'observation stop time in BJD - BJDREF')
                hdu1.header.update('LC_START',mjd[0],'observation start time in MJD')
                hdu1.header.update('LC_END',mjd[-1],'observation stop time in MJD')
                hdu1.header.update('DEADC',0.92063492,'deadtime correction')
                hdu1.header.update('TIMEPIXR',0.5,'bin time beginning=0 middle=0.5 end=1')
                hdu1.header.update('IERRELA',5.78E-07,'[d] relative time error')
                hdu1.header.update('TIERABSO',1.0E-07,'[d] absolute time error')
                hdu1.header.update('INT_TIME',6.019802903270,'[s] photon accumulation time per frame')
                hdu1.header.update('READTIME',0.518948526144,'[s] readout time per frame')
                hdu1.header.update('FRAMETIM',6.538751429414,'[s] frame time (INT_TIME + READTIME)')
                hdu1.header.update('NUM_FRM',270,'number of frames per time stamp')
                hdu1.header.update('TIMEDEL',0.02043359821692,'[d] time resolution of data')
                hdu1.header.update('DATE-OBS',dateobs+':'+timeobs+'Z','TSTART as UTC calendar date')
                hdu1.header.update('DATE-END','','TSTOP as UTC calendar date')
                hdu1.header.update('BACKAPP',True,'background is subtracted')
                hdu1.header.update('DEADAPP',True,'deadtime applied')
                hdu1.header.update('VIGNAPP',True,'vignetting or collimator correction applied')
                hdu1.header.update('GAIN',gain,'[electrons/count] channel gain')
                hdu1.header.update('READNOIS',readnois,'[electrons] read noise')
                hdu1.header.update('NREADOUT',270,'number of read per cadence')
                hdu1.header.update('TIMSLICE',timslice,'time-slice readout sequence section')
                hdu1.header.update('MEANBLCK',meanblck,'[count] FSW mean black level')
                hdu1.header.update('LCFXDOFF',lcfxdoff,'long cadence fixed offset')
                hdu1.header.update('SCFXDOFF',scfxdoff,'short cadence fixed offset')
                hdulist.append(hdu1)

# create the outfile image extension

                hdu2 = ImageHDU(mask[kepid])
                hdu2.header.update('INHERIT',True,'inherit primary keywords')
                hdu2.header.update('EXTNAME','APERTURE','extension name')
                hdu2.header.update('EXTVER',1,'extension version number')
                hdu2.header.update('TELESCOP','Kepler','telescope')
                hdu2.header.update('INSTRUME','Kepler photometer','detector type')
                hdu2.header.update('OBJECT','KIC %s' % kepid,'string version of kepID')
                hdu2.header.update('KEPLERID',str(kepid),'unique Kepler target identifier')
                hdu2.header.update('RADESYS','ICRS','reference frame of celestial coordinates')
                hdu2.header.update('RA_OBJ',0.0,'[deg] right ascension from KIC')
                hdu2.header.update('DEC_OBJ',0.0,'[deg] declination from KIC')
                hdu2.header.update('EQUINOX',2000.0,'equinox of celestial coordinate system')
                hdu2.header.update('WCSAXES',2,'number of WCS axes')
                hdu2.header.update('CTYPE1','RA---TAN','right ascension coordinate type')
                hdu2.header.update('CTYPE2','DEC--TAN','declination coordinate type')
                hdu2.header.update('CRPIX1',crpix1,'[pixel] reference pixel along image axis 1')
                hdu2.header.update('CRPIX2',crpix2,'[pixel] reference pixel along image axis 2')
                hdu2.header.update('CRVAL1',crval1,'[deg] right ascension at reference pixel')
                hdu2.header.update('CRVAL2',crval2,'[deg] delination at reference pixel')
                hdu2.header.update('CUNIT1','deg','physical unit in column dimension')
                hdu2.header.update('CUNIT2','deg','physical unit in row dimension')
                hdu2.header.update('CDELT1',-0.001107921039197,'[deg] pixel scale in RA dimension')
                hdu2.header.update('CDELT2',0.001107921039197,'[deg] pixel scale in Dec dimension')
                hdu2.header.update('PC1_1',-0.48060068512850346,'linear transformation matrix element cos(th)')
                hdu2.header.update('PC1_2',0.8815429040054924,'linear transformation matrix element -sin(th)')
                hdu2.header.update('PC2_1',-0.8748562715817224,'linear transformation matrix element sin(th)')
                hdu2.header.update('PC2_2',-0.4760223379649577,'linear transformation matrix element cos(th)')
                hdu2.header.update('WCSNAMEP','PHYSICAL','name of world coordinate system alternate P')
                hdu2.header.update('WCSAXESP',2,'table column physical WCS dimensions')
                hdu2.header.update('CTYPE1P','RAWX','table column physical WCS axis 1 type, CCD col')
                hdu2.header.update('CTYPE2P','RAWY','table column physical WCS axis 2 type, CCD row')
                hdu2.header.update('CUNIT1P','PIXEL','table column physical WCS axis 1 unit')
                hdu2.header.update('CUNIT2P','PIXEL','table column physical WCS axis 2 unit')
                hdu2.header.update('CRVAL1P',crval1p[kepid],'table column physical WCS axis 1 ref value')
                hdu2.header.update('CRVAL2P',crval2p[kepid],'table column physical WCS axis 2 ref value')
                hdu2.header.update('CDELT1P',1.0,'table column physical WCS axis 1 step')
                hdu2.header.update('CDELT2P',1.0,'table column physical WCS axis 2 step')
                hdu2.header.update('CRPIX1P',1,'table column physical WCS axis 1 reference')
                hdu2.header.update('CRPIX2P',1,'table column physical WCS axis 2 reference')
                hdulist.append(hdu2)            

# write output files
                        
                outname = str(kepid)
                if kepid < 10: outname = '0' + outname
                if kepid < 100: outname = '0' + outname
                if kepid < 1000: outname = '0' + outname
                if kepid < 10000: outname = '0' + outname
                if kepid < 100000: outname = '0' + outname
                if kepid < 1000000: outname = '0' + outname
                if kepid < 10000000: outname = '0' + outname
                if kepid < 100000000: outname = '0' + outname
                outname = 'mod%s.%s_kplr%9s-%13s_lpd-targ.fits' % (module[kepid],output[kepid],outname,lastcad)
                try:
                    os.remove(outname)
                except:
                    pass
                hdulist.writeto(outname,checksum=True)

# runtime

        print 'Channel #%2d - no of targets: %5d - time taken: %4dm' % (ext,len(channel.keys()),(ti.time() - start) / 60)

# close the Pixel Mapping Reference File

    try:
        instr.close()
    except:
        message =  'ERROR -- Cannot close the pixel reference file'
        status = PyMSG.err(None,message)


#-------------------------------
def ccd2mask(crpix1,crpix2,crval1,crval2,cdelt1,cdelt2,ccd_column,ccd_row):
    
    mask_column = (ccd_column - crval1) * cdelt1 + crpix1 - 1
    mask_row = (ccd_row - crval2) * cdelt2 + crpix2 - 1

    return mask_column, mask_row



#-------------------------------
def usage():

    print ' -----------------------------------------------------------------------'
    print ' Martin Still NASA ARC Oct 01, 2013'
    print ' '
    print ' Extract Target Pixel Files from Cadence Pixel Files'
    print ' '
    print ' Typical usage:'
    print ' CadencePixels.py --cpflist=infiles.lis --pmrf=/Volumes/data/Kepler/PMRF/q5-lc.fits'
    print ' '
    print '   --cpflist    ASCII list of input cadence pixel filenames'
    print '   --pmrf       Pixel *** Reference File (i.e. which pixels are collected)'
    print ' -----------------------------------------------------------------------'
    PyMSG.exit(None)

#-------------------------------
if __name__ == "__main__":
    main()

# Copyright (C) 2015 Associated Universities, Inc. Washington DC, USA.
# 
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2 of the License, or
# (at your option) any later version.

# This program is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
# General Public License for more details.
# 
# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
# 
# Correspondence concerning GBT software should be addressed as follows:
#       GBT Operations
#       National Radio Astronomy Observatory
#       P. O. Box 2
#       Green Bank, WV 24944-0002 USA

from get_data import get_data

import pyfits
import numpy
from scipy import constants
import astropy.time as apTime
import argparse
import sys
import os

sdfToAipsVersion = "0.0"

def read_command_line(argv):
    """Read options from the command line."""
    # if no options are set, print help
    if len(argv) == 1:
        argv.append("-h")

    parser = argparse.ArgumentParser(epilog="sdfToAips version: %s" % sdfToAipsVersion)
    parser.add_argument("-c","--channels", type=str,
                        help="Optional channel range to use.  "
                             "'<start>:<end>' counting from 0.")
    parser.add_argument("-a","--average", type=int,
                        help="Optionally average channels, keeping only nchan/naverage channels")
    parser.add_argument("-s","--scans", type=str, 
                        help="Only use data from these scans.  comma separated list or <start>:<end> range syntax or combination of both")
    parser.add_argument("-m","--maxtsys", type=float,
                        help="max Tsys value to use")
    parser.add_argument("-z","--mintsys", type=float,
                        help="min Tsys value to use")
    parser.add_argument("SDFITSfiles", type=str, nargs="+",
                        help="The calibated SDFITS files to use.")
    parser.add_argument("--clobber", default=False, action="store_true",
                        help="Overwrites existing output file if set.")

    # what does idlToSdfits use for output when not specified?  Eventually make this optional.
    parser.add_argument("-o","--output",type=str,required=True,
                        help="output name")
    parser.add_argument("-v","--verbose", type=int, default=4,
                        help="set the verbosity level-- 0-1: none, "
                        "2:errors only, 3:+warnings, "
                        "4(default):+user info, 5:+debug")
    parser.add_argument("-V","--version",action="version",version="sdfToAips version: %s" % sdfToAipsVersion)
    
    args = parser.parse_args()

    return args


# parse_channels, parse_scans, and format_scans are all as they are in gbtgridder
# find a way to share this code

def parse_channels(channelString,verbose=4):
    "Turn a valid channel range into start and end channels"
    start = None
    end = None
    if channelString is not None:
        print "channelString (%s)" % channelString
        # there must be a ":"
        items = channelString.split(":")
        if len(items) != 2:
            if verbose > 1:
                print "Unexpected channels argument, must contain exactly one ':' - %s" % channelString
            return (-1,1)
        if len(items[0]) > 0:
            try:
                # subtract by 1 to go from FITS convention to python
                start = int(items[0])-1
                # fixes for start < 0 happen when used
            except(ValueError):
                if verbose > 1:
                    print repr(':'.join(items[0])), 'not convertable to integer'
                raise
        if len(items[1]) > 0:
            try:
                # subtract by 1 to go from FITS to python convention
                end = int(items[1])-1
            except(ValueError):
                if verbose > 1:
                    print repr(':'.join(items[1])), 'not convertable to integer'
                raise
    return (start,end)

def parse_scans(scanlist):
    """Given a range string, produce a list of integers

    Inclusive and exclusive integers are both possible.

    The range string 1:4,6:8,10 becomes 1,2,3,4,6,7,8,10
    The range string 1:4,-2 becomes 1,3,4
        
    Keywords:
       rangelist -- a range string with inclusive ranges and
                    exclusive integers

    Returns:
       a (list) of integers

    >>> cl = CommandLine()
    >>> cl._parse_range('1:4,6:8,10')
    [1, 2, 3, 4, 6, 7, 8, 10]
    >>> cl._parse_range('1:4,-2')
    [1, 3, 4]
    """
    # copied from _parse_range in gbtpipeline

    oklist = set([])
    excludelist = set([])

    scanlist = scanlist.replace(' ', '')
    scanlist = scanlist.split(',')

    # item is single value or range
    for item in scanlist:
        item = item.split(':')

        # change to ints
        try:
            int_item = [int(ii) for ii in item]
        except(ValueError):
            print repr(':'.join(item)), 'not convertable to integer'
            raise

        if 1 == len(int_item):
            # single inclusive or exclusive item
            if int_item[0] < 0:
                excludelist.add(abs(int_item[0]))
            else:
                oklist.add(int_item[0])

        elif 2 == len(int_item):
            # range
            if int_item[0] <= int_item[1]:
                if int_item[0] < 0:
                    print item[0], ',', item[1], 'must start with a '
                    'non-negative number'
                    return []

                if int_item[0] == int_item[1]:
                    thisrange = [int_item[0]]
                else:
                    thisrange = range(int_item[0], int_item[1]+1)

                for ii in thisrange:
                    oklist.add(ii)
            else:
                print item[0], ',', item[1], 'needs to be in increasing '
                'order'
                raise
        else:
            print item, 'has more than 2 values'

    for exitem in excludelist:
        try:
            oklist.remove(exitem)
        except(KeyError):
            oklist = [str(item) for item in oklist]
            print 'ERROR: excluded item', exitem, 'does not exist in '
            'inclusive range'
            raise

    return sorted(list(oklist))

def format_scans(scanlist):
    "Turn a list of scans into a string using range syntax where appropriate"
    result = None
    rangeCount = 0
    lastScan = -1
    for scan in sorted(scanlist):
        if result is None:
            result = "%s" % scan
        else:
            if (scan-lastScan) == 1:
                # it's a range
                rangeCount += 1
            else:
                if rangeCount > 0:
                    # a previous range has ended
                    result += ":%s" % lastScan
                    rangeCount = 0
                # either way, this one is printed - new item
                result += ",%s" % scan
        lastScan = scan
    if rangeCount > 0:
        # a final range needs to be terminated
        result += ":%s" % lastScan
    return result

def sdfToAips(args):

    chanStart, chanStop = parse_channels(args.channels,verbose=args.verbose)
    if (chanStart is not None and chanStart < 0) or (chanStop is not None and chanStop < 0):
        return

    if chanStart is None:
        chanStart = 0

    average = args.average

    scanlist = args.scans
    if args.scans is not None:
        scanlist = parse_scans(scanlist)    

    # this is really similar to what gbtgridder does.
    # try and use gbtgridder get_data ...
    # more code could probably be shared

    # extract everything from the SDFITS files
    # eventually only a single SDFITS file should be opened at a time
    # probably also just providing smaller chunks within that file

    xsky = None
    ysky = None
    wt = None
    data = None
    nchan = None
    frest = None
    faxis = None
    stokes = None
    source = None
    dataUnits = None
    calibType = None
    veldef = None
    specsys = None
    coordType = (None,None)
    radesys = None
    equinox = None
    observer = None
    telescop = None
    fontend = None
    dateobs = None
    jdobs = None
    projid = None
    sitelong = None
    sitelat = None
    siteelev = None
    bw = None

    # need to get the entire dateobs column back
    # converted to MJD
    
    outputFiles = {}
    if args.verbose > 3:
        print "Loading data ... "

    for thisFile in args.SDFITSfiles:
        try:
            if args.verbose > 3:
                print "   ",thisFile
            dataRecord = get_data(thisFile,nchan,chanStart,chanStop,average,scanlist,
                                  args.mintsys,args.maxtsys,
                                  verbose=args.verbose)
            if dataRecord is None:
                # there was a problem that should not be recovered from
                # reported by get_data, no additional reporting necessary here
                return

            if len(dataRecord) == 0:
                # empty file, skipping
                continue

            if xsky is None:
                xsky = dataRecord["xsky"]
                ysky = dataRecord["ysky"]
                wt = dataRecord["wt"]
                data = dataRecord["data"]
                nchan = dataRecord["nchan"]
                chanStart = dataRecord["chanStart"]
                chanStop = dataRecord["chanStop"]
                frest = dataRecord["restfreq"]
                faxis = dataRecord["freq"]
                stokes = dataRecord["stokes"]
                source = dataRecord["source"]
                dataUnits = dataRecord["units"]
                calibType = dataRecord["calibtype"]
                veldef = dataRecord["veldef"]
                specsys = dataRecord["specsys"]
                coordType = (dataRecord["xctype"],dataRecord["yctype"])
                radesys = dataRecord["radesys"]
                equinox = dataRecord["equinox"]
                telescop = dataRecord["telescop"]
                frontend = dataRecord["frontend"]
                observer = dataRecord["observer"]
                dateobs = dataRecord["date-obs"]
                jdobs = dataRecord["jdobs"]
                projid = dataRecord["projid"]
                bw = dataRecord["bandwidth"]
                sitelong = dataRecord["sitelong"]
                sitelat = dataRecord["sitelat"]
                siteelev = dataRecord["siteelev"]
                uniqueScans = numpy.unique(dataRecord["scans"])

            else:
                xsky = numpy.append(xsky,dataRecord["xsky"])
                ysky = numpy.append(ysky,dataRecord["ysky"])
                stokes = numpy.append(stokes,dataRecord["stokes"])
                jdobs = numpy.append(jdobs,dataRecord["jdobs"])
                wt = numpy.append(wt,dataRecord["wt"])
                data = numpy.append(data,dataRecord["data"],axis=0)
                uniqueScans = numpy.unique(numpy.append(uniqueScans,dataRecord["scans"]))

        except(AssertionError):
            if args.verbose > 1:
                print "There was an unexpected problem processing %s" % thisFile
            raise

    if xsky is None:
        print "No valid data found, nothing to convert"
        return

    # set the reference sky position using the mean x and y positions
    # need to worry about points clearly off the grid
    #   e.g. a reference position incorrectly included in the data to be gridded.

    # this masks out antenna positions exactly equal to 0.0 - unlikely to happen
    # except when there is no valid antenna pointing for that scan.
    nonZeroXY = (xsky!=0.0) & (ysky!=0.0)

    # watch for the pathological case where there is no good antenna data
    # which can not be gridded at all
    if numpy.all(nonZeroXY == False):
        # always print this out, independent of verbosity level
        print "All antenna pointings are exactly equal to 0.0, can not grid this data"
        return

    if args.verbose > 3 and numpy.any(nonZeroXY == False):
        print "%d spectra will be excluded because the antenna pointing is exactly equal to 0.0 on both axes - unlikely to be a valid position" % (nonZeroXY == False).sum()
        
    # the reference sky position is just the mean, excluding zero antenna values
    refXsky = numpy.mean(xsky[nonZeroXY])
    refYsky = numpy.mean(ysky[nonZeroXY])

    f_delta = faxis[1]-faxis[0]
    f_crpix = len(faxis)/2 + 1  # FITS is 1-relative
    if int(f_crpix) == f_crpix:
        f_crval = faxis[int(f_crpix)-1]
    else:
        crpix_lower = numpy.floor(f_crpix)
        f_crval = (f_crpix-crpix_lower)*f_delta + faxis[crpix_lower-1]        

    # This is what idlToSdfits does (next 2 lines of code)
    # telescop diameter, in meters
    # !CHANGE GRIDDER to USE reference freq value not mean(faxis)
    diam = 100.0
    beam_fwhm = 1.2 * constants.c * (180.0/constants.pi) / (diam * f_crval)
    # the 747.6 and 763.8 values above are equivalent to diam of 99.3 and 97.2 m in this equation, respectively

    velref = None
    velrefOffset = None
    altrpix = None
    if veldef == "OPTI":
        velrefOffset = 0
    elif veldef == "RADI":
        velrefOffset = 256
    # no other veldef values are recognized by AIPS
    if velrefOffset is not None:
        if specsys == "LSRK":
            velref = 1 + velrefOffset
        elif (specsys == "HELIOCEN") or (specsys == "BARYCENT"):
            # not exactly the same, close enough, AIPS only recognize HEL
            velref = 2 + velrefOffset
        elif specsys == "TOPOCENT":
            velref = 3 + velrefOffset
        else:
            if args.verbose > 2:
                print "The AIPS VELREF convention does not recognize %s, VELREF and ALT* values will not be set" % specsys
    # no other specsys values are recognized by AIPS
    if velref is not None:         
        # ALT axis
        # altrval is always 0.0 here so ...
        # altrpix is channel where v == 0.0, when frest==f_frame independent of velocity defn
        altrpix = f_crpix + (frest-f_crval)/f_delta

    # image parameters
    # 6 cells per beam, integer number of ", back to degrees
    cellSize = numpy.ceil(beam_fwhm*3600.0/6.0)/3600.0  

    # image,
    # add 10% to 1/2 range and truncate
    # then add in 45" in units of cells, rounded up
    # multiple by 2 to ensure an even number
    # need to watch for wrap in coordinates, what does idlToSdfits do?
    # set a minimum size?
    # watch for a maximum size?
    xRange = xsky.max()-xsky.min()
    yRange = ysky.max()-ysky.min()
    xsize = (int((1.1*xRange/2.0)/cellSize) + int(numpy.ceil(45.0/3600.0/cellSize)))*2
    ysize = (int((1.1*yRange/2.0)/cellSize) + int(numpy.ceil(45.0/3600.0/cellSize)))*2

    dateobs = dateobs[0:10]
    jd0 = apTime.Time(dateobs,format='isot',scale='utc').jd


    print "nchan : ", nchan
    print "chanStart : ", chanStart
    print "chanStop : ", chanStop
    print "refXsky : ", refXsky
    print "refYsky : ", refYsky
    print "beam_fwhm : ", beam_fwhm
    print "f_delta : ", f_delta
    print "f_crpix : ", f_crpix
    print "f_crval : ", f_crval
    print "altrpix : ", altrpix
    print "velref  : ", velref
    print "cellSize: ", cellSize
    print "xsize   : ", xsize
    print "ysize   : ", ysize
    print "date-obs : ", dateobs
    print "jd0    : ", jd0
    print "jdobs[0]: ", jdobs[0]-jd0
    print "wt[0]  : ", wt[0]

    # assemble the output file
    npol = 1
    uvdata = numpy.zeros(shape=(len(data),1,1,len(data[0]),npol,3),dtype=numpy.float32)
    uvdata[:,0,0,:,0,0] = data
    # is there a faster way to do this?  Eventually the wt will be a vector an this will be unnecessary.
    for i in range(len(data)):
        uvdata[i,0,0,:,0,2] = wt[i]

    # ra param is the xsky values
    # dec param is the ysky values
    # beam param needs to be created : just set it all to 1 for now
    # DATE param is the jdobs+jd0 value
    # Is scaling even necessary?
    xsky /= numpy.degrees(1.0)
    ysky /= numpy.degrees(1.0)
    jdobs -= jd0
    gd = pyfits.GroupData(uvdata,parnames=['RA','DEC','BEAM','DATE'],
                          pardata=[xsky,ysky,numpy.zeros(len(data)),jdobs], 
                          bitpix=-32)
    hdu = pyfits.GroupsHDU(gd)

    # placeholder for correct values
    # idlToSdfits has LPC* values - appear to never be set, not available in sdfits
    hdu.header['ORIGIN'] = ('NRAO, Green Bank','WV 24944  304-456-2011')
    hdu.header['TELESCOP'] = telescop
    hdu.header['SITELONG'] = (sitelong,'Telescope location')
    hdu.header['SITELAT'] = (sitelat,'Telescope location')
    hdu.header['SITEELEV'] = (siteelev,'Telescope location')
    # physical units for the value, calibration type ("Ta","Tmb","Jy") in the comments
    hdu.header['BUNIT'] = (dataUnits,calibType)

    # add in axis description keywords
    # the basic data + weights
    hdu.header['CTYPE2'] = ('COMPLEX','Signal, Baseline, Weight')
    hdu.header['CRVAL2'] = 1.
    hdu.header['CDELT2'] = 1.
    hdu.header['CRPIX2'] = 1.
    hdu.header['CROTA2'] = 0.

    # The STOKES axis
    hdu.header['CTYPE3'] = ('STOKES','X,Y or R,L as -5,-6, or -1,-2')
    hdu.header['CRVAL3'] = -1.   # X, stand in here for real value when done properly
    hdu.header['CDELT3'] = -1.   # to get to Y if necessary
    hdu.header['CRPIX3'] = 1.
    hdu.header['CROTA3'] = 0.

    # Frequency axis
    hdu.header['CTYPE4'] = ('FREQ','frequency axis in %s' % specsys) 
    hdu.header['CRVAL4'] = f_crval
    hdu.header['CDELT4'] = f_delta
    hdu.header['CRPIX4'] = f_crpix
    hdu.header['CROTA4'] = 0.
    # See if AIPS does anything with this
    hdu.header['SPECSYS'] = specsys

    # ALT axis for velocity for use by AIPS
    # aide in human readability, replace this comment with one specific for this VELREF value
    hdu.header['VELREF'] = (velref,'1 LSR, 2 HEL, 3 OBS + 256 Radio')
    hdu.header['ALTRVAL'] = 0.
    hdu.header['ALTRPIX'] = (altrpix,'Channel where V=0 in VELREF frame')

    # rest frequency
    hdu.header['RESTFREQ'] = frest
    hdu.header['BANDWIDT'] = (bw,"bandwidth in telescope rest frame")

    # RA,DEC - center of image
    hdu.header['CTYPE5'] = ('RA','Right Ascension')
    hdu.header['CRVAL5'] = (refXsky,'Right Ascension of Center')
    hdu.header['CDELT5'] = (cellSize,'Desired Image Cellsize (degree)')
    hdu.header['CRPIX5'] = (xsize,'Desired Image X Size (pixels)')
    hdu.header['CROTA5'] = 0.

    hdu.header['CTYPE6'] = ('DEC','Declination')
    hdu.header['CRVAL6'] = (refYsky,'Declination of Center')
    hdu.header['CDELT6'] = (cellSize,'Desired Image Cellsize (degree)')
    hdu.header['CRPIX6'] = (ysize,'Desired Image Y Size (pixels)')
    hdu.header['CROTA6'] = 0.

    # EPOCH for equatorial
    # try other coordinate systems to see how AIPS behaves without the fudge we do now
    hdu.header['EPOCH'] = 2000.0   # use real value here eventually
    # same here - is this really the WCS appropriate way to do this?
    hdu.header['COORDTYP'] = ('J2000','Coordinate Type')

    hdu.header['OBJECT'] = source

    # are these necessary
    # idlToSdfits appears to use the center sky frequency here
    hdu.header['IMCLASS'] = ('%s' % int(numpy.round(f_crval/1.e6)),'Class is Center Freq (MHz)')
    hdu.header['IMNAME'] = (source,'Name is source name')

    # actually check whether SORTORD is 'TB' or ''
    hdu.header['SORTORD'] = ('','Data is NOT time sorted')
    hdu.header['DATE-OBS'] = (dateobs,'UTC Date of observation')
    # use actual first scan seen
    hdu.header['SCAN'] = (uniqueScans[0],'GBT M&C Scan Number for first spectrum') 
    hdu.header['OBSERVER'] = observer
    # use actual projid
    hdu.header['PROJID'] = (projid,"Project ID")
    # this appears to NOT be picking up the receiver - idlToSdfits does
    hdu.header['INSTRUME'] = (frontend,'Front End')
    # ned to find a way to record the back end
    # idlToSdfits includes REFRACT, DAXMODEL, DELMODEL, not needed

    hdu.header['BMAJ'] = (beam_fwhm,'Angular Resolution Estimate (FWHM Degrees)')
    hdu.header['BMIN'] = (beam_fwhm,'Angular Resolution Estimate (FWHM Degrees)')
    hdu.header['BPA'] = 0.
    hdu.header['CELLSIZE'] = (round(cellSize*3600.),'Desired Image Cellsize (arcsec)')

    # HISTORY cards, including receiver and samplers
    hdu.header['HISTORY'] = ('sdfToAips v%s' % sdfToAipsVersion,'conversion tool and version number')

    hdu.writeto(args.output)

    del hdu

    fout = pyfits.open(args.output,mode='update')
    fout[0].header.insert('PTYPE1',('PSCAL1',numpy.degrees(1.0),'Degrees per radian'),after=True)
    fout[0].header.insert('PTYPE2',('PSCAL2',numpy.degrees(1.0),'Degrees per radian'),after=True)
    fout[0].header.insert('PTYPE4',('PZERO4',jd0,'JD of DATE-OBS'),after=True)
    fout.flush()
    fout.close()
    del fout

if __name__ == "__main__":
    
    args = read_command_line(sys.argv)

    # argument checking - perhaps this should be a separate function
    if args.mintsys is not None and args.mintsys < 0:
        print "mintsys must be > 0"
        sys.exit(1)

    if args.maxtsys is not None and args.maxtsys < 0:
        print "maxtsys must be > 0"
        sys.exit(1)

    if args.maxtsys is not None and args.mintsys is not None and args.maxtsys <= args.mintsys:
        print "maxtsys must be > mintsys"
        sys.exit(1)

    # set default output file here - eventually

    # if output function exist, remove if clobber else an error
    if os.path.exists(args.output):
        if not args.clobber:
            if args.verbose > 1:
                print args.output + ' exists, remove or use --clobber to overwrite'
            sys.exit(1)
        os.remove(args.output)
        if args.verbose > 3:
            print "existing " + args.output + " removed"

    # all of the sdfits files must exist
    for sdf in args.SDFITSfiles:
        if not os.path.exists(sdf):
            if arg.verbose > 1:
                print sdf + ' does not exist'
            sys.exit(1)

    try:
        sdfToAips(args)
    except ValueError, msg:
        if args.verbose > 1:
            print 'ERROR : ', msg
        sys.exit(-1)
    

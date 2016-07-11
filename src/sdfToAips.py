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

from gbt_stokes_iter import gbt_stokes_iter
from get_metadata import get_metadata
from resize_groups import resize_groups
from get_faxis import get_faxis

import astropy.io.fits as pyfits
import numpy
from scipy import constants
import astropy.time as apTime
import argparse
import sys
import os
import math

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
                        help="Optionally average channels, keeping only nchan/average channels")
    parser.add_argument("-s","--scans", type=str, 
                        help="Only use data from these scans.  comma separated list or <start>:<end> range syntax or combination of both")
    parser.add_argument("-w","--window", type=int, default=-1,
                        help="Data is always selected from one spectral window (IFNUM), this allows you to specify which one.  Defaults to smallest IFNUM found.")
    parser.add_argument("-m","--maxtsys", type=float,
                        help="max Tsys value to use")
    parser.add_argument("-z","--mintsys", type=float,
                        help="min Tsys value to use")
    parser.add_argument("--hdu",type=int,
                        help="Only find data at this HDU in all FITS files.  Counts from 1",default=None)
    parser.add_argument("SDFITSfiles", type=str, nargs="+",
                        help="The calibrated SDFITS files to use.")
    parser.add_argument("--clobber", default=False, action="store_true",
                        help="Overwrites existing output file if set.")
    parser.add_argument("-1",action='store_true',default=False,
                        dest="singlestokes",
                        help="Force single stokes structure even when more than one polarization is seen")

    # what does idlToSdfits use for output when not specified?  Eventually make this optional.
    parser.add_argument("-o","--output",type=str,required=True,
                        help="output name")
    parser.add_argument("-v","--verbose", type=int, default=4,
                        help="set the verbosity level-- 0-1: none, "
                        "2:errors only, 3:+warnings, "
                        "4(default):+user info, 5:+debug")
    parser.add_argument("-V","--version",action="version",version="sdfToAips version: %s" % sdfToAipsVersion)

    # this argument will likely be removed in production.  Used here during testing
    parser.add_argument("-g","--groupbuff",type=int,help="Group buffer size on output", default=200)
    
    args = parser.parse_args()

    return args


# parse_channels, parse_scans, and format_scans are all as they are in gbtgridder
# find a way to share this code

def parse_channels(channelString,verbose=4):
    "Turn a valid channel range into start and end channels"
    start = None
    end = None
    if channelString is not None:
        if verbose > 4:
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

    groupBuff = args.groupbuff
    if groupBuff <= 0:
        print "Usage error: -g,--groupbuff must be > 0"
        sys.exit(1)

    if args.verbose > 3:
        print "Loading data ... "

    # this does NOT pull in all of the data, just the few columns necessary to procede
    metaData = get_metadata(args.SDFITSfiles,scanlist=scanlist,ifnum=args.window,hdu=args.hdu,verbose=args.verbose)
    if len(metaData) == 0:
        print "There was a problem with the data, nothing to convert"
        sys.exit(-1)

    if args.verbose > 3:
        print "Processing ..."

    if chanStop is None or chanStop >= metaData['nchan']:
        chanStop = metaData['nchan']-1

    uniqueScans = numpy.unique(metaData['scan'])
    uniqueScans.sort()

    xsky = metaData['crval2']
    ysky = metaData['crval3']

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
        print "ERROR: All antenna pointings are exactly equal to 0.0, can not grid this data"
        return

    if args.verbose > 3 and numpy.any(nonZeroXY == False):
        print "%d spectra will be excluded because the antenna pointing is exactly equal to 0.0 on both axes - unlikely to be a valid position" % (nonZeroXY == False).sum()
        
    # the reference sky position is just the mean, excluding zero antenna values
    # this will not work if the X positions are near the coordinate wrap
    refXsky = numpy.mean(xsky[nonZeroXY])
    refYsky = numpy.mean(ysky[nonZeroXY])

    # frequency axis using the first row values.
    # code farther down warns the user when things are "too far" from this axis
    # definitio of "too far" needs to be refined

    # NOTE: get_faxis wants chanStart and chanStop in FITS terms, 1-relative
    # but they are 0-relative here, hence the +1 in the call below.
    nchan = metaData['nchan']
    (f_crval, f_delta, f_crpix, nout, bw) = get_faxis(metaData['crval1'][0],
                                                      metaData['crpix1'][0],
                                                      metaData['cdelt1'][0],
                                                      metaData['vframe'][0],
                                                      nchan,
                                                      average,
                                                      chanStart+1,
                                                      chanStop+1)

    # This is the formula that idlToSdfits uses for beam_fwhm
    # diam is telescop diameter, in meters

    # idlToSdfits uses the doppler shifted frequency, but sky freq should be used
    # not significantly different for how it's intended to be used - but when
    # comparing the results with idlToSdfits files that might be noticed
    diam = 100.0
    beam_fwhm = 1.2 * constants.c * (180.0/constants.pi) / (diam * metaData['crval1'][0])

    velref = None
    velrefOffset = None
    altrpix = None
    if metaData['veldef'] == "OPTI":
        velrefOffset = 0
    elif metaData['veldef'] == "RADI":
        velrefOffset = 256
    # no other veldef values are recognized by AIPS
    if velrefOffset is not None:
        if metaData['specsys'] == "LSRK":
            velref = 1 + velrefOffset
        elif (metaData['specsys'] == "HELIOCEN") or (metaData['specsys'] == "BARYCENT"):
            # not exactly the same, close enough, AIPS only recognize HEL
            velref = 2 + velrefOffset
        elif specsys == "TOPOCENT":
            velref = 3 + velrefOffset
        else:
            if args.verbose > 2:
                print "The AIPS VELREF convention does not recognize %s, VELREF and ALT* values will not be set" % metaData['specsys']
    # no other specsys values are recognized by AIPS
    if velref is not None:         
        # ALT axis
        # altrval is always 0.0 here so ...
        # altrpix is channel where v == 0.0, when restfreq==f_frame independent of velocity defn
        altrpix = f_crpix + (metaData['restfreq'][0]-f_crval)/f_delta

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

    dateobs = metaData['date-obs'][0]
    dateobs = dateobs[0:10]
    jd0 = apTime.Time(dateobs,format='isot',scale='utc').jd

    # STOKES axis
    stokesVals = metaData['crval4']
    uniqStokes = numpy.unique(stokesVals)
    uniqStokes.sort()
    nstokes = len(uniqStokes)
    crvalStokes = -1
    cdeltStokes = -1

    if args.singlestokes:
        # ignore actual stokes values and force into single stokes array
        nstokes = 1
        iStokes = numpy.zeros_like(stokesVals)
        # this should use the largest stokes value, which given that these are likely < 0
        # will be the smallest absolute stokes value, so XX if XX, YY are used or LL or LL, RR are used.
        crvalStokes = uniqStokes[-1]
        if len(uniqStokes) > 1 and args.verbose > 2:
            print "WARNING: Multi-stokes data was found.  Will be forced to a single single STOKES axis, labelled as : ", crvalStokes
    else:
        # if these are sequential values, good to go
        if abs(uniqStokes[0]-uniqStokes[-1]) == nstokes-1:
            # this should be the norm
            if uniqStokes[0] < 0:
                crvalStokes = uniqStokes[-1]
                cdeltStokes = -1
            else:
                crvalStokes = uniqStokes[0]
                cdeltStokes = 1
            # appropriate 0-relative index into stokes axis by data row number
            iStokes = (stokesVals - crvalStokes)/cdeltStokes
        else:
            # always print out this warning
            if args.verbose > 2:
                print "WARN: Unexpected polarization values, ignoring polarization information in output UVFITS. Labelling as all XX"
            nstokes = 1
            iStokes = numpy.zeros_like(stokesVals)

    # set up the data iterator using these stokes values
    data_iter = gbt_stokes_iter(args.SDFITSfiles,uniqStokes,scanlist=scanlist,ifnum=metaData['ifnum'],hdu=args.hdu,verbose=args.verbose)


    if args.verbose > 4:
        # very useful during debugging, not so much most of the time
        print "nchan : ", metaData['nchan']
        print "chanStart : ", chanStart
        print "chanStop : ", chanStop
        print "actual nchan (after selection and averaging): ", nout
        print "refXsky : ", refXsky
        print "refYsky : ", refYsky
        print "beam_fwhm : ", beam_fwhm
        print "f_delta : ", f_delta
        print "f_crpix : ", f_crpix
        print "f_crval : ", f_crval
        print "altrpix : ", altrpix
        print "velref  : ", velref
        print "stokes crval : ", crvalStokes
        print "stokes increment : ", cdeltStokes
        print "n stokes : ", nstokes
        print "cellSize: ", cellSize
        print "xsize   : ", xsize
        print "ysize   : ", ysize
        print "date-obs : ", metaData['date-obs'][0]
        print "jd0    : ", jd0
        print "jdobs[0]: ", metaData['jdobs'][0]-jd0

    # assemble the output file

    # guess at the number of output rows.  This isn't perfect.
    # for now, will warn and skip if this isn't sufficient eventually will beed to extend
    if nstokes == 1:
        nOutrows = len(stokesVals)
    else:
        nOutrows = 0
        for pol in uniqStokes:
            nOutrows = max(nOutrows,(stokesVals==pol).sum())

    # set up the initial size and empty contents of the file
    # work out the most appropriate groupBuff value given nOutrows
    groupBuff = min(groupBuff,nOutrows)
    nGroups = round(float(nOutrows)/groupBuff)
    groupBuff = int(math.ceil(float(nOutrows)/nGroups))

    if args.verbose > 4:
        print "Using groupBuff size : ", groupBuff
        print "nOutrows : ", nOutrows
        print "nGroups : ", nGroups

    # the data array
    # shape : (#records, 1 (field RA center), 1 (field DEC center), 
    #          nchan, nstokes, nparam [3 = data, baseline==0, wt])
    # associated group fields are RA, DEC, FEED, JD, same number of output rows
    uvdata = numpy.zeros(shape=(groupBuff,1,1,nout,nstokes,3),dtype=numpy.float32)

    # vector used to hold the weights for each row
    vectorWts = numpy.ones(shape=(nout),dtype=numpy.float32)

    # arrays to put the group values, one per row.
    outXsky = numpy.zeros(groupBuff)
    outYsky = numpy.zeros(groupBuff)
    outBeam = numpy.zeros(groupBuff)
    outRelJdObs = numpy.zeros(groupBuff)

    # write it out first, actually fill in values later.  Will be used when buffering
    # output shortly.
    gd = pyfits.GroupData(uvdata,parnames=['RA','DEC','BEAM','DATE'],
                          pardata=[outXsky,outYsky,outBeam,outRelJdObs], 
                          bitpix=-32)
    hdu = pyfits.GroupsHDU(gd)

    currentSize = groupBuff

    # Missing from original idlToSdfits has
    #    LPC* values - appear to never be set by idlToSdfits
    #                  not available in sdfits, omit here
    #    REFRACT, DAXMODEL, DELMODEL, not needed

    # Origin and telescope location
    hdu.header['ORIGIN'] = ('NRAO, Green Bank','WV 24944  304-456-2011')
    hdu.header['TELESCOP'] = metaData['telescop'][0]
    hdu.header['SITELONG'] = (metaData['sitelong'],'Telescope location')
    hdu.header['SITELAT'] = (metaData['sitelat'],'Telescope location')
    hdu.header['SITEELEV'] = (metaData['siteelev'],'Telescope location')

    # physical units for the value, calibration type ("Ta","Tmb","Jy") in the comments
    hdu.header['BUNIT'] = (metaData['units'],metaData['calibtype'])

    # axis description keywords

    # the basic data + weights
    hdu.header['CTYPE2'] = ('COMPLEX','Signal, Baseline, Weight')
    hdu.header['CRVAL2'] = 1.
    hdu.header['CDELT2'] = 1.
    hdu.header['CRPIX2'] = 1.
    hdu.header['CROTA2'] = 0.

    # The STOKES axis
    hdu.header['CTYPE3'] = ('STOKES','X,Y or R,L as -5,-6, or -1,-2')
    # some readers complain if these values aren't floats
    hdu.header['CRVAL3'] = float(crvalStokes)
    hdu.header['CDELT3'] = float(cdeltStokes)
    hdu.header['CRPIX3'] = 1.
    hdu.header['CROTA3'] = 0.

    # Frequency axis
    hdu.header['CTYPE4'] = ('FREQ','frequency axis in %s' % metaData['specsys']) 
    hdu.header['CRVAL4'] = (f_crval,'Frequency at reference pixel (Hz)')
    hdu.header['CDELT4'] = (f_delta,'Channel spacing (Hz)')
    hdu.header['CRPIX4'] = (f_crpix,'Reference pixel')
    hdu.header['CROTA4'] = 0.
    # See if AIPS does anything with this
    hdu.header['SPECSYS'] = metaData['specsys']
    # ALT axis for velocity for use by AIPS
    # aid for human readability
    #   Todo: replace this comment with one specific for this VELREF value
    hdu.header['VELREF'] = (velref,'1 LSR, 2 HEL, 3 OBS + 256 Radio')
    hdu.header['ALTRVAL'] = 0.
    hdu.header['ALTRPIX'] = (altrpix,'Channel where V=0 in VELREF frame')

    # rest frequency and total bandwidth
    hdu.header['RESTFREQ'] = (metaData['restfreq'][0], "line rest frequency (Hz)")
    hdu.header['BANDWIDT'] = (bw,"bandwidth in telescope rest frame (Hz)")

    # Desired output map - center RA, DEC, image size, cell size
    #   RA/DEC used here even though actual coordinate system may differ
    #   ToDo: work out how best to convey original sky coordinate system
    hdu.header['CTYPE5'] = ('RA','Right Ascension')
    hdu.header['CRVAL5'] = (refXsky,'Right Ascension of Center (deg)')
    hdu.header['CDELT5'] = (cellSize,'Desired Image Cellsize (deg)')
    hdu.header['CRPIX5'] = (xsize,'Desired Image X Size (pixels)')
    hdu.header['CROTA5'] = 0.

    hdu.header['CTYPE6'] = ('DEC','Declination')
    hdu.header['CRVAL6'] = (refYsky,'Declination of Center (deg)')
    hdu.header['CDELT6'] = (cellSize,'Desired Image Cellsize (deg)')
    hdu.header['CRPIX6'] = (ysize,'Desired Image Y Size (pixels)')
    hdu.header['CROTA6'] = 0.

    # EPOCH for equatorial
    # try other coordinate systems to see how AIPS behaves without the fudge we do now
    #   ToDo:  need real EPOCH and COORDTYP here when appropriate
    hdu.header['EPOCH'] = 2000.0   # use real value here eventually
    #   Q: same here - is this really the WCS appropriate way to do this?
    hdu.header['COORDTYP'] = ('J2000','Coordinate Type')

    # source name
    hdu.header['OBJECT'] = metaData['object'][0]

    # idlToSdfits appears to use the center sky frequency here
    #   Q: how does AIPS actually use this value?
    hdu.header['IMCLASS'] = ('%s' % int(numpy.round(f_crval/1.e6)),'Class is Center Freq (MHz)')
    #   Q: how does AIPS use this value?
    hdu.header['IMNAME'] = (metaData['object'][0],'Name is source name')

    # ToDo: actually check whether SORTORD is 'TB' or '' and use 'TB' here when
    #       appropriate. It might help AIPS be a bit faster.
    #       it's probably only likely to be sorted by time if there's just a single input SDFITS file
    hdu.header['SORTORD'] = ('','Data is NOT time sorted')

    # informational
    hdu.header['DATE-OBS'] = (dateobs,'UTC Date of observation')
    # use actual first scan seen
    hdu.header['SCAN'] = (uniqueScans[0],'GBT M&C Scan Number for first spectrum') 
    hdu.header['OBSERVER'] = metaData['observer'][0]
    # use actual projid
    hdu.header['PROJID'] = (metaData['projid'][0],"Project ID")
    # this appears to NOT be picking up the receiver - idlToSdfits does
    hdu.header['INSTRUME'] = (metaData['frontend'][0],'Front End')
    # ToDo: need to find a way to record the backend name

    # approximate beam information
    # assumes circular
    hdu.header['BMAJ'] = (beam_fwhm,'Angular Resolution Estimate (FWHM Degrees)')
    hdu.header['BMIN'] = (beam_fwhm,'Angular Resolution Estimate (FWHM Degrees)')
    hdu.header['BPA'] = 0.
    # Q: why this AND the CDELT5 and CDELT6 values above?
    hdu.header['CELLSIZE'] = (round(cellSize*3600.),'Desired Image Cellsize (arcsec)')

    # HISTORY cards, including receiver and samplers
    hdu.header.add_history(('sdfToAips v%s' % sdfToAipsVersion) + ' conversion tool and version number')
    if args.channels is None:
        hdu.header.add_history('sdfToAips no channel selection')
    else:
        hdu.header.add_history('sdfToAips channels: '+args.channels)
    if args.average is not None and average > 1:
        hdu.header.add_history('sdfToAips average: %s channels' % average)
    else:
        hdu.header.add_history('sdfToAips no channel averaging')
    if args.scans is None:
        hdu.header.add_history('sdfToAips no scan selection')
    else:
        hdu.header.add_history('sdfToAips scans ; '+args.scans)
    hdu.header.add_history('sdfToAips used ifnum=%d' % metaData['ifnum'])
    if args.mintsys is None and args.maxtsys is None:
        hdu.header.add_history('sdfToAips no tsys selection')
    else:
        if args.mintsys is not None:
            hdu.header.add_history('sdfToAips mintsys : %f' % args.mintsys)
        if args.maxtsys is not None:
            hdu.header.add_history('sdfToAips maxtsys : %f' % args.maxtsys)
    if args.singlestokes:
        hdu.header.add_history('sdfToAips single STOKES axis forced')
    hdu.header.add_history('sdfToAips sdfits files ...')
    for thisFile in args.SDFITSfiles:
        # protect against long file names - don't use more than one HISTORY row
        # to document this.  80 chars total, 8 for "HISTORY ", 11 for "sdfToAips: "
        # leaving 61 for the file name
        if len(thisFile) > 61:
            thisFile = "*"+thisFile[-60:]
        hdu.header.add_history("sdfToAips: " + thisFile)
    # FITS defininition boilerplate reference
    hdu.header.add_comment("  FITS (Flexible Image Transport System) format is defined in 'Astronomy")
    hdu.header.add_comment("  and Astrophysics', volume 376, page 359; bibcode: 2001A&A...376..359H")
    
    hdu.writeto(args.output)

    del hdu

    # creating then closing and reopening it is how we avoid having to create the entire
    # group buffer in full at file creation time

    # dictionary mapping tuple : (relJdObs,feed) to outRow, remembering also which crvals have already been used
    # so (outRow,usedStokes) = rowMap[(relJdObs[i],feed[i])]
    # where rowList is a list of (outRow,usedStokes) for the given key.  There will normally be one element
    # in that list except when a there are multiple instances of the same stokes for a given time and feed.
    # In that case elements are added as needed.  Elements are popped off when all usedStokes are true.
    # The entire key is removed if there are no longer an elements in that list.
    # if (relJdObs[i],feed[i]) not found or iStokes[i] is True in usedStokes, add a new row for that tuple

    rowMap = {}

    usedStokesTemplate = numpy.zeros(nstokes,dtype=bool)
    outRowCount = 0
    nTsysFlagged = 0.0

    fout = pyfits.open(args.output,mode='update',memmap=True)
    phdu = fout[0]

    # set up boxcar convolution as appropriat
    box = None
    if average is not None and average > 1:
        box = numpy.ones(average,'float32')/average

    chanSpaceWarning = False
    chanCenterWarning = False

    for row in data_iter:
        data = row['DATA'][chanStart:(chanStop+1)]
        if box is not None:
            smoothedData = numpy.convolve(box,data,mode='valid')
            # decimate and save
            data = smoothedData[0:-1:average]

        # compare frequency axes - only if not already warned
        # and only care about them if the verbosity level is > 3, else just plow ahead
        if args.verbose > 3 and (not chanSpaceWarning or not chanCenterWarning):
            (thisCrval, thisDelta, thisCrpix, thisNout, thisBW) = \
                get_faxis(row['CRVAL1'],row['CRPIX1'],row['CDELT1'],row['VFRAME'],nchan,average,chanStart+1,chanStop+1)
            deltaDiff = (thisDelta-f_delta)/f_delta

            chan1diff = (thisCrval-f_crval)/f_delta + deltaDiff - thisCrpix*thisDelta/f_delta + f_crpix
            chanNdiff = (thisCrval-f_crval)/f_delta + nchan*deltaDiff - thisCrpix*thisDelta/f_delta + f_crpix

            if not chanSpaceWarning and abs(deltaDiff) > 0.050:
                print "WARNING: channel spacings differ by more than 0.5%."
                chanSpaceWarning = True

            if not chanCenterWarning and (abs(chan1diff) > 0.1 or abs(chanNdiff) > 0.1):
                print "WARNING: channel centers differ by more than 10%"
                chanCenterWarning = True
        
        # replace NaNs with zeros - mask also used to set weights
        # should translation of NaNs be done before or after averaging?
        nanDataMask = numpy.isnan(data)
        data[nanDataMask] = 0.0

        # other fields needed from row
        # scaling of the group values is done here and the header updated later
        # scale degrees to radians
        thisXsky = row['CRVAL2'] / numpy.degrees(1.0)
        thisYsky = row['CRVAL3'] / numpy.degrees(1.0)
        # use relative JD by subtracting jd0 value
        thisDateobs = row['DATE-OBS']
        thisJDrel = apTime.Time(thisDateobs,format='isot',scale='utc').jd - jd0
        thisBeam = row['FEED']
        thisStokes = row['CRVAL4']
        # convert to pixel along output STOKES axis - this is zero relative, which is appropriate here
        if nstokes == 1:
            thisIstokes = 0
        else:
            thisIstokes = (thisStokes - crvalStokes)/cdeltStokes

        # work out which row this goes into
        outRow = None
        if nstokes > 1:
            key = (thisJDrel,thisBeam)
            keyKnown = False
            if key in rowMap:
                # this row already exists, output file must be big enough already
                rowMapList = rowMap[key]
                keyKnown = True
                toRemove = []
                # loop to find first list that this can be added to - stokes not already seen
                for listIndex in range(len(rowMapList)):
                    thisOutRow, usedStokes = rowMapList[listIndex]
                    if not usedStokes[thisIstokes]:
                        outRow = thisOutRow
                        usedStokes[thisIstokes] = True
                        if numpy.all(usedStokes):
                            # this item can be deleted
                            toRemove.append(listIndex)
                        else:
                            # update this list
                            rowMapList[listIndex] = (outRow,usedStokes)
                            rowMapList[key] = rowMapList
                        # either way, we're done here
                        break
                    # otherwise this is already used - move to next item to consider
                # remove any items, reverse order of toRemove
                toRemove.reverse()
                for listIndex in toRemove:
                    del rowMapList[listIndex]

                # remove the list associated with key if it's empty
                if len(rowMapList) == 0:
                    del rowMap[key]

        if outRow is None:
            # this row does not yet exist
            outRow =  outRowCount
            outRowCount += 1
            # is the output file big enough
            if outRowCount > currentSize:
                # no - close the file first
                del phdu
                fout.close()
                if outRowCount > nOutrows:
                    # it's requesting a row past where it estimated it would finish
                    # just keep growing by groupBuff, will need to trim at the end
                    currentSize += groupBuff
                else:
                    # add on groupBuff rows to not exceed nOutrows
                    currentSize = min(currentSize+groupBuff,nOutrows)

                if args.verbose > 4:
                    print "expanding output file to %d, next row is %d" % (currentSize, outRow)
                resize_groups(args.output,currentSize)

                # and reopen the file and reset phdu
                fout = pyfits.open(args.output,mode='update',memmap=True)
                phdu = fout[0]

            # it's only necessary to set the group parameters once for each outRow
            phdu.data['RA'][outRow] = thisXsky
            phdu.data['DEC'][outRow] = thisYsky
            phdu.data['BEAM'][outRow] = thisBeam
            phdu.data['DATE'][outRow] = thisJDrel

            # remember this in any case
            if nstokes > 1:
                # only need to do this for the multi-stokes case
                usedStokes = usedStokesTemplate.copy()
                usedStokes[thisIstokes] = True
                if keyKnown:
                    # append it to the end of the current list
                    rowMapList = rowMap[key]
                    rowMapList.append((outRow,usedStokes))
                else:
                    # new list
                    rowMap[key] = [(outRow,usedStokes)]

        # weight
        thisTsys = row['TSYS']
        # default is to have wt be 0.0
        scalarWt = 0.0
        # do not use Tsys == 0.0 values
        if thisTsys != 0.0:
            # tsys flagging
            if (args.mintsys is None or thisTsys >= args.mintsys) and (args.maxtsys is None or thisTsys <= args.maxtsys):
                # no Tsys flagging or the Tsys is in the specified range, OK to use it
                # normalize to Tsys = 25.0 as per idlToSdfits
                relTsys = thisTsys / 25.0
                scalarWt = row['EXPOSURE']/(relTsys*relTsys)
            else:
                # Tsys flagging, Tsys out of specified range
                # scalarWt already set to 0.0
                nTsysFlagged += 1

        # finally set the data weights
        vectorWts[:] = scalarWt
        vectorWts[nanDataMask] = 0.0

        phdu.data['DATA'][outRow,0,0,:,thisIstokes,0] = data
        phdu.data['DATA'][outRow,0,0,:,thisIstokes,2] = vectorWts

    fout.close()
    # probably overkill
    del fout
    del phdu

    # trim if necessary
    if outRowCount < currentSize:
        if args.verbose > 4:
            print "Trimming output file to final size: %d, was %d" % (outRowCount,currentSize)
        resize_groups(args.output,outRowCount)

    fout = pyfits.open(args.output,mode='update',memmap=True)
    fout[0].header.insert('PTYPE1',('PSCAL1',numpy.degrees(1.0),'Degrees per radian'),after=True)
    fout[0].header.insert('PTYPE2',('PSCAL2',numpy.degrees(1.0),'Degrees per radian'),after=True)
    fout[0].header.insert('PTYPE4',('PZERO4',jd0,'JD of DATE-OBS'),after=True)
    if args.mintsys is not None or args.maxtsys is not None:
        msg = 'sdfToAips N spectra outside tsys range: %d' % nTsysFlagged
        fout[0].header.add_history(msg)
        if args.verbose > 3:
            print msg

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
    

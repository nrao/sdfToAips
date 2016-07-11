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

import numpy
import fitsio
import astropy.time as apTime

def getSDFitsColumns(sdfitsFileList,colnames,scanlist,ifnum):
    # colnames is a list or an array of column names
    # this follows the Green Bank convention
    # if a given name is not found in the table and exists as a keyword
    # that value will be used as a constant valued column 
    # if not found at all that name will be missing in the returned dictionary
    # DATA_TUNIT is translated to the TUNIT for the DATA column
    # DATA_WIDTH is translated to the width (number of elements) in the DATA column
    result = {}
    rowCount = 0
    for file in sdfitsFileList:
        thisFIO = fitsio.FITS(file)
        for i in range(1,len(thisFIO)):
            if thisFIO[i].get_extname() == "SINGLE DISH":
                # must recheck columns for each new HDU
                # translations for use later in reconstructing result
                dataTunitName = None
                hduColnames = thisFIO[i].get_colnames()
                thisHdr = thisFIO[i].read_header()
                colList = []
                hdrValues = {}
                colValues = None
                # FITS columns count from 1
                dataColNum = hduColnames.index('DATA')+1
                nrows = thisFIO[i].get_nrows()
                for col in colnames:
                    if col == 'DATA_TUNIT':
                        dataTunitName = 'TUNIT%d' % dataColNum
                        col = dataTunitName
                    # we'll get DATA_WIDTH later
                    if col in hduColnames:
                        colList.append(col)
                    else:
                        if col in thisHdr:
                            hdrValues[col] = thisHdr[col]
                        # otherwise silently ignored
                if scanlist is not None:
                    # SCAN is required, even when not requested
                    if 'SCAN' not in colnames:
                        colList.append('SCAN')
                if len(colList) > 0:
                    colValues = thisFIO[i][colList][:]
                    # SCAN selection happens here
                    if scanlist is not None:
                        rowMask = numpy.full(colValues.shape,True,dtype=bool)
                        allScans = colValues['SCAN']
                        uniqueScans = numpy.unique(allScans)
                        for scan in uniqueScans:
                            if scan not in scanlist:
                                rowMask[allScans==scan] = False
                        colValues = colValues[rowMask]
                        nrows = len(colValues)

                    # IFNUM selection is always done
                    # unless IFNUM is not present (older data) in which case there
                    # will be no IFNUM returned in result and no selection can happen here
                    if 'IFNUM' in colList:
                        if not ifnum >= 0:
                            # no value or negative ifnum supplied
                            ifnum = min(colValues['IFNUM'])
                        colValues = colValues[colValues['IFNUM']==ifnum]
                        nrows = len(colValues)

                    for col in colList:
                        outcol = col
                        if col == dataTunitName:
                            outcol = 'DATA_TUNIT'
                        if col == 'SCAN' and 'SCAN' not in colnames:
                            # this is where we drop SCAN from the results if not requested
                            continue
                        if rowCount == 0:
                            # nothing there yet, just add to result as is
                            result[outcol] = colValues[col]
                        else:
                            # more complicated
                            if outcol in result:
                                # append
                                result[outcol] = numpy.append(result[outcol],colValues[col])
                            else:
                                # it's there in this HDU but was not seen in previous ones
                                # refuse to handle that case - ValueError
                                raise ValueError("%s column not found in all SDFITS tables" % col)

                for kw in hdrValues:
                    val = hdrValues[kw]
                    outkw = kw
                    if kw == 'DATA_TUNIT':
                        outkw = 'DATA_TUNIT'
                    if type(val) == str:
                        # support strings up to 32 characters here 
                        colVal = numpy.empty(nrows,dtype='S32')
                    else:
                        colVal = numpy.empty(nrows,dtype=type(val))
                    colVal.fill(hdrValues[kw])
                    if rowCount == 0:
                        result[outkw] = colVal
                    else:
                        result[outkw] = numpy.append(result[kw],colVal)

                if 'DATA_WIDTH' in colnames:
                    dataWidth = int(thisHdr['TFORM%d'%dataColNum].strip()[:-1])
                    colVal = numpy.empty(nrows,dtype=type(dataWidth))
                    colVal.fill(dataWidth)
                    if rowCount == 0:
                        result['DATA_WIDTH'] = colVal
                    else:
                        result['DATA_WIDTH'] = numpy.append(result['DATA_WIDTH'],colVal)

                # check to make sure all fields in result have the same length == rowCount
                rowCount += nrows
                for field in result:
                    if len(result[field]) != rowCount:
                        raise ValueError("%s column/keyword not found in all SDFITS tables" % field)
        thisFIO.close()
    return result        


def get_metadata(sdfitsFileList, scanlist=None, ifnum=None, verbose=4):
    """
    Given an sdfits file, return the associated meta information
    for the list of SDFITS files for the given scans and for data
    where the TSYS values are between the supplied min and max.
    verbose is desired output level: 0-1:none, 2:errors, 3:+warnings, 4(default):+user, 5:+debug
    """
    result = {}

    # no checks yet on CTYPE1 value (must be frequency, assumes OBS frame for now)
    # and CTYPE4 value (which must be 'STOKES')
    # DATA_TUNIT is translated in getSDFITSColumns to appropriate TUNIT<colnum>
    # string given the location of the DATA column.
    # DATA_WIDTH contains the number of elements in the DATA column for each row.

    columns = ['SCAN','DATE-OBS','TSYS','EXPOSURE','DATA_TUNIT','DATA_WIDTH',
               'CRVAL1','CDELT1','CRPIX1','CTYPE1','VFRAME','RESTFREQ',
               'CRVAL2','CTYPE2','CRVAL3','CTYPE3','CRVAL4','CTYPE4',
               'RADESYS','EQUINOX','OBJECT',
               'VELDEF','TELESCOP','FRONTEND','OBSERVER','PROJID',
               'SITELONG','SITELAT','SITEELEV','FEED','IFNUM']

    # need to add ability to ensure that the NCHAN value is the same across all tables
    colValues = getSDFitsColumns(sdfitsFileList,columns,scanlist,ifnum)

    if len(colValues) == 0:
        if verbose > 1:
            print "No single dish FITS rows were found in those files."
        return result

    if len(colValues['DATA_WIDTH']) == 0:
        if verbose > 3:
            print "No data found for the list of selected scan numbers and spectral window."
        return result

    # copy columns to result
    for col in colValues:
        # sanity check to ensure that columns expected to be constant are constant
        if col in ['DATA_WIDTH','DATA_TUNIT','CTYPE1','CTYPE2','CTYPE3','CTYPE4',
                   'RADESYS','EQUINOX','VELDEF','IFNUM',
                   'SITELONG','SITELAT','SITEELEV']:
            # OBJECT, TELESCOP, FRONTEND, OBSERVER, and PROJID are only taken from the 
            # first row, but they're only used informationally
            # and so it's not an error if they vary.  The constant SITE*
            # keywords should ensure that TELESCOP is the same, but we don't check that.
            if not(numpy.all(colValues[col]==colValues[col][0]) or numpy.all(numpy.isnan(colValues[col]))):
                # it varies and not because they're all NaNs
                if verbose > 1:
                    if col == "DATA_TUNIT":
                        print "The DATA units vary across multiple SDFITS tables.  Can not continue."
                    elif col == "DATA_WIDTH":
                        print "The DATA width (number of channels) varies across multiple SDFITS tables.  Can not continue."
                    else:
                        print "An SDFITS column expected to be constant was found to not be constant.  Column: %s" % col
                return {}
            # their actual values are used later in setting fields in result
        else:
            # save all other columns to result, using lower-case name
            if type(colValues[col][0]) is numpy.string_:
                # strip off trailing spaces
                result[col.lower()] = numpy.core.defchararray.rstrip(colValues[col])
            else:
                result[col.lower()] = colValues[col]

    # CTYPE4 must be STOKES
    if colValues['CTYPE4'][0].strip() != "STOKES":
        if verbose > 1:
            print "4th axis of DATA array must be STOKES (CRVAL4)"
        return {}

    # the number of channels - already known to be constant
    result['nchan'] = colValues['DATA_WIDTH'][0]

    # the IFNUM used
    if 'IFNUM' in colValues:
        # already known to be constant
        result['ifnum'] = colValues['IFNUM'][0]
    else:
        # older data, not available, set ifnum to -1 to signal that
        result['ifnum'] = -1

    # try and work out calibration type from units and set units to 
    # actual physical units, start with TUNIT<data_column> value
    # already known to be constant
    result['units'] = colValues['DATA_TUNIT'][0].rstrip()
    # currently the pipeline sets this to the argument value of the
    # requested calibration units, e.g. Ta, Tmb, Jy.  This should
    # really be physical units with the calibration type indicated
    # separately.  Attempt to detangle that here ... anything that
    # isn't "Jy" is assumed to be "K" and the original unit
    # string is copied over to "calibtype".
    result['calibtype'] = result['units']
    if result['units'] != 'Jy':
        result['units'] = 'K'

    # decompose VELDEF into velocity definition
    # (which is still called VELDEF) and SPECSYS - appropriate for WCS
    # spectral coordinate convention.
    # This will throw an exception if there is hyphen in VELDEF.
    # A proper SDFITS file should always have that hyphen.
    veldef = colValues['VELDEF'][0].strip()
    veldef, dopframe = veldef.split('-')
    result['veldef'] = veldef
    # translate dopframe into specsys
    # I'm not 100% sure that COB from the GBT is the same as CMDBIPOL from Greisen et al
    specSysDict = {'OBS':'TOPOCENT',
                   'GEO':'GEOCENTR',
                   'BAR':'BARYCENT',
                   'HEL':'HELIOCEN',
                   'GAL':'GALACTOC',
                   'LSD':'LSRD',
                   'LSR':'LSRK',
                   'LGR':'LOCALGRP',
                   'COB':'CMBDIPOL'}
    if dopframe in specSysDict:
        result['specsys'] = specSysDict[dopframe]
    else:
        print "WARN: unrecognized frequency reference frame %s ... using OBS" % dopframe
        result['specsys'] = specSysDict['OBS']

    # convert DATE-OBS strings to JD
    result['jdobs'] = apTime.Time(colValues['DATE-OBS'],format='isot',scale='utc').jd
    
    # remaining constant column values just copied over as is, stripping out
    # trailing characters from any strings
    for col in ['CTYPE1','CTYPE2','CTYPE3',
                'RADESYS','EQUINOX',
                'SITELONG','SITELAT','SITEELEV']:
        if type(colValues[col][0]) is numpy.string_:
            # strip off trailing spaces
            result[col.lower()] = numpy.core.defchararray.rstrip(colValues[col][0])
        else:
            result[col.lower()] = colValues[col][0]

    return result



    

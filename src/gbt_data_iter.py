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

import fitsio
import numpy

class gbt_data_iter:
    # an iterator for use with GBT calibrated SDFITS files
    # Uses fitsio and allows for optional SCAN number, polarization (CRVAL4) and spectral window (IFNUM) selection.
    
    def __init__(self, sdfitsFileList, scanlist=None, crval4=None, ifnum=None, iter_row_buffer=100, verbose=4):

        self.files = sdfitsFileList

        self.verbose = verbose

        # selection criteria
        self.scanlist = scanlist
        self.crval4 = crval4
        self.ifnum = ifnum

        # the current opened file and pointer to location in self.files
        self.fio = None
        self.ifile = -1

        # indications as to where things are/what comes next
        self.hdu = 0
        self.irow = 0
        self.nrows = 0
        self.iter_row_buffer = iter_row_buffer

        # the primary header
        self.primaryHeader = None

        # the current header
        self.header = None

        # the fitsio iterator is used when no selection is in play
        self.fio_iter = None

        # otherwise we select the rows and buffer them here
        self.rows = None
        self.bufferSize = iter_row_buffer
        # points at the current location in rows and dataBuffer < bufferSize
        self.iBuff = 0
        self.dataBuffer = None

        # start at the first SINGLE DISH extension with nrows > 0
        # if none found, self.hdu will be set to len(self.fio) and self.nrows = 0
        self.openNextFile()

    def __del__(self):
        if self.fio is not None:
            self.fio.close()
            self.fio = None

    def __iter__(self):
        return self

    def nextSDFhdu(self):
        # go to the next SINGLE DISH extension with nrows > 0
        # set self.hdu to that extension number
        # set self.nrows = NAXIS2 for that extension
        # set self.header to the header for that extension
        # at end (including nore more non-empty SINGLE DISH tables)
        #   set nrows = 0, hdu = len(self.fio), header=None
        if self.hdu < len(self.fio):
            self.hdu += 1
            self.header = None
            self.nrows = 0
            self.irow = 0
            self.fio_iter = None
            self.rows = None
            self.bufferSize = self.iter_row_buffer
            self.iBuff = 0
            self.dataBuffer = None

            # otherwise those are already set to appropriate values at end

        while (self.hdu < len(self.fio)):
            if self.fio[self.hdu].get_extname() == 'SINGLE DISH':
                # HDU of appropriate type
                self.nrows = self.fio[self.hdu].get_nrows()
                if self.nrows > 0:
                    # this works
                    self.header = self.fio[self.hdu].read_header()
                    # do IFNUM selection only when requested (IFNUM >= 0) and the column exists
                    #   The IFNUM column is absent in older SDFITS files
                    do_ifnum = (self.ifnum >= 0) and ('IFNUM' in self.fio[self.hdu].get_colnames())
                    if self.scanlist is None and self.crval4 is None and (not do_ifnum):
                        # Use fitsio iterator if no selection
                        self.fio_iter = self.fio[self.hdu].__iter__()
                        break
                    else:
                        # buffer internally
                        self.fio_iter = None
                        # selection might result in an effectively empty table
                        # the columns
                        if do_ifnum:
                            cols = ['scan','crval4','ifnum']
                        else:
                            cols = ['scan','crval4']
                        scanPol = self.fio[self.hdu][cols][:]
                        rowMask = numpy.full(scanPol.shape,True,dtype=bool)

                        if self.crval4 is not None:
                            rowMask[scanPol['CRVAL4'] != self.crval4] = False
                        if self.scanlist is not None and rowMask.sum() > 0:
                            uniqueScans = numpy.unique(scanPol['SCAN'])
                            for scan in uniqueScans:
                                if scan not in self.scanlist:
                                    rowMask[scanPol['SCAN']==scan] = False
                        if do_ifnum:
                            rowMask[scanPol['IFNUM'] != self.ifnum] = False

                        if rowMask.sum() > 0:
                            # valid rows exist
                            self.rows = numpy.arange(self.nrows)[rowMask]
                            self.iBuff = 0
                            self.nrows = len(self.rows)
                            self.bufferSize = min(self.iter_row_buffer,self.nrows)
                            self.dataBuffer = self.fio[self.hdu][self.rows[0:self.bufferSize]]
                            break
                        else:
                            # effectively empty
                            self.header = None
                            self.nrows = 0

                # else this is an empty table
            # if it gets here, try the next HDU
            self.hdu += 1

    def openNextFile(self):
        # open the next file and point at the first SINGLE DISH HDU with nrows > 0
        # self.ifile will be > len(self.files) if there is no such next file
        self.ifile += 1
        if self.fio is not None:
            self.fio.close()
            self.fio = None
        self.primaryHeader = None

        while self.ifile < len(self.files):
            thisFile = self.files[self.ifile]
            if self.verbose > 3:
                print thisFile, self.ifile
            self.fio = fitsio.FITS(thisFile,iter_row_buffer=self.iter_row_buffer)
            self.primaryHeader = self.fio[0].read_header()
            self.hdu = 0
            self.nextSDFhdu()
            if self.hdu < len(self.fio):
                break
            self.ifile += 1

    def next(self):
        # avoid the StopIteration from fitsio iterator 
        # by keeping track here
        if self.irow < self.nrows:
            self.irow += 1
        else:
            # try the next SINGLE DISH extension with nrows > 0
            self.nextSDFhdu()
            if self.fio is not None and self.hdu < len(self.fio):
                self.irow += 1
            else:
                # none found here, move to next file
                self.openNextFile()
                if self.fio is not None and self.hdu < len(self.fio):
                    self.irow += 1
                else:
                    # end of the road
                    raise StopIteration

        if self.fio_iter is not None:
            return self.fio_iter.next()

        # using buffer
        # self.iBuff points to what should be returned next
        # do we need more buffer
        if self.iBuff >= self.bufferSize:
            # (self.irow-1) is the next element in self.rows to fetch
            # going past the end is already caught in the above code, so OK here 
            # there must be more to fetch
            self.bufferSize = min(self.iter_row_buffer,self.nrows-self.irow+1)
            self.dataBuffer = self.fio[self.hdu][self.rows[(self.irow-1):(self.irow-1+self.bufferSize)]]
            self.iBuff = 0
        thisBuffRow = self.iBuff
        self.iBuff += 1
        return self.dataBuffer[thisBuffRow]



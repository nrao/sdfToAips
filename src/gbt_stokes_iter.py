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

from gbt_data_iter import gbt_data_iter

class gbt_stokes_iter:
    # A data iterator for gbt data, breaks the data up into invidual gbt_data_iter
    # instances, one for each unique CRVAL4 (stokes) values.  Iterates them together.
    # the invidual unique stokes values to use must be known first.

    def __init__(self, sdfitsFileList, stokesList, scanlist=None, ifnum=None, iter_row_buffer=100, hdu=None, verbose=4):

        self.iterList = []
        self.counters = []

        for stokes in stokesList:
            self.iterList.append(gbt_data_iter(sdfitsFileList,scanlist=scanlist,crval4=stokes,ifnum=ifnum,iter_row_buffer=iter_row_buffer,hdu=hdu,verbose=verbose))
            self.counters.append(0)

        self.nIter = len(stokesList)
        # start point before the list, nextIter should also be invoked first before using this
        self.iIter = -1

        

    def __del__(self):
        for iter in self.iterList:
            if iter is not None:
                del iter
        self.iterList = []
        self.counters = []

    def __iter__(self):
        return self

    def nextIter(self):
        iCount = 0
        while iCount < self.nIter:
            iCount += 1
            self.iIter += 1
            if self.iIter >= self.nIter:
                self.iIter = 0
            if self.iterList[self.iIter] is not None:
                return
        # nothing left
        self.iIter = 0

    def next(self):
        self.nextIter()
        while self.iterList[self.iIter] is not None:
            try:
                val = self.iterList[self.iIter].next()
                self.counters[self.iIter] += 1
                return val
            except StopIteration:
                # forget about this iterator
                self.iterList[self.iIter] = None
                self.nextIter()
        # if it gets to here, nothing left to iterate over
        raise StopIteration

    def getCounters(self):
        return self.counters
        


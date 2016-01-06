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
import operator
import math
import os

def resize_groups(fitsFileName, newgcount):
    """Resize a FITS Groups primary HDU to the given new GCOUNT value"""

    # assumptions
    assert os.path.exists(fitsFileName), "%s does not exist" % fitsFileName
    assert newgcount >= 0, "GCOUNT is not >= 0 : %d" % GCOUNT

    # use fitsio to update GCOUNT and work out what the size of the
    # file should be

    fio = fitsio.FITS(fitsFileName,fitsio.READWRITE)
    hdr = fio[0].read_header()
    if hdr['GCOUNT'] == newgcount:
        # nothing to do, assume it's all good
        return

    # what should the size be
    absBitpix = abs(hdr['BITPIX'])
    pcount = hdr['PCOUNT']
    # ignore first dimension (last dimension here), always 0
    naxes = fio[0].get_dims()[:-1]
    nelms = reduce(operator.mul,naxes)
    # size in bytes
    dataSize = absBitpix * newgcount * (pcount + nelms) / 8
    # read_header_list seems to return the entire header except END as a list
    hdrSize = 80*len(fio[0].read_header_list())+1

    # 2880 bytes per block
    hdrBlocks = int(math.ceil(float(hdrSize)/2880.0))
    dataBlocks = int(math.ceil(float(dataSize)/2880.0))

    expectedSize = (hdrBlocks+dataBlocks)*2880L

    # now actually make the change
    fio[0].write_key('GCOUNT',newgcount)
    fio.close()
    del fio

    # and check the size
    actualSize = os.path.getsize(fitsFileName)

    # this should never happen
    assert actualSize >= expectedSize, "actual size is < expected size, which should never happen %d %d" % (actualSize,expectedSize)

    if actualSize > expectedSize:
        # it needs to be trimmed
        ftmp = open(fitsFileName,'a')
        ftmp.truncate(expectedSize)
        ftmp.close()
    
    return

    

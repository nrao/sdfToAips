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
# 

from scipy import constants
import math

def get_faxis(crval, crpix, cdelt, vframe, nchan, naverage, chan_start, chan_stop):
    # crval, crpix, cdelt describe the observed frequency axis in FITS channels from 1:nchan
    # vframe is the line of sight velocity (m/s) of the Doppler tracked frame relative to the observer
    # nchan is the total number channels in the spectrum
    # naverage is the number of channels to average.  An integer >= 1.
    # chan_start and chan_stop describe the channels to select, before averaging.
    #   These are FITS channels, counting from 1.
    # there are no sanity checks done on the reasonableness of any of these values

    # quantities necessary to transform to the Doppler tracked frame
    beta = vframe/constants.c
    doppler = math.sqrt((1.0+beta)/(1.0-beta))

    if naverage > 1:
        # final number of channels after selection and averaging
        # this matches what numpy.convolve does and is what is expected
        # if channel selection happens before boxcar smoothing (in which
        # case the extra channels aren't available
        # idlToSdfits appears to be greedy, using additional channels beyond
        # chan_stop so that the end channel is not always lost in the convolution
        nout = (chan_stop-chan_start+1)/naverage
        if naverage*nout == (chan_stop-chan_start+1):
            # the channel on the end is always lost, even in this case
            nout -= 1
    else:
        nout = (chan_stop-chan_start+1)
        naverage = 1

    # full bandwidth in observed frame after channel selection and averaging
    bw = abs(cdelt*naverage*nout)

    # the new reference pixels is at the center of the new frequency axis
    new_crpix = float(nout)/2. + 1.

    # channel spacing transformed to Doppler tracked frame
    new_cdelt = cdelt * doppler * naverage

    # the new crval is the value on frequency axis in the Doppler tracked frame after
    # channel selection and averaging that corresponds to new_crpix
    # which works out to this quantity
    new_crval = (crval+cdelt*(chan_start+(new_crpix-1)*naverage-crpix) + cdelt*(naverage-1.)/2.0) * doppler

    # and return things
    return (new_crval, new_cdelt, new_crpix, nout, bw)

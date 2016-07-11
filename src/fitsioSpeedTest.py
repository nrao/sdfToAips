
import sys
import fitsio
import time
import os

if len(sys.argv) != 2:
    print "Usage: python fitsioSpeedTest.py sdfits_file"
    sys.exit(1)

if not os.path.exists(sys.argv[1]):
    print "File not found : ", sys.argv[1]
    sys.exit(1)

t0 = time.time()
fio = fitsio.FITS(sys.argv[1],iter_row_buffer=1000)
colNames = fio[1].get_colnames()
scanCol = colNames.index('SCAN')
for row in fio[1]:
    # just to fish something out of it
    s = row[scanCol]

fio.close()
print "by number time (s) : ", time.time()-t0

fio == None
colNames = None
scanCol = None

t0 = time.time()
fio = fitsio.FITS(sys.argv[1],iter_row_buffer=1000)
for row in fio[1]:
    # just to fish something out of it
    s = row['SCAN']

fio.close()
print "by name time (s) : ", time.time()-t0



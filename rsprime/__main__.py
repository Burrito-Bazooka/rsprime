# encoding: UTF-8
# Copyright (c) 2014 Ryan Castellucci <code@ryanc.org>
# See LICENSE.txt for license terms

import sys
from rscoder import RSCoder

import sys
coder = RSCoder(59,58,52)
data = sys.argv[-1]
if "-d" in sys.argv:
    if coder.verify(data) is False:
        print 'WARNING: errors present, correction attempted'
    print coder.decode(data)
else:
    print coder.encode(data)

# vim: sw=4 ts=4 et ai si bg=dark

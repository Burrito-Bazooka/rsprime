# encoding: UTF-8
# Copyright (c) 2014 Ryan Castellucci <code@ryanc.org>
# See LICENSE.txt for license terms

from polynomial import Polynomial
from rscoder import RSCoder
from mapper import Mapper
from pfint import PFint

"""This module implements Reed-Solomon Encoding.
It supports arbitrary configurations for n and k, the codeword length and
message length. This can be used to adjust the error correcting power of the
code.

Warning: Because of the way I've implemented things, leading null bytes in a
message are dropped. Be careful if encoding binary data, pad the data yourself
to k bytes per block to avoid problems. Also see the nostrip option to
decode().

When called as a script, this file encodes data from standard in and outputs it
to standard out, using the standard RS code 255,223. This is suitable for
encoding text and trying it out, but don't try to encode binary data with it!

When encoding, it outputs blocks of 255 bytes, 223 of them are data (padded
with leading null bytes if necessary) and then 32 bytes of parity data.

Use the -d flag to decode data on standard in to standard out. This reads in
blocks of 255 bytes, and outputs the decoded data from them. If there are less
than 16 errors per block, your data will be recovered.
"""

# vim: sw=4 ts=4 et ai si bg=dark

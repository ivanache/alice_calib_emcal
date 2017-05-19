#! /usr/bin/env python
from os import path
import sys

import bad_channel_lhc15o_2
import bad_channel_lhc15o_3

in2not3 = []
bad_all2 = sorted(bad_channel_lhc15o_2.bad_all)
bad_all3 = sorted(bad_channel_lhc15o_3.bad_all)

for n in bad_all2:
    if not (n in bad_all3):
        print n, 'bad channel missing in run 3 map'
        if not (n in in2not3):
            in2not3.append(n)

bad_all = sorted(bad_channel_lhc15o_3.bad_all + in2not3)
print in2not3, '\n\n'
print bad_all3, '\n\n'
print bad_all


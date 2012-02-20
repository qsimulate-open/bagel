#!/usr/bin/python

import sys
import re
import os

files = os.listdir('.')

for ff in files:
    if ff[-3:] != ".cc" and ff[-2:] != ".h": continue

    filename = ff;
    fp = open(filename,"r")
    lines = fp.read().split("\n")
    
    if lines[0][0:2] != "//" or lines[1][0:2] != "//" or lines[2][0:2] != "//" or lines[3][0:2] != "//":
        print filename
        continue

    tag = "// Newint"
    if lines[1][0:len(tag)] == tag: continue
    
    year = re.search("20[0-1][0-9]", lines[2])
    yearprint = ""
    if (year):
        yearprint = year.group(0)
    else:
        yearprint = "2009"

    
    FILE = "\
//\n\
// Newint - Parallel electron correlation program.\n\
// Filename: " + filename + "\n\
// Copyright (C) " + yearprint + " Toru Shiozaki\n\
//\n\
// Author: Toru Shiozaki <shiozaki.toru@gmail.com>\n\
// Maintainer: Shiozaki group\n\
//\n\
// This file is part of the Newint package (to be renamed).\n\
//\n\
// The Newint package is free software; you can redistribute it and\/or modify\n\
// it under the terms of the GNU Library General Public License as published by\n\
// the Free Software Foundation; either version 2, or (at your option)\n\
// any later version.\n\
//\n\
// The Newint package is distributed in the hope that it will be useful,\n\
// but WITHOUT ANY WARRANTY; without even the implied warranty of\n\
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the\n\
// GNU Library General Public License for more details.\n\
//\n\
// You should have received a copy of the GNU Library General Public License\n\
// along with the Newint package; see COPYING.  If not, write to\n\
// the Free Software Foundation, 675 Mass Ave, Cambridge, MA 02139, USA.\n\
//\n\n"
    
    out = FILE + "\n".join(lines[4:])
    fp.close()

    fp = open(filename, "w")
    fp.write(out)
    fp.close()
    

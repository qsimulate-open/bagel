#!/usr/bin/python
import sys
import re


argvs = sys.argv
argc = len(argvs)

if (argc < 2):
    print "type an input"
    sys.exit();

filename = argvs[1];
fp = open(filename, "r");
lines = fp.read().split("\n")

# converts NWChem format to my format

tag = "****"
atom = "\
//\n\
// BAGEL - Parallel electron correlation program.\n\
// Filename: " + filename + ".json\n\
// Copyright (C) 2013 Toru Shiozaki\n\
//\n\
// Author: Toru Shiozaki <shiozaki@northwestern.edu>\n\
// Maintainer: Shiozaki group\n\
//\n\
// This file is part of the BAGEL package.\n\
//\n\
// The BAGEL package is free software; you can redistribute it and/or modify\n\
// it under the terms of the GNU Library General Public License as published by\n\
// the Free Software Foundation; either version 2, or (at your option)\n\
// any later version.\n\
//\n\
// The BAGEL package is distributed in the hope that it will be useful,\n\
// but WITHOUT ANY WARRANTY; without even the implied warranty of\n\
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the\n\
// GNU Library General Public License for more details.\n\
//\n\
// You should have received a copy of the GNU Library General Public License\n\
// along with the BAGEL package; see COPYING.  If not, write to\n\
// the Free Software Foundation, 675 Mass Ave, Cambridge, MA 02139, USA.\n\
//\n\
"

content=[]
nextskip = 0

atom += "\n{\n"

label = ""
numbers = []

for i in range(0,len(lines)):
    if (nextskip == 1):
        nextskip = 0
        continue
    ll = lines[i]
    # search for the tag
    if (len(ll) >= len(tag) and ll[0:len(tag)] == tag):
        # if buffer is used, flash
        if (len(numbers) > 0):
          m = len(numbers[0])
          atom += "      \"prim\" : ["
          for p in numbers:
            atom += p[0] + ", "
          atom = atom[:-2] + "],\n"

          atom += "      \"cont\" : ["
          for p in range(1,m):
            atom += "["
            for q in numbers:
              atom += q[p] + ", "
            atom = atom[:-2] + "],\n"
          atom = atom[:-2] + "]\n"
          atom += "    }\n"
          atom += "  ]"

        # then go to the next atom
        if (i != len(lines)-2):
          if (len(numbers) > 0):
            atom += ",\n"
          numbers = []
          atom += "  \"" + lines[i+1][0:2].strip(" ") + "\" : [\n" 
          nextskip = 1
    elif (ll[0:1] != " " and ll[0:1] != "*"):
        if (len(ll) == 0): continue
        # if buffer is used, flash
        if (len(numbers) > 0):
          m = len(numbers[0])
          atom += "      \"prim\" : ["
          for p in numbers:
            atom += p[0] + ", "
          atom = atom[:-2] + "],\n"

          atom += "      \"cont\" : ["
          for p in range(1,m):
            atom += "["
            for q in numbers:
              atom += q[p] + ", "
            atom = atom[:-2] + "],\n"
          atom = atom[:-2] + "]\n"
          atom += "    }, {\n"
          numbers = []
        else:
          atom += "    {\n"
        atom += "      \"angular\" : \"" + ll[0].lower() + "\",\n"
    elif (ll[0:1] == " "):
        tmp = ll[0:].split()
        numbers.append(tmp)

atom += "\n\
}"
atom = atom.replace("D-", "E-")
atom = atom.replace("D+", "E+")
atom = atom.replace("d-", "E-")
atom = atom.replace("d+", "E+")

fp2 = open(filename + ".json", "w");
fp2.write(atom)
fp2.close()
fp.close()

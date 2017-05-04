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

# converts Gaussian format to BAGEL format

tag = "****"
atom = ""

content=[]
nextskip = 0

atom += "{\n"

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

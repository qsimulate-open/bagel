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
atom = "//\n// Author : Toru Shiozaki\n// Date   : April 2012\n//\n"
content=[]
nextskip = 0

for i in range(0,len(lines)):
    if (nextskip == 1):
        nextskip = 0
        continue
    ll = lines[i]
    # search for the tag
    if (len(ll) >= len(tag) and ll[0:len(tag)] == tag):
        atom += "\nAtom:" + lines[i+1][0:2] + "\n"
        nextskip = 1
    elif (ll[0:1] != " " and ll[0:1] != "*"):
        if (len(ll) == 0): continue
        atom += ll[0].lower() + "   "
    elif (ll[0:1] == " "):
        if (atom[-1] == "\n"): atom+= "    "
        atom += "  " + ll[0:] + "\n"
atom += "\n// end of file"

fp2 = open(filename + ".basis", "w");
fp2.write(atom)
fp2.close()
fp.close()

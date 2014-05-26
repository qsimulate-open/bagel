#!/usr/bin/python

import sys

argvs = sys.argv
argc = len(argvs)

if (argc < 2):
    print "type an input"
    sys.exit();

filename = argvs[1];
fp = open(filename, "r");
lines = fp.read().split("\n")

out = []
prevshell = ""
currshell = ""
prevprims = ""

angkey = "      \"angular\" : \""
primkey = "      \"prim\" : ["

done = {}

for i in lines:
    if len(i) == 0: out.append("")
    elif i[0] == "{" or i[0] == "/": out.append(i)
    elif i[0] == "}": continue
    # atom
    elif i[0:3] == "  \"":
        prevshell = ""
        prevprims = []

        if i[3:5] != "H\"":
            out = out[:-1]
            out.append("    }")
            out.append("  ],")
            done.clear()
        out.append(i)
        out.append("    {")

    elif len(i) > len(angkey) and i[0:len(angkey)] == angkey:
        currshell = i[len(angkey)]
    elif len(i) > len(primkey) and i[0:len(primkey)] == primkey:
        prevshell = currshell
        if prevprims != i:
            prims = i[len(primkey):-2].split(", ")
            for j in prims:
                if done.has_key(currshell + j): continue
                done[currshell + j] = j
                out.append("      \"angular\" : \"" + currshell + "\",")
                out.append("      \"prim\" : [" + j + "],")
                out.append("      \"cont\" : [[ 1.0 ]]")
                out.append("    }, {")
            out = out[:-1]
            out.append("    }, {")
        prevprims == i
out = out[:-1]
out.append("    }")
out.append("  ]")
out.append("}")

print "\n".join(out)

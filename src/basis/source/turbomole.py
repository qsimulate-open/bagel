#!/usr/bin/python
import sys
import string

if len(sys.argv) < 2:
    sys.exit("specify the filename")

fp = open(sys.argv[1], "r");
lines = fp.read().split("\n")

dict = ["s", "p", "d", "f", "g", "h", "i"]
dictn = []
for i in range(0,100):
    dictn.append(str(i))

out = []
out.append("//")
out.append("// Author : Toru Shiozaki")
out.append("// Date   : Feb 2013")
out.append("//")
out.append("")

set = ""
atommark = 0
for l in lines:
    if len(l) == 0:
        continue
    elif l[0] == "*" or l[0] == ".":
        continue
    elif l.strip() == "":
        continue

    if len(l) > 6:
        s = l.strip()
        if s[3] in dict and s[0] in dictn:
            set = s[3]
            continue

    if atommark == 1:
        chars = list(l.split(" ")[0])
        str = "Atom:"+chars[0].upper()
        if len(chars) > 1: str += "".join(chars[1:])
        out.append(str)
        atommark = 0
        continue

    if len(l) > 3:
        if l == "$end":
            out.append("")
            continue
        elif l[0] == "$":
            atommark = 1
            continue

    #when does not start with number skip it
    if not (l.strip()[0] in dictn): continue

    if set != "":
        out.append(set + l)
        set = ""
    else:
        out.append(l)

out.append("")
out.append("//")
out.append("//")

for l in out:
    l = l.replace("D-", "E-")
    l = l.replace("D+", "E+")

fp.close()
fp = open(sys.argv[1]+".basis", "w")
fp.write("\n".join(out))
fp.close()

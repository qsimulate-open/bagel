#!/usr/bin/python
import sys
import string

def flash(numbers):
  m = len(numbers[0])
  atom = "      \"prim\" : ["
  for p in numbers:
    atom += p[0] + ", "
  atom = atom[:-2] + "],\n"

  atom += "      \"cont\" : ["
  for p in range(1,m):
    atom += "["
    for q in numbers:
      atom += q[p] + ", "
    atom = atom[:-2] + "],\n"
  atom = atom[:-2] + "]"
  return atom

if len(sys.argv) < 2:
    sys.exit("specify the filename")

fp = open(sys.argv[1], "r");
lines = fp.read().split("\n")

dict = ["s", "p", "d", "f", "g", "h", "i"]
dictn = []
for i in range(0,100):
    dictn.append(str(i))

out = []
out.append("{")

numbers = []

set = ""
atommark = 0
first = 1
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
        if len(numbers) > 0:
          #flash
          out.append(flash(numbers))
          out.append("    }")
          numbers = []
        if (first == 0):
          out.append("  ],")
        first = 0
        chars = list(l.split(" ")[0])
        str = "  \""+chars[0].upper()
        if len(chars) > 1: str += "".join(chars[1:])
        str += "\" : [\n"
        str += "    {" 
        out.append(str)
        atommark = 0
        continue

    if len(l) > 3:
        if l == "$end":
            continue
        elif l[0] == "$":
            atommark = 1
            continue

    #when does not start with number skip it
    if not (l.strip()[0] in dictn): continue

    if set != "":
        if len(numbers) > 0:
          #flash
          out.append(flash(numbers))
          str = "    }, {"
          numbers = []
          out.append(str)
        out.append("      \"angular\" : \"" + set + "\",")
        numbers.append(l.split())
        set = ""
    else:
        numbers.append(l.split())

#flash here too
out.append(flash(numbers))
out.append("    }")
out.append("  ]")
out.append("}")

for l in out:
    l = l.replace("D-", "E-")
    l = l.replace("D+", "E+")

fp.close()
fp = open(sys.argv[1]+".json", "w")
fp.write("\n".join(out))
fp.close()

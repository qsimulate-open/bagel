#!/bin/sh
find . -type f -exec perl -p -i -e 's/Newint/XXX/g' {} \;

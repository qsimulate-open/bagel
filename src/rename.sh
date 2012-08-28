#!/bin/sh
find . -type f -exec perl -p -i -e 's/package \(\)/package/g' {} \;

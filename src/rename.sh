#!/bin/sh
find . -type f -exec perl -p -i -e 's/ $//g' {} \;

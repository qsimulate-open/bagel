#!/bin/sh
for FILE in *.f; do cp $FILE $FILE.bak; sed -e "s/e-/d-/g" $FILE.bak > $FILE; rm $FILE.bak; done
for FILE in *.f; do cp $FILE $FILE.bak; sed -e "s/e+/d+/g" $FILE.bak > $FILE; rm $FILE.bak; done

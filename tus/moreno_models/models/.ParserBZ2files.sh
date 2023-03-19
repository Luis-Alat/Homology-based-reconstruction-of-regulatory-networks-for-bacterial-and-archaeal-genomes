#!/bin/bash

[ ! -d ".original" ] && mkdir .original

echo "Uncompressing BZ2 files"
bzip2 -dk *bz2

mv *bz2 .original

for FILE in *; do

 echo "Processing $FILE"

 grep -v "#" $FILE | 
  nl | 
  sed -r 's/^\s+//g' | 
  perl -F"\t" -lnae 'print($F[0], "\t", $F[2], "\t", $F[1])' > "${FILE}.txt"

 rm $FILE

done

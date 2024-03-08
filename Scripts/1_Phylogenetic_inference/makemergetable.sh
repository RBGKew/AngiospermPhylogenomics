#!/bin/bash

files=$@

i=1
for file in $files; do
  n=$(grep '^>' $file | wc -l)
  let j=$j+$n
  echo $(seq -s " " $i $j) "#" $file
  let i=$j+1
done

#!/bin/bash

file="housing.data"
cp "$file" "../$file.csv"
file="../$file.csv"

perl -i -pe "s/ +/ /g" "$file"
perl -i -pe "s/^ //g" "$file"
perl -i -pe "s/ /,/g" "$file"


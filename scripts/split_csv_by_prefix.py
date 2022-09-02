#!/usr/bin/env python3

import csv
import sys

file_name=sys.argv[1]
if len(sys.argv) < 3:
    output_prefix=file_name
else:
    output_prefix=sys.argv[2]

outputs={}
files=[]

with open(file_name, newline='') as src:
    for row in csv.reader(src):
        first=row[0]
        if first not in outputs:
            dst=open(f"{output_prefix}.{first}.csv","w")
            dst_csv=csv.writer(dst)
            outputs[first]=dst_csv
            files.append(dst)
        outputs[first].writerow(row[1:])

for d in files:
    d.flush()
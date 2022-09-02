#!/usr/bin/env python3

import csv
import sys
import collections

prefixes={}


max_time=-1

for src in sys.argv[1:]:
    with open(src,"rt") as xx:
        for row in csv.reader(xx):
            row=[src]+row
            
            assert len(row)==10, row
            row=[r.strip() for r in row]
            first=row[0]
            assert(row[1]=="Prop")
            time=int(row[2])
            max_time=max(max_time,time)

            times = prefixes.setdefault(first,collections.OrderedDict())
            entries=times.setdefault(time,collections.OrderedDict())
            key=tuple(row[3:7]) # Everything except last value is key
            assert key not in entries, f"{row}, {key}"
            val=row[7:]
            assert len(val)==3, (val,row)
            entries[key] = val

variants=prefixes.keys()
print(variants)

for t in range(0,max_time+1):
    all_entries=[prefixes[v].get(t,{}) for v in variants]
    keys=collections.OrderedDict()
    for s in all_entries:
        keys.update(s.items()) # The values don't matter
    for k in keys:
        print(f'{t},{k[0]:9},{k[1]:9},{k[2]:9},{k[3]:17}',end="")
        first = None
        different = False
        count = 0
        for v in variants:
            val=prefixes[v].get(t,{}).get(k)
            if val:
                count+=1
                assert(len(val)==3), val
                print(f"  ,{val[0]:10},{val[1]:10},{val[2]:10}", end="")
                if first==None:
                    first=val
                else:
                    different=val!=first
            else:
                print("  ,,,",end="")
        if(different):
            print(", DIFF", end="")
        elif count>1:
            print(", same", end="")
        else:
            print(",     ", end="")
        print()


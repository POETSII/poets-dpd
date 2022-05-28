#!/usr/bin/env python3

from PIL import Image
import sys
import glob
import pathlib

input_pattern=sys.argv[1]

min_x=10000000000
min_y=10000000000
max_x=0
max_y=0
max_w=0
max_h=0

all={}

for i in glob.glob(input_pattern):
    p=pathlib.Path(i)
    name=p.name
    x=int(name[1])
    y=int(name[0])
    print(f"({x},{y}) - {name}")
    max_x=max(max_x,x)
    max_y=max(max_y,y)
    min_x=min(min_x,x)
    min_y=min(min_y,y)
    img=Image.open(i)
    max_w=max(max_w,img.width)
    max_h=max(max_h,img.height)
    all[(x,y)]=img

res=Image.new("RGB", (max_w*(max_x-min_x+1), max_h*(max_y-min_y+1)))
for ((x,y),img) in all.items():
    res.paste(img, ((x-min_x)*max_w,(y-min_y)*max_h))

res.save("out.png")
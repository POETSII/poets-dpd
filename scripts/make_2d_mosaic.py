#!/usr/bin/env python3

from PIL import Image
import sys
import glob
import pathlib
import re

input_pattern=sys.argv[1]

min_x=10000000000
min_y=10000000000
max_x=0
max_y=0
max_w=0
max_h=0

t_values=set()

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

    mm=re.fullmatch(".+[.]([0-9]+)[.]png", i)
    if mm:
        t=int(mm.group(1))
    else:
        t=""

    t_values.add(t)
    all[(t,x,y)]=img

assert len(all)>0, "No images found"

images={}
for t in t_values:
    images[t]=Image.new("RGB", (max_w*(max_x-min_x+1), max_h*(max_y-min_y+1)))

for ((t,x,y),img) in all.items():
    images[t].paste(img, ((x-min_x)*max_w,(y-min_y)*max_h))

for (t,img) in images.items():
    sys.stderr.write(f"Writing {t}\n")
    if t=="":
        img.save("out.png")
    else:
        img.save(f"out.{t:08d}.png")
    
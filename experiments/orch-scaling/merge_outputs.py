#!/usr/bin/env python3
import glob
import re

directory="outputs"

split=re.compile(directory+"/(.+)_([0-9]+)x([0-9]+)x([0-9]+)_([0-9]+)[.]([^.]+).*")

orch_line=re.compile("STATS_fba956f3:\s+load:([^,]+), place:([^,]+), compile:([^,]+), run:([^,]+)")


print("algorithm,platform,width,height,depth,steps,volume,beads,bead_steps,status,"
    +"total_time,total_bead_steps_per_sec,total_steps_per_sec,"
    +"compile_time,compile_time_per_volume,"
    +"run_time,run_time_bead_steps_per_sec"
)

# Load orchestrator output log
for i in glob.glob(directory+"/*.orch.log"):
    m=split.match(i)
    assert m

    generator=m.group(1)
    width=int(m.group(2))
    height=int(m.group(3))
    depth=int(m.group(4))
    steps=int(m.group(5))
    platform=m.groups(6)

    status=None

    with open(i,"rt") as src:
        for line in src:
            m=orch_line.match(line)
            if m:
                load_time=float(m.group(1))
                place_time=float(m.group(2))
                compose_time=float(m.group(3))
                run_time=float(m.group(4))
                compile_time=load_time+place_time+compose_time

            if line.strip()=="Received success from devices.":
                status="Success"
            elif line.strip()=="Received explicit failure from devices.":
                status="CheckFailed"
            elif line.strip()=="Timeout or unexpected outcome while composing app.":
                status="ComposeTimeout"
            elif line.strip().startswith("Timeout while running app"):
                status="RunTimeout"

    assert run_time is not None

    if status==None:
        continue

    volume=width*height*depth
    beads=volume*3
    bead_steps=beads*steps
    total_time=load_time+place_time+compose_time+run_time
    total_bead_steps_per_sec=bead_steps/total_time
    total_steps_per_sec=steps/total_time
    compile_time_per_volume=compile_time/volume
    run_time_bead_steps_per_sec=bead_steps/run_time

    if status=="RunTimeout" or status=="ComposeTimeout":
        total_time=""
        total_bead_steps_per_sec=""
        total_steps_per_sec=""
        run_time_bead_steps_per_sec=""
        run_time=""
    if status=="ComposeTimeout":
        compile_time_per_volume=""
        compile_time=""

    print(f"basic_v5,orch,{width},{height},{depth},{steps},{volume},{beads},{bead_steps},{status},"
        +f"{total_time},{total_bead_steps_per_sec},{total_steps_per_sec},"
        +f"{compile_time},{compile_time_per_volume},"
        +f"{run_time},{run_time_bead_steps_per_sec}"
    )



# Load orchestrator output log
for i in glob.glob(directory+"/*.polite*.csv"):
    m=split.match(i)
    assert m, i

    generator=m.group(1)
    width=int(m.group(2))
    height=int(m.group(3))
    depth=int(m.group(4))
    steps=int(m.group(5))
    platform=m.group(6)
    status=None

    with open(i,"rt") as src:
        for line in src:
            if line.strip()=="Error: max SRAM partition size exceeded":
                status="FailureOutOfSRAM"
            elif line.strip()!="":
                parts=line.split(",")
                assert len(parts)==9, parts

                load_time=float(parts[5])
                run_time=float(parts[6])

                status="Success"

    if status==None:
        continue

    volume=width*height*depth
    beads=volume*3
    bead_steps=beads*steps
    total_time=load_time+run_time
    compile_time=load_time
    total_bead_steps_per_sec=bead_steps/total_time
    total_steps_per_sec=steps/total_time
    compile_time_per_volume=compile_time/volume
    run_time_bead_steps_per_sec=bead_steps/run_time

    print(f"basic_v5,{platform},{width},{height},{depth},{steps},{volume},{beads},{bead_steps},{status},"
        +f"{total_time},{total_bead_steps_per_sec},{total_steps_per_sec},"
        +f"{compile_time},{compile_time_per_volume},"
        +f"{run_time},{run_time_bead_steps_per_sec}"
    )



This repo is an adaptation of Jonny Beaumont and Shane Fleming's DPD code in
dpd-baremetal, with the main purpose being to add arbitrary degree hookean
bonds and angle bonds. 

A key interface is [`DPDEngine`](include/dpd/core/dpd_engine.hpp), which provides a base-class for exposing and
interacting with different DPD engines.

Compilation
===========

Dependencies
-----------

**TBB** : Many of the tools rely on Threaded Building Blocks (TBB) for parallelisation,
and so will not build without it.

**CMake** : This repo does not currently need cmake, but Osprey (see below) does.

**bats-core** : For internal testing of tools the [BATS-core](https://github.com/bats-core/bats-core)
framework is used. This only needs to be installed if you want to run tests of DPD-POETS.
You can install from the bats-core repo, or `apt install bats` should work on recent (18.04+?) ubuntus,
or `brew install bats-core` works with homebrew.

**Osprey-DPD** : This repo is intended to interact with [Osprey-DPD](https://github.com/Osprey-DPD/osprey-dpd.git),
but it needs some extra extensions for full inter-operability. The extensions are available at git@github.com:m8pple/osprey-dpd.git.
To install them do:
```
$ git clone https://github.com/m8pple/osprey-dpd.git
$ mkdir build
$ cd build
$ cmake ..
$ make dpd-poets
```
That should leave you with a binary `dpd-poets`, which is Osprey with the POETS DPD extensions.

**Tinsel** : this repo contains Tinsel implementations of POETS, and needs HostLink and
POLite to build them. However, it uses quite a customised version of HostLink and POLite,
so you should use the tinsel sub-module included with this repo:
```
$ git submodule update --init --recursive
```

Compile
-------

To compile, do:
```
$ make all
```
This will compile all the engines, tests, and tools.

If you are on a machine which has RISCV/Tinsel tools installed, then you
can enable them by setting `ENABLE_RISCV=1`:
```
$ ENABLE_RISCV=1 make all
```

Test
----

To test the entire system, you can do:
```
$ make test
```
or use `bats` on the individual bats files.



World state
===========

All engines work on a single shared representation of a DPD simulation
at some point in time. This contains both the static and dynamic properties
of the simulations, and is completely self-contained. The main things in
a simulation are:

- BeadTypes : each bead has a 0-based index, and a unique textual name.
   Optionally it may also be flagged as stationary (i.e. immobile).

- PolymerTypes : each bead belongs to a polymer (possibly a single-bead
   monomer), and each polymer has a type. polymer types have a 0-based
   index and a unique textual name. PolymerTypes are templates for the
   actual polymers - any change to a PolymerType affects all the polymers
   of that type. The polymer type consists of:

   - bead_types : an array of BeadType indices, describing how many beads
      there are, and what types they have.

   - bonds : an array of hookean bonds, which contains the index of the head
      and tail beads within the polymer, as well as the bond-length and
      restoring force.

   - bond_pairs : an array of angle bonds, each of which contains the index
      of the head and tail bonds in the polymer, as well as the resting angle
      and restoring force.

- polymers : an array of polymer instances. Each polymer has a 0-based polymer index,
   the index of it's polymer type, and an array of bead indices describing which
   beads are in it.

- beads : an array of bead instances describing every bead in the
   simulation. Each bead has a 0-based bead index, the index of the polymer
   that contains it, and the dynamic properties of the bead: x, v, and f.
   Each bead is in exactly 1 polymer.

   - hash: beads also have an implied "hash", which is a compact 32-bit representation
      of the bead identity, containing the bead's polymer index, offset within the
      polymer, and its bead type. Typically there are huge numbers of monomers
      (single bead polymers), and only a small number of multi-bead polymers, so 
      the hash is constructed to make best use of the available bits.

The world-state structure is represented using [include/dpcp/core/dpd_state.h](https://github.com/POETSII/poets-dpd/blob/master/include/dpd/core/dpd_state.hpp).

Engines
=======

Life-cycle
----------

An engine is an instance of the class `DPDEngine` in [`dpd/core/dpd_engine.hpp`](include/dpd/core/dpd_engine.hpp).
Each engine can work on at most once world-state at one time, and each world-state can be
attached to at most one engine at one time.

The main methods on `DPDEngine` define the life-cycle and interactions:

- `CanSupport(pState)` : Used to check if the engine can support a given WorldState. Not all engines
   support all features, so this method can determine when engines are valid for a state
   that has just been created/read.

- `Attach(pState)` : This method attaches a world-state to an engine, and lets the engine
   pre-initialise data-structures and get ready for execution. While attached, the world-state
   should not be modified but can be read.

- `Run(nSteps)` : Runs the world-state forwards for the given number of steps. While running
   the world-state is completely owned by the engine. Once Run finishes, the world-state is
   readable again.

- `Run(nIntervals,nStepsPerInterval,fIntervalCallback)` : A more advanced version of run
   which optimises doing huge numbers of runs, and takes into account the possibility of
   high-IO and startup costs (e.g. in hardware). At the end of each interval the call-back
   will execute, allowing examination of the state, but the engine may continue processing
   in parallel.

- `Attach(null)` : Detaches any world-state from the engine. After this point, the world-state
   is no longer owned, and can be read/written as normal.


Engine types
------------

There are a number of different engines, which have evolved from each other.

- NaiveDPDEngine : The simplest DPD engine that could work. This is considered the base-line for functional testing. It calculates all forces twice (once
    for each bead.

   - NaiveDPDEngineCore : A variant of NaiveDPDEngine which abstracts the maths out into a set of static functions. This is mainly used to test the static 
   maths functions so that they can be used elsewhere.

   -  NaiveDPDEngineHalfStep : This is a re-ordering of the DPD maths loop, so that
      rather than doing "position,force,velocity" as seperate steps, it is done
      as two steps "(velocity;position),force". Before the simulation loop starts the engine pre-corrects the world backwards by one position step.

      - NaiveDPDEngineHalfStepTBB : A parallel software version.

   -  NaiveDPDEngineHalfMerge : An alternative algorithm which only calculates
      the forces in one direction (i.e. more efficient).

      -  NaiveDPDEngineHalfMergeTBB : Parallel version that uses some blocking
         tricks to safely parallelise the force calculation while respecting
         dependencies.

-  BasicDPDEngine : initial software verson of an engine designed to be
    ported to hardware. Breaks the world up into unit-cells, and using
    pseudo-messages and barriers to schedule execution

File formats
============

A given WorldState object can be read or written to files in a number of
formats. Reading or writing is not guaranteed to be lossless, and
quantisation and compression is likely when saving and loading files.

Two main formats are available:

- Text : a simple textual version which is just the WorldState data-structure
   described as text, with redundancy removed.

- Binary : a much more efficient version which packs bead data down a lot
   more. It recognises that there is significant correlation between things
   like polymer and bead ids, and that absolute position is most important.
   Velocity and force need relative accuracy, and are much less important
   for capturing the state of the simulation. The main things encoded per bead are:

   - position : Absolute position encoded with 16 fractional bits and 16 integer
      bits.
   
   - velocity : Velocity vector encoded as a 4-byte floating-point format
      where there is a 5-bit exponent shared among 3 9-bit signed components
      for x, y, and z.

   - force : Force vector encoded in the same way as velocity

   Overall the beads require about 1+12+4+4=21 bytes per complete bead, versus
   about 60 for a decimal text version.

Files can be read and written using the `read_world_state` and `write_world_state`
functions in [`dpd/core/dpd_state_io.hpp`](include/dpd/core/dpd_state_io.hpp).

When reading and writing files the functions will automatically detect the
presence of the `.gz` extension on file-names, and compress/decompress with
`gzip`/`gunzip` while reading and writing. If `gzip` is not available, the functions
will fail.

Programs
========

All programs (should) have self-describing usage strings if you
run them with no arguments. So if you want to know what parameters
`bin/run_world` takes, do:
```
$ bin/run_world
```

Simulation programs
-------------------

- `bin/run_world` : Main program for DPD simulation. This reads an initial world state,
   and then steps it forwards by a specified step interval. At the end of each
   interval it extracts the state and/or dumps it as vtk/povray/png, and then carries on.

- `bin/step_world` : Alternative program for DPD simulation. A bit simpler than run_world,
   this just runs the world forwards for a certain number of time-steps, then saves the
   state.

- `bin/step_world` : Certain engines need the world state to have "settled" before they
   can execute them, particularly with respect to bond lengths. This programme will run
   the world until all bond-lengths have dropped below a certain distance, then write
   out the relaxed world.

State modification/extraction/conversion programs
-------------------------------------------------

- `bin/world_state_to_binary` : Converts any input world-state to a binary form.
   Mainly used it you want to make a textual world-state smaller.

- `bin/world_state_to_pov` : Converts a world-state to povray output form, compatible
   with the povray format used in Osprey-DPD. It will not output any bead-types with 
   the name "W".

- `bin/world_state_to_vtk` : Converts a world-state to vtk form, compatible
   with tools like ParaView. It will not output any bead-types with 
   the name "W".

- `bin/world_state_diff` : Takes two different world states and compares them, both in terms of 
  the structure of the world (number of beads, polymer topologies, ...), and then the positions
  of the beads in the world. If it finds any differences then it will return with a
  non-zero exit code.

- `bin/change_world_dt` : Reads a world and writes out an identical world with a different `dt`
   (time-step). This is quite a heavy-weight process for such a small change, but is usually
   quite fast as long as your aren't doing it a lot.

Testing
-------

- `bin/test/test_engine` : Applies a set of small functional tests to a given engine.

- `bin/test/test_engine_diff` : Applies a set of small functional tests to a given engine, while also
   diffing the output against a reference engine (naive_dpd_engine).

- `bin/test/engine_diff` : Takes a world state file and runs it through two different engines.
   It tracks the maximum and average error in things like bead position, velocity, and force.
   Error is measured over a chosen interval, and the number of consecutive intervals can be 
   chosen. This is mainly intended to look at things like round-off error due to different
   data-types and execution orders, though it can also find gross errors and logic problems.

Creating world states
---------------------

The files in `src/create_state` get compiled into binaries in `bin/create_state`, and can be
used to create worlds of different types.

Getting started
---------------

### Create a state file

We'll use `bin/create_state/create_moving_water_state` to create a
world which has:

- Dimensions : 12 x 12 x 12
- Density : 3
- Delta-t : 0.01

The output of the programme is written to stdout, so we'll save it in `water.state`:
```
$ make bin/create_state/create_state_dimers_and_water
$ bin/create_state/create_state_dimers_and_water 12 12 12 3 0.01 > dimers.state
```
If you view `dimers.state` with a text editor you'll see the complete
world description. For example, look at the first 23 lines:
```
$ head -n 30 dimers.state
```

Example:
<details>
   <summary>Click to expand</summary>

```
WorldState v0 3 2 5184 4860
T 0 0.01
Lambda 0.5
Origin 0 0 0
Box 12 12 12
Seed 1
ConservativeStrength 0 100.000000 200.000000 50.000000
ConservativeStrength 1 200.000000 100.000000 100.000000
ConservativeStrength 2 50.000000 100.000000 100.000000
DissipativeStrength 0 4.000000 4.000000 4.000000
DissipativeStrength 1 4.000000 4.000000 4.000000
DissipativeStrength 2 4.000000 4.000000 4.000000

# BeadTypes
BeadType 0 W 0.5
BeadType 1 A 0.5
BeadType 2 B 0.5

# polymerTypes
PolymerType 0 W 1 0 0
BeadTypeIndices 0

PolymerType 1 P 2 1 0
BeadTypeIndices 1 2
Bond 0 1 100.000000 0.500000


# Beads
B 0 0 0 0  4.66309 3.54688 2.54395  0 0 0  0 0 0
B 1 1 0 0  4.90527 3.83496 2.87305  0 0 0  0 0 0
```
</details>

You should see the initial size of the world, the bead types, polymers,
and then lots of beads. 

### View a state file

Even before simulating you can look at the structure by converting
state files into other representations.

#### VTK / ParaView

Assuming you have ParaView installed:
```
$ make bin/world_state_to_vtk
$ bin/world_state_to_vtk dimers.state dimers.vtk
```
You can then open dimers.vtk in ParaView, and look at the amazing dimers.
Reccommended procedure for viewing is:

1. Open the vtk file.
2. In the "pipeline browser" on the left, click the "eye" icon next to the file to make it visible.
3. In the "properties" pane below "pipeline browser", change "representation" to "Point Gaussian".
4. Use the "Reset Camera Closest" button (circle with four arrows pointing out, I guess) in "Camera Controls toolbar to re-centre. 

Note that beads with name "W" are filtered out, so we can only see the "A" and "B" beads in the polymers.

#### Povray output and rendering

You can render the non-water state to a povray file:
```
$ make bin/world_state_to_pov
$ bin/world_state_to_pov dimers.state dimers.pov
```
The resulting text-file can be rendered directly in povray, if installed:
```
$ povray -D -V dimers.pov +Odimers.png
```

### Simulating the state file

There are a number of ways of simulating. The simplest is `bin/step_engine`.
For example, we can step the simulation forwards by 1000 time-steps:
```
$ bin/step_world naive_dpd_engine dimers.state dimers-t1000.state 1000
```
The command line parameters are just the selected engine, the input
state file, the desired output state file, and the number of time-steps.

You can use the previous tools to look at `dimers-t1000.state` to see how
it has changed.

Selecting a different engine allows for faster and/or hardware accelerated
engines to be used. Most programes taking an engine ask for it as their
first input parameter - if you don't specify parameters they will print
all known engines as part of their usage message.

A decent engine is "naive_dpd_engine_half_merge_tbb", which is algorithmically
a bit better, and parallelised using TBB. Running it with this engine
should be much faster:
```
$ bin/step_world naive_dpd_engine_half_merge_tbb dimers.state dimers-t1000.state 1000
```

Engines should (!) mostly be identical except for floating-point and other
round-off and ordering errors, with random-numbers stable for any pair
of beads. If we try comparing the results of the two engines they will be
very different:
```
$ bin/world_state_diff dimers-t1000-naive.state dimers-t1000-tbb.state 
src1_file=dimers-t1000-naive.state, src2_file=dimers-t1000-tbb.state
  Bead id 0, x1=(2.22892,5.67518,0.322306), x2=(4.01138,5.1076,0.182428), dist=1.87587
```
But this difference is expected over 1000 time-steps.

Comparing over a 10 time-steps we see no errors between the two:
```
$ bin/step_world naive_dpd_engine dimers.state dimers-t10-naive.state 10
$ bin/step_world naive_dpd_engine_half_merge_tbb dimers.state dimers-t10-merge.state 10
$ bin/world_state_merge dimers-t10-naive.state dimers-t10-merge.state
```

### Running longer experiments

For long-runs with large worlds over many time-steps it is usually
much more efficient to let the simulation engine run as fast as possible,
while snapshots are stored to disk in parallel. The `bin/run_world`
program is useful for this, as it allows state and display snapshots
to be gathered at specific intervals.

The main parameters to run_world are:

1. engine-name : Which engine to use
2. src-file : the initial state file to read
3. output-base-name : a prefix which will be used as the base name for state and display snap-shot.
4. interval_count : The total number of intervals to run
5. state_interval_size : how long to wait between writing a full state file.
6. snapshot_interval_size : how long betewen writing a display snapshot.

The state and snapshot interval size must either be the same, or one must
factor the other. The total number of time-steps run will be $interval_count*min(state_interval_size,snapshot_interval_size)$.

Additional modifiers can be used to control what type of display snapshot is taken:
```
   --vtk-snapshot : Dump all non water beads into a vtk snapshot of positions.
   --povray-snapshot : Dump all non water beads into a povray snapshot of positions.
   --povray-render : Dump positions and render to a png as well (implies --povray-snapshot).
   --gzip-snapshot : Snapshots (e.g. povray and vtk) will be gzipped.
```

Let us run the dimers example for 10000 time-steps, taking a vtk snapshot every
100 time-steps, and storing the files in a directory called `dimers-out` with prefix `dimers-`.
```
$ mkdir dimers-out
$ bin/run_world naive_dpd_engine_half_merge_tbb dimers.state dimers-out/dimers- 100 100 100 --vtk-snapshot
```
This should result in the directory `dimers-out` being filled with a large number
of `dimers-.00000??00.state.gz` and `dimers-.00000??00.vtk` files. Note that
these vtk sequences can be directly recognised and loaded as an animated
sequence in paraview.
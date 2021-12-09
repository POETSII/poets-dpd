

This repo is an adaptation of Jonny Beaumont and Shane Fleming's DPD code in
dpd-baremetal, with the main purpose being to add arbitrary degree hookean
bonds and angle bonds. 

A key interface is [`DPDEngine`](include/dpd/core/dpd_engine.hpp), which provides a base-class for exposing and
interacting with different DPD engines.

Engines
-------

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

    - 

Programs
--------

There are a set of programs which are able to work with engines, of which the most important are:

- `bin/run_world` : The actual DPD simulator. This reads an initial world state,
   and then steps it forwards by a specified step interval. At the end of each
   interval it extracts the state and dumps it as vtk, then it moves to the next interval.

- `bin/test_engine_diff` : Applies a set of small functional tests to a given engine, while also
   diffing the output against a reference engine.

- `bin/world_state_diff` : Takes two different world states and compares them, both in terms of 
  the structure of the world (number of beads, polymer topologies, ...), and then the positions
  of the beads in the world.


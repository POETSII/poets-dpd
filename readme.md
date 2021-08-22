

This repo is an adaptation of Jonny Beaumont and Shane Fleming's DPD code in
dpd-baremetla, with the main purpose being to add arbitrary degree hookean
bonds and angle bonds. 

A key interface is `DPDEngine`, which provides a base-class for exposing and
interacting with different DPD engines.

There are a set of programs which are able to HERE

- `bin/run_world` : The actual DPD simulator. This reads an initial world state,
   and then steps it forwards by a specified step interval. At the end of each
   interval it extracts the state and dumps it as vtk, then it moves to the next interval.

- `bin/test_engine_diff` : Applies a set of small functional tests to a given engine.
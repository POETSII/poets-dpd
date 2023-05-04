DMPCI grid tools
================

The scripts starting with the `dmpci_grid_' prefix are a set of scripts for
managing large batches of DPD runs in a task management system such as 
[SLURM](https://slurm.schedmd.com/documentation.html). It was developed in an
ad-hoc way for the Imperial HPC system, then tidied up and refined for the
Southampton [Iridis](southampton.ac.uk/isolutions/staff/iridis.page) system.
In principle it could work for any SLURM based system, though there are
inevitably going to be assumptions baked in that need to be fixed.

Pre-requisites
--------------

The dmpci grid scripts rely on the POETS-DPD binaries mentioned in the [main readme](../readme.md), and
a few other standard tools:

-   The POETS-DPD binaries - These are expected to be located _relative to the dmpci grid scripts_. So if you
    are running 'POETS_DPD_DIR/scripts/dmpci_grid_init.sh', then it is expected that all POETS-DPD binaries will be located at
    `POETS_DPD_DIR/bin`. This should work regardless of how the script is called, as it will work out where
    it is located before using relative paths to binaries.

-   A customised version of Osprey DPD - This should be the version mentioned in the main readme.md which
    can export POETS-DPD format world states. The dmpci scripts will search for `dpd-poets` in the following places:
    1.  In the current PATH
    2.  $HOME/projects/osprey-dpd/build
    3.  $HOME/osprey-dpd/build
    4.  $HOME/POETS/osprey-dpd/build

- (Optionally) povray - This can be used during execution to render povray files into png files.
    It is assumed that the executable `povray` is on the path, both in the submission node and on
    the execuction nodes.

- (Optionally) PIL (python library) - This is used during packaging to create mosaics. Currently
    this is only needed on the submission node, or wherever the post-processing happens.

- (Optionally) ffmpeg - This will be used in packaging scripts to convert individual png files into videos.
    This is also only needed on the submission or post-processing node.

Workflow
--------

The top-level work-flow is as follows:

1. Create grid : (Manual) User prepares a zip file containing a "grid", which is one or more dmpci files.

2. Upload grid : (Manual) User copies grid to the HPC system, e.g. using `scp` or `rsync` 

3. Init grid : (Script : `dmpci_grid_init.sh`) The zip is used to initialise a working aread for the grid. This working area
   contains all the directories and meta-data needed to track the individual DPD tasks within the grid.

4. Submit grid : (Script : `dmpci_grid_submit.sh`) The initialised grid is submitted to SLURM, with one task for each grid.

5. Execution : (Automatic) The tasks within the grid are executed using SLURM as CPUs become available.
    - Check grid status : (Script : `dmpci_grid_status.sh`) During execution we can check the status of each
        task in a grid, and see whether they are queued, executing, finished, or failed.

6. Packaging: (Script) : `dmpci_grid_package.sh` : The outputs of the run are post-processed and packaged. Steps involved are:
    - Renaming : various files are renamed to match Osprey output formats
    - Mosaicing : the individual runs within the grid are combined into mosaics using python+PIL
    - Video conversion : individual pngs and mosaic pngs are convered into videos
    - Archiving : a single zip file is created containing all outputs.

7. Download : (Manual) User download the massive zip file, e.g. using `scp` or `rsync`

Walkthrough
-----------

### 1. Creating a grid

The input to the flow is a zip file containing one or more dmpci files - see the Osprey DPD
manual for a description of the file format. Typically we would expect the zip file to
be a "grid" of files, so there are a set of related dmpci files that explore some
2D parameter space. The position in the 2D grid is indicated by the 2 digit decimal
suffix, for example the six files `dmpci.{11,12,21,22,31,32}` represent a 3x3 grid.
It is not required (probably...) that files have integer grid suffixes, but certain features
like auto-mosaicing will fail/complain if they don't.

An example of an input grid is in `2B6a.zip`, which is a 4x5 grid.

### 2. Uploading a grid

Usually the grid will be prepared locally, then uploaded to the HPC system for execution.
Any method for copying files onto the HPC system can be used, including `ssh`, `rsync`,
drag-and-drop through VSCode...

### 3. Initialising the grid task

Before submitting the grid of simulations for execution we need to do some prep work, including
creating working and output directories, generating task scripts, and tagging status files.
Each initialised grid (i.e. input zip) has a single directory containing information about
all sub tasks, and this is generated using `dmpci_grid_init.sh`. These working directories
can get very large (10GB+ per grid), so you usually _don't_ want them in your home directory. Most HPC
systems have some kind of big scratch disk which is not backed up, and this is
the best place for it. In iridis this would be `/scratch/{USER}'.

Often one might want to submit/manage multiple grids at once, so a suggested method for
tracking and managing jobs this is to create a master directory for grids started on the
same date. So we could use `/scratch/{USER}/dpd_dmpci_grids/YYYY-MM-DD` as the directory for
all grids started on that day.

Let's assume we want to run `~/2B6a.zip` on `2023-05-01`. We could do the following:
```
$ mkdir -p /scratch/$USER/dpd_dmpci_grids/2023-05-01
$ cd /scratch/$USER/dpd_dmpci_grids/2023-05-01
$ cp ~/2B6a.zip .
$ dmpci_grid_init.sh 2B6a.zip
```

This should create a directory called `/scratch/$USER/dpd_dmpci_grids/2023-05-01/2B6a` containing
the initialised working areas for the grid. The internal structure of the grid is described
later in this document. We can interrogate the status of the grid using `dmpci_grid_status.sh`:
```
$ dmpci_grid_status.sh 2B6a
```
This should list the status of all the sub-tasks in the grid, and show them all as `Ready`.
As a convenience, you can use the status script with no parameters, and it will show the
status of all grids in the current directory:
```
$ dmpci_grid_status.sh
```
This is useful for monitoring the status of all grids submitted on a particular day.

### 4. Submitting the grid

Currently all the tasks in the grid are `Ready`, so now we need to submit them:
```
$ dmpci_grid_enqueue.sh 2B6a
```
This will use `sbatch` to submit them all to SLURM, and mark each task in the grid as
`Queued`. If we run:
```
$ dmpci_grid_status.sh 2B6a
```
Then the status will have changed. Usually the tasks will stay `Queued` until
the scheduler catches up and nodes are free. If you are lucky it may already have changed to
`Running` for some jobs. If your grid tasks are super-quick, or you are super-slow, they may even
have transitioned to `Finished`.

You can also look at tasks from SLURM's point of view using `squeue`, which shows
all current queued and executing tasks. If you only want to see your tasks,
using `squeue --user=$USER`.

### 5. Execution

The grid should now be executing, and it is just a waiting game. How long it takes
is very variable, and can depend on:

- The number of other users
- The specified job execution time limits in the grid
- The number of jobs you have recently submitted
- Configuration/priority changes within the HPC system

Execution can be monitored with `dmpci_grid_status.sh` and/or `squeue`.

### 5. Packaging outputs

Once all the jobs within the grid have finished, we can package the results:
```
$ dmpci_grid_package_output.sh 2B6a
```
This will kick off lots of copying and post-processing, but eventually
results in a single zip file containing all the outputs. The zip file
will be of the form `2B6a-DATETIME.zip` where DATETIME is the
date and time at the point of packaging.

For large grids this step may take a long time (5-60 minutes), particularly on the iridis submission
nodes as they limit CPUs per user. It can make sense to run this under
`screen` or `tmux` in case the connection drops or you close your laptop.


### 6. Downloading outputs

Use `scp`, `rsync`, or whatever to get the giant zip file back. Check the
file size first, and be prepared to wait if you are not on campus.

Customisation
-------------





Internal structure for a DMPCI grid
-----------------------------------

Input:
- Some/Directory/{BASE}.zip
The zip file must contain dmpci files, all at the same level

Working directory:
./${BASE}/dmpci_grid.name    # The name of the grid. Should match folder
./${BASE}/input/dmpci.AAAA, dmpci.BBBB, ...
./${BASE}/scripts/
    AAAA.sh   # bash/slurm script for each task
    BBBB.sh
./${BASE}/working/
    AAAA/
        dmpci.AAAA           # Modified version
        AAAA.begin.state.gz  # Output from Osprey after 1 step
        AAAA.init.state.gz   # Relaxed version
./${BASE}/output/
    status.txt      # Status of each the jobs (updated only when packaging)
    AAAA/
        dmi???.AAAA
        AAAA.status.txt
            - Empty/missing : job not started
            - "RUNNING,PID,START_TIME"
            - "FINISHED"
            - "FAILED"
            - "TIMEOUT"
        AAAA.log            # Logs from stderr of various programs (e.g. from slurm)
        AAAA.NNNN.rst   # solvent free osprey output
        AAAA.NNNN.vtk.gz
        AAAA.NNNN.png
        AAAA.NNNN.state

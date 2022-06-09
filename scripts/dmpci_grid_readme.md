Structure for dmpci grids:

Input:
- Some/Directory/{BASE}.zip
The zip file must contain dmpci files, all at the same level

Working:
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

#!/bin/bash
HERE="$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )"
LD_LIBRARY_PATH=${HERE} ${HERE}/run_world $@

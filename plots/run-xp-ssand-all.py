#!/usr/bin/env python3
from expTools import *

easypapOptions = {
    "-k ": ["ssandPile"],
    "-i ": [100],
    "-s ": [512],
    "-ts ": [8, 16, 32, 64]
}

# OMP Internal Control Variable
ompICV = {
    "OMP_SCHEDULE=": ["static", "static,1", "dynamic"],
    "OMP_NUM_THREADS=": [1] + list(range(2, 49, 4)) + [47, 48]
}

nbrun = 3

easypapOptions["-v "] = ["seq", "omp", "tiled", "omp_tiled", "omp_taskloop"]
easypapOptions["-of "] = ["xp/ssand-xp-all.csv"]

# Lancement des experiences
execute('./run ', ompICV, easypapOptions, nbrun, verbose=False, easyPath=".")


easypapOptions["-v "] = ["tiled", "lazy"]
easypapOptions["-of "] = ["xp/ssand-xp-TiledVsLazy.csv"]

# Lancement des experiences
execute('./run ', ompICV, easypapOptions, nbrun, verbose=False, easyPath=".")


easypapOptions["-v "] = ["omp_tiled", "omp_lazy"]
easypapOptions["-of "] = ["xp/ssand-xp-OmpTiledVsOmpLazy.csv"]

# Lancement des experiences
execute('./run ', ompICV, easypapOptions, nbrun, verbose=False, easyPath=".")
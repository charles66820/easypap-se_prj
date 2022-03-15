#!/usr/bin/env python3
from expTools import *

easypapOptions = {
    "-k ": ["ssandPile"],
    "-i ": [100],
    "-v ": ["seq", "omp", "tiled", "omp_tiled", "omp_taskloop", "lazy", "omp_lazy"],
    "-s ": [512],
    "-ts ": [32, 64],
    "-of ": ["ssand-xp-all.csv"]
}

# OMP Internal Control Variable
ompICV = {
    "OMP_SCHEDULE=": ["dynamic"],
    "OMP_NUM_THREADS=": [1] + list(range(2, 49, 4)) + [47, 48]
}

nbrun = 3
# Lancement des experiences
execute('./run ', ompICV, easypapOptions, nbrun, verbose=False, easyPath=".")

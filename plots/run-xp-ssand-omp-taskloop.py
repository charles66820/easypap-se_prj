#!/usr/bin/env python3
from expTools import *

easypapOptions = {
    "-k ": ["ssandPile"],
    "-i ": [10],
    "-v ": ["omp"],
    "-s ": [512],
    "-ts ": [8, 16, 32],
    "-of ": ["ssand-omp-taskloop.csv"]
}

# OMP Internal Control Variable
ompICV = {
    "OMP_SCHEDULE=": ["static", "static,1", "dynamic"],
    "OMP_NUM_THREADS=": [1] + list(range(2, 49, 4)) + [47, 48]
}

nbrun = 3
# Lancement des experiences
execute('./run ', ompICV, easypapOptions, nbrun, verbose=False, easyPath=".")

del easypapOptions["-ts "]
easypapOptions["-th "] = [1]
easypapOptions["-tw "] = [32, 64, 128, 256, 512]

#execute('./run ', ompICV, easypapOptions, nbrun, verbose=False, easyPath=".")

# Lancement de la version seq avec le nombre de thread impose a 1

easypapOptions = {
    "-k ": ["ssandPile"],
    "-i ": [10],
    "-v ": ["omp_taskloop"],
    "-s ": [512],
    "-ts ": [8, 16, 32],
    "-of ": ["ssand-omp-taskloop.csv"]
}
ompICV = {
    "OMP_SCHEDULE=": ["static", "static,1", "dynamic"],
    "OMP_NUM_THREADS=": [1] + list(range(2, 49, 4)) + [47, 48]
}
execute('./run ', ompICV, easypapOptions, nbrun, verbose=False, easyPath=".")

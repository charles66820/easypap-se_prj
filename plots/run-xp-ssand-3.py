#!/usr/bin/env python3
from expTools import *

easypapOptions = {
    "-k ": ["ssandPile"],
    "-i ": [100],
    "-s ": [512],
		"-wt ": ["opt", "avx"],
    "-th ": [1, 32, 64],
    "-tw ": [32, 64, 512]
}

# OMP Internal Control Variable
ompICV = {
    "OMP_NUM_THREADS=": [1]
}

nbrun = 3


easypapOptions["-v "] = ["omp_tiled", "omp_lazy"]
easypapOptions["-of "] = ["xp/ssand-xp-avx.csv"]

# Lancement des experiences
execute('./run ', ompICV, easypapOptions, nbrun, verbose=False, easyPath=".")


easypapOptions = {
    "-k ": ["ssandPile"],
    "-i ": [4000],
    "-s ": [512],
		"-th ": [2 ** i for i in range(0, 10)],
    "-tw ": [2 ** i for i in range(0, 10)],
		"-o ": ["-o "]
}
easypapOptions["-of "] = ["xp/ssand-xp-ocl-pef.csv"]

# Lancement des experiences
execute('./run ', ompICV, easypapOptions, nbrun, verbose=False, easyPath=".")


# easypapOptions["-v "] = ["omp_tiled", "omp_lazy"]
# easypapOptions["-of "] = ["xp/ssand-xp-OmpTiledVsOmpLazy.csv"]

# # Lancement des experiences
# execute('./run ', ompICV, easypapOptions, nbrun, verbose=False, easyPath=".")
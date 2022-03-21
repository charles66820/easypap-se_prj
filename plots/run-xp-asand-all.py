#!/usr/bin/env python3
from expTools import *

easypapOptions = {
    "-k ": ["asandPile"],
    "-i ": [100],
    "-s ": [512],
    "-ts ": [8, 16, 32, 64]
}

# OMP Internal Control Variable
ompICV = {
    "OMP_SCHEDULE=": ["dynamic"],
    "OMP_NUM_THREADS=": [1] + list(range(2, 49, 4)) + [47, 48]
}

nbrun = 3

easypapOptions["-v "] = ["seq", "omp_tiled"]
easypapOptions["-of "] = ["xp/asand-xp-seqVsOmpTiled.csv"]

# Lancement des experiences
execute('./run ', ompICV, easypapOptions, nbrun, verbose=False, easyPath=".")


easypapOptions["-v "] = ["tiled", "omp_tiled"]
easypapOptions["-of "] = ["xp/asand-xp-tiledVsOmpTiled.csv"]

# Lancement des experiences
execute('./run ', ompICV, easypapOptions, nbrun, verbose=False, easyPath=".")


easypapOptions["-v "] = ["omp_tiled", "omp_task"]
easypapOptions["-of "] = ["xp/asand-xp-OmpTiledVsOmpTask.csv"]

# Lancement des experiences
execute('./run ', ompICV, easypapOptions, nbrun, verbose=False, easyPath=".")


easypapOptions["-v "] = ["tiled", "lazy"]
easypapOptions["-of "] = ["xp/asand-xp-TiledVsLazy.csv"]

# Lancement des experiences
execute('./run ', ompICV, easypapOptions, nbrun, verbose=False, easyPath=".")


easypapOptions["-v "] = ["omp_tiled", "omp_lazy"]
easypapOptions["-of "] = ["xp/asand-xp-OmpTiledVsOmpLazy.csv"]

# Lancement des experiences
execute('./run ', ompICV, easypapOptions, nbrun, verbose=False, easyPath=".")
# SBS-offline
Reconstruction and analysis code for SuperBigBite (SBS) experiments

This code is presently independent of TreeSearch library

https://github.com/JeffersonLab/TreeSearch/

This code also requires a modification to the present analyzer
to handle the data sizes in bank decoding

Contains:
    MPDModule
    Decoder for MPD/APV25 used for GEMs

    replay.C
    Example driver script for analysis of GEM data

    SBSBigBite
    Roughing out the BigBite spectrometer for the
    SBS suite of experiments

    db_*
    Example databases for classes

Build prerequisites:

ROOT version 6
cmake version 3.9 or higher
Podd version 1.6 and above

How to build:

Assuming working ROOT build and environment setup, and that environment variable ANALYZER points to top-level installation directory for Podd:

```shell
mkdir build
cd build
cmake -DCMAKE_INSTALL_PREFIX=/path/to/desired/installation/directory ../SBS-offline
make install
```


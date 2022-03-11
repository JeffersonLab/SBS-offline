# SBS-offline
Reconstruction and analysis code for SuperBigBite (SBS) experiments

Build prerequisites:

ROOT version 6
cmake version 3.X or higher
Podd version 1.6 and above

**Recommended build procedure**:

Prior to building SBS-offline, it is optional, but highly recommended, that https://github.com/JeffersonLab/SBS-replay already be installed and that the environment variable SBS_REPLAY already be set to point to the top-level SBS-replay directory. This enables the convenience features for adding $SBS_REPLAY/replay and $SBS_REPLAY/scripts to ROOT's macro path, among other things.

Assuming you have a working ROOT and analyzer build and appropriately configured environment for these packages, and assuming that the environment variable ANALYZER points to the top-level installation directory for Podd, the recommended procedure for installing SBS-offline is to create a build directory parallel to the SBS-offline source directory, and an installation directory that is outside of both the build and source directories in a convenient location.

Starting from the directory containing the SBS-offline source directory, do: 

```shell
mkdir build
cd build
cmake -DCMAKE_INSTALL_PREFIX=/path/to/desired/installation/directory ../SBS-offline
make -jN install
```
where N is the number of cores available to use for multithread compilation. 

In the top-level install directory, you will find the folders: 

```shell 
bin etc include lib run_replay_here
```

You should add the line 
```
source /path/to/desired/installation/directory/bin/sbsenv.(c)sh 
```
to your login script, where the choice to source sbsenv.sh or sbsenv.csh depends on your preferred shell. This script sets the environment variable SBS to the top-level installation directory for SBS-offline and adds $SBS/lib to $LD_LIBRARY_PATH (or $DYLD_LIBRARY_PATH on Mac). 

The folder $SBS/etc contains a rootlogon.C file that automates the loading of the shared sbs library (libsbs.so) when you start the analyzer, and adds $SBS/include to ROOT's include path. 

The folder $SBS/include contains all the SBS-offline headers and ROOT dictionary files. 

The folder $SBS/lib contains the shared library that must be loaded (or linked to any external program) in order to use SBS-offline classes and methods. 

The folder $SBS/run_replay_here contains the "magic" .rootrc file that performs several functions. 

If you run ```analyzer``` from $SBS/run_replay_here (or any other folder to which $SBS/run_replay_here/.rootrc has been copied), then the following steps will happen when you start analyzer/root: 
* **$SBS/etc/rootlogon.C** will be executed, automating the loading of libsbs.so and addition of $SBS/include to ROOT's include path.
* If $SBS_REPLAY was defined prior to the cmake build of SBS-offline, it will add **$SBS_REPLAY/replay**, **$SBS_REPLAY/scripts**, and **$SBS_REPLAY/onlineGUIconfig** to ROOT's macro path. This means that if macro.C exists in any of these folders, then you can do .x macro.C or .L macro.C from any directory containing a copy of $SBS/run_replay_here/.rootrc


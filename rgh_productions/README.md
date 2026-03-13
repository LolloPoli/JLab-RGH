# CLAS12 RGH SIDIS Simulation Studies
Basic repository for RGH simulation work.

# Prerequisites
You will need the CLAS12 simulation software.
* [clasdis](https://github.com/JeffersonLab/clasdis/tree/master)
* gemc==5.10
* coatjava==11.0.1

These should already be installed on ifarm.  For a full list see the [CLAS12 Software Wiki](https://clasweb.jlab.org/wiki/index.php/How_to_run_the_full_simulation_train_single-threaded).

# Installation

Clone the directory
```bash
git clone https://github.com/mfmceneaney/rgh_simulation.git
```
and add the following to your startup script (`~/.cshrc` or `~/.bashrc` depending on your shell) on ifarm:
```bash
module load sqlite/dev # Copies a local snapshot of the CCDB so you don't overload the servers!
module load xrootd
```

You will also need to setup your enviroment.  Change the `RGH_SIM_WORK_DIR` and `RGH_SIM_VOL_DIR` directories as needed in `env.sh` or `env.csh`.  The `RGH_SIM_WORK_DIR` directory will be used to store and submit jobs scripts, and the `RGH_SIM_VOL_DIR` directory will be used to store the actual simulation outputs.  Then add the following to your startup script:
```bash
cd /path/to/rgh_simulation
source env.sh # env.csh if you are using tcsh/csh
cd -
```

# Getting started

## Generating events
Submit a job to generate events with `clasdis`:
```bash
touch jobs.txt
sbatch submit_clasdis.sh >> jobs.txt
```
This will create the directory `$RGH_SIM_VOL_DIR/lund` and store the lund files produced by clasdis there.  You can modify the arguments for clasdis to generate more or less events and change the kinematic range, target type, etc.

## Simulating events
Once that finishes, you need to simulate and reconstruct these events.  Modify `setup.sh` so that the total number of events, the number of events per file, and the prefix for naming files all match what you have in `submit_clasdis.sh`.  You will also want to make sure that `submit.sh` has the sbatch options you want.  Currently they are set to be optimal for 10k events per file.  Then, submit your simulation jobs with:
```bash
touch jobs.txt
./setup.sh >> jobs.txt
```

These jobs will create the following directory structure:
```bash
$RGH_SIM_VOL_DIR
├── cooked
├── dst
├── gemc
└── lund
```

### Detector simulation
Each job you submitted should now have its own copy of `submit.sh` and `job.sh` in `$RGH_SIM_WORK_DIR`.  The `job.sh` files will each check that the correct `clasdis` output file exists in `$RGH_SIM_VOL_DIR/lund` and then run the detector simulation on that file with `gemc` storing the output HIPO file in `$RGH_SIM_VOL_DIR/gemc`.  `gemc` will use whatever GCard you specify in `job.sh`

### Detector reconstruction
If that step finishes successfully, then the detector reconstruction will be run on the `gemc` output using `recon-util` from your coatjava installation.  The yaml file specifying which reconstruction services to run can be set in `job.sh` as well.  The output HIPO file will be stored in `$RGH_SIM_VOL_DIR/cooked`.

### Bank skimming
Finally, the cooked files will be filtered for the relevant banks using `hipo-utils` from your coatjava installation and placed in `$RGH_SIM_VOL_DIR/dst`.  The intermediate files from `gemc` and `recon-util` will be deleted after this step unless you change this in `job.sh`.

#

Contact: matthew.mceneaney@duke.edu

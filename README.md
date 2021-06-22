# localpattern
Local pattern description generation &amp; analysis

## Descripiton
Generate local pattern description around a hydrogen bond in a protein. 

## Software Requirements

- GSL (Gnu Scientific Library) >=1.16. 
- Perl >=5.14.2
- r
- GNU parallel (optional)


The above programs may already be present, depending on the platform. GNU parallel can be used to run computing in parallel, making use of multi-core platforms. 

Following perl modules are required. The easiest way is by cpanminus (on Ubuntu; `sudo apt-get install cpanminus`). Run `cpanm File::Slurp` and `cpanm Math::GSL`

- File::Slurp
- Math::GSL

Python-related requirements are as follows. Conda (https://docs.conda.io/en/latest/miniconda.html) makes it easy to switch between python versions.

*Note:* Conda also includes packages listed above (GSL, perl, r, cpanminus), but trying to install perl modules with conda-cpanm seems to cause problems. It is therefore recommended to build non-python requirements outside conda.

- Python: 2.7
- Pyyaml
- NumPy

If running prediction by alignment, alignment package (in GDT codebase, version >1.1.7) is required. Use pip to install.

## Data Requirement

The protein Hbond files may need to be edited to add some extra columns, which contain entries in the rotation matrix corresponding to the H-bond rotation. This can be done by running `addcols.py` script.

## Procedure

The following files must be present in the same directory as the code files;

- `config.py`
- `config.yaml`
- `Option.py`
- `Pattern.py`
- `Protein.py`
- `Rasmus2.r`
- `improve_modebox.pl`
- `StartClustSubsetVers4.txt`
- `eval.pl`

The analysis is done in the following order:

1. `gen_patrot.py`

   This generate the hbond pattern & rotation information for the training data. Each row in the output files `{proteins_id}_{step #}.patrot` contain pattern, rotation. Parameters for pattern generation (whether to use tertiary bonds, directory containing _opts files, etc.) are stored in `config.yaml`.

   Example: `python gen_patrot.py sampledata/hbonds/1ONWA00.txt 336-345 sampledata/pat`

   Average throughput is ~3.5 files per hour per cpu core.

2. `filter.py`

   Remove patterns from the previous step which have fewer than n corresponding rotation values, where n is the minimum threshold specified in `config.yaml`. Output files are called `step##_rotations.pkl`. 

   Example: `python filter.py sampledata/pat sampledata/rotations 336`

3. `clstr2.py`
   Run clustering on the results from the previous step. Clustering is done by R script (`Rasmus2.r`), and post-processing by perl script (`improve_modebox.pl`). Output files are called `step##_summary.pkl`.

   Example: `python clstr2.py sampledata/rotations/step336_rotations.pkl`

4. `eval2.py`

   Evaluate clustering results in `step##_summary.pkl` files. The actual work is done by `eval.pl` perl script. Results are `step##_assess` files.

   Example: `python eval2.py sampledata/summaries/step336_summary.pkl`

5. `predict.py`
   Prediction for a single protein using the information from `_assess` files. The results are saved in `[protein id]_pred.txt` files. 

   Example: `python predict.py sampledata/hbonds/3F4SA00.txt sampledata/assess sampledata/predictions/3F4SA00_pred.txt`

6. `alignlocpats.py`

   Prediction using alignment. Found in `scripts` under `alignment` package install directory.

To see the required (and optional) arguments for python scripts, run `python [script].py -h`. Note all scripts above can be run in parallel, either by proteins or by step numbers. GNU parallel is very good at running these "embarrasingly parallel" tasks. See https://www.gnu.org/software/parallel/parallel_tutorial.html for more info in the use of parallel.
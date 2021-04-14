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

   This generate the hbond pattern & rotation information for the training data. Each row in the output files `{proteins_id}_{step #}.patrot` contain pattern, rotation.

   Average throughput is ~3.5 files per hour per cpu core.

2. `filter.py`

   Remove patterns from the previous step which have fewer than n corresponding rotation values, where n is the minimum threshold specified in `config.yaml`. Output files are called `step##_rotations.pkl`. 

3. `clstr2.py`
   Run clustering on the results from the previous step. Clustering is done by R script (`Rasmus2.r`), and post-processing by perl script (`improve_modebox.pl`). Output files are called `step##_summary.pkl`.

4. `eval2.py`
   
   Evaluate clustering results in `step##_summary.pkl` files. The actual work is done by `eval.pl` perl script. Results are `step##_assess` files.
   
5. `predict.py`
   Prediction for a single protein using the information from `_assess` files. The results are saved in `[protein id]_pred.txt` files. 

To see the required (and optional) arguments for python scripts, run `python [script].py -h`. Note all scripts above can be run in parallel, either by proteins or by step numbers. GNU parallel is very good at running these "embarrasingly parallel" tasks. See https://www.gnu.org/software/parallel/parallel_tutorial.html for more info in the use of parallel.
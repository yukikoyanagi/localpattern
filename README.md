# localpattern
Local pattern description generation &amp; analysis

## Descripiton
Generate local pattern description around a hydrogen bond in a protein. 

## Procedure
The analysis is done in the following order:
1. patrot_prll.sh
   This generate rotations.pkl files containing {pattern: [list of rotaions]} dict in /work/austmathjea/cdp/step* directories. 
   Directory containing input protein files, tertiary files, _opts files, etc are specified in config.yaml file. Pattern is the
   string representation of pattern object. This step takes approx. 3 hours on 32 slim nodes.
2. filter.sh
   Filter the pattern-rotation data generated in step 1, so that only the patterns with minimum number of rotations or more are 
   used in clustering analysis. The minimum number is specified in the config.yaml file as min_for_cluster. Generates rotations.pkl
   files in /work/austmathjea/cdp/step*/n{min_for_cluster} directories. This step is fast; takes <5 min on 1 slim node.
3. clstr.sh
   Run clustering analysis for each step. Generates summary.pkl file in /work/austmathjea/cdp/step*/n{min_for_cluster} directories.
   One can control the steps to run the analysis by specifying -J option in sbatch; {job name}.{start step}.{end step}
   This step takes approx. 3 hours on 32 slim nodes.
4. eval.sh
   Evaluate clustering results to produce assessment files. Generates step*_assess file in /work/austmathjea/cdp/step*/n* 
   directories. The cutoff value can be specified using -J option in sbatch; {job name}-{cutoff}
   Takes < 5 min on 1 slim node.
5. pred_prll.sh
   This is the script that predicts rotations. Output prediction file is specified in the file as an argument to .py file. 
   Takes < 10 min on 16 slim nodes.

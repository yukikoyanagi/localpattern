#!/usr/bin/env bash
#
# File: prepare_eval.sh 
#
# Description: Prepare ucloud instance for pattern-rotation generation
#
#

echo Update apt cache
sudo apt-get update

echo Install GSL, cpanm
sudo apt-get install -y libgsl-dev
sudo apt-get install -y cpanminus

echo Install perl modules
cpanm --sudo File::Slurp
cpanm --sudo --force Shell::Config::Generate
cpanm --sudo Math::GSL

#echo Install R
#sudo apt-get install -y r-base

echo Install GNU Parallel
sudo apt-get install -y parallel

echo Enable conda
echo eval "$(/work/miniconda3/bin/conda shell.bash hook)" >> /home/ucloud/.bashrc
bash -i


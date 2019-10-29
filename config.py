#!/usr/bin/env python

import yaml
import os

with open("config.yaml", "r") as f:
    cfg = yaml.full_load(f)
if cfg['abacus']:
    workdir = os.environ['WORK']
    cfg['protdir'] = os.path.join(workdir, cfg['protdir'])
    cfg['optsdir'] = os.path.join(workdir, cfg['optsdir'])
    cfg['outdir'] = os.path.join(workdir, cfg['outdir'])
    cfg['tertdir'] = os.path.join(workdir, cfg['tertdir'])

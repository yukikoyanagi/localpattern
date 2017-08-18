#!/usr/bin/env python
#
# File: clstr.py
#
# Time-stamp: <2017-02-24 10:26:49 au447708>
#
# Description: Run cluster analysis on Abacus. Step no's are obtained
#  from $SLURM_JOB_NAME, which should be in the form {name}.start.end
#  or {name}-group#-#ofgroups. To overwrite existing summary.pkl
#  files, start $SLURM_JOB_NAME with 'rp'. Default is to skip the
#  step if summary.pkl file exists.
#  Change minn parameter as required.
#
# Author: Yuki Koyanagi
# History:
#  2017-01-16: Created
#  2017-01-18: Submit job for sets of patterns to avoid starting
#   too many threads.
#  2017-01-20: Use subprocess module instead of os.system() to run
#   improve_modebox script.
#  2017-02-02: Set default to skip the step if summary.pkl file exits.
#   Output summary.pkl file are written as the steps are finished.
#  2017-02-14: Do not run clustering when only one observation is
#   available. Instead save them in a separate file, singles.pkl and
#   append it later during evaluation.
#  2017-08-15: Copied from cdp to locpat repository.
#
import os
import cPickle as pickle
import subprocess
import pp
from itertools import islice, tee

from config import cfg

minn = cfg['min_for_cluster']
maxstep = cfg['max_step']

# Set up servers
ppservers = open('/tmp/nodelist').read().strip().split()
ppservers = tuple(pp + ':2048' for pp in ppservers)
job_server = pp.Server(0, ppservers=ppservers)


def filtersingles(data, outdir):
    '''
    Removes items in data which has len==1 and writes them to
    singles.pkl file in outdir directory. The remaining data is
    returned as is.
    '''
    singles = {}
    for key in data.keys():
        if len(data[key]) == 1:
            singles[key] = data.pop(key)[0]
    outf = os.path.join(outdir, 'singles.pkl')
    with open(outf, 'wb') as f:
        pickle.dump(singles, f)

    return data


def chunk(data, n=24):
    '''Splits input data (dict) into n equal chunks.'''
    it = iter(data)
    for i, iterator in enumerate(tee(it, n)):
        yield {k: data[k] for k in islice(iterator, i, None, n)}


def runclstrs(data, step):
    d = {}
    for pattern in data:
        # If R script has error runclstr returns None
        summ = runclstr(pattern, data[pattern], step)
        if summ:
            d[pattern] = summ
    return step, d


def runclstr(pattern, data, step):
    import os
    import subprocess

    # Files and dirs
    rscript = os.path.join(os.path.expandvars('$SLURM_SUBMIT_DIR'),
                           'Rasmus2.r')
    imprv = os.path.join(os.path.expandvars('$SLURM_SUBMIT_DIR'),
                         'improve_modebox.pl')
    stepd = os.path.join(os.path.expandvars('$LOCALSCRATCH'),
                         'step{}'.format(step))
    rotf = os.path.join(stepd, pattern)

    if not os.path.exists(stepd):
        try:
            os.mkdir(stepd)
        except OSError:
            # stepd now exists -prob. created by another worker proc.
            pass

    with open(rotf, 'w') as f:
        f.write('\n'.join(['\t'.join(map(str, d)) for d in data]))

    try:
        subprocess.check_output(['Rscript', rscript, rotf],
                                stderr=subprocess.STDOUT)
    except subprocess.CalledProcessError as cpe:
        print 'Error while processing step{}/{}'.format(step, pattern)
        print cpe.output
        return None

    # sumf = rotf + '_Summary.txt'

    cmd = ' '.join([imprv, rotf])
    os.system(cmd)

    sumf = rotf + '_Summary2.txt'

    with open(sumf) as f:
        res = f.read()

    return res


# Parse jobname to get start and end step no's
jobname = os.path.expandvars('$SLURM_JOB_NAME')
if '.' in jobname:
    _, startstep, endstep = jobname.split('.')
    steps = range(int(startstep), int(endstep)+1)
elif '-' in jobname:
    _, grp, total = jobname.split('-')
    steps = range(int(grp), maxstep+1, int(total))
else:
    raise SyntaxError('SLURM_JOB_NAME not set correctly.')

overwrite = jobname.startswith('rp')

cdpdir = os.path.join(os.path.expandvars('$WORK'), 'cdp')

# Submit jobs to server
clsts = []
submitted = []

for step in steps:
    # Load the rotation data
    workdir = os.path.join(cdpdir,
                           'step{}'.format(step),
                           'n{}'.format(minn))
    if (
            (not overwrite) and
            os.path.exists(os.path.join(workdir, 'summary.pkl'))):
        print 'Skipping step{}. Summary.pkl file exists.'.format(step)
        continue

    rotfile = os.path.join(workdir, 'rotations.pkl')
    with open(rotfile, 'rb') as r:
        rotations = pickle.load(r)

    rots = filtersingles(rotations, workdir)
    # Submit jobs for subset of data to avoid too many threads
    for subset in chunk(rots):
        clstr = job_server.submit(runclstrs,
                                  args=(subset, step),
                                  depfuncs=(runclstr,),
                                  group='step{}'.format(step))
        clsts.append(clstr)
    submitted.append(step)

# Get the results in dict; {pattern: content}
for step in submitted:
    job_server.wait('step{}'.format(step))
    summary = {}
    for clst in [task for task in clsts
                 if task.group == 'step{}'.format(step)]:
        step, d = clst()
        summary.update(d)
    outf = os.path.join(cdpdir,
                        'step{}'.format(step),
                        'n{}'.format(minn),
                        'summary.pkl')
    with open(outf, 'wb') as o:
        pickle.dump(summary, o)

# job_server.destroy()

#1/usr/bin/env python2

MAXSTEP = 791
RSCRIPT = 'Rasmus2.r'
IMPSCRIPT = 'improve_modebox.pl'

import argparse, re, os.path, uuid, subprocess

import Protein, Option, conv

def loadassess(ad):
    # Read assess file into dict
    adic = {}
    expr = (r'(?P<pat>^[^#]\S+)'
            r'\t\S+\t(?P<rot>\S+)'
            r'\t(?P<score>\d+.\d+)\s+')
    p = re.compile(expr, re.M)

    ad = ad.rstrip('/')
    for step in range(MAXSTEP+1):
        af = '{}/step{}_assess'.format(ad, step)
        with open(af) as fh:
            for line in fh:
                m = p.search(line)
                if m:
                    rot = map(float, re.split(r'[;,]', m.group('rot')))
                    adic[(step, m.group('pat'))] = (rot,
                                                    float(m.group('score')))
    return adic

def main(pf, ad, optd, od):
    # Files and dirs
    tmpd = os.path.expandvars('$LOCALSCRATCH')
    if tmpd == '$LOCALSCRATCH':
        tmpd = '/tmp'
    
    # Load protein
    pid = os.path.splitext(os.path.basename(pf))[0]
    prot = Protein.Protein(pid)
    prot.fromfiles(pf, None)

    # Load assess files into dict {(step, pat): (rot, #obs, score)}
    adic = loadassess(ad)

    optd = optd.rstrip('/')
    opts = []
    for step in range(MAXSTEP+1):
        of = '{}/step{}_opts'.format(optd, step)
        opt = Option.Option(of)
        opts.append(opt)

    res = []
    for bond in prot.hbonds:
        preds = []
        for step in range(MAXSTEP+1):
            pat = str(prot.findpattern2(bond, opts[step]))
            try:
                rot, score = adic[(step, pat)]
            except KeyError:
                continue
            # We don't want ambiguous clusters
            if score > 10:
                # Add bias towards larger window size
                for i in range(opts[step].window//3 + 1):
                    preds.append((rot[-1], rot[0], rot[1], rot[2]))
        #print('Collected rotations for {}/{}'.format(pid, bond.lineno))

        rotf = os.path.join(tmpd, uuid.uuid4().hex)
        with open(rotf, 'w') as fh:
            fh.write('\n'.join(
                ['\t'.join(
                    map(str, d)) for d in preds]))

        try:
            subprocess.check_output(['Rscript', RSCRIPT, rotf],
                                    stderr=subprocess.STDOUT)
        except subprocess.CalledProcessError as cpe:
            print('Error while clustering {}/{}'.format(pid, bond.lineno))
            print(cpe.output)
            return None

        try:
            subprocess.check_output(['perl', IMPSCRIPT, rotf],
                                    stderr=subprocess.STDOUT)
        except subprocess.CalledProcessError as cpe:
            print('Error while processing {}/{}'.format(pid, bond.lineno))
            print(cpe.output)
            return None

        # Now process the output from clustering
        sumf = rotf + '_Summary2.txt'
        with open(sumf) as fh:
            line = fh.readline()
        cols = line.split()
        try:
            pred = conv.box2ev(map(int, cols[1:4]))
        except ValueError:
            print('Error while processing {}'.format(rotf))
        p_rot = conv.ev2mat(*pred)
        brot = bond.rotation
        t_rot = conv.ev2mat(brot[-1], brot[0:3])
        diff = conv.so3dist(p_rot, t_rot)
        ol = '{}\t{}\t{}'.format(pid, bond.lineno, diff)
        if not od:
            print(ol)
        else:
            res.append(ol)
    if od:
        of = os.path.join(od, 'pred.{}.txt'.format(pid))
        with open(of, 'w') as fh:
            fh.write('\n'.join(res))


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('pdir',
                        help='Directory containing protein files, '
                        'or a protein H-bond file. If directory '
                        'assume we are running under slurm.')
    parser.add_argument('assessd',
                        help='Directory containing *_assess files.')
    parser.add_argument('optd',
                        help='Directory containing _opts files.')
    parser.add_argument('-o', '--outd',
                        help='Output directory.')
    args = parser.parse_args()
    main(args.pdir, args.assessd, args.optd, args.outd)

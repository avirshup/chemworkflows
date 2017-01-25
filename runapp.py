import os
import pyccc
from pyccc.workflow.runner import SerialCCCRunner as Runner
import pickle


def runapp(modname, pdbid):
    engine = pyccc.Docker()

    if modname == 'MMminimize':
        from apps.MMminimize import app
    #elif modname == 'vde':
    #    from apps.vde import vde as app
    else:
        raise ValueError('No app named %s')

    idir = 0
    while os.path.exists('outputs.%d' % idir):
        idir += 1

    outdir = 'outputs.%d' % idir

    print 'Running workflow "%s" with input "%s." Outputs will be written to "./%s".' % (
        modname, pdbid, outdir)

    runner = Runner(app,
                    engine=engine,
                    molecule_json={'pdb':pdbid})
    runner.run()

    # writing outputs ...
    os.mkdir(outdir)
    prmtop = runner.outputs['prmtop']
    prmtop.put(os.path.join(outdir, 'prmtop'))

    inpcrd = runner.outputs['inpcrd']
    inpcrd.put(os.path.join(outdir, 'inpcrd'))

    pdbfile = runner.outputs['finalpdb']
    with open(os.path.join(outdir, 'final_structure.pdb'), 'w') as outfile:
        print >> outfile, pdbfile


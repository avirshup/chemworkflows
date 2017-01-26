import os
import pyccc
from pyccc.workflow.runner import SerialCCCRunner, SerialRuntimeRunner


def get_engine():
    server = os.environ.get('CCC', None)
    if not server:
        print 'WARNING: No CCC server set in environment variable "$CCC", using default.'
        server = 'ccc.bionano.autodesk.com:9000'

    print 'Using CCC server: %s' % server
    return pyccc.engines.CloudComputeCannon(server)


def runapp(args):
    runner = SerialCCCRunner
    if args.localdocker:
        assert not args.here
        engine = pyccc.Docker()
    elif args.here:
        runner = SerialRuntimeRunner
        engine = None
    else:
        engine = get_engine()

    if args.outputdir is None:
        outdir = _get_output_dir(args.appname)
    else:
        outdir = args.outdir

    print 'Running workflow "%s" with input "%s."\nOutputs will be written to "%s".' % (
        args.appname, args.inputfile, outdir)

    APPNAMES[args.appname](engine, runner, outdir, args.inputfile)
    print 'DONE. Output is in %s' % outdir


def _get_output_dir(modname):
    idir = 0
    while os.path.exists('%s.out.%d' % (modname,idir)):
        idir += 1

    outdir = 'outputs.%d' % idir
    return outdir


def run_vde(engine, Runner, outdir, inputfile):
    import json, yaml
    from .apps.vde import vde as vdeapp

    ext = inputfile.split('.')[-1]
    if ext in ('js', 'json', 'yml', 'yaml'):
        with open(inputfile, 'r') as infile:
            inputjson = yaml.load(infile)
    else:
        with open(inputfile, 'r') as infile:
            inputjson = {'filename': inputfile,
                         'content': infile.read()}

    runner = Runner(vdeapp,
                    engine=engine,
                    molecule_json=inputjson)

    runner.run()

    if not os.path.exists(outdir):
        os.mkdir(outdir)
    with open(os.path.join(outdir, 'result.pdb'), 'w') as outfile:
        print >> outfile, runner.outputs['pdbstring']

    with open(os.path.join(outdir, 'result.json'), 'w') as outfile:
        json.dump(runner.outputs['results'], outfile)


def run_ligand_finder(engine, Runner, outdir, inputfile):
    import json
    from .apps.MMminimize import prep_active_site as ligand_app

    with open(inputfile, 'r') as infile:
        content = infile.read()
    runner = Runner(ligand_app,
                    engine=engine,
                    molecule_json={'filename': inputfile,
                                   'content': content})

    runner.run()
    if not os.path.exists(outdir):
        os.mkdir(outdir)

    with open(os.path.join(outdir, 'prep.pdb'), 'w') as outfile:
        print >> outfile, runner.outputs['pdbstring']

    result = {'ligands': runner.outputs['ligand_options']}
    result['success'] = runner.outputs['validates']
    result['errors'] = runner.outputs['validation_errors']

    with open(os.path.join(outdir, 'prep.json'), 'w') as outfile:
        json.dump(result, outfile)



def run_minimize(engine, Runner, outdir, pdbid):
    from .apps.MMminimize import full_minimization as minimizerapp
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


APPNAMES = {'MMminimize': run_minimize,
            'prep_ligands': run_ligand_finder,
            'vde': run_vde}

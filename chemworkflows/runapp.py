import os
import json
import pickle

import yaml
import pyccc
from pyccc.workflow.runner import SerialCCCRunner, SerialRuntimeRunner

from .apps import vde, MMminimize

APPNAMES = {'MMminimize': MMminimize.mm_minimization,
            'vde': vde.vde}


def main(args):
    engine, runner = get_execution_env(args)
    inputjson = process_input_file(args.inputfile)
    outdir = get_output_dir(args)
    app = APPNAMES[args.appname]

    runner = runner(app,
                    engine=engine,
                    molecule_json=inputjson)

    if args.preprocess:
        run_preprocessing(runner, outdir)
    else:  # run the whole thing
        run_workflow(runner, outdir)

    print 'DONE. Output directory:'
    print "    ", os.path.abspath(outdir)


def run_workflow(runner, outdir):
    runner.run()

    os.mkdir(outdir)
    for name, value in runner.outputs:
        if isinstance(value, basestring):
            with open(os.path.join(outdir, name), 'w') as outfile:
                print >> outfile, value
        else:
            value.put(os.path.join(outdir, name))


def run_preprocessing(runner, outdir):
    t = runner.preprocess()

    os.mkdir(outdir)
    with open(os.path.join(outdir, 'prep.pdb'), 'w') as outfile:
        print >> outfile, t
    resultjson = {}
    for field in t.outputs:
        if field == 'pdbstring':
            continue
        else:
            resultjson[field] = t.getoutput(field)
    with open(os.path.join(outdir, 'prep.json'), 'w') as outfile:
        json.dump(resultjson, outfile)
    with open(os.path.join(outdir, 'workflow_state.P'), 'w') as outfile:
        pickle.dump(runner, outfile)


def get_output_dir(args):
    if args.outputdir is None:
        idir = 0
        while os.path.exists('%s.out.%d'%(args.appname, idir)):
            idir += 1
        outdir = '%s.out.%d'%(args.appname, idir)
        return outdir

    else:
        outdir = args.outdir

    print 'Running workflow "%s" with input "%s."\nOutputs will be written to "%s".'%(
        args.appname, args.inputfile, outdir)
    return outdir


def get_execution_env(args):
    runner = SerialCCCRunner
    if args.localdocker:
        assert not args.here
        engine = pyccc.Docker()
    elif args.here:
        runner = SerialRuntimeRunner
        engine = None
    else:
        engine = get_engine()
    return engine, runner


def get_engine():
    server = os.environ.get('CCC', None)
    if not server:
        print 'WARNING: No CCC server set in environment variable "$CCC", using default.'
        server = 'ccc.bionano.autodesk.com:9000'

    print 'Using CCC server: %s' % server
    return pyccc.engines.CloudComputeCannon(server)


def process_input_file(inputfile):
    try:
        jsraw = _get_json(inputfile)
    except ValueError:
        pass
    else:
        inputjson = json.loads(jsraw)
        return inputjson

    ext = inputfile.split('.')[-1]
    if ext in ('js', 'json', 'yml', 'yaml'):
        with open(inputfile, 'r') as infile:
            inputjson = yaml.load(infile)
    else:
        with open(inputfile, 'r') as infile:
            inputjson = {'filename': inputfile,
                         'content': infile.read()}
    return inputjson


def _get_json(s):
    s = s.strip()
    while s and s[0] == s[-1] and s[0] in ('"', "'"):
        s = s[1:-1].strip()

    if not s or s[0] != '{' or s[-1] != '}':
        raise ValueError()
    else:
        return s


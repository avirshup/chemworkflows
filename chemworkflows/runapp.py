import os
import sys
import json

import dill

import yaml
import pyccc
from pyccc.workflow.runner import SerialCCCRunner, SerialRuntimeRunner

from .apps import vde, MMminimize

APPNAMES = {'MMminimize': MMminimize.mm_minimization,
            'vde': vde.vde}


def main(args):
    outdir = get_output_dir(args)
    os.mkdir(outdir)

    if args.restart:  # restarts branch here!!!
        restart_workflow(args, outdir)
        sys.exit(0)

    engine, RunnerClass = get_execution_env(args)
    inputjson = process_input_file(args.inputfile)
    app = APPNAMES[args.appname]

    runner = RunnerClass(app,
                         engine=engine,
                         molecule_json=inputjson)

    if args.preprocess:
        run_preprocessing(runner, outdir)
    else:  # run the whole thing
        run_workflow(runner, outdir)


def run_workflow(runner, outdir):
    runner.run()

    print 'DONE. Output directory:'
    print "    ", os.path.abspath(outdir)

    with open(os.path.join(outdir,'workflow_state.dill'), 'w') as outfile:
        dill.dump(runner, outfile)
    for name, value in runner.outputs.iteritems():
        filebase = os.path.join(outdir, name)
        if isinstance(value, basestring):
            with open(filebase, 'w') as outfile:
                print >> outfile, value
        elif isinstance(value, dict):
            with open(filebase + '.json', 'w') as outfile:
                json.dump(value, outfile)
        elif hasattr(value, 'put'):
            value.put(filebase)
        elif hasattr(value, 'read'):
            with open(filebase, 'w') as outfile:
                print >> outfile, value.read
        else:
            with open(filebase + '.dill', 'w') as outfile:
                dill.dump(value, outfile)



def restart_workflow(args, outdir):
    with open(args.inputfile, 'r') as infile:
        runner = dill.load(infile)

    engine, RunnerClass = get_execution_env(args)
    assert RunnerClass is runner.__class__

    print 'Restarting workflow %s' % runner.workflow.name

    run_workflow(runner, outdir)


def run_preprocessing(runner, outdir):
    t = runner.preprocess()

    print 'FINISHED preprocessing. Output directory:'
    print "    ", os.path.abspath(outdir)

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
    with open(os.path.join(outdir, 'workflow_state.dill'), 'w') as outfile:
        dill.dump(runner, outfile)



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


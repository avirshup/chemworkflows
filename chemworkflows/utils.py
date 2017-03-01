from __future__ import print_function
import sys
import pyccc
import pickle


def run_mdt(fn, *args, **kwargs):
    """ Runs a python command in an MDT container
    """
    pflush('Running python:%s ...' % fn.__name__, end='')
    wait = kwargs.pop('wait', True)
    inputs = kwargs.pop('inputs', {})

    pycall = pyccc.PythonCall(fn, *args, **kwargs)
    job = pyccc.PythonJob(config.ENGINE, config.MDTIMAGE,
                          pycall,
                          inputs=inputs)

    with open('job%s.pkl' % job.jobid, 'w') as jobfile:
        pickle.dump(job, jobfile)
    pflush('id:%s ...' % job.jobid, end='')

    return finish_job(job, wait)


def finish_job(job, wait):
    if wait:
        job.wait()
        pflush('done', end='\n')
    else:
        pflush('submitted', end='\n')
        return job
    try:
        s = job.stderr
    except Exception:
        print('Failed to get stderr...')
    else:
        print('STDERR:')
        if s.strip():
            print(s)
        print("STDOUT:")
        if job.stdout.strip():
            print(job.stdout)
    return job


def get_asset(filename):
    import os
    thispath = os.path.dirname(os.path.abspath(__file__))
    filepath = os.path.join(thispath, 'assets', filename)
    return open(filepath)


def submit_job(job, image=None, wait=True):
    """ Runs a job object. This is usually something created by MDT
    that drives an external tool
    """
    pflush('Running job: %s ...' % job.name, end='')
    job.engine = config.ENGINE
    if image is not None: job.image = image

    job.submit()

    with open('job%s.pkl' % job.jobid, 'w') as jobfile:
        pickle.dump(job, jobfile)

    pflush('id:%s ...' % job.jobid, end='')

    return finish_job(job, wait)



def pickle_trajectory(job):
    r = job.result
    r.write('traj.pkl')
    r.mol.write('out.pkl')


def pickle_mol(job):
    r = job.result
    r.mol.write('out.pkl')


def finish_job_callback():
    import moldesign as mdt

    job_description = mdt.read('job.description.pkl')
    finished_job = mdt.read('job.run.pkl')


def pflush(s, **kwargs):
    print(s, **kwargs)
    sys.stdout.flush()



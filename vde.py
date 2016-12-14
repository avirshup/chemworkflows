import json
import sys

import yaml

import common
from utils import run_mdt, submit_job
import config

# what would make this easier? Make ALL nodes use compatible i/o formats
# Then we don't need to translate before/after sending to NWChem
# Otherwise we're stuck with these prepare input/run job/parse output chains

# Other options (not exclusive):
# 1. Somehow turn vanilla python scripts into workflow with clever magic (probably not possible)
# 2. Some sort of workflow-compatible tail-call-recursion for the translate-run-translate patterns
# 3. wrappers to just make this all less verbose (done)
# 4. Better integration with MDT's job launching system (probably necessary if we want this to be
#     usable ... not sure what it looks like

# 2/4: When MDT generates a job object, it should be ready to run. Actually, tail-call recursion
# works well here - jobs return job objects, which can be launched. It will be a common task to
# configure them, make a function for this.

# Currently, any python functions called here may NOT return MDT objects; pyccc
# will try to depickle them (and fail)



def main(inputjson):
    print 'input:', inputjson

    # 1. read input
    build_mol = run_mdt(common.read_molecule, inputjson)

    # 2. doublet minimization (3 jobs - generate / run / parse)
    make_doublet_job = run_mdt(prep_doublet_minimization,
                               inputs={'in.pkl': build_mol.get_output('out.pkl')})
    doublet_job = make_doublet_job.result
    # note: retrieving the result here means we will probably be copying
    # the input files onto the master server - they're part of the job object
    # That's potentially bad, need to change it somehow

    submit_job(doublet_job,
               image=config.NWCHEMIMAGE)
    doublet_result = run_mdt(finish_doublet, doublet_job,
                             inputs={'in.pkl': make_doublet_job.get_output('out.pkl')})

    # 3. vertical singlet calculation
    make_singlet_job = run_mdt(prep_singlet_energy,
                               inputs={'in.pkl': doublet_result.get_output('out.pkl')})
    singlet_job = make_singlet_job.result

    submit_job(singlet_job,
               image=config.NWCHEMIMAGE)
    singlet_result = run_mdt(finish_singlet, singlet_job,
                             inputs={'in.pkl': make_singlet_job.get_output('out.pkl')})

    # 4. Write the PDB file
    pdbjob = run_mdt(to_pdb,
                     inputs={'in.pkl': doublet_result.get_output('out.pkl')})

    # 5. Parse results and write output
    resultjob = run_mdt(get_results,
                      inputs={'doublet.pkl': doublet_result.get_output('out.pkl'),
                              'singlet.pkl': singlet_result.get_output('out.pkl')
                      })
    with open('output.json', 'w') as outjson:
        json.dump(resultjob.result, outjson)

    pdbjob.get_output('out.pdb').put('out.pdb')


def prep_doublet_minimization():
    import moldesign as mdt

    mol = mdt.read('in.pkl')

    # anion doublet minimization
    mol.charge = -1 * mdt.units.q_e
    mol.set_energy_model(mdt.models.NWChemQM,
                         basis='sto-3g', theory='uks', multiplicity=2)
    mol.write('out.pkl')

    # PROBLEM: the when_finished callback can't be pickled. We need to abstract this?
    # Also: actually, we can pickle this, just need a custom __reduce__ or whatever ...
    newjob = mol.energy_model._make_minimization_job(50)
    newjob.when_finished = None
    return newjob


def finish_doublet(job):
    import moldesign as mdt

    mol = mdt.read('in.pkl')

    result = mol.energy_model.finish_min(job)
    result.write('traj.pkl')
    result.mol.write('out.pkl')


def prep_singlet_energy():
    import moldesign as mdt

    mol = mdt.read('in.pkl')

    # anion doublet minimization
    mol.charge = 0 * mdt.units.q_e
    mol.set_energy_model(mdt.models.NWChemQM,
                         basis='sto-3g', theory='rks', multiplicity=1)
    mol.write('out.pkl')

    # TODO: don't need to do this explicitly
    newjob = mol.energy_model._make_calculation_job(requests=['potential_energy'])
    newjob.when_finished = None
    return newjob


def finish_singlet(job):
    import moldesign as mdt

    mol = mdt.read('in.pkl')
    properties = mol.energy_model.finish(job)
    mol.positions = properties.positions  # they might be slightly different ...
    mol.properties = properties
    mol.write('out.pkl')


def get_results():
    import moldesign as mdt

    def json_quantity(q):
        return {'value': q.magnitude,
                'units': str(q.units)}

    singlet = mdt.read('singlet.pkl')
    doublet = mdt.read('doublet.pkl')

    return {'vde': json_quantity(singlet.potential_energy - doublet.potential_energy)}


def to_pdb():
    import moldesign as mdt
    mol = mdt.read('in.pkl')
    mol.write('out.pdb')




if __name__ == '__main__':
    with open(sys.argv[1], 'r') as jsonfile:
        inputs = yaml.load(jsonfile)

    main(inputs)

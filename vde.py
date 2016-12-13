import sys
import json
import yaml
import pyccc

import common

# what would make this easier? Make ALL nodes use compatible i/o formats
# Then we don't need to translate before/after sending to NWChem
# Otherwise we're stuck with these prepare input/run job/parse output chains

# Other options (not exclusive):
# 1. Somehow turn vanilla python scripts into workflow with clever magic (probably not possible)
# 2. Some sort of workflow-compatible tail-call-recursion for the translate-run-translate patterns
# 3. Nice wrappers to just make this all less verbose (probably the easiest)
# 4. Better integration with MDT's job launching system (probably necessary if we want this to be
#     usable ... not sure what it looks like

# 2/4: When MDT generates a job object, it should be ready to run. Actually, tail-call recursion
# works well here - jobs return job objects, which can be launched. It will be a common task to
# configure them, make a function for this.


# ENGINE = pyccc.CloudComputeCannon('http://cloudcomputecannon.bionano.autodesk.com:9000')
ENGINE = pyccc.Docker('unix://var/run/docker.sock')
MDTIMAGE = 'mdtscripts'
NWCHEMIMAGE = 'nwchem'


def main(inputjson):
    print 'input:', inputjson

    # 1. read input
    print 'reading input',
    inputjob = pyccc.PythonJob(ENGINE, MDTIMAGE,
                               pyccc.PythonCall(common.read_molecule, inputjson))
    print 'job:%s' % inputjob.jobid,
    inputjob.wait()
    print 'done'
    sys.stdout.flush()


    # 2. prep doublet minimization
    print 'prepping minimization',
    prepdoublet = pyccc.PythonJob(ENGINE, MDTIMAGE,
                                  pyccc.PythonCall(prep_doublet_minimization),
                                  inputs={'in.pkl': inputjob.get_output('out.pkl')})
    print 'job:%s' % inputjob.jobid,
    prepdoublet.wait()
    print 'done'
    sys.stdout.flush()


    # 3. run doublet minimization
    print 'running minimization',
    minjob = prepdoublet.result
    minjob.engine = ENGINE
    minjob.image = 'nwchem'
    minjob.submit()
    print 'job:%s' % inputjob.jobid,
    minjob.wait()
    print 'done'
    sys.stdout.flush()


    # 4. read minimization results
    processminjob = pyccc.PythonJob(ENGINE, MDTIMAGE,
                                    pyccc.PythonCall(process_minimization, minjob),
                                    inputs={'in.pkl': inputjob.get_output('out.pkl')})
    processminjob.wait()
    doublet_energy = processminjob.result


    # 5. prep singlet calculation
    prepsinglet = pyccc.PythonJob(ENGINE, MDTIMAGE,
                                  pyccc.PythonCall(prep_singlet_energy),
                                  inputs={'in.pkl': processminjob.get_output('out.pkl')})
    prepsinglet.wait()


    # 6. Run singlet calculation
    minjob = prepdoublet.result
    minjob.engine = ENGINE
    minjob.image = 'nwchem'
    minjob.submit()
    minjob.wait()


    # 7. read calculation results
    processenergy = pyccc.PythonJob(ENGINE, MDTIMAGE,
                                    pyccc.PythonCall(process_calculation, minjob),
                                    inputs={'in.pkl': processminjob.get_output('out.pkl')})
    processenergy.wait()
    singlet_energy = processminjob.result


    # 8. Get final PDB structure
    pdbjob = pyccc.PythonJob(ENGINE, MDTIMAGE,
                             pyccc.PythonCall(process_calculation),
                             inputs={'in.pkl': processminjob.get_output('out.pkl')})
    pdbjob.wait()


    # 8. return results
    assert singlet_energy['units'] == doublet_energy['units']
    vde = {'units': singlet_energy['units'],
           'value': singlet_energy['value'] - doublet_energy['value']}
    results = {'singlet_energy': singlet_energy,
               'doublet_energy': doublet_energy,
               'vertical detachment energy': vde}

    with open('output.json', 'w') as outjson:
        json.dump(results, outjson)

    pdbjob.get_output('out.pdb').put('out.pdb')


def prep_doublet_minimization():
    import moldesign as mdt

    mol = mdt.read('in.pkl')

    # anion doublet minimization
    mol.charge = -1 * mdt.units.q_e
    mol.set_energy_model(mdt.models.NWChemQM,
                         basis='sto-3g', theory='uks', multiplicity=2)
    mol.write('out.pkl')
    return mol.energy_model._make_minimization_job(50)


def process_minimization(oldjob):
    import moldesign as mdt
    mol = mdt.read('in.pkl')
    mdt.models.NWChemQM.finish_min(mol.energy_model, oldjob)
    mol.write('out.pkl')

    return json_quantity(mol.properties.potential_energy.to(mdt.units.eV))


def prep_singlet_energy():
    import moldesign as mdt

    mol = mdt.read('in.pkl')

    # anion doublet minimization
    mol.charge = 0 * mdt.units.q_e
    mol.set_energy_model(mdt.models.NWChemQM,
                         basis='sto-3g', theory='rks', multiplicity=1)
    mol.write('out.pkl')
    return mol.energy_model._make_calculation_job(50)


def process_calculation(job):
    import moldesign as mdt
    mol = mdt.read('in.pkl')
    mdt.models.NWChemQM.finish(mol.energy_model, job)
    mol.write('out.pkl')
    return json_quantity(mol.properties.potential_energy.to(mdt.units.eV))


def to_pdb():
    import moldesign as mdt
    mol = mdt.read('in.pkl')
    mol.write('out.pdb')


def json_quantity(q):
    return {'value': q.magnitude,
            'units': str(q.units)}



if __name__ == '__main__':
    with open(sys.argv[1], 'r') as jsonfile:
        inputs = yaml.load(jsonfile)

    main(inputs)

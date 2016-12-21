import common
import config

from pyccc import workflow


DOUBLET_PARAMS = dict(theory='uks',
                      functional='b3lyp',
                      basis='6-31g**',
                      charge=-1,
                      multiplicity=2)

SINGLET_PARAMS = dict(theory='uks',
                      functional='b3lyp',
                      basis='6-31g**',
                      charge=0,
                      multiplicity=1)


vde = workflow.Workflow('vertical detachment energy',
                        default_docker_image=config.MDTIMAGE)

read_molecule = vde.task(common.read_molecule,
                         inputjson=vde.input('molecule_json'))


@vde.task(name='doublet minimization',
          image=config.MDTNWCIMAGE,
          mol=read_molecule['mol'],
          params=DOUBLET_PARAMS,
          nsteps=50)
def minimize(mol, params, nsteps=None):
    import moldesign as mdt
    mol.set_energy_model(mdt.models.NWChemQM,
                         **params)
    traj = mol.minimize(nsteps=nsteps)

    return {'traj': traj,
            'mol': mol}


@vde.task(name='singlet single point',
          image=config.MDTNWCIMAGE,
          mol=minimize['mol'],
          params=workflow.Constant(SINGLET_PARAMS))
def single_point(mol, params, nsteps=None):
    import moldesign as mdt
    mol.set_energy_model(mdt.models.NWChemQM,
                         **params)
    mol.calculate()

    return {'mol': mol}


@vde.task(name='singlet single point',
          mol=vde.tasks['doublet minimization']['mol'])
def write_pdb(mol):
    return {'pdbfile': mol.write(format='pdb')}


@vde.task(doublet=vde.tasks['doublet minimization']['mol'],
          singlet=vde.tasks['singlet single point']['mol'])
def get_results(singlet, doublet):
    def json_quantity(q):
        return {'value': q.magnitude,
                'units': str(q.units)}

    return {'vde': json_quantity(singlet.potential_energy - doublet.potential_energy)}


vde.set_outputs(pdb=write_pdb['pdbfile'],
                results_json=get_results['vde'])
"""
Calculate the vertical detachment energy of an open shell, anionic species.


"""

import common
import config

from pyccc import workflow

vde = workflow.Workflow('Photoelectron spectrum calculator',
                        default_docker_image=config.MDTIMAGE)

read_molecule = vde.task(common.read_molecule,
                         inputjson=vde.input('molecule_json'))


@vde.task(image=config.MDTNWCIMAGE,
          mol=read_molecule['mol'],
          nsteps=50)
def minimize_doublet(mol, nsteps=None):
    import moldesign as mdt
    params = dict(theory='uks',
                  functional='b3lyp',
                  basis='6-31g**',
                  charge=-1)
    mol.set_energy_model(mdt.models.NWChemQM,
                         **params)
    traj = mol.minimize(nsteps=nsteps)
    return {'traj': traj,
            'mol': mol}


@vde.task(name='singlet single point',
          image=config.MDTNWCIMAGE,
          mol=minimize_doublet['mol'])
def single_point_singlet(mol):
    import moldesign as mdt
    from moldesign import units as u

    mol.charge = 0 * u.q_e
    params = dict(theory='rks',
                  functional='b3lyp',
                  basis='6-31g**',
                  charge=0)
    mol.set_energy_model(mdt.models.NWChemQM,
                         **params)
    mol.calculate()

    return {'mol': mol}


@vde.task(name='singlet single point',
          mol=vde.tasks['doublet minimization']['mol'])
def write_pdb(mol):
    return {'pdbfile': mol.write(format='pdb')}


@vde.task(doublet=vde.tasks['minimize_doublet']['mol'],
          singlet=vde.tasks['single_point_singlet']['mol'])
def get_results(singlet, doublet):

    resultjson = {'vde': (doublet.potential_energy - singlet.potential_energy).to_json()}

    return {'results.json': resultjson}


vde.set_outputs(pdb=write_pdb['pdbfile'],
                results_json=get_results['results.json'])

vde.metadata = workflow.Metadata(authors=["Aaron Virshup", "Marat Valiev"],
                                 affiliations=["Autodesk, Inc.",
                                               "Pacific Northwest National Laboratory"],
                                 corresponding="Marat Valiev",
                                 description=__doc__
                                 )

"""
Calculate the vertical detachment energy of an open shell, anionic species.


"""

from .. import common

from pyccc import workflow

_VERSION = "0.0.1b1"
MDTIMAGE = 'docker.io/avirshup/mst:mdt_subprocess-%s' % _VERSION
NWCHEMIMAGE = 'docker.io/avirshup/mst:mdt_nwchem-%s' % _VERSION
MDTAMBERTOOLS = 'docker.io/avirshup/mst:mdt_ambertools-%s' % _VERSION

vde = workflow.Workflow('Photoelectron spectrum calculator',
                        default_docker_image=MDTIMAGE)

read_molecule = vde.task(common.read_molecule,
                         description=vde.input('molecule_json'))


@vde.preprocessor
@vde.task(mol=read_molecule['mol'])
def validate(mol):
    from moldesign import units as u

    success = True
    errors = []

    if mol.num_atoms > 30:
        success = False
        errors.append("Molecule too large; max 30 atoms.")

    neutral_electrons = sum(atom.atnum for atom in mol.atoms)
    if neutral_electrons % 2 != 0:
        success = False
        errors.append("Molecule must have a closed shell neutral state.")

    if neutral_electrons > 200:
        success = False
        errors.append("Too many electrons (max 200).")

    return {'success': success,
            'errors': ' '.join(errors),
            'pdbstring': mol.write('pdb')}


@vde.task(image=NWCHEMIMAGE,
          mol=read_molecule['mol'],
          nsteps=50)
def minimize_doublet(mol, nsteps=None):
    import moldesign as mdt
    mol.charge = -1 * mdt.units.q_e
    params = dict(theory='uks',
                  functional='b3lyp',
                  basis='6-31g',
                  charge=-1 * mdt.units.q_e,
                  multiplicity=2)
    mol.set_energy_model(mdt.models.NWChemQM,
                         **params)
    traj = mol.minimize(nsteps=nsteps)
    return {'traj': traj,
            'mol': mol,
            'pdbstring':mol.write(format='pdb')}


@vde.task(image=NWCHEMIMAGE,
          mol=minimize_doublet['mol'])
def single_point_singlet(mol):
    import moldesign as mdt
    from moldesign import units as u

    mol.charge = 0 * mdt.units.q_e
    params = dict(theory='rks',
                  functional='b3lyp',
                  basis='6-31g')
    mol.set_energy_model(mdt.models.NWChemQM,
                         **params)
    mol.calculate()
    return {'mol': mol}


@vde.task(doublet=minimize_doublet['mol'],
          singlet=single_point_singlet['mol'])
def get_results(singlet, doublet):
    results = {'vde': (doublet.potential_energy - singlet.potential_energy).to_json(),
               'singlet_energy': singlet.potential_energy.to_json(),
               'doublet_energy': doublet.potential_energy.to_json()}
    return {'results': results}


outputs = {'final_structure.pdb': minimize_doublet['pdbstring'],
           'results': get_results['results']}

vde.set_outputs(**outputs)

#vde.metadata = workflow.Metadata(authors=["Aaron Virshup", "Marat Valiev"],
#                                 affiliations=["Autodesk, Inc.",
#                                               "Pacific Northwest National Laboratory"],
#                                 corresponding="Marat Valiev",
#                                 description=__doc__
#                                 )

"""
Calculate the vertical detachment energy of an open shell, anionic species.


"""

from .. import common
from ..utils import get_asset

from pyccc import workflow

_VERSION = "0.0.1b6"
MDTIMAGE = 'docker.io/avirshup/mst:mdt_subprocess-%s' % _VERSION
NWCHEMIMAGE = 'docker.io/avirshup/mst:mdt_nwchem-%s' % _VERSION
MDTAMBERTOOLS = 'docker.io/avirshup/mst:mdt_ambertools-%s' % _VERSION


METADATA = dict(id='0',
                version_id='1beta6',
                title='Refine a protein-ligand complex',
                selectLigands=False,
                bgIndex=3,
                bgColor='#292E60',
                color='#2FE695',
                comingSoon=False,
                bgImg=get_asset('boundligand.png'),
                description='Calculate the electron binding energy of an' +
                            ' anionic doublet species using DFT',)

vde = workflow.Workflow('Photoelectron spectrum calculator',
                        default_docker_image=MDTIMAGE)

read_molecule = vde.task(common.read_molecule,
                         description=vde.input('molecule_json'),
                         __interactive__=True)


@vde.preprocessor
@vde.task(mol=read_molecule['mol'],
          __interactive__=True)
def validate(mol):
    from moldesign import units as u

    success = True
    errors = []

    MAXATOMS = 50
    MAXELECTRONS = 200

    if mol.num_atoms > MAXATOMS:
        success = False
        errors.append("Molecule too large; max %d atoms." % MAXATOMS)

    neutral_electrons = sum(atom.atnum for atom in mol.atoms)
    if neutral_electrons % 2 != 0:
        success = False
        errors.append("Molecule must have a closed shell neutral state.")

    if neutral_electrons > MAXELECTRONS:
        success = False
        errors.append("Too many electrons (max %d)." % MAXELECTRONS)

    return {'success': success,
            'errors': ' '.join(errors),
            'pdbstring': mol.write('pdb')}


@vde.task(__image__=NWCHEMIMAGE,
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


@vde.task(__image__=NWCHEMIMAGE,
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
    from moldesign import units as u
    vde = singlet.potential_energy - doublet.potential_energy
    vde_json = vde.to(u.eV).to_json()
    vde_json['name'] = 'VDE'
    results = {'vde': (vde).to(u.eV).to_json(),
               'singlet_energy': singlet.potential_energy.to(u.eV).to_json(),
               'doublet_energy': doublet.potential_energy.to(u.eV).to_json(),
               'output_values': [vde_json]}
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

""" A workflow to automatically assign a forcefield to a biomolecule, including those with
non-standard residues.

Currently hardcoded to use AMBER14 and GAFF. User can choose the charge model for any
ligands

"""
import pickle

from .. import common, interactive
from ..common import missing_internal_residues

from pyccc import workflow

_VERSION = "0.0.alpha0"
MDTIMAGE = 'docker.io/avirshup/mst:mdt_subprocess-%s' % _VERSION
NWCHEMIMAGE = 'docker.io/avirshup/mst:mdt_nwchem-%s' % _VERSION
MDTAMBERTOOLS = 'docker.io/avirshup/mst:mdt_ambertools-%s' % _VERSION

prep_active_site = workflow.Workflow('Prep active site calculations',
                                     default_docker_image=MDTIMAGE)

read_molecule = prep_active_site.task(common.read_molecule,
                                      description=prep_active_site.input('molecule_json'))


@prep_active_site.task(mol=read_molecule['mol'])
def validate(mol):
    missing = missing_internal_residues(mol)
    all_errors = []

    for chain_name, reslist in missing.iteritems():
        all_errors.append('The following residues are not present in the PDB file: %s'
                          % ','.join(name+num for num, name in reslist.iteritems()))

    return {'success': not all_errors,
            'errors': all_errors}


@prep_active_site.task(mol=read_molecule['mol'])
def get_ligands(mol):
    """ Return a dict of possible ligands in the molecule.
    dict is of the form {ligand_name: [atom_idx1, atom_idx2, ...], ...}
    where the atom_idxN are the 0-based indicies of the atoms comprising each potential ligand
    """
    ligands = [residue for residue in mol.residues
               if residue.type == 'unknown']

    found_ligands = {}
    for ligand in ligands:
        bound_residues = list(ligand.bonded_residues)

        if len(bound_residues) == 0:
            found_ligands[ligand.name] = [atom.index for atom in ligand.atoms]
        elif len(bound_residues) == 1:
            nbr = bound_residues[0]
            if len(list(nbr.bonded_residues)) != 1:
                continue
            else:
                found_ligands['%s - %s' % (ligand.name, nbr.name)] = \
                    [atom.idx for atom in ligand.atoms+nbr.atoms]

    return {'ligand_options': found_ligands,
            'pdbfile': mol.write('pdb'),
            'molpickle': pickle.dumps(mol)}

prep_active_site.set_outputs(ligand_options=get_ligands['ligand_options'],
                             pdbstring=get_ligands['pdbfile'],
                             molpickle=get_ligands['molpickle'],
                             validates=validate['success'],
                             validation_errors=validate['errors'])


run_optimization = workflow.Workflow('Run QM/MM active site calculation',
                                     default_docker_image=MDTIMAGE)


@run_optimization.task(pickleurl=run_optimization.input('pickleurl'))
def getinputmol(pickleurl):
    import pickle
    import requests
    r = requests.get(pickleurl)
    mol = pickle.loads(r.text)
    return {'mol': mol}


@run_optimization.task(mol=getinputmol['mol'],
                       ligand_atom_ids=run_optimization.input('atom_ids'),
                       image=MDTAMBERTOOLS)
def prep_ligand(mol, ligand_atom_ids):
    """
    Create force field parameters for the chosen ligand
    """
    import moldesign as mdt

    ligand = mdt.Molecule([mol.atoms[idx] for idx in ligand_atom_ids])
    #ligh = mdt.set_hybridization_and_ph(ligand, 7.4)
    ligh = ligand
    params = mdt.interfaces.ambertools.parameterize(ligh, charges='gasteiger')
    return {'ligand_parameters': params,
            'ligand': ligh}


@run_optimization.task(mol=getinputmol['mol'],
                       ligand_atom_ids=run_optimization.input('atom_ids'),
                       ligand_params=prep_ligand['ligand_parameters'],
                       image=MDTAMBERTOOLS)
def prep_forcefield(mol, ligand_atom_ids, ligand_params):
    """
    Assign forcefield to the protein/ligand complex
    """
    import moldesign as mdt

    mdt.guess_histidine_states(mol)
    stripmol = mdt.Molecule([atom for atom in mol.atoms
                             if atom.residue.type in ('dna', 'protein')
                             or atom.index in ligand_atom_ids])
    withff = mdt.interfaces.ambertools.assign_forcefield(stripmol, parameters=ligand_params)

    return {'molecule': withff,
            'prmtop': withff.ff.amber_params.prmtop,
            'inpcrd': withff.ff.amber_params.inpcrd}


@run_optimization.task(mol=prep_forcefield['molecule'])
def mm_minimization(mol):
    import moldesign as mdt
    mol.set_energy_model(mdt.models.OpenMMPotential, implicit_solvent='obc')
    traj = mol.minimize(nsteps=2000)

    return {'trajectory': traj,
            'molecule': mol,
            'initial_energy': traj.frames[0].potential_energy.to_json(),
            'final_energy': mol.potential_energy.to_json(),
            'rmsd': traj.rmsd()[-1].to_json()}


@run_optimization.task(traj=mm_minimization['trajectory'])
def result_coordinates(traj):
    m = traj.mol
    return {'finalpdb': m.write('pdb')}


result_display = run_optimization.task(interactive.protein_minimization_display,
                                       initial_energy=mm_minimization['initial_energy'],
                                       final_energy=mm_minimization['final_energy'],
                                       rmsd=mm_minimization['rmsd'],
                                       trajectory=mm_minimization['trajectory'])

run_optimization.set_outputs(prmtop=prep_forcefield['prmtop'],
                             inpcrd=prep_forcefield['inpcrd'],
                             finalpdb=result_coordinates['finalpdb'])

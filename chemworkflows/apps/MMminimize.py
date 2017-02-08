""" A workflow to automatically assign a forcefield to a biomolecule, including those with
non-standard residues.

Currently hardcoded to use AMBER14 and GAFF. User can choose the charge model for any
ligands

"""
from .. import common, interactive
from ..common import missing_internal_residues

from pyccc import workflow

_VERSION = "0.0.1b1"
MDTIMAGE = 'docker.io/avirshup/mst:mdt_subprocess-%s' % _VERSION
NWCHEMIMAGE = 'docker.io/avirshup/mst:mdt_nwchem-%s' % _VERSION
MDTAMBERTOOLS = 'docker.io/avirshup/mst:mdt_ambertools-%s' % _VERSION

minimization = workflow.Workflow('Refine ligand binding site',
                                 default_docker_image=MDTIMAGE)

read_molecule = minimization.task(common.read_molecule,
                                  description=minimization.input('molecule_json'))


@minimization.task(mol=read_molecule['mol'])
def get_ligands(mol):
    """ Return a dict of possible ligands in the molecule.
    dict is of the form {ligand_name: [atom_idx1, atom_idx2, ...], ...}
    where the atom_idxN are the 0-based indices of the atoms comprising each potential ligand

    Ligands can span, at most, 2 different residues
    """
    ligands = [residue for residue in mol.residues
               if residue.type == 'unknown']

    found_ligands = {}
    mv_ligand_strings = {}  # temp - selection strings for molecule viewer
    for ligand in ligands:
        bound_residues = list(ligand.bonded_residues)

        if len(bound_residues) == 0:
            found_ligands[ligand.name] = [atom.index for atom in ligand.atoms]
            mv_ligand_strings[ligand.name] = ['1.%s.%s-%d'%(ligand.chain.name,
                                                            ligand.chain.name,
                                                            ligand.pdbindex)]
        elif len(bound_residues) == 1:
            nbr = bound_residues[0]
            if len(list(nbr.bonded_residues)) != 1:
                continue
            else:
                ligname = '%s - %s' % (ligand.name, nbr.name)
                found_ligands[ligname] = [atom.idx for atom in ligand.atoms+nbr.atoms]
                mv_ligand_strings[ligname] = ['1.%s.%s-%d'%(l.chain.name,
                                                            l.chain.name,
                                                            l.pdbindex)
                                              for l in (ligand, nbr)]

    return {'ligand_options': found_ligands,
            'mv_ligand_strings': mv_ligand_strings}


@minimization.preprocessor
@minimization.task(mol=read_molecule['mol'],
                   ligands=get_ligands['ligand_options'],
                   mv_ligand_strings=get_ligands['mv_ligand_strings'])
def validate(mol, ligands, mv_ligand_strings):
    missing = missing_internal_residues(mol)
    all_errors = []
    success = True

    if len(ligands) == 0:
        success = False
        all_errors.append('No ligands found in this structure.')

    for chain_name, reslist in missing.iteritems():
        success = False
        all_errors.append('The following residues are not present in the PDB file: %s'
                          % ','.join(name+num for num, name in reslist.iteritems()))

    return {'success': success,
            'errors': ' '.join(all_errors),
            'pdbstring': mol.write('pdb'),
            'ligands': ligands,
            'mv_ligand_strings': mv_ligand_strings}


atomselection = minimization.task(interactive.SelectAtomsFromOptions(),
                                  taskname='user_atom_selection',
                                  choices=get_ligands['ligand_options'])


@minimization.task(mol=read_molecule['mol'],
                   ligand_atom_ids=atomselection['atom_ids'],
                   ligandname=atomselection['ligandname'],
                   image=MDTAMBERTOOLS)
def prep_ligand(mol, ligand_atom_ids, ligandname):
    """
    Create force field parameters for the chosen ligand
    """
    import moldesign as mdt

    print 'Parameterizing ligand "%s"' % ligandname

    ligand = mdt.Molecule([mol.atoms[idx] for idx in ligand_atom_ids])
    ligh = mdt.set_hybridization_and_ph(ligand, 7.4)
    #ligh = ligand
    params = mdt.interfaces.ambertools.parameterize(ligh, charges='gasteiger')
    return {'ligand_parameters': params,
            'ligand': ligh}


@minimization.task(mol=read_molecule['mol'],
                   ligand_atom_ids=atomselection['atom_ids'],
                   ligand_params=prep_ligand['ligand_parameters'],
                   image=MDTAMBERTOOLS)
def prep_forcefield(mol, ligand_atom_ids, ligand_params):
    """
    Assign forcefield to the protein/ligand complex
    """
    import moldesign as mdt

    for residue in mol.residues:
        if residue.resname == 'HIS':
            print 'Guessing histidine protonation based on hydrogens present in file:'
            mdt.guess_histidine_states(mol)
            break

    stripmol = mdt.Molecule([atom for atom in mol.atoms
                             if atom.residue.type in ('dna', 'protein')
                             or atom.index in ligand_atom_ids])
    withff = mdt.interfaces.ambertools.assign_forcefield(stripmol, parameters=ligand_params)

    return {'molecule': withff,
            'prmtop': withff.ff.amber_params.prmtop,
            'inpcrd': withff.ff.amber_params.inpcrd}


@minimization.task(mol=prep_forcefield['molecule'])
def mm_minimization(mol):
    import moldesign as mdt
    mol.set_energy_model(mdt.models.OpenMMPotential, implicit_solvent='obc')
    traj = mol.minimize(nsteps=2000)

    results = {'initial_energy': traj.frames[0].potential_energy.to_json(),
               'final_energy': mol.potential_energy.to_json(),
               'rmsd': traj.rmsd()[-1].to_json()}

    minstep_file, minstep_filenames = _traj_to_tarred_pdbs(traj)

    return {'trajectory': traj,
            'molecule': mol,
            'pdbstring': mol.write(format='pdb'),
            'results': results,
            'minstep_frames': minstep_filenames,
            'minsteps.tar.gz': minstep_file}

minimization.set_outputs(**{'prmtop': prep_forcefield['prmtop'],
                            'inpcrd': prep_forcefield['inpcrd'],
                            'results': mm_minimization['results'],
                            'final_structure.pdb': mm_minimization['pdbstring'],
                            'minsteps.tar.gz': mm_minimization['minsteps.tar.gz'],
                            'minstep_frames': mm_minimization['minstep_frames']
                            })

def _traj_to_tarred_pdbs(traj):
    """ Quick hack to create a tarfile with frames from a trajectory
    """
    import os, pyccc

    minstep_filenames = []
    for i, frame in enumerate(traj.frames):
        fname = 'minstep.%d.pdb' % i
        frame.write(fname)
        minstep_filenames.append(fname)

    os.system('tar czvf minsteps.tar.gz %s' % ' '.join(minstep_filenames))
    return pyccc.LocalFile('minsteps.tar.gz'), minstep_filenames

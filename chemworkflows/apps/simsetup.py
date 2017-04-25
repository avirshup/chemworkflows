""" A workflow to automatically assign a forcefield to a biomolecule, including those with
non-standard residues.

Currently hardcoded to use AMBER14 and GAFF. User can choose the charge model for any
ligands

"""
from .. import common, interactive
from ..common import missing_internal_residues

from pyccc import workflow

_VERSION = "0.0.1b6"
MDTIMAGE = 'docker.io/avirshup/mst:mdt_subprocess-%s' % _VERSION
MDTAMBERTOOLS = 'docker.io/avirshup/mst:mdt_ambertools-%s' % _VERSION
MDTLAMMPS = 'docker.io/avirshup/mst:mdt_lammps-%s' % _VERSION

simsetup = workflow.Workflow('Set up a simulation',
                             default_docker_image=MDTIMAGE)

read_molecule = simsetup.task(common.read_molecule,
                              description=simsetup.input('molecule_json'),
                              __interactive__=True)


@simsetup.task(mol=read_molecule['mol'],
               __interactive__=True)
def get_ligands(mol):
    """ Return a dict of possible ligands in the molecule.
    dict is of the form {ligand_name: [atom_idx1, atom_idx2, ...], ...}
    where the atom_idxN are the 0-based indices of the atoms comprising each potential ligand

    Ligands can span, at most, 2 different residues
    """
    if mol.num_residues == 1:
        ligands = [mol.residues[0]]
    else:
        ligands = [residue for residue in mol.residues if residue.type == 'unknown']

    found_ligands = {}
    mv_ligand_strings = {}  # temp - selection strings for molecule viewer
    for ligand in ligands:
        bound_residues = list(ligand.bonded_residues)

        if len(bound_residues) == 0:
            found_ligands[ligand.name] = [atom.index for atom in ligand.atoms]
            mv_ligand_strings[ligand.name] = ['1.%s.%s-%s'%(ligand.chain.name,
                                                            ligand.chain.name,
                                                            ligand.pdbindex)]
        elif len(bound_residues) == 1:
            nbr = bound_residues[0]
            if len(list(nbr.bonded_residues)) != 1:
                continue
            else:
                ligname = '%s - %s' % (ligand.name, nbr.name)
                found_ligands[ligname] = [atom.idx for atom in ligand.atoms+nbr.atoms]
                mv_ligand_strings[ligname] = ['1.%s.%s-%s'%(l.chain.name,
                                                            l.chain.name,
                                                            l.pdbindex)
                                              for l in (ligand, nbr)]

    return {'ligand_options': found_ligands,
            'mv_ligand_strings': mv_ligand_strings}


@simsetup.preprocessor
@simsetup.task(mol=read_molecule['mol'],
               ligands=get_ligands['ligand_options'],
               __interactive__=True)
def validate(mol, ligands):
    return {'success': True,
            'errors': ''}


atomselection = simsetup.task(interactive.SelectAtomsFromOptions(),
                              __taskname__='user_atom_selection',
                              choices=get_ligands['ligand_options'],
                              __interactive__=True)


@simsetup.task(mol=read_molecule['mol'],
               ligand_atom_ids=atomselection['atom_ids'],
               ligandname=atomselection['ligandname'],
               __image__=MDTAMBERTOOLS)
def prep_ligand(mol, ligand_atom_ids, ligandname):
    """
    Create force field parameters for the chosen ligand
    """
    import moldesign as mdt

    print 'Parameterizing ligand "%s"' % ligandname

    ligand = mdt.Molecule([mol.atoms[idx] for idx in ligand_atom_ids])
    ligh = mdt.set_hybridization_and_ph(ligand, 7.4)
    params = mdt.interfaces.ambertools.parameterize(ligh, charges='gasteiger')
    return {'ligand_parameters': params,
            'ligand': ligh}


@simsetup.task(mol=read_molecule['mol'],
               ligmol=prep_ligand['ligand'],
               ligand_atom_ids=atomselection['atom_ids'],
               ligand_params=prep_ligand['ligand_parameters'],
               __image__=MDTAMBERTOOLS)
def prep_forcefield(mol, ligmol, ligand_atom_ids, ligand_params):
    """
    Assign forcefield to the protein/ligand complex
    """
    import moldesign as mdt

    # hackity hackity hackity hack
    if mol.num_residues == 1:
        return {'molecule': ligmol}

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


@simsetup.task(mol=prep_forcefield['molecule'],
               __image__=MDTLAMMPS)
def write_lammps_setup(mol):
    # note: this is REALLY REALLY HACKY because there's a version mismatch that is not worth fixing
    import moldesign as mdt
    import parmed

    class FakeModel(object):
        pass

    model = FakeModel()
    model.mol = mol
    mol.ff.parmed_obj = parmed.amber.AmberParm(mol.ff.amber_params.prmtop.open())
    lines = mdt.models.LAMMPSPotential._get_lammps_data.__func__(model)

    return {'mol.pdb': mol.write('pdb'),
            'lammps.data': lines}


simsetup.set_outputs(**{'mol.pdb': write_lammps_setup['mol.pdb'],
                        'lammps.data': write_lammps_setup['lammps.data']
                        })

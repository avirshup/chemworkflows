import json


def missing_internal_residues(mol):
    """ Return a list of missing internal residues
    """
    missing = mol.metadata.get('missing_residues', None)
    if missing is None:
        return {}

    missing_internal = {}

    for chain in mol.chains:
        missres = missing.get(chain.pdbname, {})
        if not missres: continue
        imin, imax = chain.n_terminal.pdbindex, chain.c_terminal.pdbindex
        for resindex, resname in missres.iteritems():
            if imin < resindex < imax:
                missing_internal.setdefault(chain.pdbname, {})[resindex] = resname

    return missing_internal


def read_molecule(description):
    """ All-purpose routine for initializing molecules.
    The input "description" must be a yaml or JSON file with exactly one
    of the following key-value pairs:
     - SMILES: [a smiles string]
     - IUPAC: [an IUPAC string]
     - inCHI: [an inchi identifier]
     - PDB: [a 4-letter PDB code]
     - filename: filename (indicates that the molecule should be read from the
               passed "molfile"

    If "asfile" is passed, then "molfile" should also be present. The format
    will be determined from the filename passed in the description JSON
    """
    import moldesign as mdt

    #d = json.loads(description)

    d = description
    assert len(d) == 1, "%s %s" %(d, len(d))

    if 'filename' in d:
        format, compression = mdt.fileio._get_format(d['filename'])
        m = mdt.read(description['content'], format=format)
    elif 'smiles' in d:
        m = mdt.from_smiles(d['smiles'])
    elif 'iupac' in d:
        assert len(d) == 1
        m = mdt.from_name(d['iupac'])
    elif 'inchi' in d:
        assert len(d) == 1
        m = mdt.from_inchi(d['inchi'])
    else:
        raise ValueError(description)

    m.write('out.pkl')
    m.write('out.pdb')


#def write(fmt):
#   import moldesign as mdt
#    mol = mdt.read('in.pkl')

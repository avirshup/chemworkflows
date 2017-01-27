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
    if 'filename' in d:
        format, compression = mdt.fileio._get_format(d['filename'], None)
        m = mdt.read(description['content'], format=format)
    elif 'smiles' in d:
        m = mdt.from_smiles(d['smiles'])
    elif 'iupac' in d:
        assert len(d) == 1
        m = mdt.from_name(d['iupac'])
    elif 'inchi' in d:
        assert len(d) == 1
        m = mdt.from_inchi(d['inchi'])
    elif 'pdb' in d:
        assert len(d) == 1
        m = mdt.from_pdb(d['pdb'])

    elif 'input' in d:
        assert len(d) == 1
        content = d['input']

        if len(d) == 4 and d[0].isdigit():
            try:
                m = mdt.from_pdb(content)
            except:
                pass
            else:
                print 'Reading molecule as PDB ID %s' % d
                return {'mol':m}

        try:
            m = mdt.from_smiles(content)
        except:
            pass
        else:
            print 'Reading molecule as smiles %s' % d
            return {'mol':m}

        try:
            m = mdt.from_name(content)
        except:
            pass
        else:
            print 'Reading molecule as IUPAC name %s' % d
            return {'mol': m}

        try:
            m = mdt.from_inchi(content)
        except:
            pass
        else:
            print 'Reading molecule as inchi string %s' % d
            return {'mol': m}

        raise ValueError(d)


    else:
        raise ValueError(description)

    return {'mol':m}


#def write(fmt):
#   import moldesign as mdt
#    mol = mdt.read('in.pkl')

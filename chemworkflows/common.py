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
        data = d['input']

        if len(data) == 4 and data[0].isdigit():
            try:
                m = mdt.from_pdb(data)
            except Exception as e:
                print 'Not recognized as PDB ID', e
            else:
                print 'Reading molecule as PDB ID "%s"' % data
                return {'mol':m}

        try:
            m = mdt.from_smiles(data)
        except Exception as e:
            print 'Not recognized as SMILES name', e
        else:
            print 'Reading molecule as smiles "%s"' % data
            return {'mol':m}

        try:
            m = mdt.from_name(data)
        except Exception as e:
            print 'Not recognized as IUPAC name', e
        else:
            print 'Reading molecule as IUPAC name "%s"' % data
            return {'mol': m}

        try:
            m = mdt.from_inchi(data)
        except Exception as e:
            print 'Not recognized as INCHI string', e
        else:
            print 'Reading molecule as inchi string "%s"' % data
            return {'mol': m}

        raise ValueError("Failed to parse input data '%s' as PDB id, SMILES, IUPAC, or INCHI" %
                         data)


    else:
        raise ValueError(description)

    return {'mol': m}


#def write(fmt):
#   import moldesign as mdt
#    mol = mdt.read('in.pkl')

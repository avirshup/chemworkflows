# MDT Workflow apps

These are lightweight workflow scripts that use CCC to run MDT on the cloud. The scripts don't do any computation locally - everything happens in external docker containers, these just link the all the inputs and outputs together.

Note MDT does NOT need to be installed locally. The only requirements are the CCC libraries

### Install

```
git clone https://github.com/avirshup/chemworkflows.git
pip install -r chemworkflows/requirements.txt
pip install -e ./chemworkflows
```

### Run

This example runs the VDE calculator on N<sub>2</sub>.

1. The job can be launched with a JSON description of the input structure, or with the path to an input file.
 
```bash
$ chemworkflow vde '{"iupac":"water"}' 
$ chemworkflow vde ethylene.xyz
```


1. Running it will launch a series of docker computes. (By default, the will run in the cloud. Use the `--localdocker` flag to run them in your local docker instance).
```bash
$ chemworkflow vde '{"smiles":"N#N"}' 
input: {'smiles': 'N#N'}
Running python:read_molecule ...id:ddabf6a27bdd7f05267061064bb54f1a3d7e386522537565d9ff448e27416aa5 ...done
Running python:prep_doublet_minimization ...id:5a775a76244ecb54c46b357c131c6951ba5200d5d6a3b1266b845f741d69265b ...done
Running job: nwchem/in.pkl ...id:18baf291a9d82348a425e47828f5754157378f9fbb722067c439070454869605 ...done
Running python:finish_doublet ...id:897afe9d9db1c37788997889f475881d3eb5efef34fd432d885843e9ad23402b ...done
Running python:prep_singlet_energy ...id:3f5ac8584a016e4904bebc90c70b6e279530dd85c516406d05f639bbfc21491d ...done
```

4. When the job is done, it will produce a directory with the workflow's outputs:

```bash
$ ls vde.out.0.json
final_structure.pdb
results.json
workflow_state.dill
```







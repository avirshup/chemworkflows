# MDT Workflow apps

These are lightweight workflow scripts that use CCC to run MDT on the cloud. The scripts don't do any computation locally - everything happens in external docker containers, these just link the all the inputs and outputs together.

Note MDT does NOT need to be installed locally. The only requirements are the CCC libraries

### Install

```
git clone https://github.com/avirshup/mdtscripts.git
cd mdtscripts
pip install -r requirements.txt
```

### Run

This example runs the VDE calculator on N<sub>2</sub>.

1. The basic input file, which has the same format for all workflows, describes the molecular structure input:
 ```
$ cat inputs/n2.json
{"smiles": "N#N"}
```

2. Pass the input file's path to the `vde` python script:
```
python vde.py inputs/n2.json
```

3. When the job is done, it will produce two files in the working directory:

```
$ cat result.json  # contains numerical results
{...}

$ cat out.pdb  # contains final structure
[...]
```







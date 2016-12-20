import pyccc

ENGINE = pyccc.CloudComputeCannon('http://cloudcomputecannon.bionano.autodesk.com:9000')
#ENGINE = pyccc.Docker()
MDTIMAGE = 'docker.io/avirshup/mdtapptest'
NWCHEMIMAGE = 'docker.io/avirshup/nwchem:twelve'

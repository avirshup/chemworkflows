import pyccc

ENGINE = pyccc.CloudComputeCannon('http://cloudcomputecannon.bionano.autodesk.com:9000')
#ENGINE = pyccc.Docker()

MDTIMAGE = 'docker.io/avirshup/mst:moldesign_complete'
NWCHEMIMAGE = 'docker.io/avirshup/mst:mdt_nwchem'
MDTAMBERTOOLS = 'docker.io/avirshup/mst:mdt_ambertools'


# Get backbone

import prody
 
# Get the backbone structure
# Chains are deleted by Gromacs and the recovered by prody
# WARNING: This step is performed with prody instead of Gromacs since it would delete chains data
# DANI: Hay que provar que funciona, que antes se hacía con gromacs + prody y ahora solo prody
def get_backbone (
    topology_reference,
    output_backbone_filename : str):

    topology = topology_reference.topology
    backbone = topology.select('backbone')
    prody.writePDB(output_backbone_filename, backbone)

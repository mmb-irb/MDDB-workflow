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
    # WARNING: The default 'backbone' prody selection is equal to 'protein and name N CA C O'
    # WARNING: The equivalent 'Backbone' gromacs selection would be 'name N CA C' in prody
    # Since this backbone is used with the PCA results, which come from gromacs, it must match
    # Oxygens are removed, since gromacs results do not include oxygens (e.g. all backbone oxygens)
    # Protein restriction is removed, since gromacs would consider other things (e.g. caping residues)
    backbone = topology.select('name N CA C')
    prody.writePDB(output_backbone_filename, backbone)

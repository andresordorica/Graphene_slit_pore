import mbuild as mb
import numpy as np
import os
from foyer import Forcefield



def Get_ff_path(ff_name):
    """Get the path to a force field xml file """
    """in a directory of the same name."""
    cache_dir = '{}/files/'.format(os.getcwd())
    ff_path = os.path.join(cache_dir, ff_name + '.xml')
    return ff_path

def get_ntype(name):
    atom_bond = {
                "graphene": "C-C" 
                  
            }
    case = atom_bond[name]
    return case

def get_structure(compound, compund_name = "graphene",force_field_name = "carbon", create_bonds = False):
    ff_ = Forcefield(forcefield_files=[Get_ff_path(force_field_name)])

    if create_bonds ==True:
        list_pos =[]
        for particle in block.particles():
            pos_1 = particle.pos
            l = []
            for p in block.particles():
                if particle != p:
                    pos_2 = p.pos
                    distance = np.sqrt(((pos_1[0]-pos_2[0])**2)+((pos_1[1]-pos_2[1])**2)+((pos_1[2]-pos_2[2])**2))
                    l.append(distance)
                
            alf = min(l)
            list_pos.append(alf)

        var = max(list_pos)*1.3

        case = (get_ntype(compund_name)).split("-")
        compound.generate_bonds(case[0],case[1], dmin = 0.0, dmax = var)
        
    else:
        structure = ff_.apply(compound, residues = 'gra', assert_bond_params = False, assert_dihedral_params = False, assert_angle_params = False)

        for residue in structure.residues:
            residue.name = "gra"


        for at in  structure.adjust_types:
            at.chgscale = 0.5

    return structure




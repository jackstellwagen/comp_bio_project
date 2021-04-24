#import module_manager
#module_manager.review()

import Bio
from Bio import PDB
from Bio.PDB.ResidueDepth import residue_depth, get_surface
from Bio.PDB.DSSP import DSSP

#keys = AA index
#values = [residue, residue code, phi, psi, surface depth]

result = {}
structName = "1HMP"
fileName = "/Users/kcooo/Downloads/pdb4wxv.pdb"

def readPDBFile(structName, fileName):

    #code from https://warwick.ac.uk/fac/sci/moac/people/students/peter_cock/python/ramachandran/calculate/#BioPython
    for model in Bio.PDB.PDBParser().get_structure(structName, fileName) :
        
        surface = get_surface(model) #numpy array of all surface vertices of folded protein

        for chain in model :

            polypeptides = Bio.PDB.PPBuilder().build_peptides(chain)
            dssp = DSSP(chain, fileName)


            for poly_index, poly in enumerate(polypeptides) :

                startIndex = poly[0].id[1]
                endIndex = poly[-1].id[1]

                phiPsi = poly.get_phi_psi_list()

                for residue in range(len(poly)):

                    #takes first data point for overlapping chains
                    if poly[residue].id[1] not in result:

                        res = poly[residue].resname

                        #encoded as int from 0-1
                        resCode = Bio.PDB.Polypeptide.three_to_index(poly[residue].resname)/20

                        phi, psi = phiPsi[residue]
                        #phi/psi cannot be calculated at edges
                        if phi != None:
                            phi /= 360
                        if psi != None:
                            psi /= 360

                        #average depth of all atoms in residue from surface
                        depth = residue_depth(poly[residue], surface) 

                        result[poly[residue].id[1]] = [res, resCode, phi, psi, depth]
    print(result)
    return result


readPDBFile(structName, fileName)
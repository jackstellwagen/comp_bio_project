import module_manager
module_manager.review()

import Bio
from Bio import PDB
from Bio.PDB.ResidueDepth import ResidueDepth
#from Bio.PDB.DSSP import DSSP

#keys = AA index
#values = [residue, residue code, phi, psi]
result = {}
structName = "1HMP"
fileName = "/Users/kcooo/Downloads/pdb4wxv.pdb"

def readPDBFile(structName, fileName):

    #code from https://warwick.ac.uk/fac/sci/moac/people/students/peter_cock/python/ramachandran/calculate/#BioPython
    for model in Bio.PDB.PDBParser().get_structure(structName, fileName) :
        
        depth = ResidueDepth(model) 
        
        for chain in model :
            depth = ResidueDepth(chain)
            #dssp = DSSP(chain, "/Users/kcooo/Downloads/pdb4wxv.pdb", dssp='mkdssp')
            polypeptides = Bio.PDB.PPBuilder().build_peptides(chain)

            for poly_index, poly in enumerate(polypeptides) :

                startIndex = poly[0].id[1]
                endIndex = poly[-1].id[1]

                phiPsi = poly.get_phi_psi_list()

                for residue in range(len(poly)):
                    if poly[residue].id[1] not in result:

                        res = poly[residue].resname
                        resCode = Bio.PDB.Polypeptide.three_to_index(poly[residue].resname)/20
                        phi, psi = phiPsi[residue]
                        if phi != None:
                            phi /= 360
                        if psi != None:
                            psi /= 360

                        result[poly[residue].id[1]] = [res, resCode, phi, psi]
    print(result)
    return result


readPDBFile(structName, fileName)
#import module_manager
#module_manager.review()

import Bio
from Bio import PDB
from Bio.PDB.ResidueDepth import residue_depth, get_surface
from Bio.PDB.DSSP import DSSP
from Bio.PDB import Selection

#keys = AA index
#values = [residue, residue code, x, y, z, phi, psi, surface depth, 
#          number of residues within a certain distance, num charged residues,
#          num polar residues, num non-polar residues, secondary structure, 
#          NH_O_1_relidx, NH_O_1_energy, O_NH_1_relidx, O_NH_1_energy,
#          NH_O_2_relidx, NH_O_2_energy, O_NH_2_relidx, O_NH_2_energy]

result = {}
structName = "1HMP"
fileName = "/Users/kcooo/Downloads/pdb4wxv.pdb"

def readPDBFile(structName, fileName):

    charged = {"ARG", "HIS", "LYS", "ASP", "GLU"}
    polar = {"SER", "THR", "TYR", "ASN", "GLN"}
    nonPolar = {"ALA", "VAL", "ILE", "LEU", "MET", "PHE", "PRO", "TRP", "GLY", "CYS"}

    #code from https://warwick.ac.uk/fac/sci/moac/people/students/peter_cock/python/ramachandran/calculate/#BioPython
    structure = Bio.PDB.PDBParser().get_structure(structName, fileName)
    for model in structure:
        
        surface = get_surface(model) #numpy array of all surface vertices of folded protein

        for chain in model :

            polypeptides = Bio.PDB.PPBuilder().build_peptides(chain)
            #dssp = DSSP(chain, fileName, dssp = "mkdssp")

            atoms  = Bio.PDB.Selection.unfold_entities(chain, "A")
            ns = Bio.PDB.NeighborSearch(atoms)

            for poly_index, poly in enumerate(polypeptides) :

                startIndex = poly[0].id[1]
                endIndex = poly[-1].id[1]

                phiPsi = poly.get_phi_psi_list()

                for residue in range(len(poly)):
                    #takes first data point for overlapping chains
                    if poly[residue].id[1] not in result:

                        res = poly[residue].resname
                        x,y,z = poly[residue]["CA"].coord

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

                        #return number of total atoms and grouped residues within searchRadius of resiude's alpha carbon
                        searchRadius = 5
                        closeAtoms = ns.search(poly[residue]["CA"].coord, searchRadius)
                        numCloseAtoms = len(closeAtoms)
                        residues, curCharged, curPolar, curNonPolar = set(), set(), set(), set()
                        for atom in closeAtoms:
                            currentRes = atom.get_parent()
                            curName = currentRes.resname
                            residues.add(currentRes)
                            if curName in charged:
                                curCharged.add(currentRes)
                            elif curName in polar:
                                curPolar.add(currentRes)
                            else:
                                curNonPolar.add(currentRes)

                        #if (chain, poly[residue].id) in dssp:
                        #    secondary = Bio.PDB.DSSP.ss_to_index(dssp[(chain, poly[residue].id)][2])
                        #    energyList = list(dssp[(chain, poly[residue].id)][6:])

                        result[poly[residue].id[1]] = [res, resCode, x, y, z, phi, psi, depth, 
                                                        len(residues), len(curCharged), len(curPolar), len(curNonPolar)] 
                                                        # + [secondary] + energyList
    print(result)
    return result


readPDBFile(structName, fileName)
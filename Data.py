#import module_manager
#module_manager.review()

import Bio
from Bio import PDB
from Bio.PDB.ResidueDepth import residue_depth, get_surface
from Bio.PDB.DSSP import DSSP
from Bio.PDB import Selection
from Bio.PDB.PDBList import PDBList
#outputs dictionary indexed by ("group name", index)
#from pka_processing.py import get_pka_dict

#keys = AA index
#values = [residue, residue code, x, y, z, phi, psi, surface depth, 
#          number of residues within a certain distance, num charged residues,
#          num polar residues, num non-polar residues, secondary structure, 
#          NH_O_1_relidx, NH_O_1_energy, O_NH_1_relidx, O_NH_1_energy,
#          NH_O_2_relidx, NH_O_2_energy, O_NH_2_relidx, O_NH_2_energy]

#https://www.kosbie.net/cmu/spring-15/15-112/notes/notes-functions-redux-and-web-and-file-io.html
def readFile(filename, mode="rt"):
    with open(filename, mode) as fin:
        return fin.read()

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

#readPDBFile(structName, fileName)

string = """966c    A       1.90    BS06    RS2     A       1       N180 L181 A182 V215 H218 E219 H222 H228 L235 Y237 P238 S239 Y240 T241   N73 L74 A75 V108 H111 E112 H115 H121 L128 Y130 P131 S132 Y133 T134      M236 E219;E219  M129 E112;E112  3.4.24.-        0004222,0006508,0008237,0008270,0031012         ki=23nM (RS2)   Ki=23nM (RS2)           P03956  10074939        RWEQTHLTYRIENYTPDLPRADVDHAIEKAFQLWSNVTPLTFTKVSEGQADIMISFVRGDHRDNSPFDGPGGNLAHAFQPGPGIGGDAHFDEDERWTNNFREYNLHRVAAHELGHSLGLSHSTDIGALMYPSYTFSGDVQLAQDDIDGIQAIYGRSQ"""

pdbList = PDBList()

def readAnnotations(path):
    #string = readFile(path)
    for protein in string.splitlines():
        elems = protein.split("  ")
        elems.remove("")
        elems = [i for i in elems if i != ""]
        print(elems)
        PDBId = elems[0]
        print("PDB ID", PDBId)
        PDBFile = pdbList.retrieve_pdb_file(pdb_code = PDBId, file_format = "pdb")
        features = readPDBFile(PDBId, PDBFile)
        binding = elems[7]
        catalytic = elems[9]
        resBinding = [0 for i in range(len(features))]
        resCat = [0 for i in range(len(features))]
        for res in binding:
            index = int(res.strip()[1:])
            


         
readAnnotations("HI")

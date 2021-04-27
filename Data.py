#import module_manager
#module_manager.review()

import Bio
from Bio import PDB
from Bio.PDB.ResidueDepth import residue_depth, get_surface
from Bio.PDB.DSSP import DSSP
from Bio.PDB import Selection
from Bio.PDB.PDBList import PDBList
from pka_processing.py import get_pka_dict

#keys = AA index
#values = [residue, residue code, x, y, z, phi, psi, surface depth, 
#          number of residues within a certain distance, num charged residues,
#          num polar residues, num non-polar residues, secondary structure, 
#          NH_O_1_relidx, NH_O_1_energy, O_NH_1_relidx, O_NH_1_energy,
#          NH_O_2_relidx, NH_O_2_energy, O_NH_2_relidx, O_NH_2_energy]

################################################################################
#https://www.kosbie.net/cmu/spring-15/15-112/notes/notes-functions-redux-and-web-and-file-io.html
def readFile(filename, mode="rt"):
    with open(filename, mode) as fin:
        return fin.read()

def writeFile(path, contents):
    with open(path, "wt") as f:
        f.write(contents)
###############################################################################

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

                        if (chain, poly[residue].id) in dssp:
                            secondary = (Bio.PDB.DSSP.ss_to_index(dssp[(chain, poly[residue].id)][2]))/7
                            energyList = list(dssp[(chain, poly[residue].id)][6:])

                        featuresList = [resCode, phi, psi, depth, 
                                        len(residues), len(curCharged), len(curPolar), len(curNonPolar)] + [secondary] + energyList

                        if None not in featuresList: #removes residues for which phi/psi cannot be calculated
                            result[poly[residue].id[1]] = featuresList
    return result

pdbList = PDBList()

#returns dict mapping PDB ID to (input array, binding output list, catalytic output list)
def readAnnotations(path):
    result = dict()

    string = readFile(path)
    for protein in string.splitlines():
        
        #parsing ANNOTATIONS file table
        elems = protein.split("  ")
        elems.remove("")
        elems = [i for i in elems if i != ""]
        PDBId = elems[0]
        print("extracting features for PDB ID: ", PDBId)

        #retrieve PDB file and create features list
        PDBFile = pdbList.retrieve_pdb_file(pdb_code = PDBId, file_format = "pdb")
        features = readPDBFile(PDBId, PDBFile)

        #integrate pKa 
        pkaDict = get_pka_dict(PDBFile)

        for key in features:
            features[key].append(99.999)

        for (group, res) in pKaDict:
            pka = pkaDict[(group,res)]
            if pka != None and res in features:
                features[res][-1] = pka/99.999

        lenProtein = len(features)
        resOrder = sorted(features.keys())
        inputList = []
        for key in resOrder:
            inputList.append(features[key])

        #creates 1D output lists for binding and catalytic sites
        binding = elems[7]
        catalytic = elems[9]
        resBinding = [0 for i in range(len(features))]
        resCat = [0 for i in range(len(features))]

        for site in binding.split(";"):
            for res in site.split():
                resId = int(res.strip()[1:])
                resBinding[resOrder.index(resId)] = 1
        for site in catalytic.split(";"):
            for res in site.split():
                resId = int(res.strip()[1:])
                resCat[resOrder.index(resId)] = 1

        result[PDBId] = (inputList, resBinding, resCat)
        writeToFile(PDBId, inputList, resBinding, resCat)

    normalize(result)

    return result

def normalize(proteinDict):
    maxSurfaceDepth = 0
    maxTotal = 0
    maxPolar = 0
    maxCharged = 0
    maxNonPolar = 0
    for protein in proteinDict:
        for residue in proteinDict[protein][0]:

            surfaceDepth = residue[3]
            if surfaceDepth > maxSurfaceDepth:
                maxSurfaceDepth = surfaceDepth

            numAtoms = residue[4:8]
            if numAtoms[0] > maxTotal:
                maxTotal = numAtoms[0]
            if numAtoms[1] > maxCharged:
                maxCharged = numAtoms[1]
            if numAtoms[2] > maxPolar:
                maxPolar = numAtoms[2]
            if numAtoms[3] > maxNonPolar:
                maxNonPolar = numAtoms[3] 
    for protein in proteinDict:
        for residue in proteinDict[protein][0]:
            residue[3] /= maxSurfaceDepth
            residue[4] /= maxTotal
            residue[5] /= maxCharged
            residue[6] /= maxPolar
            residue[7] /= maxNonPolar

def writeToFile(id, inputList, output1, output2):
    string = ""
    string += "---> input" + "\n"
    for line in inputList:
        for item in line:
            string += str(item) + ","
        string = string[:-1] + "\n"
    string += "---> binding residues" + "\n"
    for line in output1:
        string += str(line) + "\n"
    string += "---> catalytic residues" +"\n"
    for line in output2:
        string += str(line) + "\n"
    path = "TrainingData" + "\\" + str(id) + ".txt"
    writeFile(path, string)

def readFromFile(path):
    data = readFile(path)
    writingInput = False
    writingOutput1 = False
    writingOutput2 = False
    
    inputList, output1, output2 = [], [], []
    for line in data.splitlines():

        if "---> input" in line:
            writingInput = True
        elif "---> binding residues" in line:
            writingOutput1 = True
            writingInput = False
        elif "---> catalytic residues" in line:
            writingOutput2 = True
            writingOutput1 = False
        
        else:

            if writingInput:
                row = []
                for item in line.split(","):
                    row.append(float(item))
                inputList.append(row)
            elif writingOutput1:
                for item in line.split(","):
                    output1.append(int(item))
            else:
                for item in line.split(","):
                    output2.append(int(item))

    return (inputList, output1, output2)

readFromFile("TrainingData\966c.txt")

readAnnotations("HI")

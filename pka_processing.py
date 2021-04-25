from propka import run

###Note, modifications were made to the following PROPKA files:
    #molecular_container.py - write_pka - now returns value of internal call
        #to write_pka.
    #output.py - Commented out open/write/close on file, return the determinant
        #section.

###Feature Processing: pKa determination

#Returns a dictionary of the important pkas in the protein in pdb_fp. 
#The result is indexed into by the group and then by the residue number.
def get_pka_dict(pdb_fp):
    protein = run.single(pdb_fp, write_pka = False)
    pka_string = protein.write_pka()
    
    pka_string = pka_string.strip()
    
    #Split into entries.
    pka_entries = pka_string.split("\n")
    
    #Remove headers
    pka_entries = pka_entries[3:]
    
    #Remove extraneous whitespace and split into columns
    for i in range(len(pka_entries)):
        entry = pka_entries[i]
        entry = entry.strip()
        entry = entry.replace("\t", " ")
        while "  " in entry:
            entry = entry.replace("  ", " ")
        pka_entries[i] = entry.split(" ")
    
    #Filter out terminal pka's
    def entry_not_terminal(entry):
        if entry[0] == "N+":
            return False
        if entry[0] == "C-":
            return False
        return True
    
    pka_entries = list(filter(entry_not_terminal, pka_entries))
    pka_dict = dict()
    for entry in pka_entries:
        residue_position = int(entry[1])
        group = entry[2]
        pka = entry[3]
        if group not in pka_dict:
            pka_dict[group] = dict()
        if residue_position in pka_dict[group]:
            assert(False)
        pka_dict[group][residue_position] = float(pka)
    print(pka_dict)
    return pka_dict
    
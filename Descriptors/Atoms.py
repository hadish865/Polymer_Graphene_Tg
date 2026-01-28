class Atoms:
    def __init__(self):
        self.atomic_number = -1 
        self.bonded_atoms    = []  # this are only the atoms that are bonded and come after this atom in the formula
        self.previous_bonded_atoms = [] # this are only the atoms that are bonded and come before this atom in the formula
        self.index 		   = 0
        self.ring = False
        self.aromatic_ring = False
        self.ionic = False
        self.charge = 0

    def Identify_atom(self,index,structure):
        # this function returns the atomic number of atoms and for numbers is return 0
        self.index = index
        if structure[index].lower() == 'c': 
            if structure[index+1].lower() == 'l': return 17
            else: return 6
        if structure[index].lower() == 'n':
            if structure[index].lower() == 'a': return 11
            else: return 7
        if structure[index].lower() == 'o': return 8
        if structure[index].lower() == 'f': return 9
        if structure[index].lower() == 's': 
            if structure[index+1].lower() == 'i': return 14
            else: return 16
        if structure[index].isnumeric(): return 0
        return -1

    def Identify_end_atom(self,index,structure):
    # this function returns the atomic number of atoms 
        self.index = index
        if structure[index].lower() == 'c': return 6 
        if structure[index].lower() == 'n': return 7 
        if structure[index].lower() == 'o': return 8
        if structure[index].lower() == 'f': return 9
        if structure[index].lower() == 's': return 16 
        if structure[index].isnumeric(): return 0
        return -1

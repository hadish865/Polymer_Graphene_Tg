import re
from math import pow

class Utilities:
    def __init__(self):
        self.structure = ""
        self.atoms = None
        self.ring_indexes = None
        self.bridges = None
        self.functional_groups = {
                            "COOH": {
                                "donors": ["OH"],
                                "acceptors": ["O=C--", "O--"]
                            },
                            "COOC": {
                                "donors": [],
                                "acceptors": ["O=C--", "O--"]
                            },
                            "OH": {
                                "donors": ["OH"],
                                "acceptors": ["O--"]
                            },
                            "--NH": {
                                "donors": ["NH"],
                                "acceptors": ["N---"]
                            },
                            "=NH": {
                                "donors": ["NH"],
                                "acceptors": ["N=-"]
                            },
                            "-NH2": {
                                "donors": ["NH","NH"],
                                "acceptors": ["N---"]
                            },
                            "N": {
                                "donors": [],
                                "acceptors":["N---"]
                            },
                            "=N": {
                                "donors": [],
                                "acceptors": ["N=-"]
                            },
                            "SH": {
                                "donors": ["SH"],
                                "acceptors": ["S--"]
                            },
                            "CO": {
                                "donors": [],
                                "acceptors": ["O=C--"]
                            },
                            "Cl": {
                                "donors": [],
                                "acceptors": ["Cl"]
                            },
                            "F": {
                                "donors": ["FH"],
                                "acceptors": ["F-"]
                            },
                            "O": {
                                "donors": [],
                                "acceptors": ["O--"]
                            }
                        }
        self.H_bond_energy = {
                    ("OH","O=C--"): -4.66,
                    ("O=C--","OH"): -4.66,
                    ("OH","O--"): -5.02,
                    ("O--","OH"): -5.02,
                    ("OH","N---"): -7.65,
                    ("N---","OH"): -7.65,
                    ("N---","NH"): 6.505,
                    ("NH","N---"): 6.505,
                    ("OH","N=-"): -6.45,
                    ("N=-","OH"): -6.45,
                    ("NH","O=C--"): -3.47,
                    ("O=C--","NH"): -3.47,
                    ("NH","O--"): -2.99,
                    ("O--","NH"): -2.99,
                    ("NH","N=-"): -3.59,
                    ("N=-","NH"): -3.59,
                    ("OH","F-"): -3.70,
                    ("F-","OH"): -3.70,
                    ("OH","Cl-"): -2.2,
                    ("Cl-","OH"): -2.2,
                    ("OH","S--"): -4.18,
                    ("S--","OH"): -4.18,
                    ("NH","S--"): -3.59,
                    ("S--","NH"): -3.59,
                    ("SH","N=-"): -2.39,
                    ("N=-","SH"): -2.39,
                    ("NH","F-"): -3.7,
                    ("F-","NH"): -3.7,
                    ("SH","O--"):-2,
                    ("SH","O=C--"): -2,
                    ("O--","SH"): -2,
                    ("O=C--","SH"): -2,
                    ("FH","F-"): -6.8,
                    ("F-","FH"): -6.8
                }
        self.electronegativity = {  #based on pauling scale
                            1 : 2.20,    # H
                            5 : 2.04,   # B
                            6 : 2.55,   # C
                            7 : 3.04,   # N
                            8 : 3.44,   # O
                            9 : 3.98,   # F
                            14: 1.90,  # Si
                            15: 2.19,  # P
                            16: 2.58,  # S
                            17: 3.16,  # Cl
                            35: 2.96,  # Br
                            53: 2.66   # I
                        }

    def Create_bond_table(self):
        file = open("BondTable.txt")
        lines = file.readlines()
        bonds = [re.split(",|\n",line) for line in lines if (not line[0] == '#') and (not line == '\n')]
        self.bond_table = {(int(bond[0]),int(bond[1]),int(bond[2])): float(bond[3]) for bond in bonds}

    def Calc_bond_length(self,atom1,atom2,bond_order): # this function returns -1 if it doesn't find the specified bond length
        atomic_number1 = atom1.atomic_number if (not atom1.aromatic_ring) else atom1.atomic_number * 11
        atomic_number2 = atom2.atomic_number if (not atom2.aromatic_ring) else atom2.atomic_number * 11
        try:
            return self.bond_table[(min(atomic_number1,atomic_number2),max(atomic_number1,atomic_number2),bond_order)]
        except KeyError:
            print("please add bond length between these two atoms:" )
            print(f"atom1 atomic number\t{atomic_number1}\tindex\t{atom1.index}")
            print(f"atom2 atomic number\t{atomic_number2}\tindex\t{atom2.index}")
            print(f"bond order:\t{bond_order}")
            return -1

    def Set_structure(self,structure):
        self.structure = structure

    def Set_atoms(self,atoms):
        self.atoms = atoms

    def Set_ring_indexes(self,ring_indexes):
        self.ring_indexes = ring_indexes

    def Set_bridges(self,bridges):
        self.bridges = bridges

    def Sync(self,structure,atoms,ring_indexes,main_chain_indexes = None,bridges=None,bond_lengths=None,previous_bond_lengths=None,):
        self.Set_structure(structure)
        self.Set_atoms(atoms)
        self.Set_ring_indexes(ring_indexes)
        if not (main_chain_indexes == None): self.main_chain_indexes = main_chain_indexes
        if not (bridges == None): self.bridges = bridges
        if not (bond_lengths == None): self.bond_lengths = bond_lengths
        if not (previous_bond_lengths == None): self.previous_bond_lengths = previous_bond_lengths

    def Print_atom(self,index=0,atoms=None,structure="",end_='\n'): # if you want to use this continously please use the Set_structure method before this method to set a value for the self.structure variable which can make this method more efficient
        symbols = ['*','=','#','(',')']
        if (not self.structure == "") and (not self.atoms == None):
            if (index == len(self.structure)-1):
                print(f"atom:  {self.structure[index]}  atomic number:  {self.atoms[index].atomic_number}",end=end_)
            elif (self.atoms[index+1].atomic_number == -1) and (self.structure[index+1] not in symbols):
                print(f"atom:  {self.structure[index]}{self.structure[index+1]}  atomic number: {self.atoms[index].atomic_number}",end=end_)
            else:
                print(f"atom:  {self.structure[index]}  atomic number:  {self.atoms[index].atomic_number}",end=end_)
        else:
            try:
                if (atoms[index+1].atomic_number == -1) and (structure[index+1] not in symbols):
                    print(f"atom:  {structure[index]}{structure[index+1]} atomic number:  {atoms[index].atomic_number}",end=end_)
                else:
                    print(f"atom:  {structure[index]}  atomic number:  {atoms[index].atomic_number}",end=end_)
            except:
                CRED = "\033[91m"
                CEND = "\033[0m"
                print(CRED+"please either pass a structure to the method parameters or use the Set_structure method before using this structure"+CEND)

    
    def extract_donors_acceptors(self, input_groups, functional_groups):
        donors = []
        acceptors = []

        for group, count in input_groups.items():
            if group not in self.functional_groups:
                continue

            for donor in self.functional_groups[group]["donors"]:
                donors.extend([donor] * count)

            for acceptor in self.functional_groups[group]["acceptors"]:
                acceptors.extend([acceptor] * count)

        return donors, acceptors

    def Calc_hydrogen_bond_energy(self, groups1, groups2):
        donors1, acceptors1 = self.extract_donors_acceptors(groups1, self.functional_groups)
        donors2, acceptors2 = self.extract_donors_acceptors(groups2, self.functional_groups)

        total_energy = 0.0

        for donor in donors1:
            for acceptor in acceptors2:
                energy = self.H_bond_energy.get((donor, acceptor))
                if energy:
                    total_energy += energy

        for donor in donors2:
            for acceptor in acceptors1:
                energy = self.H_bond_energy.get((donor, acceptor))
                if energy:
                    total_energy += energy
        return total_energy
    
    def compute_dipoles(self,atoms):
        dipoles = []

        for atom_dict in atoms:
            for atom_number, bonds in atom_dict.items():
                for bonded_atom, bond_length in bonds:
                    en1 = self.electronegativity[atom_number]
                    en2 = self.electronegativity[bonded_atom]
                    diff = abs(en1 - en2)

                    if diff >= 0.45:
                        dipole = round(diff * bond_length ,4)
                        dipoles.append(dipole)
        return dipoles

    def Calc_polar_interaction(self, atoms1, atoms2):
        dipoles1 = self.compute_dipoles(atoms1)
        dipoles2 = self.compute_dipoles(atoms2)
        interaction_total = 0.0

        for d1 in dipoles1:
            for d2 in dipoles2:
                interaction_total += pow(d1, 2) * pow(d2, 2)

        return interaction_total

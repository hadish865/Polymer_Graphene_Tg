import xlsxwriter 
from Utilities import Utilities
from Atoms import Atoms

# Fajans' Rules for bond polarity: bond dipole = (amount of bar seperation) * distance 
# amount of bar seperation is approximated using the difference between electronegativities of atoms
# Linus Pauling derivation of electronegativities
# F: 4.0, O: 3.5, N: 3.0, C: 2.5, CL: 3.0, S: 2.5, Si = 1.8
# bond lengths
# benzene rings: 1.4 A, C=O: 1.23 A, C-C: 1.54
# C-F: 1.35 A, C-Cl: 1.78 A, C-O: 1.43 A, N-C: 1.47, Si-C: 1.86, S=O: 1.42
# C Atomic radius: 0.7, N Atomic radius: 0.65, O Atomic radius: 0.6, 
# F Atomic radius: 0.5, Si Atomic radius: 1.1, Cl Atomic radius: 1, S Atomic radius: 1 

STRUCTURE_LINKAGES = [] # where two structures link toghether, every entry is a list that contains the first and second atoms that link toghether
class Structure_Processor:
    def __init__(self,smiles_structure,GO_flag=False,RGO_flag=False):
        self.smiles_structure = smiles_structure
        self.smiles_structures = [[]  for _ in range(self.smiles_structure.count('.')+1)]
        self.structures = [[]  for _ in range(len(self.smiles_structures))]  
        self.structures_atoms_list = [[] for _ in range(len(self.smiles_structure))]
        self.len_structure = len(self.smiles_structure)
        self.descriptors = []
        self.GO_flag = GO_flag
        self.RGO_flag = RGO_flag

    def Pre_process_structures(self):
        structure_index = 0
        bracket_counter = 0
        index = 0
        while index < self.len_structure:
            s = self.smiles_structure[index]
            if (s == '.'):
                STRUCTURE_LINKAGES.append([index-1,index+1])
                structure_index += 1
                index += 1
                continue
            elif (s == '['):
                bracket_counter += 1
                index += 1
                self.smiles_structures[structure_index].append(s)
                continue
            elif (s == ']'):
                bracket_counter -= 1
                self.smiles_structures[structure_index].append(s)
                if (bracket_counter == 0) and (not structure_index == 0):
                    structure_index -= 1
                index += 1
                continue
            # create the ionic charges
            elif (s == '{'):
                bar = 0
                increment = 0
                if (self.smiles_structure[index+2] == '}'):
                    increment = 3
                    bar = 1
                else: 
                    bar = int(self.smiles_structure[index+2])
                    increment = 4
                if (self.smiles_structure[index+1] == '-'): bar = -1 * bar
                self.structures_atoms_list[structure_index][-1].ionic = True
                self.structures_atoms_list[structure_index][-1].charge = bar
                index += increment 
                continue
            self.structures[structure_index].append(s)
            self.smiles_structures[structure_index].append(s)
            self.structures_atoms_list[structure_index].append(Atoms())
            index += 1

    def Create_descriptors(self):
        for structure_index,structure_list in enumerate(self.structures):
            structure = "".join(structure_list)
            smiles_structure = "".join(s for i,s in enumerate(self.smiles_structures[structure_index]) if not ((i == 0 or i == len(self.smiles_structures[structure_index])-1) and (s == '[' or s == ']')))
            self.descriptors.append(Descriptors(structure,self.structures_atoms_list[structure_index],smiles_structure,GO_flag=self.GO_flag,RGO_flag=self.RGO_flag))

class Descriptors:
    def __init__(self,structure,atoms,smiles_structure,GO_flag=False,RGO_flag=False):
        self.structure = structure 
        self.smiles_structure = smiles_structure
        self.TWO_WORD_ATOMS = [14,17] # this list contains atomic number of the atoms that have two word names
        self.len_structure = len(self.structure)	
        self.bridges  = []
        self.Create_bridges()
        self.atoms : list[Atoms] = atoms
        self.Complete_atoms()
        self.bond_lengths = [{} for _ in range(self.len_structure)] # self.bond_lengths[i][j]: the i index corresponds to atom i and j index corresponds to an atom with index j that is bonded to an atom with index i 
        self.previous_bond_lengths = [{} for _ in range(self.len_structure)] # the same as bond_lengths but for previous bonds
        self.utilities = Utilities()
        self.utilities.Create_bond_table()
        self.Create_bonds()
        self.Create_previous_bonds()
        self.functional_groups = {"COOH":0, "COOC": 0, "CO" : 0, "OH": 0, "O": 0, "--NH": 0, "=NH": 0, "-NH2": 0, "N": 0, "=N":0, "SH":0, "Cl":0, "F":0}
        self.molecule = False
        self.one_side_connected = False # for structures that might be connected from only one end of the structure
        self.main_chain_indexes = []
        self.Set_main_chain_indexes()
        self.ring_indexes = []
        self.Create_rings()
        self.number_of_inner_benzene_rings, self.number_of_side_benzene_rings = self.Set_benzene_rings()
        self.Create_functional_groups()
        #self.Filter_functional_groups()
        self.GO = GO_flag
        self.RGO = RGO_flag

    def Calc_polymer_descriptors(self):
        self.number_of_atoms = self.Calc_number_of_atoms()
        self.flexibility     = self.Calc_flexibility()
        self.occupied_volume = self.Calc_occupied_volume()
        self.important_functional_groups = {"CH": 0, "C=OH":0, "C=OOH":0, "OH":0,"NH":0, "NO2":0, "B":0,"BCH": 0, "BC=OH":0, "BC=OOH":0, "BOH":0,"BNH":0, "BNO2":0} 
        self.number_of_important_functional_groups = {"CH": 0, "C=OH":0, "C=OOH":0, "OH":0,"NH":0, "NO2":0, "B":0,"BCH": 0, "BC=OH":0, "BC=OOH":0, "BOH":0,"BNH":0, "BNO2":0} 
        self.G_important_functional_groups_energy_value = {"CH": 17.2, "C=OH":24.5, "C=OOH":22.7, "OH":18.8,"NH":27.7, "NO2":27.9, "B":16.5,"BCH": 17.2, "BC=OH":25.5, "BC=OOH":27.1, "BOH":23.4,"BNH":20.8, "BNO2":26.9} 
        self.GO_important_functional_groups_energy_value = {"CH": 17.2, "C=OH":20.7, "C=OOH":34.1, "OH":25.3,"NH":27.7, "NO2":27.9, "B":14.9,"BCH": 17.2, "BC=OH":16.4, "BC=OOH":22.7, "BOH":27.3,"BNH":23.3, "BNO2":26.9} 
        self.Calc_important_groups_energy_value() # this is for the graphene project only

    def Create_bridges(self):
        ib = [] # incomplete bridges
        cb = [] # completed bridges
        for i in range(self.len_structure):
            if self.structure[i] == '(': ib.append(i)
            elif self.structure[i] == ')': 
                beginning = ib.pop()
                cb.append([beginning,i])
        len_bridges = len(cb)
        for i in range(len_bridges):
            for j in range(i+1, len_bridges):
                if cb[i][1] == cb[j][0] - 1: cb[i].append(j)
        self.bridges = cb

    def Find_connected_bridges(self,index) -> list:
        # this function only returns the bridges connected to an index
        if (self.atoms[index].atomic_number in self.TWO_WORD_ATOMS): 
            index += 1
        len_bridges = len(self.bridges)
        connected_bridges = []
        for i in range(len_bridges):
            if (index == self.bridges[i][0] - 1): connected_bridges.append(self.bridges[i]) 
        return connected_bridges

    def Find_end_of_connecting_bridges(self, connecting_bridges:list) -> int:
        # after finding the connecting bridges this finds the end of the bridges
        end = connecting_bridges[0][1]
        for bridge in connecting_bridges:
            if bridge[1] > end: end = bridge[1]
            if len(bridge) == 3: connecting_bridges.append(self.bridges[bridge[2]])
        return end

    def Complete_atoms(self):
        index = 0
        while index < self.len_structure-1:
            self.atoms[index].atomic_number = self.atoms[index].Identify_atom(index,self.structure)
            self.atoms[index].index = index
            if (self.structure[index].islower()): self.atoms[index].aromatic_ring = True
            index += 1		
        self.atoms[index].atomic_number = self.atoms[index].Identify_end_atom(index,self.structure)
        self.atoms[index].index = index
        if (self.structure[index].islower()): self.atoms[index].aromatic_ring = True

    def Find_index_connected_before_number(self,index): # this function was orginaly was used for finding which atom connects to a specific number index (before the number) but since has been modified to work with symbols as well
        while(True):
            if not (self.atoms[index].atomic_number == -1 or self.atoms[index].atomic_number == 0): return index
            else: index -= 1

    def Find_index_connected_after_number(self,index): # this is used for finding which atom connects to a specific index (after the number or for atoms that two word names like Si)
        while(index < self.len_structure-1):
            if not (self.atoms[index].atomic_number == -1 or self.atoms[index].atomic_number == 0): return index
            else: index += 1
        return index

    def Set_bond(self,index1,index2):
        if self.structure[index2] == '=': 
            if self.atoms[index2+1].atomic_number == 0: 
                self.Set_bond(index1,index2+2)
            else:
                self.atoms[index1].bonded_atoms.append((index2+1,2))
                self.bond_lengths[index1].update({index2+1:self.utilities.Calc_bond_length(self.atoms[index1],self.atoms[index2+1],2)})
        elif self.structure[index2] == '#':
            self.atoms[index1].bonded_atoms.append((index2+1,3))
            self.bond_lengths[index1].update({index2+1:self.utilities.Calc_bond_length(self.atoms[index1],self.atoms[index2+1],3)})
        elif not (self.atoms[index2].atomic_number == -1 or self.atoms[index2].atomic_number == 0): # if its an atom
            self.atoms[index1].bonded_atoms.append((index2,1))
            self.bond_lengths[index1].update({index2:self.utilities.Calc_bond_length(self.atoms[index1],self.atoms[index2],1)})
        elif self.structure[index2] == '(':
            connecting_bridges = self.Find_connected_bridges(index1)
            for bridge in connecting_bridges:
                begining_bridge_index = bridge[0]
                self.Set_bond(index1,begining_bridge_index+1)
                if len(bridge) == 3:
                    connecting_bridges.append(self.bridges[bridge[2]])
            index_of_end_bridge = self.Find_end_of_connecting_bridges(connecting_bridges)
            if not (index_of_end_bridge == self.len_structure-1 or index_of_end_bridge == 0 or index_of_end_bridge == -1):
                self.Set_bond(index1,index_of_end_bridge+1)

    def Create_bonds(self):  
        for index,atom in enumerate(self.atoms[0:self.len_structure-1]):
            if not (atom.atomic_number == 0 or atom.atomic_number == -1  or atom.atomic_number in self.TWO_WORD_ATOMS): 
                self.Set_bond(index,index+1)
            # for atoms that have either an number after them or a two word name
            # if it's Si (for now it's only used for Si but may support other atoms that have two word names, for example Cl has two words but doesn't have any bonds after it so we dont even have to care for it) if you need to add any new atoms just add them to the list below and it will handle it accordingly

            elif (atom.atomic_number in self.TWO_WORD_ATOMS) and not (index + 1 == self.len_structure-1): 
                self.Set_bond(index,index+2)
            # for atoms that have an number after them we need to do add something
            try:
                if (self.atoms[index+1].atomic_number == 0) and not(self.structure[index+2] == ')' or self.atoms[index].atomic_number == -1):
                    index_after_number = self.Find_index_connected_after_number(index+1)
                    if not (index_after_number  == self.len_structure -1):
                        self.Set_bond(index,index_after_number)
            except:
                pass
        if not len([1 for bl in self.bond_lengths if (-1 in bl.values())]) == 0: 
            CRED = "\033[91m"
            CEND = "\033[0m"
            print(CRED+"some bond lengths could not be found, please add them to the BondTable.txt file"+CEND)
            print("program will exit for now")
            exit() 
                        
    def Create_previous_bonds(self):
        index = 0
        while index < self.len_structure:
            if not len(self.atoms[index].bonded_atoms) == 0:
                for _,bond in enumerate(self.atoms[index].bonded_atoms):
                    self.atoms[bond[0]].previous_bonded_atoms.append([index,bond[1]])
                    self.previous_bond_lengths[bond[0]].update({index:self.bond_lengths[index][bond[0]]})
            index += 1

    def Is_it_inside_bridge(self,index):
        # cheks to see if the index is inside a bridge or not
        for bridge in self.bridges:
            if (index > bridge[0]) and (index < bridge[1]): return True 
        return False

    def Find_index_connected_before_bridge(self,index) -> int: # index is the index of the bridge end =  ")" or any atom that is inside a bridge, if an atom is not inside a bridge it will return the index again
        for bridge in self.bridges:
                if (index <= bridge[1]) and (index >= bridge[0]):
                    return self.Find_index_connected_before_bridge(bridge[0]-1)

        if self.atoms[index].atomic_number == 0:
            index = self.Find_index_connected_before_number(index)
        return index

    def Set_main_chain_indexes(self): # if the structure doesn't have *s the function will presume that the first index is the begining of the molecules main chain axis and the last atom that is not inside a branch (bridge) is the last atom in the main chain axis, and every atom that is not inside a branch will be considered inside the main chain
        begin_index  = self.structure.find("*") 
        end_index    = self.structure.rfind("*",begin_index,self.len_structure)
        method2 = self.Is_it_inside_bridge(end_index)
        if (begin_index == -1):
            CRED = "\033[91m"
            CEND = "\033[0m"
            print(CRED+"could not find main chain, if you have entered a polymer chain please insert '*'s at the begining and end of you're chain, if not then ignore this message"+CEND)
            self.molecule = True
            begin_index = 0
            end_index = self.Find_index_connected_before_number(self.len_structure-1)
            method2 = False
        elif (self.structure.count('*') == 1):
            CRED = "\033[91m"
            CEND = "\033[0m"
            print(CRED+"chain was only connected at one point"+CEND)
            end_index = self.Find_index_connected_before_number(self.len_structure-1)
            self.one_side_connected = True
            method2 = False
        # method 1: check to see if the *s are at the two ends of the SMILES reoresentation 
        if not (method2):
            index = begin_index
            while index <= end_index:
                if self.structure[index] == '(': index = self.Find_end_of_connecting_bridges(self.Find_connected_bridges(index-1)) + 1
                else:
                    if not (self.structure[index].isnumeric() or self.atoms[index].atomic_number == -1): self.main_chain_indexes.append(index)
                    index += 1
        # method 2: check the paths till you find the * 
        else:
            if (self.structure[end_index-1] == ')'): # end_index is the index of the *
                end_index = self.Find_index_connected_before_bridge(end_index-1)
            elif self.atoms[end_index-1].atomic_number == 0:
                end_index = self.Find_index_connected_before_number(end_index)
            else: 
                end_index -= 1

            self.main_chain_indexes.append(end_index) 
            current_index = self.atoms[end_index].previous_bonded_atoms[0][0] # move one atom back
            while (current_index > begin_index + 1): # begin_index+1 is the first atom (begin_index is the * index)
                self.main_chain_indexes.append(current_index)
                current_index = self.atoms[current_index].previous_bonded_atoms[0][0]
            self.main_chain_indexes.append(current_index)
            self.main_chain_indexes.reverse()
        #if (method2): # for checking to see which method found the main chain indexes but it's not necessary
        #    print("main chain found using method 2:")
        #else:
        #    print("main chain found using method 1:")

    def Create_rings(self):
        self.ring_starts = []
        self.ring_ends = []
        i = 1
        while True: # search for ring indexes
            index = self.structure.find(str(i))
            if not index == -1: 
                begin_index  = self.Find_index_connected_before_number(index) 
                end_index    = self.Find_index_connected_before_number(self.structure.rfind(str(i),begin_index,self.len_structure))
                bond_order = 2 if self.structure[begin_index+1] == '=' else 1 
                bond_length = self.utilities.Calc_bond_length(self.atoms[begin_index],self.atoms[end_index],bond_order)
                if [begin_index,1] not in self.atoms[end_index].bonded_atoms:
                    self.atoms[end_index].bonded_atoms.append([begin_index,bond_order])
                    self.bond_lengths[end_index].update({begin_index:bond_length})
                if [end_index,1] not in self.atoms[begin_index].previous_bonded_atoms:
                    self.atoms[begin_index].previous_bonded_atoms.append([end_index,bond_order])
                    self.previous_bond_lengths[begin_index].update({end_index:bond_length})
                self.ring_starts.append(begin_index)
                self.ring_ends.append(end_index)
                i += 1
            else: break

        for j in range(1,i):
            self.Set_ring(str(j))

    def Set_ring(self,ring_number):
        begin_index  = self.Find_index_connected_before_number(self.structure.find(ring_number)) 
        end_index    = self.Find_index_connected_before_number(self.structure.rfind(ring_number,begin_index,self.len_structure))
        ring_indexes = []
        current_index = end_index
        while (current_index > begin_index + 1): # begin_index+1 is the first atom (begin_index is the * index)
            ring_indexes.append(current_index)
            if (current_index in self.ring_ends) and (not current_index == end_index):
                current_index = self.atoms[current_index].bonded_atoms[-1][0] 
            else:
                current_index = self.atoms[current_index].previous_bonded_atoms[0][0]
        ring_indexes.append(current_index)
        ring_indexes.reverse()
        self.ring_indexes.append(ring_indexes)
        for index in ring_indexes: self.atoms[index].ring = True
        #        bond_counter = 5 # the number of bonds - 1
        #        ring_indexes_dict = {begin_index:self.atoms[begin_index].bonded_atoms.copy()} 
        #        checked_bonds = []
        #        ring_indexes = [begin_index]
        #        flag = False
        #        while not len(ring_indexes_dict) == 0: 
        #            bond = ring_indexes_dict[ring_indexes[-1]][0]
        #            if   (bond[0] == end_index and bond_counter == 1): flag = True
        #            elif (bond[0] == end_index) and (not bond_counter == 1): 
        #                ring_indexes_dict[ring_indexes[-1]].remove(bond)
        #            elif (not bond[0] == end_index) and (bond_counter == 1):
        #                ring_indexes_dict[ring_indexes[-1]].remove(bond)
        #            elif (not bond[0] == end_index) and (not bond_counter == 1): # moving to another node if that node is not already checked
        #                change_node = False
        #                if bond not in checked_bonds:
        #                    checked_bonds.append(bond)
        #                    ring_indexes.append(bond[0])
        #                    ring_indexes_dict.update({bond[0]:self.atoms[bond[0]].bonded_atoms.copy()})
        #                    bond_counter -= 1
        #                    change_node = True
        #                if not (change_node):
        #                    ring_indexes_dict[ring_indexes[-1]].remove(bond)
        #            if len(ring_indexes_dict[ring_indexes[-1]]) == 0:
        #                ring_indexes_dict.pop(ring_indexes[-1])
        #                ring_indexes.pop()
        #                bond_counter += 1
        #            if (flag): break
        #        ring_indexes_dict.update({end_index:self.atoms[end_index].bonded_atoms.copy()})
        #        ring_indexes.append(end_index)
        #        for index in ring_indexes: self.atoms[index].ring = True
        #        self.ring_indexes.append(ring_indexes)

    def is_ring_inside_main_chain(self,index): # this can be any index from the ring indexes and the function finds the ring automaticly
        for ring in self.ring_indexes:
            if index in ring:
                number_of_atoms_inside_main_chain = 0
                for i in ring:
                    if i in self.main_chain_indexes: number_of_atoms_inside_main_chain += 1
                if number_of_atoms_inside_main_chain >= 2: return True
                else: return False

    def Is_it_side_chain(self,index):
        if index in self.main_chain_indexes: return False
        # for rings, we check to see if atleast two atoms are inside the main chain or not
        if self.atoms[index].ring: return not self.is_ring_inside_main_chain(index)
        return True

    def Calc_important_groups_energy_value(self): # todo: this function does not calcaulate n=oo groups, add if necessary  
        #index = 0
        begin_index  = self.structure.find("*")+1 
        end_index    = self.structure.rfind("*",begin_index,self.len_structure)-1
        for index,atom in enumerate(self.atoms):
            number_of_bonds = 0
            number_of_bonds  = sum([i[1] for i in atom.bonded_atoms])
            number_of_bonds += sum([i[1] for i in atom.previous_bonded_atoms])
            benzene_flag = False
            for bond in atom.bonded_atoms:
                if (self.atoms[bond[0]].aromatic_ring): benzene_flag = True 
            for bond in atom.previous_bonded_atoms:
                if (self.atoms[bond[0]].aromatic_ring): benzene_flag = True 
            if (index == end_index) or (index == begin_index): number_of_bonds += 1
            if atom.atomic_number == 6 and (not atom.aromatic_ring): # carbon 
                if (number_of_bonds == 4): # because hydrogens are not counted in the structure
                    #index += 1
                    continue
                else:
                    if len([1 for bond in atom.bonded_atoms if (self.atoms[bond[0]].atomic_number == 8 and bond[1] == 2)]) > 0:# we have c=oh
                        if (self.GO):
                            if (benzene_flag):
                                self.important_functional_groups["BC=OH"] += self.GO_important_functional_groups_energy_value["BC=OH"]
                                self.number_of_important_functional_groups["BC=OH"] += 1
                            else:
                                self.important_functional_groups["C=OH"] += self.GO_important_functional_groups_energy_value["C=OH"]
                                self.number_of_important_functional_groups["C=OH"] += 1 
                        elif (self.RGO):
                            if (benzene_flag):
                                self.important_functional_groups["BC=OH"] += (self.GO_important_functional_groups_energy_value["BC=OH"] + self.G_important_functional_groups_energy_value["BC=OH"])/2 
                                self.number_of_important_functional_groups["BC=OH"] += 1 
                            else:
                                self.important_functional_groups["C=OH"] += (self.GO_important_functional_groups_energy_value["C=OH"] + self.G_important_functional_groups_energy_value["C=OH"])/2 
                                self.number_of_important_functional_groups["C=OH"] += 1 
                        else:
                            if (benzene_flag):
                                self.important_functional_groups["BC=OH"] += self.G_important_functional_groups_energy_value["BC=OH"]
                                self.number_of_important_functional_groups["BC=OH"] += 1 
                            else:
                                self.important_functional_groups["C=OH"] += self.G_important_functional_groups_energy_value["C=OH"]
                                self.number_of_important_functional_groups["C=OH"] += 1 
                    else: # we have normal ch bonds
                        if (self.GO):
                            if (benzene_flag):
                                self.important_functional_groups["BCH"] += (4-number_of_bonds)*self.GO_important_functional_groups_energy_value["BCH"]
                                self.number_of_important_functional_groups["BCH"] += 1
                            else:
                                self.important_functional_groups["CH"] += (4-number_of_bonds)*self.GO_important_functional_groups_energy_value["CH"]
                                self.number_of_important_functional_groups["CH"] += 1 
                        elif (self.RGO):
                            if (benzene_flag):
                                self.important_functional_groups["BCH"] += (4-number_of_bonds)*(self.GO_important_functional_groups_energy_value["BCH"] + self.G_important_functional_groups_energy_value["BCH"])/2 
                                self.number_of_important_functional_groups["BCH"] += 1 
                            else:
                                self.important_functional_groups["CH"] += (4-number_of_bonds)*(self.GO_important_functional_groups_energy_value["CH"] + self.G_important_functional_groups_energy_value["CH"])/2 
                                self.number_of_important_functional_groups["CH"] += 1 
                        else:
                            if (benzene_flag):
                                self.important_functional_groups["BCH"] += (4-number_of_bonds)*self.G_important_functional_groups_energy_value["BCH"]
                                self.number_of_important_functional_groups["BCH"] += 1 
                            else:
                                self.important_functional_groups["CH"] += (4-number_of_bonds)*self.G_important_functional_groups_energy_value["CH"]
                                self.number_of_important_functional_groups["CH"] += 1 
            elif atom.atomic_number == 8: # oxygen 
                if (number_of_bonds == 2):
                    #index += 1
                    continue
                else:
                    previous_bonded_atom_index = atom.previous_bonded_atoms[0][0]
                    previous_bonded_atom = self.atoms[previous_bonded_atom_index] 
                    if (previous_bonded_atom.atomic_number == 6 and len([1 for bond in previous_bonded_atom.bonded_atoms if (self.atoms[bond[0]] == 8 and bond[1] == 1)]) > 0): # we have c=ooh
                        if (self.GO):
                            if (benzene_flag):
                                self.important_functional_groups["BC=OOH"] += self.GO_important_functional_groups_energy_value["BC=OOH"]
                                self.number_of_important_functional_groups["BC=OOH"] += 1
                            else:
                                self.important_functional_groups["C=OOH"] += self.GO_important_functional_groups_energy_value["C=OH"]
                                self.number_of_important_functional_groups["C=OOH"] += 1 
                        elif (self.RGO):
                            if (benzene_flag):
                                self.important_functional_groups["BC=OOH"] += (self.GO_important_functional_groups_energy_value["BC=OOH"] + self.G_important_functional_groups_energy_value["BC=OOH"])/2 
                                self.number_of_important_functional_groups["BC=OOH"] += 1 
                            else:
                                self.important_functional_groups["C=OOH"] += (self.GO_important_functional_groups_energy_value["C=OOH"] + self.G_important_functional_groups_energy_value["C=OOH"])/2 
                                self.number_of_important_functional_groups["C=OOH"] += 1 
                        else:
                            if (benzene_flag):
                                self.important_functional_groups["BC=OOH"] += self.G_important_functional_groups_energy_value["BC=OOH"]
                                self.number_of_important_functional_groups["BC=OOH"] += 1 
                            else:
                                self.important_functional_groups["C=OOH"] += self.G_important_functional_groups_energy_value["C=OOH"]
                                self.number_of_important_functional_groups["C=OOH"] += 1 
                    else: # we have normal oh bonds
                        if (self.GO):
                            if (benzene_flag):
                                self.important_functional_groups["BOH"] += self.GO_important_functional_groups_energy_value["BOH"]
                                self.number_of_important_functional_groups["BOH"] += 1
                            else:
                                self.important_functional_groups["OH"] += self.GO_important_functional_groups_energy_value["OH"]
                                self.number_of_important_functional_groups["OH"] += 1 
                        elif (self.RGO):
                            if (benzene_flag):
                                self.important_functional_groups["BOH"] += (self.GO_important_functional_groups_energy_value["BOH"] + self.G_important_functional_groups_energy_value["BOH"])/2 
                                self.number_of_important_functional_groups["BOH"] += 1 
                            else:
                                self.important_functional_groups["OH"] += (self.GO_important_functional_groups_energy_value["OH"] + self.G_important_functional_groups_energy_value["OH"])/2 
                                self.number_of_important_functional_groups["OH"] += 1 
                        else:
                            if (benzene_flag):
                                self.important_functional_groups["BOH"] += self.G_important_functional_groups_energy_value["BOH"]
                                self.number_of_important_functional_groups["BOH"] += 1 
                            else:
                                self.important_functional_groups["OH"] += self.G_important_functional_groups_energy_value["OH"]
                                self.number_of_important_functional_groups["OH"] += 1 
            elif atom.atomic_number == 7: # nitrogen 
                if (self.GO):
                    if (benzene_flag):
                        self.important_functional_groups["BNH"] += (3-number_of_bonds)*self.GO_important_functional_groups_energy_value["BNH"]
                        self.number_of_important_functional_groups["BNH"] += 1
                    else:
                        self.important_functional_groups["NH"] += (3-number_of_bonds)*self.GO_important_functional_groups_energy_value["NH"]
                        self.number_of_important_functional_groups["NH"] += 1 
                elif (self.RGO):
                    if (benzene_flag):
                        self.important_functional_groups["BNH"] += (3-number_of_bonds)*(self.GO_important_functional_groups_energy_value["BBH"] + self.G_important_functional_groups_energy_value["BNH"])/2 
                        self.number_of_important_functional_groups["BNH"] += 1 
                    else:
                        self.important_functional_groups["NH"] += (3-number_of_bonds)*(self.GO_important_functional_groups_energy_value["NH"] + self.G_important_functional_groups_energy_value["NH"])/2 
                        self.number_of_important_functional_groups["NH"] += 1 
                else:
                    if (benzene_flag):
                        self.important_functional_groups["BNH"] += (3-number_of_bonds)*self.G_important_functional_groups_energy_value["BNH"]
                        self.number_of_important_functional_groups["BNH"] += 1 
                    else:
                        self.important_functional_groups["NH"] += (3-number_of_bonds)*self.G_important_functional_groups_energy_value["NH"]
                        self.number_of_important_functional_groups["NH"] += 1 
            #index += 1
        total_number_of_benzenes = self.number_of_inner_benzene_rings + self.number_of_side_benzene_rings 
        if (self.GO):
            self.important_functional_groups["B"] = (total_number_of_benzenes)*self.GO_important_functional_groups_energy_value["B"]
            self.number_of_important_functional_groups["B"] += total_number_of_benzenes
        elif (self.RGO):
            self.important_functional_groups["B"] = (total_number_of_benzenes)*(self.GO_important_functional_groups_energy_value["B"] + self.G_important_functional_groups_energy_value["B"])/2 
            self.number_of_important_functional_groups["B"] += total_number_of_benzenes 
        else:
            self.important_functional_groups["B"] = (total_number_of_benzenes)*self.G_important_functional_groups_energy_value["B"]
            self.number_of_important_functional_groups["B"] += total_number_of_benzenes
        # for now i'm going to add the BCH energies to CH
        self.important_functional_groups["CH"] += self.important_functional_groups["BCH"]
        self.number_of_important_functional_groups["CH"] += self.number_of_important_functional_groups["BCH"]

    def is_it_rotatable(self,index): # for now if an atom is connected to any ring it is considered not rotatable
        if (self.atoms[index].atomic_number == -1) or (self.atoms[index].atomic_number == 0) or (self.atoms[index].ring): return False

        for bonds in self.atoms[index].bonded_atoms:
            if (not bonds[1] == 1) or (self.atoms[bonds[0]].ring): return False
        for bonds in self.atoms[index].previous_bonded_atoms:
                if (not bonds[1] == 1) or (self.atoms[bonds[0]].ring): return False
        return True
            
    def Calc_number_of_atoms(self):
        number_of_atoms = 0
        for atom in self.atoms:
            if not (atom.atomic_number == -1 or atom.atomic_number == 0): number_of_atoms += 1
        return number_of_atoms

    def Set_benzene_rings(self): # this method calculates the number of benzene rings inside and outside the main chain 
        inner_benzene_rings, side_benzene_rings = (0,0)
        for ring in self.ring_indexes:
            number_of_atoms_inside = 0 # this are the number of atoms that are inside the main chain, if there are atleast two the ring is considered to be inside the main chain
            if not self.atoms[ring[0]].aromatic_ring: continue 
            for index in ring:
                if not self.Is_it_side_chain(index): number_of_atoms_inside += 1
            if number_of_atoms_inside >= 2: inner_benzene_rings += 1
            else: side_benzene_rings += 1
        return (inner_benzene_rings, side_benzene_rings) 

    def Calc_flexibility(self):
        ether_bond_coeff = 7
        normal_bond_coeff = 3
        stiff_coeff = 5
        flexibility = 0
        i = 0
        while i < self.len_structure:                    
            if (self.atoms[i].atomic_number == 0) or (self.atoms[i].atomic_number == 1) or (self.atoms[i].atomic_number == -1) or (self.Is_it_side_chain(i) or self.atoms[i].ring): i += 1
            else:
                for bonds in self.atoms[i].bonded_atoms:
                    if (not bonds[1] == 1) or (self.Is_it_side_chain(bonds[0])) or (not self.is_it_rotatable(bonds[0])): continue
                    else:
                        if (self.atoms[i].atomic_number == 8) or (self.atoms[bonds[0]].atomic_number == 8): flexibility += ether_bond_coeff
                        else: flexibility += normal_bond_coeff 
                i += 1
        # check the first and last index
        begin_index  = self.structure.find("*")+1 
        end_index    = self.structure.rfind("*",begin_index,self.len_structure)-1
        if (self.is_it_rotatable(begin_index)) and (self.is_it_rotatable(end_index)):
            if (self.atoms[begin_index].atomic_number == 8) or (self.atoms[end_index].atomic_number == 8): flexibility += ether_bond_coeff
            else: flexibility += normal_bond_coeff 
        total_number_of_benzene_rings = self.number_of_side_benzene_rings + self.number_of_inner_benzene_rings
        flexibility -= total_number_of_benzene_rings * stiff_coeff
        return flexibility/self.number_of_atoms

    def Calc_occupied_volume(self):
        # bond lengths
        # benzene rings: 1.4 a, c=o: 1.23 a, c-c: 1.54
        # c-f: 1.35 a, c-cl: 1.78 a, c-o: 1.43 a, n-c: 1.47, si-c: 1.86, s=o: 1.42, c-s: 1.6, c#n: 1.14
        # c atomic radius: 0.7, n atomic radius: 0.65, o atomic radius: 0.6, 
        # f atomic radius: 0.5, si atomic radius: 1.1 ,cl atomic radius: 1, s atomic radius: 1
        # for benzene we calculate the diameter of the structure plus 2*c atomic radius: 2 * 1.4 + 2 * 0.7 = 4.2
        # for c-c:  1.54 + 0.7 + 0.7  = 2.97
        # for c=c:  1.31 + 0.7 + 0.7  = 2.71
        # for c-n:  1.47 + 0.7 + 0.65 = 2.82
        # for c#n:  1.14 + 0.7 + 0.65 = 2.49
        # for c-o:  1.43 + 0.7 + 0.6  = 2.73
        # for c=o:  1.23 + 0.7 + 0.6  = 2.53
        # for c-f:  1.35 + 0.7 + 0.5  = 2.55
        # for c-si: 1.86 + 0.7 + 1.1  = 3.66
        # for c-s:  1.6  + 0.7 + 1.0  = 3.30
        # for c-cl: 1.78 + 0.7 + 1.0  = 3.48
        # for s=o:  1.42 + 1   + 0.6  = 3.02
        def indexes_part_of_the_same_ring(index1,index2):
            for ring in self.ring_indexes:
                if (index1 in ring) and (index2 in ring): return True
            return False

        occupied_volume = 0
        index = 0
        for index,atom in enumerate(self.atoms):
            if (not atom.atomic_number == -1) and (not atom.atomic_number == 0) and (self.Is_it_side_chain(index)):
                for bond in atom.previous_bonded_atoms:
                    if (not indexes_part_of_the_same_ring(index,bond[0])) and (not self.Is_it_side_chain(bond[0])):
                        occupied_volume += self.previous_bond_lengths[index][bond[0]]
                for bond in self.atoms[index].bonded_atoms:
                    if indexes_part_of_the_same_ring(index,bond[0]): continue
                    occupied_volume += self.bond_lengths[index][bond[0]]
        occupied_volume += self.number_of_side_benzene_rings * 4.2
        return occupied_volume

    def Calc_shortest_distance_from_main_chain(self,index): # make sure index does not have the atomic number of -1 or 0
        # this function returns the distance and also the main chain index that it lands on
        distance = 0
        current_index = index
        while current_index not in self.main_chain_indexes:
            if (not self.atoms[current_index].ring) or (self.atoms[current_index].ring and current_index in self.ring_starts):
                next_index = self.atoms[current_index].previous_bonded_atoms[0][0]
                distance += self.previous_bond_lengths[current_index][next_index]
                current_index = next_index
            elif current_index in self.ring_ends: # were at the end of the ring and if we move one step forward we get to the begining of the ring
                next_index = self.atoms[current_index].bonded_atoms[-1][0]
                distance += self.bond_lengths[current_index][next_index] # the last added bond the ring bond
                current_index = next_index
            else: # were in a ring and have to find the closest path to the first atom of the ring
                current_ring_indexes = []
                for ring_index in self.ring_indexes: # find the ring we are in right now
                    if current_index in ring_index:
                        current_ring_indexes = ring_index
                        break
                distance_from_the_first_index = current_ring_indexes.index(current_index) # find out where are we in the ring relative to the first index of the ring 
                distance_from_the_last_index = (len(current_ring_indexes) - 1) - distance_from_the_first_index # our relative distance from the last index of the chain
                if distance_from_the_first_index < distance_from_the_last_index: # its easier to get to the first index 
                    next_ring_index = current_ring_indexes[distance_from_the_first_index-1] 
                    while True:
                        distance += self.previous_bond_lengths[current_index][next_ring_index]
                        current_index = next_ring_index
                        distance_from_the_first_index -= 1
                        if (distance_from_the_first_index == 0): break
                        else:
                            next_ring_index = current_ring_indexes[distance_from_the_first_index-1] 
                else: # its easier to get to the last index and go to the first index from there
                    next_ring_index = current_ring_indexes[distance_from_the_first_index+1] 
                    while True:
                        distance += self.bond_lengths[current_index][next_ring_index]
                        current_index = next_ring_index
                        distance_from_the_first_index += 1
                        if (distance_from_the_first_index == len(current_ring_indexes)-1): break
                        else:
                            next_ring_index = current_ring_indexes[distance_from_the_first_index+1] 
                    distance += self.bond_lengths[current_index][current_ring_indexes[0]]
                    current_index = current_ring_indexes[0]
        return (distance,current_index) 

    def Calc_shortest_distance_from_first_atom(self,index): # make sure index does not have the atomic number of -1 or 0
        distance,distance_from_main_chain = (0,0) 
        current_index = index
        first_index = self.main_chain_indexes[0] 

        if self.Is_it_side_chain(current_index): # get to the main chain
            (distance_from_main_chain,current_index) = self.Calc_shortest_distance_from_main_chain(current_index)
        if current_index == first_index: return (distance,distance_from_main_chain) 
        while True:
            if (current_index in self.main_chain_indexes): 
                i = self.main_chain_indexes.index(current_index)
                next_index = self.main_chain_indexes[i-1] 
                distance += self.previous_bond_lengths[current_index][next_index]
                current_index = next_index
                if (current_index == first_index): break
                else:
                    next_index = self.main_chain_indexes[i-1] 
            elif current_index in self.ring_ends: # were at the end of the ring and if we move one step forward we get to the begining of the ring
                next_index = self.atoms[current_index].bonded_atoms[-1][0]
                distance += self.bond_lengths[current_index][next_index] # the last added bond the ring bond
                current_index = next_index
            else: # the index is part of a ring and we need to move to the first atom of the ring
                current_ring_indexes = []
                for ring_index in self.ring_indexes: # find the ring we are in right now
                    if current_index in ring_index:
                        current_ring_indexes = ring_index
                        break
                distance_from_the_first_index = current_ring_indexes.index(current_index) # find out where are we in the ring relative to the first index of the ring 
                distance_from_the_last_index = (len(current_ring_indexes) - 1) - distance_from_the_first_index # our relative distance from the last index of the chain
                if distance_from_the_first_index < distance_from_the_last_index: # its easier to get to the first index 
                    next_ring_index = current_ring_indexes[distance_from_the_first_index-1] 
                    while True:
                        distance += self.previous_bond_lengths[current_index][next_ring_index]
                        current_index = next_ring_index
                        distance_from_the_first_index -= 1
                        if (distance_from_the_first_index == 0): break
                        else:
                            next_ring_index = current_ring_indexes[distance_from_the_first_index-1] 
                else: # its easier to get to the last index and go to the first index from there
                    next_ring_index = current_ring_indexes[distance_from_the_first_index+1] 
                    while True:
                        distance += self.bond_lengths[current_index][next_ring_index]
                        current_index = next_ring_index
                        distance_from_the_first_index += 1
                        if (distance_from_the_first_index == len(current_ring_indexes)-1): break
                        else:
                            next_ring_index = current_ring_indexes[distance_from_the_first_index+1] 
                    distance += self.bond_lengths[current_index][current_ring_indexes[0]]
                    current_index = current_ring_indexes[0]
        return (distance+distance_from_main_chain,distance_from_main_chain) 
    
    def Calc_first_or_last_atom_penalty(self, index):
        begining_or_end_atom_penalty = 0
        if not (self.molecule or self.one_side_connected):
                begining_or_end_atom_penalty = 1 if (index == self.main_chain_indexes[0] or index == self.main_chain_indexes[-1]) else 0
        elif self.one_side_connected:
                begining_or_end_atom_penalty = 1 if index == self.main_chain_indexes[0] else 0
        return begining_or_end_atom_penalty
    
    def Create_functional_groups(self): # functional groups for calculating hydrogen bond power
        # groups: COOH, COOC, CO, OH, O, --NH, =NH, -NH2, N, =N, SH, Cl, F
        counted_indexes = [] # this is for structures like COOC, when we encounter one oxygen we add the whole structure to the self.functional_groups dictionary and add the index of both oxygen groups to [counted_groups] list
        for index,atom in enumerate(self.atoms):
            if index in counted_indexes: continue
            begining_or_end_atom_penalty = self.Calc_first_or_last_atom_penalty(index)
            number_of_bonded_atoms  = len(atom.bonded_atoms) + len(atom.previous_bonded_atoms) + begining_or_end_atom_penalty 
            double_bond_flag = len([1 for bond in atom.previous_bonded_atoms if bond[1] == 2]) == 1 or len([1 for bond in atom.bonded_atoms if bond[1] == 2]) == 1  
            if (atom.atomic_number == 6 and double_bond_flag):
                oxygen_indexes = [bond[0] for bond in atom.bonded_atoms if self.atoms[bond[0]].atomic_number == 8]
                match len(oxygen_indexes):
                    case 1: 
                        self.functional_groups["CO"] += 1
                    case 2:
                        if len(self.atoms[oxygen_indexes[1]].bonded_atoms) == 1: self.functional_groups["COOC"] += 1
                        else: self.functional_groups["COOH"] += 1
                    case _:
                        print("could not identify carbon, ","index:",atom.index)
                counted_indexes.extend(oxygen_indexes)

            elif (atom.atomic_number == 8):
                match number_of_bonded_atoms:
                    case 2: # ether
                        self.functional_groups["O"] += 1
                    case 1: # OH
                        if not double_bond_flag:
                            self.functional_groups["OH"] += 1
                    case _:
                        if self.molecule: # H2O -> O
                            self.functional_groups["OH"] += 2
                        else:
                            print("could not identify oxygen, ","index:",atom.index)

            elif (atom.atomic_number == 7):
                double_bond_flag = "1" if ((2 in [bond[1] for bond in atom.bonded_atoms]) or (2 in [bond[1] for bond in atom.previous_bonded_atoms])) else "0" # 1 means true and 0 means false in this contex
                triple_bond_flag = "1" if (3 in [bond[1] for bond in atom.previous_bonded_atoms]) else "0"
                code = str(number_of_bonded_atoms) + double_bond_flag + triple_bond_flag
                match code: # think about this part more
                    case "300": # three single bonds
                        self.functional_groups["N"] += 1
                    case "200": # connected to two atoms and no double or triple bonds
                        self.functional_groups["--NH"] += 1
                    case "210": # connected to two atoms and one double bond bond
                        self.functional_groups["=N"] += 1
                    case "100":
                        self.functional_groups["-NH2"] += 1
                    case "110":
                        self.functional_groups["=NH"] += 1
                    case "101":
                        self.functional_groups["#N"] += 1
                    case _:
                        print("could not identify nitrogen, ","index:",atom.index)
            
            elif (atom.atomic_number == 17):
                self.functional_groups["Cl"] += 1
            
            elif (atom.atomic_number == 9):
                self.functional_groups["F"] += 1
            
            elif (atom.atomic_number == 16):
                if len(atom.bonded_atoms) == 0: self.functional_groups["SH"] += 1

    def Filter_functional_groups(self):
        return  {k:v for k,v in self.functional_groups.items() if v != 0}

    def Calc_polar_bonds(self):
        # {6: [(9, 1.0), (17, 1.5)]},
        polar_bonds = []
        for atom in self.atoms:
            if (atom.atomic_number == -1) or (atom.atomic_number == 0) or (atom.ring):
                continue
            if (atom.atomic_number == 7) or (atom.atomic_number == 8):
                first_or_end_atom_penalty = self.Calc_first_or_last_atom_penalty(atom.index) # add one bond if its the first or last atom of main chain
                number_of_bonds = first_or_end_atom_penalty + sum([bond[1] for bond in (atom.bonded_atoms)]) + sum([bond[1] for bond in (atom.previous_bonded_atoms)])  
                if (atom.atomic_number == 7) and not (number_of_bonds == 3):
                    polar_bonds.append({7:[(1, 1.01 + 0.65 + 0.25)]*(3 - number_of_bonds)})
                elif (atom.atomic_number == 8) and not (number_of_bonds == 2):
                    polar_bonds.append({8:[(1, 0.9575 + 0.6 + 0.25)]*(2- number_of_bonds)})
            side_chain_flag = self.Is_it_side_chain(atom.index) or self.molecule 
            bonded_side_chain_indexes = [] 
            if side_chain_flag :
                bonded_side_chain_indexes = [bond[0] for bond in atom.bonded_atoms] # all atoms are side chain
            else:
                bonded_side_chain_indexes = [bond[0] for bond in atom.bonded_atoms if (self.Is_it_side_chain(bond[0]) and (not self.atoms[bond[0]].ring))]

            if bonded_side_chain_indexes:
                polar_bonds.append({atom.atomic_number:[(self.atoms[index].atomic_number, self.bond_lengths[atom.index][index]) for index in bonded_side_chain_indexes]})
        return polar_bonds
    
    def Save_descriptors(self,structure_name="test",name="test"):
        # setting up workbook
        workbook = xlsxwriter.Workbook(name+".xlsx")  
        worksheet = workbook.add_worksheet("important values")
        # setting up formats
        titles_format = workbook.add_format()
        titles_format.set_bold()
        titles_format.set_align("center")
        titles_format.set_align("vcenter")
        # format for other cells 
        normal_format = workbook.add_format()
        normal_format.set_align("center")
        normal_format.set_align("vcenter")

        worksheet.write(0,0,"structure name",titles_format)
        worksheet.write(0,1,"structure",titles_format)
        worksheet.write(0,2,"flexibility",titles_format)
        worksheet.write(0,3,"occupied volume",titles_format)
        worksheet.write(0,4,"CH",titles_format)
        worksheet.write(0,5,"number of CH",titles_format)
        worksheet.write(0,6,"C=OH",titles_format)
        worksheet.write(0,7,"number of C=OH",titles_format)
        worksheet.write(0,8,"C=OOH",titles_format)
        worksheet.write(0,9,"number of C=OOH",titles_format)
        worksheet.write(0,10,"OH",titles_format)
        worksheet.write(0,11,"number of OH",titles_format)
        worksheet.write(0,12,"NH",titles_format)
        worksheet.write(0,13,"number of NH",titles_format)
        worksheet.write(0,14,"B",titles_format)
        worksheet.write(0,15,"number of B",titles_format)
        worksheet.write(0,16,"BCH",titles_format)
        worksheet.write(0,17,"number of BCH",titles_format)
        worksheet.write(0,18,"BC=OH",titles_format)
        worksheet.write(0,19,"number of BC=OH",titles_format)
        worksheet.write(0,20,"BC=OOH",titles_format)
        worksheet.write(0,21,"number of BC=OOH",titles_format)
        worksheet.write(0,22,"BOH",titles_format)
        worksheet.write(0,23,"number of BOH",titles_format)
        worksheet.write(0,24,"BNH",titles_format)
        worksheet.write(0,25,"number of BNH",titles_format)
        
        worksheet.write(1,0,structure_name,normal_format)
        worksheet.write(1,1,self.structure,normal_format)
        worksheet.write(1,2,self.flexibility,normal_format)
        worksheet.write(1,3,self.occupied_volume,normal_format)
        worksheet.write(1,4,self.important_functional_groups["CH"],normal_format)
        worksheet.write(1,5,self.number_of_important_functional_groups["CH"],normal_format)
        worksheet.write(1,6,self.important_functional_groups["C=OH"],normal_format)
        worksheet.write(1,7,self.number_of_important_functional_groups["C=OH"],normal_format)
        worksheet.write(1,8,self.important_functional_groups["C=OOH"],normal_format)
        worksheet.write(1,9,self.number_of_important_functional_groups["C=OOH"],normal_format)
        worksheet.write(1,10,self.important_functional_groups["OH"],normal_format)
        worksheet.write(1,11,self.number_of_important_functional_groups["OH"],normal_format)
        worksheet.write(1,12,self.important_functional_groups["NH"],normal_format)
        worksheet.write(1,13,self.number_of_important_functional_groups["NH"],normal_format)
        worksheet.write(1,14,self.important_functional_groups["B"],normal_format)
        worksheet.write(1,15,self.number_of_important_functional_groups["B"],normal_format)
        worksheet.write(1,16,self.important_functional_groups["BCH"],normal_format)
        worksheet.write(1,17,self.number_of_important_functional_groups["BCH"],normal_format)
        worksheet.write(1,18,self.important_functional_groups["BC=OH"],normal_format)
        worksheet.write(1,19,self.number_of_important_functional_groups["BC=OH"],normal_format)
        worksheet.write(1,20,self.important_functional_groups["BC=OOH"],normal_format)
        worksheet.write(1,21,self.number_of_important_functional_groups["BC=OOH"],normal_format)
        worksheet.write(1,22,self.important_functional_groups["BOH"],normal_format)
        worksheet.write(1,23,self.number_of_important_functional_groups["BOH"],normal_format)
        worksheet.write(1,24,self.important_functional_groups["BNH"],normal_format)
        worksheet.write(1,25,self.number_of_important_functional_groups["BNH"],normal_format)
        worksheet.autofit()
        workbook.close()

    def Main_chain_test(self):
        for index in self.main_chain_indexes:
            if not (self.atoms[index].atomic_number == -1 or self.atoms[index].atomic_number == 0):
                self.utilities.Print_atom(index)

    def Atom_recognition_test(self):
        for i in range(self.len_structure):
            if not (self.atoms[i].atomic_number == -1 or self.atoms[i].atomic_number == 0):
                self.utilities.Print_atom(i)

                if not len(self.atoms[i].previous_bonded_atoms) == 0:
                    print("previously bonded atoms:")
                    for bond in self.atoms[i].previous_bonded_atoms:
                        self.utilities.Print_atom(bond[0],end_="  ")
                        print(f"bond order: {bond[1]} bond length: {self.previous_bond_lengths[i][bond[0]]}")

                if not len(self.atoms[i].bonded_atoms) == 0:
                    print("bonded atoms after specified atom:")
                    for bond in self.atoms[i].bonded_atoms:
                        self.utilities.Print_atom(bond[0],end_="  ")
                        print(f"bond order: {bond[1]}  bond length: {self.bond_lengths[i][bond[0]]}")
                if not (i == self.len_structure-1): print("=================================================")

    def Bridge_recognition_test(self):
        for bridge in self.bridges:
            for j in range(bridge[0],bridge[1]+1):
                print(self.structure[j],end="")
            if len(bridge) == 3: 
                print(" linked: ",end = "")
                for k in range(self.bridges[bridge[2]][0],self.bridges[bridge[2]][1]+1):
                    print(self.structure[k],end="")
            print()

    def Connecting_bridges_test(self):
        for i in range(self.len_structure):
            bridges = self.Find_connected_bridges(i)
            if not len(bridges) == 0: 
                if (self.atoms[i-1].atomic_number in [14,17]):
                    print("atom:", self.structure[self.atoms[i-1].index]+self.structure[self.atoms[i].index],end=" ")
                else:
                    print("atom:", self.structure[self.atoms[i].index],end=" ")
            for bridge in bridges:
                for j in range(bridge[0],bridge[1]+1):
                    print(self.structure[j],end="")
                if len(bridge) == 3: 
                    bridges.append(self.bridges[bridge[2]])
                    print("\nlinked bridge:", end=" ")
                print()

    def Ring_recognition_test(self):
        for ring_indexes in self.ring_indexes:
            for index in ring_indexes:
                print(self.structure[index],end=",")
            print()
        if len(self.ring_indexes) == 0: print("no rings where found in this structure")

    #def Shortest_distance_from_main_chain_test(self):
    #    for index,atom in enumerate(self.atoms):
    #        if not (atom.atomic_number == -1 or atom.atomic_number == 0):
    #            self.utilities.Print_atom(index)
    #            if self.Is_it_side_chain(index):
    #                print(f"shortest distance from main chain: {self.Calc_shortest_distance_from_main_chain(index)[0]}")
    #            else:
    #                print(f"shortest distance from main chain: {0}")
    #            print("================================================================")

    #def Shortest_distance_from_first_atom_test(self):
    #    for index,atom in enumerate(self.atoms):
    #        if not (atom.atomic_number == -1 or atom.atomic_number == 0):
    #            self.utilities.Print_atom(index)
    #            print(f"shortest distance from the first atom of the chain: {self.Calc_shortest_distance_from_first_atom(index)[0]}")
    #            print("================================================================")
                
    # def Hydrogen_bondin_groups_test(self):
    #     self.Calc_hydrogen_DA_groups()

    def Draw_structure_RDKit(self,smiles_structure):
        from rdkit import Chem
        from rdkit.Chem import Draw
        #from rdkit.Chem import Descriptors
        from PIL import Image

        # rdkit
        m = Chem.MolFromSmiles(smiles_structure)
        Draw.ShowMol(m)
        print("test complete!")
        print("================================================================")
        # to draw and save image of molecule
        #image = Draw.MolToImage(m)
        #image.show()
        #image.save("toloene.png")
        #exit()

        # to see all the descriptors for moluclue
        #results = Descriptors.CalcMolDescriptors(m)
        #file = open("descriptors.txt","w")
        #for key,value in results.items():
        #    file.write(str(key)+"    "+str(value)+"\n")
        #file.close()
    def functional_groups_recognition_test(self):
        for key, value in self.functional_groups.items():
            print(f"funcitonal group: {key} -> {value}")

    def Hydrogen_bond_between_same_structure_test(self):
        print(f"the hydrogen bond energy between the same structures is: {self.utilities.Calc_hydrogen_bond_energy(self.functional_groups,self.functional_groups)}")
    
    def Calc_polar_interaction_between_same_structure_test(self):
        polar_bonds = self.Calc_polar_bonds()
        return self.utilities.Calc_polar_interaction(polar_bonds, polar_bonds)

    def Run_all_tests(self):
        self.utilities.Set_structure(self.structure)
        self.utilities.Set_atoms(self.atoms)
        print(self.smiles_structure)
        print("===========================main chain test=====================================")
        self.Main_chain_test()
        print("===========================atom recognition test=====================================")
        self.Atom_recognition_test()
        print("===========================bridge recognition test=====================================")
        self.Bridge_recognition_test()
        print("===========================connecting bridges test=====================================")
        self.Connecting_bridges_test()
        print("===========================ring recognition test=====================================")
        self.Ring_recognition_test()
        print("===========================benzene rings recognition test=====================================")
        print("number of inner benzene rings",self.number_of_inner_benzene_rings)
        print("number of side benzene rings",self.number_of_side_benzene_rings)
        print("===========================functional groups test=====================================")
        self.functional_groups_recognition_test()
        print("===========================hydrogen bonding energy test=====================================")
        self.Hydrogen_bond_between_same_structure_test()
        print("===========================polar interactions test=====================================")
        print(self.Calc_polar_interaction_between_same_structure_test())
        # print("===========================Saving descriptors test=====================================")
        # self.Save_descriptors()
        # print("Save complete")
        print("===========================testing structure using rdkit=====================================")
        self.Draw_structure_RDKit(self.smiles_structure)
#Run_all_tests()

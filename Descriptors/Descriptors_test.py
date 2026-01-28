import Descriptors

def Test():
    #smiles_structure = '*C(=O)c1cc2cc(ccc2cc1)C(=O)[Si](c3cc(C)c(C(=O)O{-})cc3)([Cl])C*'
    smiles_structure = '*C(=O)c1cc2cc(ccc2cc1)C(=O)[Si](c3cc(C)c(C(=O)O)cc3)([Cl])C*'
    #smiles_structure = "*NC=1CC(CCC1)C*"
    structure_proccesor = Descriptors.Structure_Processor(smiles_structure)
    structure_proccesor.Pre_process_structures()
    structure_proccesor.Create_descriptors()
    for descriptor in structure_proccesor.descriptors:
        if not descriptor.molecule:
            descriptor.Calc_polymer_descriptors()
        descriptor.Run_all_tests()

Test()

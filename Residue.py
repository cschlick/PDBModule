class Residue():
    def __init__(self,residue_name,atoms):
        self.type = "Residue"
        self.name = residue_name
        self.residue_name = self.name
        self.atoms = atoms
        self.res_seq_num = atoms[0].res_seq_num
        self.chain_id = None
        self.first_atom_number = atoms[0].serialnumber
        #self.validate()
        self.chain_id = atoms[0].chain_id

    def update(self):
        for atom in self.atoms:
            atom.res_seq_num = self.res_seq_num
            atom.chain_id = self.chain_id

    def validate(self):
        self.atom_nums = []
        self.atom_names = []
        for atom in self.atoms:
            self.atom_nums.append(atom.res_seq_num)
            self.atom_names.append(atom.residue_name)

    def print(self):

        for atom in self.atoms:
            atom.print()
        print("\n")

    def add_atom(self,atom):
        self.atoms.append(atom)
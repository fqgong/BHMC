import random
import sys
import re
import math
import numpy as np
from ase import Atoms,Atom
from pyxtal import pyxtal
from pyxtal.symmetry import Group



class Gen_Structure(Atoms):
    def __init__(self,*args,**kw):
        super().__init__(*args,**kw)
    def wrap_atoms(self):
        rel_pos = []
        diss = []
        self.chemical_symbols=self.get_chemical_symbols()
        for i,atom in enumerate(self):
            vec = self.get_distance(0,i,mic=True,vector=True)
            rel_pos.append(vec)
            diss.append(self.get_distance(0,i,mic=True))
        new_pos = []
        cell = max(diss)+10
        for pos in rel_pos:
            new_pos.append([coord+cell/2 for coord in pos])
        self = Gen_Structure(symbols=self.chemical_symbols, positions=new_pos, cell=[cell,cell,cell],pbc=True)
        return self

class Gen_Initial_Structure(Gen_Structure):
    def __init__(self,*args,**kw):
        super().__init__(*args, **kw)
    def bySymmetry(self,dimension=0):
        symbols_dict = {}
        for atom in self.get_chemical_symbols():
            symbols_dict[atom] = self.get_chemical_symbols().count(atom)
        chemical_symbols, chemical_nums = list(symbols_dict.keys()), list(symbols_dict.values())
        max_diss = 5
        while max_diss > 4:
            structure = pyxtal()
            sg=random.randint(1,56)
            if Group(sg,dim=dimension).check_compatible(numIons=chemical_nums)[0]:
                structure.from_random(dim=dimension, group=sg, species=chemical_symbols, numIons=chemical_nums)
                stc = Gen_Structure(structure.to_ase())
                self = stc.wrap_atoms()
                diss = []
                for i in range(len(self)):
                    dis = []
                    for j in range(len(self)):
                        dis.append(self.get_distance(i, j))
                    diss.append(sorted(dis)[1])
                max_diss = max(diss)
        return self
    def byRandom(self,size_dict):
        xMin,xMax = re.findall("[-+]?\d*\.\d+|\d+", size_dict['x_range'])
        yMin,yMax = re.findall("[-+]?\d*\.\d+|\d+", size_dict['y_range'])
        zMin,zMax = re.findall("[-+]?\d*\.\d+|\d+", size_dict['z_range'])
        vacuumMin,vacuumMax = re.findall("[-+]?\d*\.\d+|\d+", size_dict['vacuum_range'])
        cell_size_x = float(xMax) - float(xMin)
        cell_size_y = float(yMax) - float(yMin)
        cell_size_z = float(vacuumMax) - float(vacuumMin)
        
        while True:
            atoms_list = []
            for i in range(len(self)):
                condition = 0
                if not atoms_list:
                    atom = Atom('X',position=(0,0,0))
                else:
                    atom_temp = None
                    while condition == 0:
                        theta = random.uniform(0, 2*math.pi)
                        phi = random.uniform(0, 2*math.pi)
                        r = random.uniform(1.8,3.0)
                        
                        random_list_index = random.randint(0,len(atoms_list)-1)
                        x = atoms_list[random_list_index].position[0]+(r*math.cos(theta)*math.sin(phi))
                        y = atoms_list[random_list_index].position[1]+(r*math.sin(theta)*math.sin(phi))
                        z = atoms_list[random_list_index].position[2]+(r*math.cos(phi))
                        atom_temp = Atom('X',position=(x,y,z))
                        for atm in atoms_list:
                            dist = math.sqrt(math.pow(atom_temp.position[0] - atm.position[0], 2) +
                                             math.pow(atom_temp.position[1] - atm.position[1], 2) +
                                             math.pow(atom_temp.position[2] - atm.position[2], 2))
                            if dist < 1.5 or dist > (len(self)+1) * random.uniform(0.2,0.4)*1.8:
                                condition = 0
                                break
                            else:
                                condition = 1
                    atom = atom_temp
                atoms_list.append(atom)
            stc = Atoms(atoms_list)
            
            sys_x_size = max(stc.get_positions()[:,0]) - min(stc.get_positions()[:,0])
            sys_y_size = max(stc.get_positions()[:,1]) - min(stc.get_positions()[:,1])
            sys_z_size = max(stc.get_positions()[:,2]) - min(stc.get_positions()[:,2])
            if sys_x_size <= cell_size_x and sys_y_size <= cell_size_y and sys_z_size <= cell_size_z:
                break
            
        displacement_to_starting_point_x = float(xMin) - min(stc.get_positions()[:,0])
        displacement_to_starting_point_y = float(yMin) - min(stc.get_positions()[:,1])
        displacement_to_starting_point_z = float(zMin) - min(stc.get_positions()[:,2])
        x_random_top = cell_size_x - sys_x_size
        y_random_top = cell_size_y - sys_y_size
        z_random_top = float(zMax) - float(zMin)
        random_displacement_x = random.uniform(0, x_random_top)
        random_displacement_y = random.uniform(0, y_random_top)
        random_displacement_z = random.uniform(0, z_random_top)
             
        fitted2cell_poss = stc.get_positions()+(displacement_to_starting_point_x + random_displacement_x,
                                                displacement_to_starting_point_y + random_displacement_y,
                                                displacement_to_starting_point_z + random_displacement_z)
        fitted_stc = Gen_Structure(self.get_chemical_formula(),positions=fitted2cell_poss)
        self = fitted_stc.wrap_atoms()
            
        return self
        

class Perturb_Structure(Gen_Structure):
    def __init__(self,*args,**kw):
        super().__init__(*args,**kw)
        self.chemical_symbols = self.get_chemical_symbols()
        self.coords = self.get_positions()
    def random_move(self,stepsize=0.1):
        new_coords = []
        for coord in self.coords.reshape(-1):
            random.seed()
            new_coords.append(coord + 2*stepsize*random.uniform(-1,1))
        new_pos = np.array(new_coords).reshape(int(len(self.coords)),3)
        new_stc = Gen_Structure(Atoms(self.chemical_symbols, positions=new_pos, pbc=True))
        self = new_stc.wrap_atoms()
        return self

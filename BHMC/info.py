import os, re
from ase.io import read,write




class Info_Distill:
    def __init__(self,num_iter,num_task):
        path = 'BHMC-workflow/iter.{}/task.{}/'.format(num_iter,num_task)
        self.path = os.path.abspath(path)

class CP2K_Info(Info_Distill):
    def __init__(self,num_iter,num_task):
        super().__init__(num_iter,num_task)
        self.final_pos = read(os.path.join(self.path,'BH_GEO-pos-1.xyz'),'-1')
        self.energy = self.final_pos.info['E'] * 2.72113838565563E+01
        final_frc = read(os.path.join(self.path,'BH_GEO-frc-1.xyz'),'-1')
        self.force = max(abs(final_frc.get_positions().reshape(-1)*2.72113838565563E+01/5.29177208590000E-01))

class LAMMPS_Info(Info_Distill):
    def __init__(self,num_iter,num_task):
        super().__init__(num_iter,num_task)
        with open(os.path.join(self.path,'log.lammps'),'r') as f:
            lines = f.readlines()
        for i,line in enumerate(lines):
            if re.findall(r'Energy initial, next-to-last, final',line)!=[]:
                self.energy = float(re.split(r'\s+',lines[i+1])[3])
            if re.findall(r'Force max component initial, final',line)!=[]:
                self.force = float(re.findall(r'[-+]?(?:(?:\d*\.\d+)|(?:\d+\.?))(?:[Ee][+-]?\d+)?',line)[1])
        self.final_pos = read(os.path.join(self.path,'BH_GEO-pos-1.xyz'),'-1')

class REFINE_Info(Info_Distill):
    def __init__(self,num_iter,num_task=1,with_fp=False):
        super().__init__(num_iter,num_task=1)
        self.mlp_final_pos = read(os.path.join(self.path,'BH_GEO-pos-1.xyz'),'-1')
        with open(os.path.join(self.path,'log.lammps'),'r') as f:
            lines = f.readlines()
        for i,line in enumerate(lines):
            if re.findall(r'Force max component initial, final',line)!=[]:
                self.mlp_force = float(re.findall(r'[-+]?(?:(?:\d*\.\d+)|(?:\d+\.?))(?:[Ee][+-]?\d+)?',line)[1])
        if with_fp:
            self.fp_final_pos = read(os.path.join(self.path,'fp/BH_GEO-pos-1.xyz'),'-1')
            self.energy = self.fp_final_pos.info['E'] * 2.72113838565563E+01
            fp_final_frc = read(os.path.join(self.path,'fp/BH_GEO-frc-1.xyz'),'-1')
            self.fp_force = max(abs(final_frc.get_positions().reshape(-1)*2.72113838565563E+01/5.29177208590000E-01))

class OTHER_Info(Info_Distill):
    def __init__(self,num_iter,num_task):
        super().__init__(num_iter,num_task)
        pass

class Info_Factory:
        def style(self,style,num_iter,num_task):
            get_info = {
                'cp2k': CP2K_Info,
                'lammps': LAMMPS_Info,
                'refine': REFINE_Info,
            }.get(style,OTHER_Info)
            return get_info(num_iter,num_task) 

import math
import random
from run.run import Run,CP2KRun,LAMMPSRun




class Estimation:
    def __init__(self,inp_dict):
        self.force_criterion = inp_dict['max_force']
        self.style = inp_dict['geo_type']
        self.machine_dict = inp_dict['machine_dict']
        self.kbT = 8.617332478e-05 * inp_dict['temperature']

    def frc_est(self,force):
        if force > self.force_criterion:
            run_type, machine_dict ={
                'cp2k':(CP2KRun,self.machine_dict),
                'lammps':(LAMMPSRun,self.machine_dict),
            }.get(self.style,None)
            return run_type(machine_dict) 
        else: 
            return 0
    def mc_est(self,e1,e2):
        if  math.exp((e1-e2)/self.kbT) > random.uniform(0,1):
            return 0
        else: 
            run_type, machine_dict = {
                'cp2k':(CP2KRun,self.machine_dict),
                'lammps':(LAMMPSRun,self.machine_dict),
            }.get(self.style,None)
            return run_type(machine_dict)

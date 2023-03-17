from dpdispatcher import Resources,Task,Submission,Job,Machine
from dpdispatcher.ssh_context import SSHContext
from dpdispatcher.lsf import LSF
from dpdispatcher.slurm import Slurm




class Run:
    def __init__(self,inp_dict):
        self.inp_dict = inp_dict
        self.machine = Machine.load_from_dict(self.inp_dict['machine'])
        self.resources = Resources.load_from_dict(self.inp_dict['resources'])
#        self.machine=self.inp_dict['machine']
#        self.resources=self.inp_dict['resources']
    def sub_task(self,task_lst):
        submission = Submission(
            machine=self.machine,
            work_base='BHMC-workflow/',
            resources = self.resources,
            forward_common_files=[],
            backward_common_files=[],
            task_list=task_lst
        )
        return submission.run_submission()

class CP2KRun(Run):
    def __init__(self,inp_dict):
        super().__init__(inp_dict)
    def gen_task(self,num_iter,num_task):
        task = Task(
            command=self.inp_dict['command'],
            task_work_path='iter.{}/task.{}/'.format(num_iter,num_task),
            forward_files=['input.inp','coord.xyz'],
            backward_files=['*'],
            outlog='cp2k.log',
            errlog='cp2k.err'
        )
        return task

class LAMMPSRun(Run):
    def __init__(self,inp_dict):
        super().__init__(inp_dict)
    def gen_task(self,num_iter,num_task):
        task = Task(
            command=self.inp_dict['command'],
            task_work_path='iter.{}/task.{}/'.format(num_iter,num_task),
            forward_files=['input.lammps','conf.lmp','frozen_model.pb'],
            backward_files=['BH_GEO-pos-1.xyz','log.lammps'],
            outlog='lammps.log',
            errlog='lammps.err'
        )
        return task

import os
import sys
import json
import random
import numpy as np
from glob import glob
from make_files import Make_Folders, Make_CP2K_Files, Make_LAMMPS_Files
from gen_structure import Gen_Structure, Gen_Initial_Structure, Perturb_Structure
from run import Run, CP2KRun, LAMMPSRun
from info import Info_Distill, CP2K_Info, LAMMPS_Info, Info_Factory
from esti import Estimation

input_file = str(sys.argv[1])

ori_path = os.getcwd()
with open(os.path.join(ori_path, input_file),'r') as f:
    jdata = json.load(f)
    
chem_formula = jdata['chem_formula']
dimension = jdata['dimension']
ini_style = jdata['ini_style']
cell_size_dict = jdata.get('cell_size_dict',None)
geo_type = jdata['geo_type']
num_replicas = jdata.get('num_replicas', 1)
geo_params = jdata['geo_params']
machine_dict = jdata['machine_dict']
total_iteration = jdata['iteration']
step_size = jdata['step_size']

if geo_type == 'cp2k':
    print("Basin Hopping Monte Carlo (BHMC) workflow conducted by DFT-CP2K ")
if geo_type == 'lammps':
    print("Basin Hopping Monte Carlo (BHMC) workflow conducted by MLP-LAMMPS")

for num_task in range(num_replicas):
    folder = {
        'cp2k': Make_CP2K_Files(num_iter=0, num_task=num_task),
        'lammps': Make_LAMMPS_Files(num_iter=0, num_task=num_task)
    }.get(geo_type, None)
    folder.make_folders()
    stc = Gen_Initial_Structure(chem_formula)
    while np.all(stc.get_positions()==0):
        if ini_style == 'symmetry': 
            stc = stc.bySymmetry(dimension=dimension)
        elif ini_style == 'random':
            stc = stc.byRandom(cell_size_dict)
    folder.make_files(stc=stc, params_dict=geo_params)
    
task_lst = []
for num_task in range(num_replicas):
    run_type = {
        'cp2k': CP2KRun,
        'lammps': LAMMPSRun,
    }.get(geo_type, None)
    task_lst.append(run_type(machine_dict).gen_task(num_iter=0,num_task=num_task))

geo_work = Run(inp_dict = machine_dict)
geo_work.sub_task(task_lst=task_lst)

estis = 1
infos = {}
while not np.all(np.array(estis) == 0):
    estis = []
    for num_task in range(num_replicas):
        infos[str(num_task)] = Info_Factory().style(geo_type, num_iter=0, num_task=num_task)
        esti = Estimation(jdata)
        estis.append(esti.frc_est(infos[str(num_task)].force))
    task_lst = []
    for num_task,esti in enumerate(estis):
        if esti != 0:
            folder = {
                'cp2k': Make_CP2K_Files(num_iter=0, num_task=num_task),
                'lammps': Make_LAMMPS_Files(num_iter=0, num_task=num_task)
            }.get(geo_type, None)
            folder.make_folders()
            stc = Gen_Initial_Structure(chem_formula)
            while np.all(stc.get_positions()==0):
                stc = stc.bySymmetry(dimension=dimension)
            folder.make_files(stc=stc, params_dict = geo_params)
            print('The Geometry Optimization of Iteration # {} Task {} has failed!'.format(0, num_task))
            print('Start again with random initial structure!')
            run_type = {
                'cp2k': CP2KRun,
                'lammps': LAMMPSRun,
            }.get(geo_type, None)
            task_lst.append(run_type(machine_dict).gen_task(num_iter=0,num_task=num_task))
    if task_lst != []:
        geo_work = Run(inp_dict = machine_dict)
        geo_work.sub_task(task_lst=task_lst)
print(" --> Geometry Optimization of Initial Structure: DONE! ")

print("BHMC cycle starts here! ")

num_iter = 1
while num_iter < total_iteration:
    # 生成结构
    previous_infos = infos
    infos = {}
    for num_task in range(num_replicas):
        folder = {
            'cp2k': Make_CP2K_Files(num_iter=num_iter, num_task=num_task),
            'lammps': Make_LAMMPS_Files(num_iter=num_iter, num_task=num_task)
        }.get(geo_type, None)
        folder.make_folders()
        stc = Perturb_Structure(previous_infos[str(num_task)].final_pos).random_move(stepsize=step_size)
        folder.make_files(stc=stc, params_dict=geo_params)
    
    # 提交计算
    task_lst = []
    for num_task in range(num_replicas):
        run_type = {
            'cp2k': CP2KRun,
            'lammps': LAMMPSRun,
        }.get(geo_type, None)
        task_lst.append(run_type(machine_dict).gen_task(num_iter=num_iter,num_task=num_task))
    geo_work = Run(inp_dict = machine_dict)
    geo_work.sub_task(task_lst=task_lst)
    
    # 判断收敛
    estis = 1
    while not np.all(np.array(estis) == 0):
        estis = []
        for num_task in range(num_replicas):
            infos[str(num_task)] = Info_Factory().style(geo_type, num_iter=num_iter, num_task=num_task)
            esti = Estimation(jdata)
            estis.append(esti.frc_est(infos[str(num_task)].force))
        task_lst = []
        for num_task, esti in enumerate(estis):
            if esti != 0:
                folder = {
                    'cp2k': Make_CP2K_Files(num_iter=num_iter, num_task=num_task),
                    'lammps': Make_LAMMPS_Files(num_iter=num_iter, num_task=num_task)
                }.get(geo_type, None)
                folder.make_folders()
                stc = Perturb_Structure(previous_infos[str(num_task)].final_pos).random_move(stepsize=step_size)
                folder.make_files(stc=stc, params_dict = geo_params)
                print('The Geometry Optimization of Iteration # {} Task {} has failed!'.format(num_iter, num_task))
                print('Start again from previous perturbed structure!')
                run_type = {
                    'cp2k': CP2KRun,
                    'lammps': LAMMPSRun,
                }.get(geo_type, None)
                task_lst.append(run_type(machine_dict).gen_task(num_iter=num_iter,num_task=num_task))
        if task_lst != []:
            geo_work = Run(inp_dict = machine_dict)
            geo_work.sub_task(task_lst=task_lst)
        
    # Metroplis Criterion
    estis = 1
    while not np.all(np.array(estis) == 0):
        estis = []
        for num_task in range(num_replicas):
            infos[str(num_task)] = Info_Factory().style(geo_type, num_iter=num_iter, num_task=num_task)
            esti = Estimation(jdata)
            estis.append(esti.mc_est(previous_infos[str(num_task)].energy,infos[str(num_task)].energy))
        task_lst = []
        for num_task,esti in enumerate(estis):
            if esti != 0:
                folder = {
                    'cp2k': Make_CP2K_Files(num_iter=num_iter, num_task=num_task),
                    'lammps': Make_LAMMPS_Files(num_iter=num_iter, num_task=num_task)
                }.get(geo_type, None)
                folder.make_folders()
                stc = Perturb_Structure(previous_infos[str(num_task)].final_pos).random_move(stepsize=step_size)
                folder.make_files(stc=stc, params_dict=geo_params)
                print('The Energy of Iteration # {} Task {} is rejected according to MC criterion'.format(num_iter,num_task))
                print('Start again from previous perturbed structure!')
                run_type = {
                    'cp2k': CP2KRun,
                    'lammps': LAMMPSRun,
                }.get(geo_type, None)
                task_lst.append(run_type(machine_dict).gen_task(num_iter=num_iter,num_task=num_task))
        if task_lst != []:
            geo_work = Run(inp_dict=machine_dict)
            geo_work.sub_task(task_lst=task_lst)
        
                        
    print(" --> Finished Iteration # " + str(num_iter))
    num_iter += 1

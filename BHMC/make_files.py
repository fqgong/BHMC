import os, shutil
from glob import glob
import numpy as np
from ase.io import read,write




class Make_Folders:
    def __init__(self,num_iter,num_task):
        self.num_iter = num_iter
        self.num_task = num_task
        self.path = os.getcwd()
        self.folder_path = self.path+'/BHMC-workflow/iter.{}'.format(self.num_iter)+'/task.{}'.format(self.num_task)
    def make_folders(self):
        if os.path.exists(self.folder_path):
            num_folders = len(glob(self.folder_path+'*'))
            shutil.move(self.folder_path,self.folder_path+'-{}'.format(num_folders))
        os.makedirs(self.folder_path)

class Make_CP2K_Files(Make_Folders):

    def __init__(self,num_iter,num_task):
        super().__init__(num_iter,num_task)
        self.default_dict={
            "GLOBAL": {
              "PROJECT": "BH_GEO",
              "RUN_TYPE": "GEO_OPT"
            },
            "FORCE_EVAL": {
              "METHOD": "QS",
              "STRESS_TENSOR": "ANALYTICAL",
              "DFT": {
                "BASIS_SET_FILE_NAME": "BASIS_MOLOPT",
                "POTENTIAL_FILE_NAME": "GTH_POTENTIALS",
                "CHARGE": 0,
                "UKS": "T",
                "MGRID": {
                  "CUTOFF": 400,
                  "REL_CUTOFF": 50,
                  "NGRIDS": 4
                },
                "QS": {
                  "EPS_DEFAULT": "1.0E-12"
                },
                "SCF": {
                  "SCF_GUESS": "RESTART",
                  "EPS_SCF": "1.0E-6",
                  "MAX_SCF": 50
                },
                "XC": {
                  "XC_FUNCTIONAL": {
                    "_": "PBE"
                  }
                }
              },
              "SUBSYS": {
                "CELL":{
                  "A": "10 .0 .0",
                  "B": ".0 10 .0",
                  "C": ".0 .0 10"
                },
                "COORD": {"@include": "coord.xyz"},
                "KIND": {
                  "_": ["H","C","N"],
                  "POTENTIAL": ["GTH-PBE-q1","GTH-PBE-q4", "GTH-PBE-q5"],
                  "BASIS_SET": ["DZVP-MOLOPT-GTH","DZVP-MOLOPT-GTH","DZVP-MOLOPT-GTH"]
                }
              }
            },
            "MOTION":{
              "GEO_OPT":{
                "TYPE": "MINIMIZATION"
              },
              "PRINT":{
                "TRAJECTORY":{
                  "_": "ON"
                },
                "FORCES":{
                  "_": "ON"
                },
                "VELOCITIES":{
                  "_": "ON"
                },
                "RESTART":{
                  "BACKUP_COPIES": 3,
                  "EACH":{ 
                    "MD": 1  
                  }
                }
              }
            }
        }


    def update_dict(self, old_d, update_d):
        """
        a method to recursive update dict
        :old_d: old dictionary
        :update_d: some update value written in dictionary form
        """
        import collections
        for k, v in update_d.items():
            if (k in old_d and isinstance(old_d[k], dict)
                and isinstance(update_d[k], collections.Mapping)):
                self.update_dict(old_d[k], update_d[k])
            else:
                old_d[k] = update_d[k]

    def iterdict(self, d, out_list, flag=None):
        """
        :doc: a recursive expansion of dictionary into cp2k input
        :k: current key
        :v: current value
        :d: current dictionary under expansion
        :flag: used to record dictionary state. if flag is None,
        it means we are in top level dict. flag is a string.
        """
        for k,v in d.items():
            k=str(k) # cast key into string
            #if value is dictionary
            if isinstance(v, dict):
                # flag == None, it is now in top level section of cp2k
                if flag==None :
                    out_list.append("&"+k)
                    out_list.append("&END "+k)
                    self.iterdict(v, out_list, k)
                # flag is not None, now it has name of section
                else:
                    index = out_list.index("&END " + flag)
                    out_list.insert(index, "&"+k)
                    out_list.insert(index+1,"&END "+k )
                    self.iterdict(v, out_list, k)
            elif isinstance(v, list):
    #            print("we have encountered the repeat section!")
                index = out_list.index("&"+flag)
                # delete the current constructed repeat section
                del out_list[index:index+2]
                # do a loop over key and corresponding list
                k_tmp_list = []
                v_list_tmp_list = []
                for k_tmp, v_tmp in d.items():
                    k_tmp_list.append(str(k_tmp))
                    v_list_tmp_list.append(v_tmp)
                for repeat_keyword in zip(*v_list_tmp_list):
                    out_list.insert(index,"&" + flag)
                    out_list.insert(index+1, "&END " + flag)
                    for idx, k_tmp in enumerate(k_tmp_list):
                        if k_tmp == "_":
                            out_list[index] = "&" + flag + " " + repeat_keyword[idx]
                        else:
                            out_list.insert(index+1, k_tmp+" "+repeat_keyword[idx])

                break
            else:
                v=str(v)
                if flag==None :
                    out_list.append(k+" "+v)
                    print (k,":",v)
                else:
                    if k == "_":
                        index = out_list.index("&" + flag)
                        out_list[index] = ("&" + flag + " " + v)

                    else:
                        index = out_list.index("&END "+flag)
                        out_list.insert(index, k+" "+v)

    def make_cp2k_input(self, sys_data, fp_params):
        #covert cell to cell string
        cell = sys_data.todict()['cell']
        cell_a = np.array2string(cell[0,:])
        cell_a = cell_a[1:-1]
        cell_b = np.array2string(cell[1,:])
        cell_b = cell_b[1:-1]
        cell_c = np.array2string(cell[2,:])
        cell_c = cell_c[1:-1]

        #get update from user
        user_config=fp_params
        #get update from cell
        cell_config={"FORCE_EVAL":{
            "SUBSYS":{
                "CELL":{
                    "A": cell_a,
                    "B": cell_b,
                    "C": cell_c
                    }
                }
            }
                }
        self.update_dict(self.default_dict, user_config)
        self.update_dict(self.default_dict, cell_config)
        #output list
        input_str = []
        self.iterdict(self.default_dict, input_str)
        string="\n".join(input_str)
        return string

    def make_cp2k_xyz(self, sys_data):
        #get atomic species
        atom_list = sys_data.get_chemical_symbols()

        #write coordinate to xyz file used by cp2k input
        coord_list = sys_data.todict()['positions']
        x = '\n'
        for kind, coord in zip(atom_list, coord_list) :
            x += str(kind) + ' ' + str(coord[:])[1:-1] + '\n'
        return x

    def make_cp2k_input_from_external(self, sys_data, exinput_path):
        # read the input content as string
        with open(exinput_path, 'r') as f:
            exinput = f.readlines()

        # find the ABC cell string
        for line_idx, line in enumerate(exinput):
            if 'ABC' in line:
                delete_cell_idx = line_idx
                delete_cell_line = line

        # remove the useless CELL line
        exinput.remove(delete_cell_line)

        # insert the cell information
        # covert cell to cell string
        cell = sys_data['cells'][0]
        cell = np.reshape(cell, [3,3])
        cell_a = np.array2string(cell[0,:])
        cell_a = cell_a[1:-1]
        cell_b = np.array2string(cell[1,:])
        cell_b = cell_b[1:-1]
        cell_c = np.array2string(cell[2,:])
        cell_c = cell_c[1:-1]

        exinput.insert(delete_cell_idx, 'A  ' + cell_a + '\n')
        exinput.insert(delete_cell_idx+1, 'B  ' + cell_b + '\n')
        exinput.insert(delete_cell_idx+2, 'C  ' + cell_c + '\n')

        return ''.join(exinput)

    def make_files(self, stc, params_dict, input_path=None):
        if params_dict:
            cp2k_input = self.make_cp2k_input(stc,params_dict)
        else:
            cp2k_input = self.make_cp2k_input_from_external(stc,input_path)
        with open(self.folder_path+'/input.inp','w') as f:
            f.write(cp2k_input)
        f.close()
        cp2k_coord = self.make_cp2k_xyz(stc)
        with open(self.folder_path+'/coord.xyz','w') as f:
            f.write(cp2k_coord)
        f.close()

class Make_LAMMPS_Files(Make_Folders):
    def __init__(self,num_iter,num_task):
        super().__init__(num_iter,num_task)

    def make_lammps_input(self, stc, lammps_dict):
        nstep = lammps_dict.get('nstep',1)
        trj_freq = lammps_dict.get('trj_freq',1)
        ener_tol = lammps_dict.get('ener_tol',0)
        force_tol = lammps_dict.get('force_tol',1.0e-3)
        max_iter = lammps_dict.get('max_iter',1000)
        max_eval = lammps_dict.get('max_eval',100000)
        eles = lammps_dict.get('eles',['Cu','C','O'])
        mass_map = lammps_dict.get('mass_map',[63.546,12.011,15.999])
        nopbc = lammps_dict.get('nopbc',False)
        ele_list = ""
        for i in eles:
            ele_list += i + " "
        ret = "variable      NSTEPS      equal      %d\n" % nstep
        ret += "variable      THERMO_FREQ      equal      %d\n" % trj_freq
        ret += "variable      DUMP_FREQ      equal      %d\n" % trj_freq
        ret += "variable      ENERGY_TOL      equal      %e\n" % ener_tol
        ret += "variable      FORCE_TOL      equal      %e\n" % force_tol
        ret += "variable      MAX_ITER      equal      %d\n" % max_iter
        ret += "variable      MAX_EVAL      equal      %d\n" % max_eval
        ret += "\n"
        ret += "#初始化模型设置\n"
        ret += "units      metal\n"
        ret += "atom_style      atomic\n"
        ret += "atom_modify      map yes\n"
        if nopbc:
            ret += "boundary      f f f\n"
        else:
            ret += "boundary      p p p\n"
        ret += "neighbor      1.0 bin\n"
        ret += "\n"
        ret += "#读取模型数据\n"
        ret += "read_data      conf.lmp\n"
        for jj in range(len(mass_map)):
            ret += "mass      %d %f\n" % (jj+1,mass_map[jj])
        ret += "\n"
        ret += "#定义原子间相互作用势\n"
        ret += "pair_style      deepmd frozen_model.pb\n"
        ret += "pair_coeff      * *\n"
        ret += "#定义输出热力学信息\n"
        ret += "thermo_style      custom step temp pe ke etotal\n"
        ret += "thermo      ${THERMO_FREQ}\n"
        ret += "\n"
        ret += "#定义输出的轨迹\n"
        ret += "dump      1 all xyz ${DUMP_FREQ} BH_GEO-pos-1.xyz\n"
        ret += "dump_modify      1 sort id\n"
        ret += "dump_modify      1 element %s\n" % ele_list
        ret += "\n"
        ret += "#结构优化\n"
        ret += "minimize      ${ENERGY_TOL} ${FORCE_TOL} ${MAX_ITER} ${MAX_EVAL}\n"
        ret += "min_style      cg\n"
        return ret

    def make_files(self, stc, params_dict, input_path=None):
        write(os.path.join(self.folder_path,'conf.lmp'),stc,format='lammps-data')
        os.symlink(params_dict['model_path'],os.path.join(self.folder_path,'frozen_model.pb'))
        if params_dict:
            lammps_input = self.make_lammps_input(stc,params_dict)
            with open(os.path.join(self.folder_path,'input.lammps'),'w') as f:
                f.write(lammps_input)
        else:
            shutil.copyfile(os.path.abspath(input_path),os.path.join(self.folder_path,'input.lammps'))

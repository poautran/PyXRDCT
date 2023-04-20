#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#    Project: PyXRDCT
#             https://github.com/poautran/PyXRDCT
#
#    Copyright (C) 2022-2023 European Synchrotron Radiation Facility, Grenoble,
#             France
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in
# all copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NON INFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
# THE SOFTWARE.

import os
import subprocess

class Slurm:
    """
    Init checks if slurm is available
    """
    def __init__(self):
        self.slurm = True
        try:
            subprocess.call(["squeue"])
            print("[INFO] Integration will run on SLURM")
        except:
            print("[WARNING] System is not compatible with SLURM allocation. Integration will run locally")
            self.slurm = False
            
    def write():
        with open('slurm.job',"w") as fh:
            fh.writelines("#!/bin/bash\n")
            fh.writelines("#SBATCH --exclude=asdhpc6-[01-08],gpbm18-[01-06],gpid16a-1802,hib3-[3002-3004,3101-3104],hpc3-[2101-2102,2603-2604,2701-2703,2801-2804,2901-2902]\n")
            fh.writelines("#SBATCH --partition=nice\n")
            fh.writelines("#SBATCH --nodes=1\n")
            fh.writelines("#SBATCH --ntasks=1\n")
            fh.writelines("#SBATCH --mem=300GB\n") 
            fh.writelines("#SBATCH --cpus-per-task=40\n")                                  
            fh.writelines("#SBATCH --job-name=NabuRec\n")
            fh.writelines('echo "[INFO] Job started at $(date) on $(hostname)"\n')
            fh.writelines("# Job steps\n")
            fh.writelines('echo "[INFO] Job step: Computing pyFAI-integration"\n')
            fh.writelines("srun -l nabu python3 integrate.py \n") # --force_use_grouped_pipeline 1
            fh.writelines('echo "[INFO] Job ended at $(date)"\n')
            
    def run():
        os.system("sbatch slurm.job")
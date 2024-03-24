import os, sys
import numpy as np

script_base = os.path.dirname(os.path.realpath(__file__))
tasks = ["TIMIT", "MOCHA"]
queue  = "mind-batch"
for id_start in np.arange(1, 311, 10):
	id_end = id_start + 10
	job = f"/userdata/sjain/myenv/bin/python3 {script_base}/process_ecog_data.py --id_start {id_start} --id_end {id_end} --tasks {' '.join(tasks)}"
	job = f"submit_job -q {queue} -e shailee.jain@ucsf.edu -c 8 -m 80 -o {script_base}/job_scripts/{id_start}_{id_end}_{'-'.join(tasks)}.txt -n ep{id_start}_{id_end} -x {job}"
	os.system(job)
import os, sys
import numpy as np

script_base = os.path.dirname(os.path.realpath(__file__))
tasks = ["interview", "BURSC"]
id_starts_queues = [
	[np.arange(1, 201, 10), "mind-batch"],
	[np.arange(201, 241, 10), "spirit-batch"],
	[np.arange(241, 271, 10), "skull-batch.q"],
	[np.arange(271, 311, 10), "pia-batch.q"],
]
for id_starts, queue in id_starts_queues:
	for id_start in id_starts:
		id_end = id_start + 10
		job = f"/userdata/sjain/myenv/bin/python3 {script_base}/process_ecog_data.py --id_start {id_start} --id_end {id_end} --tasks {' '.join(tasks)}"
		job = f"submit_job -q {queue} -e shailee.jain@ucsf.edu -c 8 -m 80 -o {script_base}/job_scripts/{id_start}_{id_end}_{'-'.join(tasks)}.txt -n ep{id_start}_{id_end} -x {job}"
		os.system(job)
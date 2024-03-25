import os, sys
import csv
from collections import defaultdict
from datetime import datetime
import json
import argparse

sys.path.insert(0, "/userdata/sjain/ECoG_data_collection/")
from ecog_preprocessing.transformData import transform as transformData
sys.path.insert(0, "/userdata/sjain/ECoG_data_collection/stimulus_align/src")
import preanalysis

if __name__ == "__main__":
	parser = argparse.ArgumentParser()
	parser.add_argument('--id_start', type=int, required=True,
						help="Subject ID with which preprocessing should start.")
	parser.add_argument('--id_end', type=int, required=True,
						help="Subject ID before which preprocessing should end.")
	parser.add_argument('--tasks', nargs="+",
						help="Tasks that should be preprocessed.")
	args = parser.parse_args()
	globals().update(args.__dict__)
	
	datapath = "/data_store1/human/prcsd_data/"
	cur_dir = os.path.dirname(os.path.realpath(__file__))
	savepath = f"{cur_dir}/data/preprocessed"
	print(tasks)

	with open(f"{cur_dir}/task_database_logging/task_database_Info_20240313.csv", newline="") as csvfile:
		reader = csv.DictReader(csvfile)
		rows = [row for row in reader]
	print("Read CSV.")

	# Read each entry from our task database.
	for row in rows:
		subj, block, task, valid = row["SubjectID"], row["Block"], row["Task"], row["Valid"]
		subj_id = int(subj[2:])
		blockname = f"{subj}_B{block}"

		if valid.lower()=="bad" or task not in tasks or subj_id<id_start or subj_id>=id_end:
			print(f"Bad params! Skipping blockname {blockname}.\n")
			continue

		blockpath = os.path.join(datapath, subj, blockname)
		storepath = os.path.join(savepath, subj, blockname)

		d1_path = os.path.join(storepath, f"{blockname}_nocar_Hilb.h5")
		d2_path = os.path.join(storepath, "ecog400/ecog.mat")
		if os.path.exists(d1_path) and os.path.exists(d2_path):
			print(f"{blockname} was processed already!")
			if not os.path.exists(d1_path.replace("_nocar", "_100Hz_nocar")):
				allband = preanalysis.preanalyze_data_h5(subj, block, 100, return_ds=True)
			continue
		try:
			print(f"Processing blockname {blockname}.\n")
			transformData(blockpath, rate=400, car=False, line_noise_notch=True, force=True, storepath=storepath)
			preanalysis.preanalyze_data_h5(subj, block, 100, storepath)
			# If able to preprocess, record the entire `row` as good.
			row["prcsd_date"] = datetime.today().strftime('%Y-%m-%d %H:%M:%S')
			row["prcsd_loc"] = storepath
			with open(f"{storepath}/preproc_log.json", "w") as f:
				json.dump(row, f)
		except FileNotFoundError:
			print(f"Didn't find data. Skipping blockname {blockname}.\n")
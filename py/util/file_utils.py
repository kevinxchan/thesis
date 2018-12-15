import os

def list_dir_abs(dir):
	for f in os.listdir(dir):
		yield os.path.abspath(os.path.join(dir, f))
from oceantracker.post_processing.read_output_files import load_output_files 
import matplotlib.pyplot as plt
from datetime import datetime
import numpy as np

path_to_dir = '/scratch/local1/output/22_07_27_comp_cost_particle_scaling_v02'
cases = load_output_files.get_case_info_files_from_dir(path_to_dir)

df = np.zeros((len(cases),2))

for ii,item in enumerate(cases):
    case_info = load_output_files.read_case_info_file(item)
    name = case_info['output_files']['output_file_base'][-4:]

    n_particles = case_info['particle_release_groups'][0]['pulse_size']
    case_run_time = case_info['run_info']['model_run_duration']
    case_run_time = datetime.strptime(case_run_time,"%H:%M:%S.%f").timestamp() - datetime(1900,1,1).timestamp()

    df[ii,0] = n_particles
    df[ii,1] = case_run_time

print(df)

plt.scatter(df[:,0],df[:,1]-min(df[:,1]))
plt.savefig(path_to_dir.split('/')[-1]+'.png')
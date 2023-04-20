from oceantracker import main
from oceantracker.util import json_util

if __name__ == "__main__":
    ncase=0

    if ncase ==0:
        fn = r'E:\OneDrive - Cawthron\H_Local_drive\ParticleTracking\bug_hunting\paramsBWDFY20120.json'
        params: object = json_util.read_JSON(fn)
        params['reader']['input_dir'] = r'G:\Hindcasts_large\OceanNumNZ-2022-06-20\final_version\2012'
        params['shared_params']['root_output_dir'] = r'F:\OceanTrackerOuput\bug_hunting'
        #params['case_list'] =  params['case_list'][23:26]

        # short run
        params['case_list'] = [params['case_list'][25]]
        params['shared_params']['processors'] = 1
        #params['base_case_params']['solver']['n_substeps'] = 10
        #params['case_list'][0]['particle_release_groups'][0]['pulse_size'] = 100


    params['shared_params']['max_duration'] = 5*24*3600

    main.run(params)
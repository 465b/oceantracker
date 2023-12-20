from oceantracker.util.ncdf_util import NetCDFhandler
from glob import  glob
from os import path
from oceantracker.common_info_default_param_dict_templates import known_readers
from pathlib import Path as pathlib_Path
from oceantracker.common_info_default_param_dict_templates import known_readers
from oceantracker.util.parameter_util import make_class_instance_from_params
from copy import deepcopy

def find_file_format_and_file_list(reader_params, msg_logger):
    found=False
    # first see if it matches known formats
    for r_name, r in known_readers.items():
        params= deepcopy(reader_params)
        params['class_name'] = r
        reader = make_class_instance_from_params('reader', params, msg_logger, default_classID='reader',check_for_unknown_keys=False)
        file_list = reader.get_file_list()
        if len(file_list) > 0 and reader.is_file_format(file_list[0]):
            # found format
            if 'class_name' not in reader_params: reader_params['class_name'] = r # dont overwrite user given class name
            found = True
            break
    if found:
        msg_logger.progress_marker('found hydro-model files of type ' + r_name.upper())
    else:
        msg_logger.msg(f'Could not set up reader, no files found matching mask = "{reader_params["file_mask"]}"  or "out2d*.nc" for schism v5, or files do no match known format',
                                   fatal_error=True,  exit_now=True)
    return reader_params, file_list





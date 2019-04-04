#!/usr/bin/env python
"""Run the modification pipeline for r9 reads using deepmod"""
########################################################################
# File: run_deep_mod_r9.py
#
# Authors: Andrew Bailey
#
# History: 04/03/19
########################################################################

import shutil
import subprocess
from argparse import ArgumentParser
from py3helpers.multiprocess import *
from py3helpers.utils import *
from py3helpers.aws import *
from run_mod_detect_pipeline import call_deep_mod_detect


def parse_args():
    parser = ArgumentParser(description=__doc__)
    # parsers for running the full pipeline
    # required arguments
    parser.add_argument('--config', '-c', action='store',
                        dest='config', required=True, type=str, default=None,
                        help="Path to config file")

    args = parser.parse_args()
    return args


class DeepModPipeline2(object):
    """Run full DeepmodPipeline for one file"""
    def __init__(self, s3_file, n_threads, config_args):
        """Make directory structure and run analysis:

        working_dir
            - self.temp_dir
                - s3 downloaded_file
                - self.fast5_output_dir
            - self.deep_mod_output_dir

        """
        self.config_args = config_args
        self.s3_handler = AwsS3()
        self.s3_file = s3_file
        self.n_threads = n_threads
        self.working_dir = self.config_args.working_dir
        self.s3_save_bucket = self.config_args.s3_save_bucket
        self.file_id = os.path.basename(os.path.splitext(self.s3_file)[0])

        # setup working directories
        if not os.path.exists(self.working_dir):
            os.mkdir(self.working_dir)

        self.deep_mod_output_dir = os.path.join(self.working_dir, "deep_mod_"+self.file_id)
        os.mkdir(self.deep_mod_output_dir)

        self.temp_dir = os.path.join(self.working_dir, self.file_id)
        os.mkdir(self.temp_dir)
        self.fast5_output_dir = os.path.join(self.temp_dir, "edited_fast5s")
        os.mkdir(self.fast5_output_dir)

    def get_multi_fast5(self):
        return self.s3_handler.download_object(self.s3_file, self.temp_dir)

    def get_and_untar_gz_fast5(self, remove_tar=True):
        tar_gz_file = self.s3_handler.download_object(self.s3_file, self.temp_dir)
        out_path = untar_gz(tar_gz_file, self.temp_dir)
        fast5_dir = find_fast5_dir(out_path)
        if remove_tar:
            os.remove(tar_gz_file)
        return fast5_dir

    def run_r9(self, delete=True, debug=True):
        fast5_dir = self.get_and_untar_gz_fast5()

        if debug:
            sucess, time_in_s = time_it(self._run_r9, fast5_dir)
            print("{} took {} seconds".format(self.file_id, time_in_s))
        else:
            try:
                sucess, time_in_s = time_it(self._run_r9, fast5_dir)
                print("{} took {} seconds".format(self.file_id, time_in_s))
            except Exception as e:
                print("[DeepModPipeline] ERROR exception ({}): {}".format(type(e), e))
        if delete:
            shutil.rmtree(self.temp_dir)

    def _run_r9(self, fast5_dir):
        """Run the pipeline"""
        start = timer()
        success = call_deep_mod_detect(python_path=self.config_args.python_path,
                                       executable=self.config_args.deep_mod_executable,
                                       fast5_dir=fast5_dir,
                                       reference=self.config_args.reference,
                                       output_dir=self.deep_mod_output_dir,
                                       base=self.config_args.base,
                                       modfile=self.config_args.modfile,
                                       file_id=self.file_id,
                                       threads=self.n_threads)
        stop = timer()
        print("[run_deepmod_on_multifast5] call_deep_mod_detect running time = {} seconds".format(stop - start))
        return success


def find_fast5_dir(top_directory, ext="fast5"):
    """Find the subdirectory with the fast5 files. Error if more than one directory has fast5 files

    :param top_directory: path to top directory
    """
    return_dir = None
    n_dirs_with_files = 0
    for root, dirs, files in os.walk(top_directory):
        if len([x for x in files if x.endswith(ext)]) > 0:
            n_dirs_with_files += 1
            return_dir = root

    assert n_dirs_with_files == 1, "Need to check directory structure. More than one dir of {} files".format(ext)
    return return_dir


def deep_mod_pipeline_wrapper(fast5_file, arguments):
    """wrap DeepModPipeline into a single function call for ease of multiprocessing"""
    dmp = DeepModPipeline2(fast5_file, arguments["n_threads"], create_dot_dict(arguments))
    dmp.run_r9(delete=True)
    return True


def multiprocess_r9_deepmod_pipeline(args):
    """Multiprocess deepmod pipeline"""
    service = BasicService(deep_mod_pipeline_wrapper)
    s3_handler = AwsS3()
    fast5_files = s3_handler.list_objects(args.s3_fast5s)
    arguments = merge_dicts([args, {"n_threads": max(1, args.n_threads // 4)}])

    total, failure, messages, output = run_service(service.run, fast5_files, {"arguments": arguments}, ["fast5_file"], 4)

    return output


def multiprocess_deepmod_pipeline(args):
    """Multiprocess deepmod pipeline"""
    service = BasicService(deep_mod_pipeline_wrapper)
    s3_handler = AwsS3()
    fast5_files = s3_handler.list_objects(args.s3_tars)
    arguments = merge_dicts([args, {"n_threads": max(1, args.n_threads // 4)}])

    total, failure, messages, output = run_service(service.run, fast5_files, {"arguments": arguments}, ["fast5_file"], 4)

    return output


def main():
    # parse args
    start = timer()

    args = parse_args()
    args = create_dot_dict(load_json(args.config))
    multiprocess_deepmod_pipeline(args)

    stop = timer()
    print("Running Time = {} seconds".format(stop - start))


if __name__ == '__main__':
    main()
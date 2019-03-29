#!/usr/bin/env python
"""Run the modification pipeline using deepmod"""
########################################################################
# File: run_mod_detect_pipeline.py
#
# Authors: Andrew Bailey
#
# History: 03/27/19
########################################################################

import shutil
import subprocess
from argparse import ArgumentParser
from py3helpers.multiprocess import *
from py3helpers.utils import *
from py3helpers.aws import *


def parse_args():
    parser = ArgumentParser(description=__doc__)
    # parsers for running the full pipeline
    # required arguments
    parser.add_argument('--config', '-c', action='store',
                        dest='config', required=True, type=str, default=None,
                        help="Path to config file")

    args = parser.parse_args()
    return args


def subprocess_call(command, name):
    try:
        print("[{}] running :'{}'".format(name, command))

        command = command.split()
        proc = subprocess.Popen(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE)

        output, errors = proc.communicate()
        errors = errors.decode().splitlines()
        for x in errors:
            print(x)
        output = output.decode().splitlines()
        for x in output:
            print(x)

    except Exception as e:
        print("[{}] ERROR exception ({}): {}".format(name, type(e), e))
        raise e

    return True


def call_run_embed_fast5(embed_fast5_executable, fast5, main_cpp_dir, nanopolish_dir, output_dir, threads=1):
    """Call run_embed_fast5 on a multi_fast5 file
    :param embed_fast5_executable: path to run_embed_fast5.py
    :param fast5: multi_fast5 to process
    :param main_cpp_dir: path to main_cpp directory
    :param nanopolish_dir: path to nanopolish directory
    :param output_dir: path to write all new fast5s
    :param threads: number of jobs to run while processing
    :return: bool
    """
    assert os.path.exists(embed_fast5_executable), "embed_fast5_executable does not exist: {}".format(
        embed_fast5_executable)
    assert os.path.exists(main_cpp_dir), "main_cpp_dir does not exist: {}".format(main_cpp_dir)
    assert os.path.exists(nanopolish_dir), "nanopolish_dir does not exist: {}".format(nanopolish_dir)
    assert os.path.exists(output_dir), "output_dir does not exist: {}".format(output_dir)
    assert os.path.exists(fast5), "fast5 does not exist: {}".format(fast5)

    embed_fast5_command = embed_fast5_executable + " -f {fast5} -o {output_dir} -m {main_cpp_dir}" \
                                                   " -n {nanopolish_dir} -j {jobs}".format(fast5=fast5,
                                                                                           main_cpp_dir=main_cpp_dir,
                                                                                           nanopolish_dir=nanopolish_dir,
                                                                                           output_dir=output_dir,
                                                                                           jobs=threads)
    success = subprocess_call(embed_fast5_command, "call_run_embed_fast5")
    return success


def call_deep_mod_detect(python_path, executable, fast5_dir, reference, output_dir, base, modfile, file_id, threads):
    """
    Call DeepMod detect 
    :param python_path: path to correct python
    :param executable: path to DeepMod.py
    :param fast5_dir: path to fast5 directory 
    :param reference: path to reference
    :param output_dir: path to output folder
    :param base: base to call
    :param modfile: path to model file
    :param file_id: some id
    :param threads: number of threads to run
    """
    assert os.path.exists(python_path), "deepmod path does not exist: {}".format(python_path)
    assert os.path.exists(executable), "deepmod path does not exist: {}".format(executable)
    assert os.path.exists(fast5_dir), "fast5_dir path does not exist: {}".format(fast5_dir)
    assert os.path.exists(reference), "reference path does not exist: {}".format(reference)
    assert os.path.exists(output_dir), "output_dir path does not exist: {}".format(output_dir)
    assert os.path.exists(modfile+".meta"), "modfile path does not exist: {}".format(modfile+".meta")
    assert os.path.exists(modfile+".index"), "modfile path does not exist: {}".format(modfile+".index")

    deep_mod_command = python_path + " " + executable + " detect --wrkBase {fast5_dir} --Ref {reference} " \
                                                        "--outFolder {output_dir}" \
                                                        " --Base {base} --modfile {modfile} " \
                                                        "--FileID {file_id} " \
                                                        "--threads {threads}".format(fast5_dir=fast5_dir,
                                                                                     reference=reference,
                                                                                     output_dir=output_dir,
                                                                                     modfile=modfile,
                                                                                     base=base,
                                                                                     threads=threads,
                                                                                     file_id=file_id)
    success = subprocess_call(deep_mod_command, "call_deep_mod_detect")
    return success


def run_deepmod_on_multifast5(multi_fast5, embed_fast5_executable, main_cpp_dir, nanopolish_dir, fast5_output_dir,
                              python_path, deep_mod_executable, deep_mod_output, reference, base, modfile, file_id,
                              threads):
    """Running deepmod with a multi_fast5 input file
    :param deep_mod_output: path to where deepmod output will be placed
    :param multi_fast5: path to multi fast5
    :param embed_fast5_executable: path to run_embed_fast5.py
    :param main_cpp_dir: path to embed build directory
    :param nanopolish_dir: path to nanopolish directory
    :param fast5_output_dir: path for all fast5s to be placed
    :param python_path: path to the python to use for deepmod
    :param deep_mod_executable: path to DeepMod.py
    :param reference: path to reference file
    :param base: nucleotide to call via DeepMod
    :param modfile: model file
    :param file_id: file id for deepmod
    :param threads: n threads to run for deepmod
    """
    start = timer()

    if not os.path.exists(fast5_output_dir):
        os.mkdir(fast5_output_dir)
    if not os.path.exists(deep_mod_output):
        os.mkdir(deep_mod_output)
    # Extract reads and embed using nanopolish
    success = call_run_embed_fast5(embed_fast5_executable, multi_fast5, main_cpp_dir, nanopolish_dir,
                                   fast5_output_dir, threads)
    stop = timer()
    print("[run_deepmod_on_multifast5] call_run_embed_fast5 running time = {} seconds".format(stop - start))

    assert success, "[run_deepmod_on_multifast5] Something Failed in call_run_embed_fast5"
    success = call_deep_mod_detect(python_path=python_path, executable=deep_mod_executable,
                                   fast5_dir=fast5_output_dir, reference=reference,
                                   output_dir=deep_mod_output, base=base,
                                   modfile=modfile, file_id=file_id, threads=threads)
    stop = timer()
    print("[run_deepmod_on_multifast5] call_deep_mod_detect running time = {} seconds".format(stop - start))
    assert success, "[run_deepmod_on_multifast5] Something Failed in call_deep_mod_detect"
    return True


class DeepModPipeline(object):
    """Run full DeepmodPipeline for one file"""
    def __init__(self, s3_file, s3_handler, working_dir, s3_save_bucket, config_args):
        """Make directory structure and run analysis:

        working_dir
            - self.temp_dir
                - s3 downloaded_file
                - self.fast5_output_dir
            - self.deep_mod_output_dir

        """
        self.config_args = config_args
        self.s3_save_bucket = s3_save_bucket
        self.s3_handler = s3_handler
        self.s3_file = s3_file
        self.file_id = os.path.basename(os.path.splitext(self.s3_file)[0])

        # setup working directories
        if not os.path.exists(working_dir):
            os.mkdir(working_dir)

        self.deep_mod_output_dir = os.path.join(working_dir, "deep_mod_"+self.file_id)
        os.mkdir(self.deep_mod_output_dir)

        self.temp_dir = os.path.join(working_dir, self.file_id)
        os.mkdir(self.temp_dir)
        self.fast5_output_dir = os.path.join(self.temp_dir, "edited_fast5s")
        os.mkdir(self.fast5_output_dir)

    def get_file(self):
        return self.s3_handler.download_object(self.s3_file, self.temp_dir)

    def _run(self):
        """Run the pipeline"""
        multi_fast5 = self.get_file()
        run_deepmod_on_multifast5(multi_fast5=multi_fast5,
                                  embed_fast5_executable=self.config_args.embed_fast5_executable,
                                  main_cpp_dir=self.config_args.main_cpp_dir,
                                  nanopolish_dir=self.config_args.nanopolish_dir,
                                  fast5_output_dir=self.fast5_output_dir,
                                  python_path=self.config_args.python_path,
                                  deep_mod_executable=self.config_args.deep_mod_executable,
                                  deep_mod_output=self.deep_mod_output_dir, reference=self.config_args.reference,
                                  base=self.config_args.base, modfile=self.config_args.modfile, file_id=self.file_id,
                                  threads=self.config_args.n_threads)

        tar_path = os.path.join(self.temp_dir, self.file_id+"_{}.tar.gz".format(len(list_dir(self.fast5_output_dir))-5))
        tar_file = tar_gz(self.fast5_output_dir, tar_path)
        self.s3_handler.upload_object(tar_file, self.s3_save_bucket)
        return True

    def run(self, delete=False, debug=True):
        if debug:
            sucess, time_in_s = time_it(self._run)
            print("{} took {} seconds".format(self.file_id, time_in_s))
        else:
            try:
                sucess, time_in_s = time_it(self._run)
                print("{} took {} seconds".format(self.file_id, time_in_s))
            except Exception as e:
                print("[DeepModPipeline] ERROR exception ({}): {}".format(type(e), e))
        if delete:
            shutil.rmtree(self.temp_dir)


def main():
    # parse args
    start = timer()

    args = parse_args()
    args = create_dot_dict(load_json(args.config))

    s3_handler = AwsS3()
    fast5_files = s3_handler.list_objects(args.s3_fast5s)
    print(fast5_files)

    dmp = DeepModPipeline(fast5_files[0], s3_handler, args.working_dir, args.s3_save_bucket, args)
    dmp.run()
    stop = timer()
    print("Running Time = {} seconds".format(stop - start))


if __name__ == '__main__':
    main()
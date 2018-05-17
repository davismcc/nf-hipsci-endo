# IMPORT GENERAL STUFF
import os
import scipy as SP
import h5py
import sys
import pdb
from optparse import OptionParser

def create_group(group_name, njobs=8000):
    print ""
    command = "bgadd /%s" % group_name
    print command
    os.system(command)
    command = "bgmod -L %d /%s" % (njobs, group_name)
    print command
    os.system(command)
    print ""

if __name__ == "__main__":
    ## By doing the main check, you can have this code only execute when you
    ## want to run the module as a program and not have it execute when someone
    ## just wants to import your module and call your functions themselves.
    parser = OptionParser()
    parser.add_option("--n_jobs", dest='n_jobs', type=int, default=100)
    parser.add_option("--from_job", dest='from_job', type=int, default=0)
    parser.add_option("--to_job", dest='to_job', type=int, default=None)
    parser.add_option("--peer", action="store_true", dest='peer', default=False)
    parser.add_option("--perm", action="store_true", dest='perm', default=False)
    (opt, args) = parser.parse_args()
    opt_dict = vars(opt)

    if opt.to_job is None:
        opt.to_job = opt.n_jobs

    # create group
    group_name = 'hipsci_eqtl_trans'
    if opt.peer:    group_name += '_peer'
    if opt.perm:    group_name += '_perm'
    create_group(group_name, opt.n_jobs)

    #Create temp dir
    temp_folder_base = './../temp/%s' % group_name

    pdb.set_trace()

    # GO! GO! GO!
    for j in range(opt.from_job, opt.to_job):
        # Create a temp folder for jobs output, 1000 jobs per temp folder
        temp_folder = os.path.join(temp_folder_base, str(int(SP.ceil(j/1000))))
        if not os.path.exists(temp_folder): os.makedirs(temp_folder)
        stdout_file = os.path.join(temp_folder,
                                   'stdout_%d_%d.txt' % (opt.n_jobs, j))
        stderr_file = os.path.join(temp_folder,
                                   'stderr_%d_%d.txt' % (opt.n_jobs, j))
        command  = "bsub -g /%s " % group_name
        command += "-o %s " % stdout_file
        command += "-e %s " % stderr_file
        command += "python eqtl_cis.py "
        command += "--n_jobs %d --job_i %d" % (opt.n_jobs, j)
        if opt.peer:    command += ' --peer'
        if opt.perm:    command += ' --perm'
        print command
        os.system(command)


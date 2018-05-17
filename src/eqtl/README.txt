This folder contains all scripts used for the eqtl mapping

eqtl_cis.py
    This script contains the cis eqtl analysis.
    It can be run in debug mode using:
    ## python eqtl_cis.py --debug
    When run in normal mode, one can split the whole
    analysis in n parts and considers part i by typing:
    # python eqtl_cis.py --n_jobs n --job_i i
    You can run with permuted data by setting
    # python eqtl_cis.py --n_jobs n --job_i i --perm

eqtl_cis_run.py
    Run the eqtl_cis.py analysis on the cluster.
    One needs to specify the number of jobs the analysis
    should be split into and which one to run form:
    # python eqtl_cis_run.py --n_jobs n_jobs --from_job job_i --to_job job_j
    (default values of job_i and job_j are 0 and n_jobs respectively)
    To run permuted data one can use
    # python eqtl_cis_run.py --n_jobs n_jobs --from_job job_i --to_job job_j --perm

eqtl_cis_summary.py
    Puts together all the results from the single runs
    and produce a single summary file

eqtl_cis_plot.py
    makes some quality control plot
    from the summary files



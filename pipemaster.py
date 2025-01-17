#!/usr/bin/env python
import os, sys
import gzip
import yaml
import time
import pandas as pd
import subprocess
from optparse import OptionParser

def main():
    parser = OptionParser(usage="usage: %prog [options]")
    parser.add_option("-s","--samplesheet", dest="samplesheet",default=None, help="Samplesheet to be used for the pipeline.")
    parser.add_option("-w","--workdir",dest="workdir",default=None, help="Directory to run the pipeline.")
    parser.add_option("-c","--cachedir",dest="cachedir",default=None,help="Directory of the cache to be used. Default: /Path/to/pipeline/cache")
    parser.add_option("-C","--cacheinit", dest="cacheinit", action = "store_true", default=False, help="Run the cache initiation.")
    parser.add_option("-u","--unlock", dest="unlock", action = "store_true", default=False, help="Unlock the workdir from previous abnormal pipeline.")
    parser.add_option("-S","--silent", dest="silent", action = "store_true", default=False, help="Keep silent, stop sending Slurm email notice.")
    parser.add_option("-D","--dry", dest="dryrun", action = "store_true", default=False, help="Perform snakemake dryrun.")
        
    (options, args) = parser.parse_args()
    
    #######################################################################################
    # Cache mode
    #######################################################################################
    ## Set default cachedir inside the pipeline folder
    pipedir = os.path.dirname(__file__)
    if not options.cachedir:
        options.cachedir = os.path.join(pipedir,'cache')
    elif options.cachedir[0] != '/':
        print('Provided cachedir is not an absolute path. Placing it inside the pipeline dir:')
        print(os.path.join(pipedir, options.cachedir))
        options.cachedir = os.path.join(pipedir, options.cachedir)
        
    ## Do cache if cacheinit provided
    if options.cacheinit:
        print('Cache initiation mode, ignore all options except -c (--cachedir)')
        print('Init cache to %s...'%(options.cachedir))
        from scripts.cache import container,reference
        
        # Init container
        with open(os.path.join(pipedir, 'config', 'container.yaml'), 'r') as infile:
            config = yaml.safe_load(infile)
        print(config)
        config['pipelinedir'] = pipedir
        config['cachedir'] = options.cachedir
        status_cache_1 = container(config)
        if not status_cache_1:
            raise Exception("Container cache failed.")
        
        # Init reference
        with open(os.path.join(pipedir, 'config', 'reference.yaml'), 'r') as infile:
            config = yaml.safe_load(infile)
        config['pipelinedir'] = pipedir
        config['cachedir'] = options.cachedir
        status_cache_2 = reference(config)
        if not status_cache_2:
            raise Exception("Reference cache failed.")
        return
    
    #######################################################################################
    # Unlock mode
    #######################################################################################
    if options.unlock:
        print('Pipeline unlock mode, ignore all options except -w (--workdir)')
        
        if os.path.exists(os.path.join(options.workdir,'.snakemake','locks','0.input.lock')):
            os.remove(os.path.join(options.workdir,'.snakemake','locks','0.input.lock'))
        if os.path.exists(os.path.join(options.workdir,'.snakemake','locks','0.output.lock')):
            os.remove(os.path.join(options.workdir,'.snakemake','locks','0.output.lock'))
            
        print('Pipeline directory:',pipedir, 'unlocked!')
        return
    
    #######################################################################################
    # Prepare workdir
    #######################################################################################
    snapshot = 'run_%s'%(time.strftime("%Y%m%d%H%M%S"))
    from scripts.load import samplesheet
    if not options.samplesheet:
        parser.error("Samplesheet (-s) is not specified.")
    if not options.workdir:
        parser.error("Workdir (-w) is not specified.")
    if not os.path.exists(os.path.join(options.cachedir, 'reference', 'ref.ok')):
        parser.error("Cache (%s) is not ready for use, check parameters!"%(options.cachedir))
    if not os.path.exists(os.path.join(os.path.realpath(options.workdir), 'Pipe_runtime', snapshot, 'logs', 'slurm')):
        os.system('mkdir -p %s'%(os.path.join(os.path.realpath(options.workdir), 'Pipe_runtime', snapshot, 'logs', 'slurm')))
    # Check and write samplesheet to the workdir
    try:
        sampledic, libdic, rundic = samplesheet(options.samplesheet)
    except:
        parser.error("Samplesheet failed to load. Check %s"%(options.samplesheet))
    
    ssdf = pd.read_csv(options.samplesheet, sep='\t')
    for i in ssdf.index:
        for read in ['Read1','Read2']:
            tmpreads = []
            for r in ssdf.loc[i,read].split(','):
                if r[0] == '/':
                    if not os.path.exists(r):
                        parser.error("Couldn't find fastq file: RG:%s,PATH:%s."%(i, r))
                    else:
                        tmpreads.append(r)
                else:
                    try:
                        realpath = os.path.realpath(r)
                        if os.path.exists(realpath):
                            tmpreads.append(realpath)
                    except:
                        raise Exception("Fastq file unaccessible: %s"%ssdf.loc[i,read])
            ssdf.loc[i,read] = ','.join(tmpreads)
    ssdf = ssdf[['Sample','Lib','RunId','Flowcell','Lane','Platform','Read1','Read2']]
    ssdf.to_csv(os.path.join(options.workdir, 'Pipe_runtime', snapshot, 'samplesheet.tsv'), sep='\t', index=None)
    
    #######################################################################################
    # Prepare config file
    #######################################################################################
    config = {}
    config['samplesheet'] = os.path.realpath(os.path.join(options.workdir, 'Pipe_runtime', snapshot, 'samplesheet.tsv'))
    config['workdir'] = os.path.realpath(options.workdir)
    config['pipelinedir'] = os.path.join(config['workdir'], 'Pipe_runtime', snapshot)
    config['cachedir'] = os.path.realpath(options.cachedir)
    config['snapshot'] = snapshot
    config['options'] = {}

    with open(os.path.join(pipedir, 'config', 'config.yaml'), 'r') as infile:
        cfg = yaml.safe_load(infile)
        for i in cfg:
            if i not in config:
                config[i] = cfg[i]

    with open(os.path.join(pipedir, 'config', 'reference.yaml'), 'r') as infile:
        refcfg = yaml.safe_load(infile)
        if 'references' not in config:
            config['references'] = {}
        for i in refcfg['references']:
            config['references'][i] = os.path.join(config['cachedir'], refcfg['references'][i]['path'])

    with open(os.path.join(pipedir, 'config', 'container.yaml'), 'r') as infile:
        tmpcfg = yaml.safe_load(infile)
        config['container'] = {}
        if not tmpcfg['docker']:
            tmpcfg['docker'] = {}
        if not tmpcfg['singularity']:
            tmpcfg['singularity'] = {}
        if not tmpcfg['simg']:
            tmpcfg['simg'] = {}
        for c in tmpcfg['docker']:
            config['container'][c] = os.path.join(config['cachedir'], 'container', '%s.simg'%(c))
        for c in tmpcfg['singularity']:
            config['container'][c] = os.path.join(config['cachedir'], 'container', '%s.simg'%(c))
        for c in tmpcfg['simg']:
            config['container'][c] = os.path.join(config['cachedir'], 'container', '%s.simg'%(c))
        
    with open(os.path.join(config['workdir'], 'Pipe_runtime', snapshot, 'config.yaml'), 'w') as savefile:
        savefile.write(yaml.dump(config))
    
    # copy cluster.yaml
    # if os.path.exists(os.path.join(config['workdir'], 'Pipe_runtime', snapshot, 'cluster.yaml')):
    #     os.remove(os.path.join(config['workdir'], 'Pipe_runtime', snapshot, 'cluster.yaml'))
    # os.system('cp %s %s'%(os.path.join(pipedir, 'config', 'cluster.yaml'), os.path.join(config['workdir'], 'Pipe_runtime', snapshot, 'cluster.yaml')))
    
    #######################################################################################
    # Prepare other files
    #######################################################################################
    
    # copy scripts
    if os.path.exists(os.path.join(config['workdir'], 'Pipe_runtime', snapshot, 'scripts')):
        os.system('rm -rf %s'%(os.path.join(config['workdir'], 'Pipe_runtime', snapshot, 'scripts')))
    os.system('cp -r %s %s'%( os.path.join(pipedir, 'scripts'), os.path.join(config['workdir'], 'Pipe_runtime', snapshot, 'scripts') ))
    
    # copy rules
    if os.path.exists(os.path.join(config['workdir'], 'Pipe_runtime', snapshot, 'rules')):
        os.system('rm -rf %s'%(os.path.join(config['workdir'], 'Pipe_runtime', snapshot, 'rules')))
    os.system('cp -r %s %s'%( os.path.join(pipedir, 'rules'), os.path.join(config['workdir'], 'Pipe_runtime', snapshot, 'rules') ))
    
    # copy slurm
    if os.path.exists(os.path.join(config['workdir'], 'Pipe_runtime', snapshot, 'slurm')):
        os.system('rm -rf %s'%(os.path.join(config['workdir'], 'Pipe_runtime', snapshot, 'slurm')))
    os.system('cp -r %s %s'%( os.path.join(pipedir, 'slurm'), os.path.join(config['workdir'], 'Pipe_runtime', snapshot, 'slurm') ))

    # copy snakefile
    if os.path.exists(os.path.join(config['workdir'], 'Pipe_runtime', snapshot, 'Snakefile')):
        os.system('rm %s'%(os.path.join(config['workdir'], 'Pipe_runtime', snapshot, 'Snakefile')))
    os.system('cp -r %s %s'%( os.path.join(pipedir, 'Snakefile'), os.path.join(config['workdir'], 'Pipe_runtime', snapshot, 'Snakefile') ))
    
    # write pipeline submission bash script: captain.sh
    with open(os.path.join(config['workdir'], 'Pipe_runtime', snapshot, 'captain.sh'), 'w') as savefile:
        with open(os.path.join(pipedir, 'captain.sh')) as infile:
            captain = infile.read()
            if options.silent:
                captain = captain.replace('#SBATCH --mail-type=BEGIN,END,FAIL\n','')
            captain = captain.replace('[[PIPENICKNAME]]','FASTQ_TO_CRAM_PIPE')
            captain = captain.replace('[[WORKDIR]]', config['workdir'])
            captain = captain.replace('[[BINDPATH]]', config['bindpath'])
            captain = captain.replace('[[SNAPSHOT]]', config['snapshot'])
            savefile.write(captain)
    os.chmod(os.path.join(config['workdir'], 'Pipe_runtime', snapshot, 'captain.sh'), 0o755)
    
    #######################################################################################
    # submit captain.sh to slurm
    #######################################################################################
    if os.path.exists(os.path.join(config['workdir'], '.snakemake', 'locks', '0.input.lock')) or os.path.exists(os.path.join(config['workdir'], '.snakemake', 'locks', '0.output.lock')):
        raise Exception("Workdir locked! Run the pipemaster with -u/--unlock to unlock the workdir first.")

    print(config)
    if options.dryrun:
        os.system('cd %s && snakemake -n'%(config['pipelinedir']))
    else:
        cmd = 'cd %s && '%config['pipelinedir']
        cmd += 'sbatch captain.sh'
        try:
            subprocess.check_call([cmd], shell=True, executable="/bin/bash")
            print('Pipeline submitted!')
        except:
            raise Exception("Pipeline submission failed.")
    
if __name__ == '__main__':
    main()

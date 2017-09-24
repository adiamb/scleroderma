import os
import glob

##### Paths
VIRTUALENV=             '/path/to/virtual/environment'
CHANGEO_PATH=           '/path/to/changeo-0.3.3/bin'
R_PATH=                 '/path/to/R-3.1.2'
IMGT_GERMLINES=         '/path/to/tigger_germlines'
SCRIPTS=                '/path/to/script/folder'
DATA=                   '/path/to/data/folder'
IMGT_RESULTS=           DATA+'/IMGT_output' 
WDIR=                   DATA+'/Snakemake_wdir'
LOGS=                   WDIR+'/logs'

if not os.path.exists(LOGS):
    os.makedirs(LOGS)
    
##### Sample lists
SAMPLE_LIST = [os.path.basename(FILE).split('_minCONSCOUNT2')[0] for FILE in glob.glob(WDIR+'/*_minCONSCOUNT2.fasta')]
PARTICIPANTS = list(set([SAMPLE.split("_")[0] for SAMPLE in SAMPLE_LIST]))
SAMPLES = {PARTICIPANT:[os.path.basename(s).split('_minCONSCOUNT2')[0] for s in glob.glob(WDIR+'/'+PARTICIPANT+'*_minCONSCOUNT2.fasta')] for PARTICIPANT in PARTICIPANTS}
def participant_from_sample(sample):
    for participant in SAMPLES.keys():
      if sample in SAMPLES[participant]:
        return participant

##### Rules
rule all:
    input:
        expand(WDIR+'/{participant}_minCONSCOUNT2_all-chunks-db-pass_FUNCTIONAL-TRUE_visits-pooled_distinct-alleles-pass_germ-pass_nomask_germ-pass.tab', participant=PARTICIPANTS)
    
rule create_germlines_nomask:
    input:
        WDIR+'/{participant}_minCONSCOUNT2_all-chunks-db-pass_FUNCTIONAL-TRUE_visits-pooled_distinct-alleles-pass_germ-pass.tab',
        WDIR+'/{participant}_personal_hIGHV.fasta'
    output:
        WDIR+'/{participant}_minCONSCOUNT2_all-chunks-db-pass_FUNCTIONAL-TRUE_visits-pooled_distinct-alleles-pass_germ-pass_nomask_germ-pass.tab',
        WDIR+'/{participant}_minCONSCOUNT2_all-chunks-db-pass_FUNCTIONAL-TRUE_visits-pooled_distinct-alleles-pass_germ-pass_nomask_germ-fail.tab'
    params: name='create_germlines_nomask', partition='general', cpus_per_task='12', mem_per_cpu='5300', logfile=LOGS+'/{participant}_create_germlines_nomask.log'
    run:
        input_concat = ' '.join(input)
        input_scratch = [os.path.basename(s) for s in input]
        input_scratch_concat = ' '.join(input_scratch)
        output_scratch = [os.path.basename(s) for s in output]
        output_scratch_concat = ' '.join(output_scratch)
        prefix = os.path.splitext(input_scratch[0])[0]+'_nomask'
        shell("cp {input_concat} $LOCAL_SATA &&\
               cd $LOCAL_SATA &&\
               mkdir germline_repo &&\
               cp {input_scratch[1]} germline_repo &&\
               cp {IMGT_GERMLINES}/IMGT_hIGHD.fasta germline_repo &&\
               cp {IMGT_GERMLINES}/IMGT_hIGHJ.fasta germline_repo &&\
               set +o nounset &&\
               source {VIRTUALENV}/bin/activate {VIRTUALENV} &&\
               set -o nounset &&\
               python {CHANGEO_PATH}/CreateGermlines.py -d {input_scratch[0]} --failed -r germline_repo -g full --outname {prefix} --vf V_CALL_GENOTYPED --sf SEQUENCE_IMGT &&\
               cp {output_scratch_concat} {WDIR}/")

rule create_germlines:
    input:
        WDIR+'/{participant}_minCONSCOUNT2_all-chunks-db-pass_FUNCTIONAL-TRUE_visits-pooled_distinct-alleles-pass.tab',
        WDIR+'/{participant}_personal_hIGHV.fasta'
    output:
        WDIR+'/{participant}_minCONSCOUNT2_all-chunks-db-pass_FUNCTIONAL-TRUE_visits-pooled_distinct-alleles-pass_germ-pass.tab',
        WDIR+'/{participant}_minCONSCOUNT2_all-chunks-db-pass_FUNCTIONAL-TRUE_visits-pooled_distinct-alleles-pass_germ-fail.tab'
    params: name='create_germlines', partition='general', cpus_per_task='12', mem_per_cpu='5300', logfile=LOGS+'/{participant}_create_germlines.log'
    run:
        input_concat = ' '.join(input)
        input_scratch = [os.path.basename(s) for s in input]
        input_scratch_concat = ' '.join(input_scratch)
        output_scratch = [os.path.basename(s) for s in output]
        output_scratch_concat = ' '.join(output_scratch)
        shell("cp {input_concat} $LOCAL_SATA &&\
               cd $LOCAL_SATA &&\
               mkdir germline_repo &&\
               cp {input_scratch[1]} germline_repo &&\
               cp {IMGT_GERMLINES}/IMGT_hIGHD.fasta germline_repo &&\
               cp {IMGT_GERMLINES}/IMGT_hIGHJ.fasta germline_repo &&\
               set +o nounset &&\
               source {VIRTUALENV}/bin/activate {VIRTUALENV} &&\
               set -o nounset &&\
               python {CHANGEO_PATH}/CreateGermlines.py -d {input_scratch[0]} --failed -r germline_repo --vf V_CALL_GENOTYPED --sf SEQUENCE_IMGT &&\
               cp {output_scratch_concat} {WDIR}/")
               
rule tigger:
    input:
        WDIR+'/{participant}_minCONSCOUNT2_all-chunks-db-pass_FUNCTIONAL-TRUE_visits-pooled.tab',
        IMGT_GERMLINES+'/IMGT_hIGHV.fasta'
    output:
        WDIR+'/{participant}_minCONSCOUNT2_all-chunks-db-pass_FUNCTIONAL-TRUE_visits-pooled_distinct-alleles-pass.tab',
        WDIR+'/{participant}_personal_hIGHV.fasta'
    params: name='tigger', partition='long', nproc='4', cpus_per_task='4', mem_per_cpu='5300', logfile=LOGS+'/{participant}_tigger.log'
    run:
        input_concat = ' '.join(input)
        input_scratch = [os.path.basename(s) for s in input]
        input_scratch_concat = ' '.join(input_scratch)
        output_scratch = [os.path.basename(s) for s in output]
        output_scratch_concat = ' '.join(output_scratch)
        shell("cp {input_concat} $LOCAL_SATA &&\
               cd $LOCAL_SATA &&\
               {R_PATH}/bin/Rscript {SCRIPTS}/apply_tigger.R {input_scratch[0]} {input_scratch[1]} {params.nproc} &&\
               cp {output_scratch_concat} {WDIR}/")        
           
rule pool_visits:
    input:
        lambda wildcards: expand(WDIR+'/{sample}_minCONSCOUNT2_all-chunks-db-pass_FUNCTIONAL-TRUE.tab', sample=SAMPLES[wildcards.participant])
    output:
        WDIR+'/{participant}_minCONSCOUNT2_all-chunks-db-pass_FUNCTIONAL-TRUE_visits-pooled.tab'
    params: name='pool_visits', partition='general', cpus_per_task='12', mem_per_cpu='5300', logfile=LOGS+'/{participant}_pool_visits.log'
    run:
        input_concat = ' '.join(input)
        input_scratch = [os.path.basename(s) for s in input]
        input_scratch_concat = ' '.join(input_scratch)
        output_scratch = os.path.basename(output[0])
        shell("cp {input_concat} $LOCAL_SATA &&\
               cd $LOCAL_SATA &&\
               {R_PATH}/bin/Rscript {SCRIPTS}/pool_visits.R {input_scratch_concat} &&\
               cp {output_scratch} {WDIR}/")         
        
rule filter_functional:
    input:
        WDIR+'/{sample}_minCONSCOUNT2_all-chunks-db-pass.tab'
    output:
        WDIR+'/{sample}_minCONSCOUNT2_all-chunks-db-pass_FUNCTIONAL-TRUE.tab',
        WDIR+'/{sample}_minCONSCOUNT2_all-chunks-db-pass_FUNCTIONAL-FALSE.tab'
    params: name='filter_functional', partition='general', cpus_per_task='4', mem_per_cpu='5300', logfile=LOGS+'/{sample}_filter_functional.log'
    run:
        input_scratch = os.path.basename(input[0])
        output_scratch = [os.path.basename(s) for s in output]
        output_scratch_concat = ' '.join(output_scratch)
        shell("cp {input} $LOCAL_SATA &&\
                cd $LOCAL_SATA &&\
                set +o nounset &&\
                source {VIRTUALENV}/bin/activate {VIRTUALENV} &&\
                set -o nounset &&\
                python {CHANGEO_PATH}/ParseDb.py split -d {input_scratch} -f FUNCTIONAL &&\
                mv {wildcards.sample}_minCONSCOUNT2_all-chunks-db-pass_FUNCTIONAL-T.tab {wildcards.sample}_minCONSCOUNT2_all-chunks-db-pass_FUNCTIONAL-TRUE.tab &&\
                touch {wildcards.sample}_minCONSCOUNT2_all-chunks-db-pass_FUNCTIONAL-F.tab &&\
                mv {wildcards.sample}_minCONSCOUNT2_all-chunks-db-pass_FUNCTIONAL-F.tab {wildcards.sample}_minCONSCOUNT2_all-chunks-db-pass_FUNCTIONAL-FALSE.tab &&\
                cp {output_scratch_concat} {WDIR}/")
           
rule make_db:
    input:
        WDIR+'/{sample}_minCONSCOUNT2.fasta'
    output:
        WDIR+'/{sample}_minCONSCOUNT2_all-chunks-db-pass.tab'
    params: name='make_db', partition='general', cpus_per_task='4', mem_per_cpu='5300', logfile=LOGS+'/{sample}_make_db.log'
    run:
        input_scratch = [os.path.basename(s) for s in input]
        output_scratch = os.path.basename(output[0])
        shell("cp {input[0]} $LOCAL_SATA &&\
               cp -r {IMGT_RESULTS}/{wildcards.sample}_minCONSCOUNT2 $LOCAL_SATA &&\
               cd $LOCAL_SATA &&\
               set +o nounset &&\
               source {VIRTUALENV}/bin/activate {VIRTUALENV} &&\
               set -o nounset &&\
               python {CHANGEO_PATH}/MakeDb.py imgt -i {wildcards.sample}_minCONSCOUNT2 -s {input_scratch[0]} &&\
               mv {wildcards.sample}_minCONSCOUNT2_db-pass.tab {wildcards.sample}_minCONSCOUNT2_all-chunks-db-pass.tab &&\
               cp {output_scratch} {WDIR}/")
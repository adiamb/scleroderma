import os
import glob

##### Paths
R_PATH=                 '/path/to/R-3.1.2'
SCRIPTS=                '/path/to/script/folder'
DATA=                   '/path/to/data/folder/'
WDIR=                   DATA+'/Snakemake_wdir'
LOGS=                   WDIR+'/logs'
if not os.path.exists(LOGS):
    os.makedirs(LOGS)
    
##### Sample lists
PARTICIPANTS = [os.path.basename(FILE).split('_')[0] for FILE in glob.glob(WDIR+'/participant*_minCONSCOUNT2.fasta')]
BIN_SIZE = 7000

##### Rules
rule all:
    input:
        WDIR+'/shared-CDR3s-1mismatch_all-sclero-data_combined.RDS'

rule list_CDR3s:
    input:
        expand(WDIR+'/{participant}_minCONSCOUNT2_all-chunks-db-pass_FUNCTIONAL-TRUE_distinct-alleles-pass_germ-pass_nomask_germ-pass_lineage-pass_mutations-pass-corrected_remove-error-clouds-pass_attach-AA-sequences-pass-corrected_attach-AA-mutationtypes-pass-corr.tab', participant=PARTICIPANTS)
    output:
        WDIR+'/AA_CDR3s.RDS',
        WDIR+'/all_AA_CDR3s.RDS'
    params: name='list_CDR3s', partition='general', cpus_per_task='4', mem_per_cpu='5300', logfile=LOGS+'/Fisher_tests.log'
    run:
        input_concat = ' '.join(input)
        input_scratch = [os.path.basename(s) for s in input]
        input_scratch_concat = ' '.join(input_scratch)
        output_scratch = [os.path.basename(s) for s in output]
        output_scratch_concat = ' '.join(output_scratch)
        shell("cp {input} $LOCAL_SATA &&\
                cd $LOCAL_SATA &&\
                {R_PATH}/bin/Rscript {SCRIPTS}/list_CDR3s_all-sclero-data.R {input_scratch_concat} &&\
                cp {output_scratch_concat} {WDIR}/")
        
rule Fisher_tests_1mismatch:
    input:
        WDIR+'/AA_CDR3s.RDS',
        WDIR+'/all_AA_CDR3s.RDS'
    output:
        WDIR+'/shared-CDR3s-1mismatch_all-sclero-data_i{START}.RDS'
    params: name='Fisher_tests_1mismatch', partition='general', cpus_per_task='4', mem_per_cpu='5300', logfile=LOGS+'/{START}_Fisher_tests_1mismatch.log'
    run:
        input_concat = ' '.join(input)
        input_scratch = [os.path.basename(s) for s in input]
        input_scratch_concat = ' '.join(input_scratch)
        output_scratch = os.path.basename(output[0])
        shell("cp {input} $LOCAL_SATA &&\
                cd $LOCAL_SATA &&\
                {R_PATH}/bin/Rscript {SCRIPTS}/Fisher_tests_1mismatch.R {wildcards.START} {BIN_SIZE} {output_scratch} {input_scratch[0]} {input_scratch[1]} &&\
                cp {output_scratch} {WDIR}/")
                
rule combine_Fisher_results:
    input:
        expand(WDIR+'/shared-CDR3s-1mismatch_all-sclero-data_i{START}.RDS', START=list(range(1,2410035,BIN_SIZE)))
    output:
        WDIR+'/shared-CDR3s-1mismatch_all-sclero-data_combined.RDS'
    params: name='combine_Fisher_results', partition='general', cpus_per_task='4', mem_per_cpu='5300', logfile=LOGS+'/combine_Fisher_results.log'
    run:
        input_concat = ' '.join(input)
        input_scratch = [os.path.basename(s) for s in input]
        input_scratch_concat = ' '.join(input_scratch)
        output_scratch = os.path.basename(output[0])
        shell("cp {input} $LOCAL_SATA &&\
                cd $LOCAL_SATA &&\
                {R_PATH}/bin/Rscript {SCRIPTS}/combine_Fisher_results.R {input_scratch_concat} {output_scratch} &&\
                cp {output_scratch} {WDIR}/")
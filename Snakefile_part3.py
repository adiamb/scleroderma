import os
import glob

##### Paths
R_PATH=                 '/path/to/R-3.1.2'
SCRIPTS=                '/path/to/script/folder'
DATA=                   '/path/to/data/folder'
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
        
RUBL_Nseq = [str(int(1e4))]
GERMLINE_FIELD = ['GERMLINE_IMGT']


##### Rules
rule all:
    input:
        expand(WDIR+'/{participant}_minCONSCOUNT2_all-chunks-db-pass_FUNCTIONAL-TRUE_visits-pooled_distinct-alleles-pass_germ-pass_nomask_germ-pass_lineage-pass_mutations-pass-corrected_remove-error-clouds-pass_attach-AA-sequences-pass-corrected_attach-AA-mutationtypes-pass-corr.tab', participant=PARTICIPANTS),
        expand(WDIR+'/{participant}_minCONSCOUNT2_all-chunks-db-pass_FUNCTIONAL-TRUE_visits-pooled_distinct-alleles-pass_germ-pass_nomask_germ-pass_lineage-pass_mutations-pass-corrected_remove-error-clouds-pass_attach-AA-sequences-pass-corrected_attach-change-table.tab', participant=PARTICIPANTS),
        expand(WDIR+'/RUBL.Nseq-{N}.{germ_field}.tab', N=RUBL_Nseq, germ_field=GERMLINE_FIELD)

rule define_lineages:
    input:
        WDIR+'/{participant}_minCONSCOUNT2_all-chunks-db-pass_FUNCTIONAL-TRUE_visits-pooled_distinct-alleles-pass_germ-pass_nomask_germ-pass.tab'
    output:
        WDIR+'/{participant}_minCONSCOUNT2_all-chunks-db-pass_FUNCTIONAL-TRUE_visits-pooled_distinct-alleles-pass_germ-pass_nomask_germ-pass_lineage-pass.tab'
    params: name='define_lineages', partition='long', cpus_per_task='12', mem_per_cpu='5300', logfile=LOGS+'/{participant}_define_lineages.log'
    run:
        input_scratch = os.path.basename(input[0])
        output_scratch = os.path.basename(output[0])
        shell("cp {input} $LOCAL_SATA &&\
                cd $LOCAL_SATA &&\
                {R_PATH}/bin/Rscript {SCRIPTS}/define_lineages.R {input_scratch} &&\
                cp {output_scratch} {WDIR}/")
                
rule compute_RUBL:
    input:
        expand(WDIR+'/{participant}_minCONSCOUNT2_all-chunks-db-pass_FUNCTIONAL-TRUE_visits-pooled_distinct-alleles-pass_germ-pass_nomask_germ-pass.tab', participant=PARTICIPANTS),
        expand(WDIR+'/{participant}_personal_hIGHV.fasta', participant=PARTICIPANTS)
    output:
        WDIR+'/RUBL.Nseq-{N}.{germ_field}.tab'
    params: name='compute_RUBL', partition='unrestricted', cpus_per_task='12', mem_per_cpu='5300', logfile=LOGS+'/compute_RUBL_Nseq-{N}_{germ_field}.log'
    run:
        input_concat = ' '.join(input)
        input_scratch = [os.path.basename(s) for s in input]
        input_scratch_concat = ' '.join(input_scratch)
        output_scratch = os.path.basename(output[0])
        shell("cp {input_concat} $LOCAL_SATA &&\
               cd $LOCAL_SATA &&\
               {R_PATH}/bin/Rscript {SCRIPTS}/compute_RUBL.R {wildcards.N} {wildcards.germ_field} {input_scratch_concat} &&\
               cp {output_scratch} {WDIR}/")     
               
rule count_mutations:
    input:
        WDIR+'/{participant}_minCONSCOUNT2_all-chunks-db-pass_FUNCTIONAL-TRUE_visits-pooled_distinct-alleles-pass_germ-pass_nomask_germ-pass_lineage-pass.tab'
    output:
        WDIR+'/{participant}_minCONSCOUNT2_all-chunks-db-pass_FUNCTIONAL-TRUE_visits-pooled_distinct-alleles-pass_germ-pass_nomask_germ-pass_lineage-pass_mutations-pass-corrected.tab'
    params: name='count_mutations', partition='long', cpus_per_task='6', mem_per_cpu='5300', logfile=LOGS+'/{participant}_count_mutations.log'
    run:
        input_scratch = os.path.basename(input[0])
        output_scratch = os.path.basename(output[0])
        shell("cp {input} $LOCAL_SATA &&\
                cd $LOCAL_SATA &&\
                {R_PATH}/bin/Rscript {SCRIPTS}/count_mutations_corrected.R {WDIR} {input_scratch} &&\
                cp {output_scratch} {WDIR}/")     

rule count_mutations_nonfunctional:
    input:
        WDIR+'/{sample}_minCONSCOUNT2_all-chunks-db-pass_FUNCTIONAL-FALSE_V-CALL-GENOTYPED-pass_germ-pass_nomask_germ-pass.tab'
    output:
        WDIR+'/{sample}_minCONSCOUNT2_all-chunks-db-pass_FUNCTIONAL-FALSE_V-CALL-GENOTYPED-pass_germ-pass_nomask_germ-pass_mutations-pass-corrected.tab'
    params: name='count_mutations_nonfunctional', partition='general', cpus_per_task='6', mem_per_cpu='5300', logfile=LOGS+'/{sample}_count_mutations_nonfunctional.log'
    run:
        input_scratch = os.path.basename(input[0])
        output_scratch = os.path.basename(output[0])
        shell("cp {input} $LOCAL_SATA &&\
                cd $LOCAL_SATA &&\
                {R_PATH}/bin/Rscript {SCRIPTS}/count_mutations_nonfunctional_corrected.R {WDIR} {input_scratch} &&\
                cp {output_scratch} {WDIR}/")     
                
rule remove_error_clouds:
    input:
        WDIR+'/{participant}_minCONSCOUNT2_all-chunks-db-pass_FUNCTIONAL-TRUE_visits-pooled_distinct-alleles-pass_germ-pass_nomask_germ-pass_lineage-pass_mutations-pass-corrected.tab'
    output:
        WDIR+'/{participant}_minCONSCOUNT2_all-chunks-db-pass_FUNCTIONAL-TRUE_visits-pooled_distinct-alleles-pass_germ-pass_nomask_germ-pass_lineage-pass_mutations-pass-corrected_remove-error-clouds-pass.tab'
    params: name='remove_error_clouds', partition='long', cpus_per_task='6', mem_per_cpu='5300', logfile=LOGS+'/{participant}_remove_error_clouds.log'
    run:
        input_scratch = os.path.basename(input[0])
        output_scratch = os.path.basename(output[0])
        shell("cp {input} $LOCAL_SATA &&\
                cd $LOCAL_SATA &&\
                {R_PATH}/bin/Rscript {SCRIPTS}/remove_error_clouds.R {input_scratch} &&\
                cp {output_scratch} {WDIR}/")     
                
rule attach_AA_sequences:
    input:
        WDIR+'/{participant}_minCONSCOUNT2_all-chunks-db-pass_FUNCTIONAL-TRUE_visits-pooled_distinct-alleles-pass_germ-pass_nomask_germ-pass_lineage-pass_mutations-pass-corrected_remove-error-clouds-pass.tab'
    output:
        WDIR+'/{participant}_minCONSCOUNT2_all-chunks-db-pass_FUNCTIONAL-TRUE_visits-pooled_distinct-alleles-pass_germ-pass_nomask_germ-pass_lineage-pass_mutations-pass-corrected_remove-error-clouds-pass_attach-AA-sequences-pass-corrected.tab'
    params: name='attach_AA_sequences', partition='general', cpus_per_task='6', mem_per_cpu='5300', logfile=LOGS+'/{participant}_attach_AA_sequences.log'
    run:
        input_scratch = os.path.basename(input[0])
        output_scratch = os.path.basename(output[0])
        shell("cp {input} $LOCAL_SATA &&\
                cd $LOCAL_SATA &&\
                {R_PATH}/bin/Rscript {SCRIPTS}/attach_AA_sequences_corrected.R {WDIR} {input_scratch} &&\
                cp {output_scratch} {WDIR}/")     
                
rule attach_AA_mutation_types:
    input:
        WDIR+'/{participant}_minCONSCOUNT2_all-chunks-db-pass_FUNCTIONAL-TRUE_visits-pooled_distinct-alleles-pass_germ-pass_nomask_germ-pass_lineage-pass_mutations-pass-corrected_remove-error-clouds-pass_attach-AA-sequences-pass-corrected.tab'
    output:
        WDIR+'/{participant}_minCONSCOUNT2_all-chunks-db-pass_FUNCTIONAL-TRUE_visits-pooled_distinct-alleles-pass_germ-pass_nomask_germ-pass_lineage-pass_mutations-pass-corrected_remove-error-clouds-pass_attach-AA-sequences-pass-corrected_attach-AA-mutationtypes-pass-corr.tab'
    params: name='attach_AA_mutation_types', partition='general', cpus_per_task='6', mem_per_cpu='5300', logfile=LOGS+'/{participant}_attach_AA_mutation_types.log'
    run:
        input_scratch = os.path.basename(input[0])
        output_scratch = os.path.basename(output[0])
        shell("cp {input} $LOCAL_SATA &&\
                cd $LOCAL_SATA &&\
                {R_PATH}/bin/Rscript {SCRIPTS}/attach_AA_mutation_types_corrected.R {WDIR} {input_scratch} &&\
                cp {output_scratch} {WDIR}/")     
                
rule attach_change_table:
    input:
        WDIR+'/{sample}_minCONSCOUNT2_all-chunks-db-pass_FUNCTIONAL-TRUE_visits-pooled_distinct-alleles-pass_germ-pass_nomask_germ-pass_lineage-pass_mutations-pass-corrected_remove-error-clouds-pass_attach-AA-sequences-pass-corrected_attach-AA-mutationtypes-pass-corr.tab'
    output:
        WDIR+'/{sample}_minCONSCOUNT2_all-chunks-db-pass_FUNCTIONAL-TRUE_visits-pooled_distinct-alleles-pass_germ-pass_nomask_germ-pass_lineage-pass_mutations-pass-corrected_remove-error-clouds-pass_attach-AA-sequences-pass-corrected_attach-change-table.tab'
    params: name='attach_change_table', partition='general', cpus_per_task='6', mem_per_cpu='5300', logfile=LOGS+'/{sample}_attach_change_table.log'
    run:
        input_scratch = os.path.basename(input[0])
        output_scratch = os.path.basename(output[0])
        shell("cp {input} $LOCAL_SATA &&\
                cd $LOCAL_SATA &&\
                {R_PATH}/bin/Rscript {SCRIPTS}/attach_change_table.R {WDIR} {input_scratch} &&\
                cp {output_scratch} {WDIR}/")                     
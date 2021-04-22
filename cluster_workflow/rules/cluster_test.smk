rule cluster_test:
    ''' Test snakemake on cluster '''
  output:
    'testfile'
  shell:
    'touch {output}'

'''
snakemake --profile ../snakemake_profile testfile
'''

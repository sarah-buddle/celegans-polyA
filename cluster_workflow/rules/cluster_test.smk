rule cluster_test:
  output:
    'testfile'
  shell:
    'touch {output}'

'''
snakemake --profile ../snakemake_profile testfile
'''

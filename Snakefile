#How to run: snakemake
#SAMPLES = ["20190124_CFMPB559run1_CF30207_01", "20190124_CFMPB559run1_CF30212_06", "20190124_CFMPB559run1_CF30209_03", "20190124_CFMPB559run1_CF30213_07", "20190124_CFMPB559run1_CF30211_05", "20190124_CFMPB559run1_CF30214_08"]
import os
path = os.getcwd()
SAMPLES = [f for f in os.listdir(path) if f.endswith('.RCC')]


rule all:
	input:
		expand('{sample}.pdf', sample=SAMPLES)


rule QC:
	input:
		'/home/b.moumbeini/snakemaker/{sample}.RCC'


	output:
		'{sample}.pdf'
		'quality_control.csv'

	shell:
		"python3 QC.py"


rule GROUP_COMPARISON:
	input:
		expand('{sample}.RCC', sample=SAMPLES)
		metadata = 'metadata.txt'
		GENE = "gene_name"

	output:
               	'{GENE}.pdf'
                'Normalized_data.csv'
		'heatmap.html'

	shell:
		"python3 group_comparison.py --input {input.GENE}"

rule GENE_EXPRESSION_VOLCANO:
	input:
		metadata = 'metadata.txt'
		normal = 'Normalized_data.csv'

	output:
		'DGE.csv'
		'NanoString_Volcano.pdf'

	shell:
		"python3 GE_volcano.py"

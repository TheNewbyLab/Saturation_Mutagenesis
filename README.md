# Saturation_Mutagenesis_MPRA
 GNE project led by Noah Workman 

## Dependencies 
Besides the initial alignment (bwa, bowtie etc.) and samtools commands everything is done in Python. 
The Python packages that are required are numpy, pandas, Biopython (link: https://biopython.org/) and Pysam (link: https://pysam.readthedocs.io/en/latest/api.html) 

## How to run 
At this moment the script is set up based on Noah's initial MRPA assay design (i.e. U6 promoter flanked by barcodes of different length). If you are using a different ddesign you'll have to update the script so it knows how long of barcodes should be expected and where the will be in the sequence. 

As of now the Process_BAM.py python script contains all the commands needed to run the script.  

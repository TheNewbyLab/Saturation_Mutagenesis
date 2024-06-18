# Saturation_Mutagenesis_MPRA

This project is led by Noah Workman and is focused on performing Saturation Mutagenesis using MPRA (Massively Parallel Reporter Assay).

## Dependencies

The majority of the analysis is performed using Python, except for the initial alignment and samtools commands. The required Python packages are:

- **numpy**
- **pandas**
- **Biopython**: [Biopython website](https://biopython.org/)
- **Pysam**: [Pysam documentation](https://pysam.readthedocs.io/en/latest/api.html)

Ensure you have these packages installed before running the script.

## Installation

You can install the required Python packages using pip:

```bash
pip install numpy pandas biopython pysam
```

## How to Run

### Initial Setup

1. **Alignment**: Perform the initial alignment using tools such as `bwa` or `bowtie`, and process the aligned reads using `samtools`.

2. **Script Configuration**: The script is currently configured for Noah's initial MRPA assay design, which uses a U6 promoter flanked by barcodes of different lengths. If your design differs, you'll need to modify the script to handle the expected barcode lengths and their positions in the sequence.

### Running the Script

1. **Edit the Script**: If necessary, update the `SM_MPRA_extract_data.py` script to match your design parameters.

2. **Execute the Script**: Run the script using the following command:

```bash
python SM_MPRA_extract_data.py path_to_bam output_file.csv
```
- `path_to_bam`: The path to your input BAM file.
- `output_file.csv`: The name of the output CSV file where the results will be saved.

### Example Command

```bash
python SM_MPRA_extract_data.py ./data/aligned_reads.bam ./results/mpra_output.csv
```

## Notes

- Ensure your input BAM file is properly formatted and contains the necessary read alignments.
- Review the script to confirm it is tailored to your specific MPRA design before running it.

## License

This project is licensed under the MIT License. See the LICENSE file for details.

## Contact

For any questions or further information, please contact Noah Workman or Daniel Rabizadeh. 


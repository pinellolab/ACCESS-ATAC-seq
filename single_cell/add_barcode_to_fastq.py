# Add cell barcodes to a bam file from a fastq file.

import argparse
import logging
import warnings
import gzip

warnings.filterwarnings("ignore")

logging.basicConfig(
    format="%(asctime)s %(levelname)-8s %(message)s",
    level=logging.INFO,
    datefmt="%Y-%m-%d %H:%M:%S",
)

def parse_args():
    parser = argparse.ArgumentParser(
        description="This script adds barcode to fastq file",
        formatter_class=argparse.RawTextHelpFormatter,
    )

    # Required parameters
    parser.add_argument("--reads_fastq", type=str, default=None)
    parser.add_argument("--barcodes_fastq", type=str, default=None)
    parser.add_argument("--out_dir", type=str, default=None)
    parser.add_argument("--out_name", type=str, default=None)
    return parser.parse_args()

def main():
    args = parse_args()

    output_fastq = f"{args.out_dir}/{args.out_name}.fastq.gz"

    # read barcodes from the fastq file
    with gzip.open(args.reads_fastq, 'r') as r_fq, gzip.open(args.barcodes_fastq, 'r') as bc_fq, gzip.open(output_fastq, 'w') as out_fq:
        while True:
            # Read 4 lines from each file (FASTQ format)
            r_lines = [r_fq.readline() for _ in range(4)]
            bc_lines = [bc_fq.readline() for _ in range(4)]
            
            # convert bytes to string
            r_lines = [line.decode("utf-8") for line in r_lines]
            bc_lines = [line.decode("utf-8") for line in bc_lines]

            # Stop at end of file
            if not r_lines[0] or not bc_lines[0]:
                break
            
            # Extract read name from both files
            read_name = r_lines[0].strip().split()[0][1:]  # remove '@' and get the read name
            barcode_name = bc_lines[0].strip().split()[0][1:]  # remove '@' and get the barcode name
            
            # Check if read name matches barcode name
            if read_name != barcode_name:
                logging.warning(f"Read name {read_name} does not match barcode name {barcode_name}. Skipping this pair.")
                continue
            
            # Extract barcode sequence
            barcode = bc_lines[1].strip()
            
            # Modify the read header (line 0)
            header = r_lines[0].strip().split(" ")[0]
            new_header = f"{header.strip()}_{barcode}\n"
            
            new_header = new_header.replace(" ", "_")
    
            # Write modified record
            out_fq.write(new_header.encode("utf-8"))
            out_fq.write(r_lines[1].encode("utf-8"))  # sequence
            out_fq.write(r_lines[2].encode("utf-8"))  # '+'
            out_fq.write(r_lines[3].encode("utf-8"))  # quality

    logging.info("Done!")

if __name__ == "__main__":
    main()
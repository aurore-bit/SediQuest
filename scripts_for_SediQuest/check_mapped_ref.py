import subprocess
import sys

def check_bam_reference(bam_file, expected_reference):
    command = f"samtools view -H {bam_file}"
    result = subprocess.run(command, shell=True, capture_output=True, text=True)
    if result.returncode != 0:
        print(f"Error running samtools: {result.stderr}")
        return False
    header_lines = result.stdout.split('\n')
    for line in header_lines:
        if line.startswith('@PG') and 'CL:' in line:
            cl_field = line.split('CL:')[1].strip()
            if '-g' in cl_field:
                ref_path = cl_field.split('-g')[1].split()[0]
                if expected_reference in ref_path:
                    return True
    return False

bam_file = sys.argv[1]
expected_reference = sys.argv[2]

if check_bam_reference(bam_file, expected_reference):
    print(f"The BAM file is mapped to the expected reference: {expected_reference}")
else:
    print(f"The BAM file is not mapped to the expected reference: {expected_reference}. Run map_all again!")
import argparse
import csv
import subprocess
import os

def run_subprocess(cmd):
    p = subprocess.run([cmd], shell=True)#, executable='bash')
    return 0

# de reads en barcode file worden meegegeven.
# van deze barcode bestand wordt de header gesplit en opgeslagen
# 

parser = argparse.ArgumentParser(description='Process input files')
parser.add_argument("--r1_in", metavar="reads1", action="store",
                    dest="reads1", help="left-hand fastq file")
parser.add_argument("--r2_in", metavar="reads2", action="store",
                    dest="reads2",
                    help="right-hand fastq file")
parser.add_argument("-b", "--barcodes", metavar="input", action="store",
                    dest="barcode", default="barcodes.tsv",
                    help="input tab separated barcode file")
parser.add_argument("--output-dir", metavar="outputdir", action="store",
                    dest="outputdir", default="",
                    help="Specify output directory, only for galaxy")
parser.add_argument('--tmpdir',
                    help='temporary directory')
args = parser.parse_args()

bc_file = open(args.barcode, 'r')
header = bc_file.readline().rstrip("\n")
bc_dict = {}
line2 = (bc_file.readlines()[1])
head = header.split(sep="\t")
line2split = line2.split(sep="\t")
nCol_bc = range(len(head))

for i in nCol_bc:
    bc_dict[head[i]] = [i]

with open(os.path.join(args.outputdir, "barcode_stacks.tsv"), 'w', newline='') as out_file:
    tsv_writer = csv.writer(out_file, delimiter='\t')
    bc_file = open(args.barcode, 'r')
    lines = bc_file.readlines()[1:]
    for line in lines:
        lineout = line.split()
        valuesout = []
        valuesout.append(lineout[int(bc_dict["Barcode_R1"][0])])
        valuesout.append(lineout[int(bc_dict["Barcode_R2"][0])])
        valuesout.append(lineout[int(bc_dict["Sample"][0])])
        tsv_writer.writerow(valuesout)
out_file.close()

if int(line2split[bc_dict["Wobble_R1"][0]]) != 0:
    cmd = f"clone_filter -1 {args.reads1} -2 {args.reads2} -o {os.path.join(args.tmpdir)} --inline_inline -igzfastq --oligo_len_1 {line2split[bc_dict['Wobble_R1'][0]]} --oligo_len_2 {line2split[bc_dict['Wobble_R2'][0]]}"
    run_subprocess(cmd)

    cloneList = os.listdir(os.path.join(args.tmpdir, "Preprocessing"))
    for clonefile in cloneList: 
        if ".1.fq.gz" in clonefile:
            read1 = os.path.join(args.tmpdir, "Preprocessing", clonefile)
        if ".2.fq.gz" in clonefile:
            read2 = os.path.join(args.tmpdir, "Preprocessing", clonefile)

cmd = f"process_radtags -1 {args.reads1} -2 {args.reads2} -b {os.path.join(args.outputdir, 'barcode_stacks.tsv')} -o {os.path.join(args.tmpdir, 'clone-stacks')} -r -D --inline_inline --renz_1 {line2split[bc_dict['ENZ_R1'][0]]} --renz_2 {line2split[bc_dict['ENZ_R2'][0]]} --retain_header --disable_rad_check"
run_subprocess(cmd)


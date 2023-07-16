import sys

if len(sys.argv) < 2:
    print("Please provide the gtf file ")
    print("Usage : python3 diploidgtfbuilder.py NameOfEnsembleGTFfile")
    sys.exit(1)

gtf=sys.argv[1] 
# gtf="examplegtf.gtf"

def parse_gtf(gtf_file,num,filepath):
    file_path = filepath
    with open(gtf_file, 'r') as file:
        with open(file_path, 'w') as outfile:
            for line in file:
                if line.startswith("#"):
                    continue
                data = line.strip().split("\t")
                temp=""
                chrom = data[0] + "_"+str(num)
                attributes = data[8].split(";")
                for item in attributes:
                    if item == "":
                        continue 
                    if item.startswith(" exon_number"):
                        temp = temp + item + ";"
                    else:
                        temp = temp + item[:-1] + "_"+str(num)+"\";" 

                yo= [chrom] + data[1:8] + [temp]

                outline = '\t'.join(yo)
                outfile.write(outline + '\n')

parse_gtf(gtf,1,"diploid1.gtf")
parse_gtf(gtf,2,"diploid2.gtf")

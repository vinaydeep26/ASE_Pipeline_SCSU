import time
import sys

start_end_time = time.time()

if len(sys.argv) < 3:
    print("Please provide the 2 file names in this order")
    print(" python3 pythonDiploidGenerator.py Vcffile.vcf Gennome.file" )
    sys.exit(1)


vcf_file= sys.argv[1] 
Genome_file = sys.argv[2]

def parse_vcf_file(vcf_file): 
    print("Parsing VCF file...")
    start_time = time.time()
    variants = []

    with open(vcf_file, 'r') as file:
        nonphasedPS=-1
        for line in file:
            line = line.strip()
            row = line.split('\t')
            if not row[0].startswith('#'):  # Skip header lines
                chromosome = row[0]
                position = int(row[1])
                ref= row[3]
                alt=row[4]
                genotype = row[9][:3] 

                variants.append([chromosome] + [position] + [ref] + [alt] + [genotype])

    print("Found " + str(len(variants)) + " Variants in the VCF file")
    end_time = time.time()
    print("Vcf file parsed. Took: " + str(end_time - start_time )+ "seconds")
    return variants

ParsedVariants=parse_vcf_file(vcf_file)
# print(ParsedVariants)

def parse_genome_file(genome_file): 
    print("Parsing Genome file...")
    start_time = time.time()
    genome_string=""

    with open(genome_file, 'r') as file:

        for line in file:

            if not line.startswith('>'):
                line=line.strip()
                genome_string+= line
            else:
                global CHROMO_NAME
                CHROMO_NAME = line

    end_time = time.time()
    print("Genome file parsed. Took: " + str(end_time - start_time )+ "seconds")
    return genome_string

G_string= parse_genome_file(Genome_file)
# print(G_string)



allele1=[]
allele2=[]
for sublist in ParsedVariants:
    pos= sublist[1]
    ref= sublist[2]
    alt=sublist[3].split(",")
    allgenotypesatthisposition = [ref]+alt
    # print(allgenotypesatthisposition)
    GT=sublist[4]
    indexGT1= int(GT[0])
    indexGT2= int(GT[2])

    getbase1= allgenotypesatthisposition[indexGT1]    
    getbase2= allgenotypesatthisposition[indexGT2]

    allele1.append([pos] + [len(ref)] + [getbase1])
    allele2.append([pos] + [len(ref)] + [getbase2])
# print(allele1[:3])


index=0
newstring=""


for item in allele1:

    pos= item[0]
    lengthtoskip= item[1]
    basetoadd=item[2]

    if pos <= len(G_string):
        newstring+= G_string[index:pos-1]
        newstring+=basetoadd
        # newstring+= " " + basetoadd + " "
        index= pos-1 + lengthtoskip
        
    else:
        newstring+= G_string[index:]
        break
newstring+= G_string[index:]

index2=0
newstring2=""

for item in allele2:
    
    pos= item[0]
    lengthtoskip= item[1]
    basetoadd=item[2]

    if pos <= len(G_string):
        newstring2+= G_string[index2:pos-1]
        newstring2+=basetoadd
        # newstring2+= " " + basetoadd + " "
        index2= pos-1 + lengthtoskip
        
    else:
        newstring2+= G_string[index2:]
        break
newstring2+= G_string[index2:]

fileone= str(Genome_file.split(".")[0]) + "_1.fa"
filetwo = str(Genome_file.split(".")[0]) + "_2.fa"

print("The outputs are in file:" + fileone + " " + filetwo)

header= CHROMO_NAME.rstrip('\n').split(" ")
header1= header[0]+"_1"
header2= header[0]+"_2"
modifiedheader1= " ".join([header1] + header[1:])
modifiedheader2= " ".join([header2] + header[1:])




with open(fileone, 'w') as file:
    file.write(modifiedheader1 + '\n')
    start = 0
    end = 100
    while start < len(newstring):
        file.write(newstring[start:end] + '\n')
        start = end
        end += 100


with open(filetwo, 'w') as file:
    file.write(modifiedheader2 + '\n')
    start1 = 0
    end1 = 100
    while start1 < len(newstring2):
        file.write(newstring2[start1:end1] + '\n')
        start1 = end1
        end1 += 100

final_end_time = time.time()
print("Whole builder Took: " + str(final_end_time - start_end_time )+ "seconds")
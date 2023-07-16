import time
import re
import sys
total_time = time.time()
# # Check if the correct number of arguments is provided
if len(sys.argv) < 4:
    print("Please provide the three file names in this order")
    print(" python coveragePhasing.py Vcffile.vcf CountsFile.txt( generatedbycountsgenerator) gene_file.gtf outputfilename" )
    sys.exit(1)


# #vcf_file = 'example.vcf'  
# # gtf_file = "examplegtf.gtf"
# # counts_file = 'base_counts2.txt'

vcf_file= sys.argv[1] 
counts_file = sys.argv[2]
gene_file = sys.argv[3]
outfile = sys.argv[4]

# vcf_file="example.vcf"
# counts_file = "base_countwithindel.txt"
# gene_file = "sample.bed"


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
                info = row[8].split(':')[-1]
                genotype = row[9][:3] 
                if info=="PS":
                    PStagval = int(row[9].split(':')[-1])
                    variants.append(row[:-1]+[genotype] + [PStagval])

                else:
                    variants.append(row[:-1]+[genotype] + [nonphasedPS])
                    nonphasedPS-=1
    print("Found " + str(len(variants)) + " Variants in the VCF file")
    end_time = time.time()
    print("Vcf file parsed. Took: " + str(end_time - start_time )+ "seconds")
    return variants

ParsedVariants=parse_vcf_file(vcf_file)



# for i in range(len(ParsedVariants)):    #just to print parsed vcf file
#     print(ParsedVariants[i])


def parse_counts_file(counts_file):
    print("Parsing Counts file...")
    start_time = time.time()
    counts = []
    with open(counts_file, 'r') as file:
        for line in file:
            line = line.strip()
            row = line.split('\t')
            chromosome = row[0]
            position = int(row[1])
            refbase=row[2]
            countdict = {}
            countdict["A"]= int(row[3])
            countdict["C"]= int(row[4])
            countdict["G"]= int(row[5])
            countdict["T"]= int(row[6])
            countdict["Indel"]= int(row[7])
            countdict["Ref"]= int(row[8])

            counts.append([chromosome,position,refbase,countdict])
    end_time = time.time()
    print("Counts file parsed. Took: " + str(end_time - start_time )+ "seconds")
    print("Found "+ str(len(counts)) + " counts")
    return counts


cont=parse_counts_file(counts_file)

start_time = time.time()
gene_data = {}
with open(gene_file, 'r') as file:
    no_of_genes=0
    for line in file:
        gene_info = line.split('\t')
        gene_id = gene_info[0]
        chrom = gene_info[1]
        start = int(gene_info[2])
        end = int(gene_info[3])
        if chrom not in gene_data:
            gene_data[chrom] = {}
        gene_data[chrom][gene_id] = {'start': start, 'end': end, 'lines': []}
        no_of_genes +=1
    # print(gene_data)
    # print(list(gene_data.items())[:4])
    # print(len(gene_data["Y"].items()))
    print("Total genes in file: " + str(no_of_genes))

#traverse the new file and group by variants found

with open(vcf_file, 'r') as file:
    for line_number, line in enumerate(file, start=-1):
        if line.startswith("#"):
                continue
        line_data = line.split('\t')
        chrom = line_data[0]
        pos = int(line_data[1])

        # Step 4: Find the gene for the current position using the gene data dictionary
        if chrom in gene_data:
            for gene_id, gene_info in gene_data[chrom].items():
                start = gene_info['start']
                end = gene_info['end']
                if start <= pos <= end:
                    gene_info['lines'].append(line_number)
                    break
        # line_number += 1
GENE_GROUPING_LIST=[]

for chrom, genes in gene_data.items():
    for gene_id, gene_info in genes.items():
        # print(f"Chromosome {chrom}, Gene {gene_id}:")
        if len(gene_info['lines'])!=0:
            GENE_GROUPING_LIST.append(gene_info['lines'])
end_time = time.time()
print("Gene Grpupings done. Took: " + str(end_time - start_time )+ "seconds")
print(GENE_GROUPING_LIST)

def phaseblock(startindex,endindex,movingsumleft,movingsumright):
    for x in range(startindex,endindex):
        # print(ParsedVariants[x])
        GT1= ParsedVariants[x][-2]
        GT=re.split(r'[/|]',GT1)

        if GT[0]==GT[1]:
            genotype= GT[0] +"/"+GT[1]
        else:
            if movingsumright > movingsumleft:
                genotype= GT[0] +"|"+GT[1]
            else:
                genotype= GT[1] +"|"+GT[0]
        ParsedVariants[x][-2]= genotype



# gene_line_numbers = group_variants_by_gene(vcf_file, genes)
# gene_line_numbers = group_variants_by_gene_try3(ParsedVariants, genes)

# Step 4: Group variants by gene
#create list of indices
indeces = GENE_GROUPING_LIST

phasedGenotypes=[]
for gene in indeces:                # to iterate over genes
    movingsumone = 0
    movingsumtwo=0
    i=0    #to calculate first occuring heterozygous genotype
    genelength =len(gene)

    while i < genelength:   #looping from the start of gene until we find a heterozygous genotype
        positions= gene[i]
        GT1= ParsedVariants[positions][-2]  # will have 0,1 1,2 or whatever the genotype is. [0,1]
        if GT1[0]!= GT1[2]:
            break
        i+=1
    if (i==len(gene)):
        continue
    firsthetvariantLocalindex = i      #THE FIRst heterozygous variant encountered
    firstpos= gene[firsthetvariantLocalindex]  #example 12
    # print("first heterozygous:" + str(firstpos))
    GT1= ParsedVariants[firstpos][-2]  # GT of the first het variant
    GT=re.split(r'[/|]',GT1)
    allgenotypesatthisposition = [ParsedVariants[firstpos][3]] + ParsedVariants[firstpos][4].split(",")
    isphased=GT1[1]   #should be / or |
    indexGT1= int(GT[0])
    indexGT2= int(GT[1])

    getbase1= allgenotypesatthisposition[indexGT1]    
    getbase2= allgenotypesatthisposition[indexGT2]        


    #if INDEL and phsaed
    
    if len(getbase1)>1 or len(getbase2)>1: 
        #can be 0/1 or 1|2 0|1                                                                     \
        if isphased=="|":
            countbase1=0
            countbase2=0
        else:
            countbase1 = cont[firstpos][3]["Ref"] 
            countbase2= cont[firstpos][3]["Indel"]

    #NOT INDEL
    else:
        countbase1 = cont[firstpos][3][getbase1]    
        countbase2= cont[firstpos][3][getbase2]                  



    CurrentPS = ParsedVariants[firstpos][-1] #currentPS      there's 2 scenarios for this 0/1 or 0|1
    # print("this starting current ps " + str(CurrentPS))
    movingsumone += countbase1    #either way we add the sums to the moving sums
    movingsumtwo += countbase2
    markedpos= gene[i]

    for j in range(firsthetvariantLocalindex+1,genelength):
        positions= gene[j]    #the actual positions of the lines corresponding to the variants in the vcf
        GT1= ParsedVariants[positions][-2]  # will have 0,1 1,2 or whatever the genotype is. [0,1]
        GT=re.split(r'[/|]',GT1)
        allgenotypesatthisposition = [ParsedVariants[positions][3]] + ParsedVariants[positions][4].split(",")  #to support multiallelic sites
        indexGT1= int(GT[0])  #would be 0 or 1 or 2
        indexGT2= int(GT[1]) #would be the second genotype

        if indexGT1==indexGT2:
            continue

        VariantPS = ParsedVariants[positions][-1]
        if VariantPS == CurrentPS:

            getbase1= allgenotypesatthisposition[indexGT1]       #getting the base of first genotype. for eg: A
            getbase2= allgenotypesatthisposition[indexGT2]       #getting the base of second genotype. for eg: T

            if len(getbase1)>1 or len(getbase2)>1:
                if isphased=="|":
                    countbase1=0
                    countbase2=0
                else:
                    #could only be 0/1
                    countbase1 = cont[firstpos][3]["Ref"] 
                    countbase2= cont[firstpos][3]["Indel"]

    #NOT INDEL
            else:
                countbase1 = cont[firstpos][3][getbase1]    
                countbase2= cont[firstpos][3][getbase2] 

            movingsumone += countbase1
            movingsumtwo += countbase2

        else:  #start of new PS
            #have a function here that takes in a start and end index and phases them according to counts
            
            phaseblock(markedpos,positions,movingsumone,movingsumtwo)   #phase from variant number x to y
            CurrentPS = ParsedVariants[positions][-1]
            markedpos= positions
            getbase1= allgenotypesatthisposition[indexGT1]       #getting the base of first genotype. for eg: A
            getbase2= allgenotypesatthisposition[indexGT2]       #getting the base of second genotype. for eg: T

            if len(getbase1)>1 or len(getbase2)>1:
                #can only be 0,1 
                if isphased=="|":
                    countbase1=0
                    countbase2=0
                else:
                    #could only be 0/1
                    countbase1 = cont[firstpos][3]["Ref"] 
                    countbase2= cont[firstpos][3]["Indel"]

    #NOT INDEL
            else:
                countbase1 = cont[firstpos][3][getbase1]    
                countbase2= cont[firstpos][3][getbase2] 


            movingsumone = countbase1
            movingsumtwo = countbase2
            #Mark the start of the ne variant block
            #Set movingSum1 and MovingSum2 according to the current variant data
    phaseblock(markedpos,gene[-1]+1,movingsumone,movingsumtwo)
    CurrentPS=0 
    #out of the last for loop, call the fucntion.
    # phasedGenotypes.append(tempphasearray)

# Open the file in write mode
with open(outfile, 'w') as file:
    # Iterate over the sublists
    for sublist in indeces:     # for [0,1,2] in [[0,1,2],[10,11,12]] 
        # Join the elements of the sublist with a tab space
        for index in sublist:


            line = '\t'.join(ParsedVariants[index][0:-1])
            # Write the line to the file
            file.write(line + '\n')

end_time = time.time()
print("Coverage BasedPhasing done. Took: " + str(end_time - total_time )+ "seconds")
print("VCF File that Contains Variants after CoverageBased Phasing : CoveragePhasedGenes.vcf ")
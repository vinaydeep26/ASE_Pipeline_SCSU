import time
import sys

if len(sys.argv) < 2:
    print("Please Provide VCF to filter")
    print("usage python3 Mu;tiallelicIndelFilter NameofVcfFile")
    sys.exit(1)

vcf_file= sys.argv[1] 

def MultiallelicVariantCoverageeliminator(vcf_file):
    print("Parsing VCF file...")
    variants=0
    start_time = time.time()
    with open(vcf_file, 'r') as file:
        with open("out.vcf", 'w') as outfile:
            for line in file:


                line = line.strip()
                row = line.split('\t')
                if not row[0].startswith('#'):
                    flag=0
                    ref= row[3]
                    alt=row[4].split(",")
                    genotype = row[9][:3]
                    G1=genotype[0]
                    G2=genotype[2]
                    Phaseindentifier = genotype[1]

                    #find if insertion or deletion
                    if len(ref)>1:
                        flag=1
                    for item in alt:
                        if len(item)>1:
                            flag=1

                    #if multiallelic and insertion/deletion

                    if flag and len(alt)!=1:

                        #   if heterozygous and non phased
                        if ((G1!=G2) and Phaseindentifier=='/'):
                            variants+=1
                            continue
                    
                    outfile.write(line + '\n')




                else: outfile.write(line + '\n')    #prints header lines
    end_time = time.time()
    print("Vcf file parsed. Took: " + str(end_time - start_time )+ "seconds")
    print("Skipped Variants :" + str(variants))

MultiallelicVariantCoverageeliminator(vcf_file)
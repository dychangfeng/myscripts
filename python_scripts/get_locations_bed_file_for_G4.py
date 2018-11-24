hg38_all = pd.read_table('/uufs/chpc.utah.edu/common/home/u0674686/ncbi_genomes/bed_files_for_all_genomes/all_G4/hg38_all.fa.all_G4.bed', header = None)

def reverse_comp(seq):
    bases_dict = {'A':'T', 'T':'A', 'G':'C', 'C':'G'}
    return ''.join(bases_dict.get(base,base) for base in seq[::-1])

def get_specific_G4(seq):
    """get the locations of the specific G4"""
    G4_loc=[]
    for j in range(len(hg38_all.iloc[:,6])):
        if seq==hg38_all.iloc[:,6].iloc[j] or reverse_comp(seq)==hg38_all.iloc[:,6].iloc[j]: ## add the reverse_complementary too
            G4_loc.append(hg38_all.iloc[:,3].iloc[j])
    return G4_loc
            

G4_34=get_specific_G4('GGGGACTGTTGTGGGGTGGGGGGAGGGGGGAGGG') ## find the locations of for this G4

G4_45_1=get_specific_G4('GGGGACTGTTGTGGGGTGGGGGGAGGGGGGAGGGATAGCATTGGG')

G4_18=get_specific_G4('GGGAGGGAGGTGGGGGGG')

G4_19=get_specific_G4('GGGAGGGAGGTGGGGGGGG')

def get_bed_specific_G4(seq_list):
    """get the bed format fo specific G4 sequence"""
    loc=[]

    for seq in seq_list:
        sl=seq.split(' ')
        chrom=sl[0]
        chrom=''.join(['chr',chrom])
        start=sl[3].split('_')[1]
        end=sl[3].split('_')[2]
        if sl[3].split('_')[3]=='rev':
            direction='-'
        else: direction='+'
        loc.append((chrom,start,end, direction))
    return loc
## save as bedfile        
pd.DataFrame.from_records(get_bed_specific_G4(G4_34)).to_csv("~/Desktop/G4_bed/G4_34.txt", sep= '\t', header=False, index=False )
pd.DataFrame.from_records(get_bed_specific_G4(G4_45_1)).to_csv("~/Desktop/G4_bed/G4_45_1.txt", sep= '\t', header=False, index=False )
pd.DataFrame.from_records(get_bed_specific_G4(G4_18)).to_csv("~/Desktop/G4_bed/G4_18.txt", sep= '\t', header=False, index=False )


pd.DataFrame.from_records(get_bed_specific_G4(G4_19)).to_csv("/uufs/chpc.utah.edu/common/home/u0674686/Desktop/G4_bed/G4_19.txt",
                                                            sep="\t", header=False, index=False)


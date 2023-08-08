'''
This script accept as input sorted sam file, a output sam file, and a mismatch rate,
then it will remove duplicates based on the barcode + UMI (edit distance <= 1),
and chromatin and start site, at the same time, it will output the duplication
number for each read.

Diego has refactored this code.
'''
from Levenshtein import distance
import sys

def rm_dup_samfile(samfile, output_file, mismatch):
    
    # initialize state
    seen_barcodes, dup_barcodes, chrom0, site0, dup_num = set(), set(), 0, 0, 0
    
    # new strategy of writing unique elements and saving only the barcode UMI,
    # then not writing if seen before or within the edit distance. 
    with open(samfile, 'r') as f1, open(output_file, 'w') as f2, open(output_file+'.csv', 'w') as f3:
        for line in f1:
            if line[0]=='@':
                f2.write(line)
            else:
                line2 = line.split()
                name = (line2[0]).split(',')
                barcode_UMI = name[0] + name[1]
                chrom1, site1 = line2[2], line2[3]
                
                # if the new site is the same as the past site
                if ((site1 == site0) and (chrom1 == chrom0)):
                    # and we've seen the exact barcode already or a barcode similar to others
                    if (barcode_UMI in seen_barcodes or barcode_UMI in dup_barcodes):
                        dup_num += 1
                    else:
                        # compute min dist from new barcode to existing barcodes
                        min_distance=min([distance(barcode_UMI, bc, score_cutoff=mismatch) \
                            for bc in seen_barcodes])
                        
                        # it's a dup if few edits away from existing barcode
                        if min_distance <= mismatch:
                            dup_barcodes.add(barcode_UMI)
                            dup_num += 1
                        else:
                            # it's not a dup so add to list of seen barcodes
                            # and write the read to the output file
                            seen_barcodes.add(barcode_UMI)
                            f2.write(line)
                        
                else:
                    # different site from past site then write read and dup count
                    f2.write(line)
                    if dup_num != 0:
                        f3.write("%d\n" % (dup_num))
                    
                    # reset state and pre-populate first read info
                    dup_num, chrom0, site0 = 1, chrom1, site1
                    seen_barcodes, dup_barcodes = set([barcode_UMI]), set([])
    f1.close()
    f2.close()
    f3.close()                

if __name__ == "__main__":
    samfile, output_file, mismatch = sys.argv[1], sys.argv[2], int(sys.argv[3])
    rm_dup_samfile(samfile, output_file, mismatch)

'''
validation test for HMM model output

    result:
        HMM model works 100 on the smallgene dataset
        would be same for larger gene
'''

import os
import argparse
import tarfile
import tempfile

import isoform
import isohint

parser = argparse.ArgumentParser(
    description='Integrated test for hmm')
parser.add_argument('file', type=str, metavar='<fastas dir>',
    help='dir for the small gene set')
parser.add_argument('hmm', type=str, metavar='<hmm path>',
    help='input hmm exe path')
parser.add_argument('--min_intron', required=False, type=int, default=35,
    metavar='<int>', help='minimum length of intron [%(default)i]')
parser.add_argument('--min_exon', required=False, type=int, default=25,
    metavar='<int>', help='minimum length exon [%(default)i]')
parser.add_argument('--flank', required=False, type=int, default=99,
    metavar='<int>', help='genomic flank on each side [%(default)i]')
arg = parser.parse_args()

def run(archive_path):
    """Process FASTA files directly from tar.gz archive"""
    
    total_files = 0
    passed_files = 0
    
    with tarfile.open(archive_path, 'r:gz') as tar:
        for member in tar.getmembers():
            if not member.name.endswith('.fa'):
                continue
                
            total_files += 1

            # Extract to memory
            fileobj = tar.extractfile(member)
            if fileobj:
                # Read content
                content = fileobj.read().decode('utf-8')
                
                with tempfile.NamedTemporaryFile(mode='w', suffix='.fa', delete=False) as tmp:
                    tmp.write(content)
                    tmp_path = tmp.name
                
                try:
                    name, seq = next(isoform.read_fasta(tmp_path))
                    dons1, accs1 = isohint.mgtag_sites(seq, arg.flank, arg.min_exon+1, 39+1)
                    result = isohint.run_hmm(arg.hmm, tmp_path)
                    dons2, accs2 = isohint.parse(result)
                    dons = [pos for pos, _ in dons2]
                    accs = [pos for pos, _ in accs2]
                    
                    if not (len(dons1) == len(dons) and
                            len(accs1) == len(accs) and
                            dons1 == dons and
                            accs1 == accs):
                        print(f"\nCheck hmm output. Splice site not match at {name}")
                        print("Dons hmm:", dons)
                        print("Dons ref:", dons1)
                        print("Accs hmm:", accs)
                        print("Accs ref:", accs1)
                        print(f"{seq}")
                    else:
                        passed_files += 1
                        
                finally:
                    # Clean up temp file
                    os.unlink(tmp_path)
    
    # Print summary
    print(f"\n{'='*50}")
    print(f"SUMMARY:")
    print(f"Total number of fasta processed: {total_files}")
    print(f"Files passed: {passed_files}")
    print(f"Files failed: {total_files - passed_files}")
    
    if total_files > 0:
        pass_percentage = (passed_files / total_files) * 100
        print(f"Pass rate: {pass_percentage:.2f}%")
    else:
        print("No files were processed.")
    print(f"{'='*50}")

run(arg.file)
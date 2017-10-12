import pandas as pd
import numpy as np
from argparse import ArgumentParser
from Bio import SeqIO, Entrez
from os import path, makedirs
def fetch_entrez(ifile, opath):
    name = ifile.split('/')[-1].split('.')[0]
    opath = path.abspath(opath) + '/' + name
    ofilefas = '%s/%s.fasta' % (opath, name)
    if not path.exists(opath):
        makedirs(opath)
    else:
        pass
        #raise IOError('Path existed.')

    dfRec = pd.read_csv(ifile, sep='\t', header=1)
    num = 20 #dfRec.shape[0]
    step = 5
    starts = np.arange(0, num, step)
    Entrez.email = "cc3423@drexel.edu"
    dat = []
    with open(ofilefas, 'w') as hands:
        for s in starts:
            if s == starts[-1]:
                end=num
            else:
                end=s + step
            handle = Entrez.efetch(db="nucleotide", rettype="fasta",
                                    retmode="xml",
                                   id=tuple(dfRec[s:end]['Accession']))
            records = Entrez.parse(handle)
            for rec in records:
                dat.append({'GI number': int(rec['TSeq_gi']),
                            'Sequence': rec['TSeq_sequence']})
                hands.write('>%s\n' % rec['TSeq_gi'])
                hands.write('%s\n' % rec['TSeq_sequence'])
    dfSeq = pd.DataFrame(dat)
    dfRec = pd.merge(dfRec, dfSeq, on='GI number', how='outer')
    ofile = '%s/%s.csv' % (opath, name)
    dfRec.to_csv(ofile, sep='\t', header=True)




if __name__ == '__main__':
    parser = ArgumentParser()
    parser.add_argument('ifile', nargs='+',
                         help = 'Input file has to be a list of accession numbers from LANL database')
    parser.add_argument('--out', dest='out')
    args = parser.parse_args()
    for f in args.ifile:
        fetch_entrez(f, args.out)
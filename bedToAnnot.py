import pandas as pd
import numpy as np
import os
import argparse

def make_df_annot(bfile,bedfile):
    df_bim = pd.read_csv(bfile + '.bim',delim_whitespace=True,usecols = [0,1,2,3],
        names = ['CHR','SNP','CM','BP'])
    df_bed = df_bim[['BP']].copy()
    chr = ['chr'+str(x) for x in df_bim['CHR']]
    df_bed['CHR'] = chr
    df_bed = df_bed[['CHR','BP','BP']]
    ucscbed_file = 'temp'+str(np.random.randint(1000000))+'.out'
    df_bed.to_csv(ucscbed_file, sep='\t', header=False, index=False)
    outname = 'temp'+str(np.random.randint(1000000))+'.out'
    os.system('~/bin/bedtools/bedtools intersect -a ' + ucscbed_file + ' -b '+bedfile+ ' > '+outname)
    if open(outname).readline() == '':
        df_int = pd.DataFrame(columns = ['BP', 'ANNOT'])
    else:
        df_int = pd.read_csv(outname, delim_whitespace=True, usecols = [1], names=['BP'])
        df_int['ANNOT'] = 1
    df_annot = pd.merge(df_bim,df_int,how='left',on='BP')
    df_annot.fillna(0, inplace=True, downcast='infer')
    df_annot['OTHER'] = 1 - df_annot['ANNOT']
    df_annot = df_annot[['CHR','BP','SNP','CM','ANNOT','OTHER']]
    os.system('rm '+outname)
    os.system('rm '+ucscbed_file)
    return df_annot

def make_annot(bfile,bedfile,annotfile):
    df_annot = make_df_annot(bfile,bedfile)
    df_annot.to_csv(annotfile,sep='\t',index=False)

def make_baseline(bfile,annotfile,bedname_file):

    bedname_list = [x.rstrip('\n') for x in open(bedname_file,'r').readlines()]
    if bedname_list[0] == 'base':
        df = make_df_annot(bfile,bedname_list[1])
    else:
        df = make_df_annot(bfile,bedname_list[0])
    df = df.iloc[:,0:-2]

    for bedfile in bedname_list:
        print bedfile
        if bedfile == 'base':
            df['base'] = 1
        else:
            df_single = make_df_annot(bfile,bedfile)
            print np.mean(df_single.iloc[:,-2])
            new_name = bedfile.split('/')[-1]
            df[new_name] = df_single.iloc[:,-2]
    df.to_csv(annotfile,sep = '\t',index=False)

def get_comp(annot):
    if annot.endswith('.gz'):
        return 'gzip'
    else:
        return None

def remove_extra(df):
    unnecessary = set(['CHR','SNP','BP','CM','MAF']) & set(df.columns)
    df.drop(unnecessary,axis=1,inplace=True)

def combine(annot1,drop1,annot2,drop2,combined_name):
    compression1 = get_comp(annot1)
    compression2 = get_comp(annot2)
    df1 = pd.read_csv(annot1,delim_whitespace=True,compression=compression1)
    if drop1 is not None:
        to_drop = df1.columns[drop1]
        del df1[to_drop]
    df2 = pd.read_csv(annot2,delim_whitespace=True,compression=compression2)
    if drop2 is not None:
        to_drop = df2.columns[drop2]
        del df2[to_drop]
    remove_extra(df2)
    df_new = pd.concat([df1,df2],axis=1)
    df_new.to_csv(combined_name,sep='\t',index=False)

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--bedfile')
    parser.add_argument('--bedname-file')
    parser.add_argument('--name')
    parser.add_argument('--dataset')
    parser.add_argument('--combine1')
    parser.add_argument('--drop1',type=int)
    parser.add_argument('--combine2')
    parser.add_argument('--annotfile')
    parser.add_argument('--drop2',type=int)

    parser.add_argument('--bedfile-single')
    parser.add_argument('--annotfile-slim')
    parser.add_argument('--bfile')
    args = parser.parse_args()

    if args.bedname_file:
        make_baseline(args.bfile,args.annotfile,args.bedname_file)
    elif args.combine1:
        combine(args.combine1,args.drop1,args.combine2,args.drop2,args.annotfile)
    elif args.name:
        for chrom in range(1,23):
            print chrom
            bfile = '../data/'+args.dataset+'.'+str(chrom)
            annotfile = '../sim/ref/'+args.dataset+'.'+args.name+'.'+str(chrom)+'.annot'
            bedname_file = '../sim/ref/'+args.name+'.txt'
            make_baseline(bfile, annotfile,bedname_file)
    elif args.bedfile_single:
        df_annot = make_df_annot(args.bfile,args.bedfile_single)
        if args.annotfile is not None:
            df_annot = df_annot.ix[:,0:-1]
            out = args.annotfile
        elif args.annotfile_slim is not None:
            df_annot = df_annot.ix[:, [-2]]
            out = args.annotfile_slim
        if out.endswith('.gz'):
            df_annot.to_csv(out[0:-3], sep='\t', index=False)
            os.system('gzip -f '+out[0:-3])
        else:
            df_annot.to_csv(out, sep='\t', index=False)
    else:
        make_annot(args.bfile,args.bedfile,args.annotfile)



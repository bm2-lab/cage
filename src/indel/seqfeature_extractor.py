
from __future__ import division
import os
from pyfasta import Fasta



def __SeqRc(str_seq):
    lst_rcseq = [c for c in str_seq[::-1]]
    for i in xrange(len(lst_rcseq)):
        if lst_rcseq[i] == 'A':
            lst_rcseq[i] = 'T'
        elif lst_rcseq[i] == 'T':
            lst_rcseq[i] = 'A'
        elif lst_rcseq[i] == 'C':
            lst_rcseq[i] = 'G'
        elif lst_rcseq[i] == 'G':
            lst_rcseq[i] = 'C'
        elif lst_rcseq[i] == 'N':
            lst_rcseq[i] = 'N'
    return ''.join(lst_rcseq)

def __CalculateGC(str_seq):
    int_gc = 0
    for c in str_seq:
        if c in ('G', 'C'):
            int_gc += 1
    flt_gc = int_gc / len(str_seq)
    return flt_gc

def __EncodeSeq(str_seq, dict_code):
    lst_code = []
    for c in str_seq:
        lst_code.extend(dict_code[c])
    return lst_code

def __FormulateFeature(lst_sg, fa_ref, dict_code, int_ups, int_dws):
    idx_sgid = 0
    idx_chr = 1
    idx_strand = 2
    idx_beg = 4
    idx_end = 5
    idx_seq = 6

    lst_sgtab = [lst_sg[idx_sgid]]
    str_strand = '1' if lst_sg[idx_strand] == '+' else '0'
    #lst_sgtab.append(str_strand)
    str_ups =''
    str_dws = ''

    if lst_sg[idx_strand] == '-':
        int_bups = int(lst_sg[idx_beg]) + int_ups
        int_eups = int(lst_sg[idx_beg]) + 1
        str_ups = fa_ref[lst_sg[idx_chr]][int_eups - 1:int_bups].upper()
        str_ups = __SeqRc(str_ups)
        int_bdws = int(lst_sg[idx_end]) - 1
        int_edws = int(lst_sg[idx_end]) - int_dws
        str_dws = fa_ref[lst_sg[idx_chr]][int_edws - 1:int_bdws].upper()
        str_dws = __SeqRc(str_dws)
    else:
        int_bups = int(lst_sg[idx_beg]) - int_ups
        int_eups = int(lst_sg[idx_beg]) - 1
        str_ups = fa_ref[lst_sg[idx_chr]][int_bups - 1:int_eups].upper()
        int_bdws = int(lst_sg[idx_end]) + 1
        int_edws = int(lst_sg[idx_end]) + int_dws
        str_dws = fa_ref[lst_sg[idx_chr]][int_bdws - 1:int_edws].upper()
    lst_sgtab.extend(__EncodeSeq(str_ups, dict_code))
    lst_sgtab.extend(__EncodeSeq(lst_sg[idx_seq], dict_code))
    lst_sgtab.extend(__EncodeSeq(str_dws, dict_code))
    return lst_sgtab
    
    

def ExtractSeqFeature(str_f_sg, str_refgem, str_of_seq, int_ups=30, int_dws=27):
    f_sg = open(str_f_sg, 'r')
    gn_sg = (str_sg.strip().split('\t') for str_sg in f_sg if str_sg.strip() != '')
    f_sgtab = open(str_of_seq, 'w')

    str_refpath = '%s/%s.fa'% (os.environ['FASTADB'], str_refgem)
    fa_ref = Fasta(str_refpath)
    dict_code = dict(N=['0','0','0','0'], A=['1','0','0','0'], C=['0','1','0','0'], G=['0','0','1','0'], T=['0','0','0','1'])
    lst_code = ['A', 'C', 'G', 'T']
    lst_sgtab = []
    #lst_header = ['sgID', 'strand']
    lst_header = ['sgID']
    lst_header.extend(['ups_%d_%s'% (i, j) for i in range(int_ups, 0, -1) for j in lst_code])
    lst_header.extend(['spa_%d_%s'% (i, j) for i in range(1, 21) for j in lst_code])
    lst_header.extend(['pam_%d_%s'% (i, j) for i in range(1, 4) for j in lst_code])
    lst_header.extend(['dws_%d_%s'% (i, j) for i in range(1, int_dws+1) for j in lst_code])
    f_sgtab.write('\t'.join(lst_header) + '\n')
    for lst_sg in gn_sg:
        lst_sgtab = __FormulateFeature(lst_sg, fa_ref, dict_code, int_ups, int_dws)
        f_sgtab.write('\t'.join(lst_sgtab) + '\n')
    f_sgtab.close()

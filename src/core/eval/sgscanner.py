import os
import re
from pyfasta import Fasta
    
    

def __ParseSeq(str_chr, int_beg, int_end, str_strand, fa_ref, re_pat):
    str_seq = fa_ref.sequence({'chr':str_chr, 'start':int_beg, 'stop':int_end, 'strand':str_strand}).upper()
    lst_pam = [range(re_obj.start(), re_obj.end()-1) for re_obj in re_pat.finditer(str_seq)]
    lst_pamcord = reduce(lambda x,y: x+y, lst_pam)
    lst_sg = []
    int_sgbeg = 0
    int_sgend = 0
    for idx in lst_pamcord:
        if str_strand == '+':
            int_sgbeg = int_beg + idx - 21
            int_sgend = int_beg + idx + 1
        elif str_strand == '-':
            int_sgbeg = int_end - idx - 1
            int_sgend = int_end - idx + 21
            
        str_sg = fa_ref.sequence({'chr':str_chr, 'start':int_sgbeg, 'stop':int_sgend, 'strand':str_strand}).upper()
        yield str(str_sg)
               
def ExtractSg(str_refgem, str_chr, int_beg, int_end, str_drct):
    str_refpath = '%s/%s.fa'% (os.environ['FASTADB'], str_refgem)
    fa_ref = Fasta(str_refpath)
    re_pat = re.compile('G+')
    gn_f_sg = __ParseSeq(str_chr, int_beg, int_end, '+', fa_ref, re_pat)
    gn_r_sg = __ParseSeq(str_chr, int_beg, int_end, '-', fa_ref, re_pat)
    if str_drct == 'f':
        for str_sg in gn_f_sg:
            yield str_sg
    elif str_drct == 'r':
        for str_sg in gn_r_sg:
            yield str_sg
    elif str_drct == 'b':
        for str_sg in gn_f_sg:
            yield str_sg
        for str_sg in gn_r_sg:
            yield str_sg


def FaToFq(gn_sg, str_of_fq):
    f_fq = open(str_of_fq, 'w')
    for i, str_sg in enumerate(gn_sg):
        str_id = '@sg%d'% (i+1)
        str_seq = str_sg
        str_sign = '+'
        str_qua = '<' * len(str_seq)
        f_fq.write(str_id + '\n')
        f_fq.write(str_seq + '\n')
        f_fq.write(str_sign + '\n')
        f_fq.write(str_qua + '\n')

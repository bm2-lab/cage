
import subprocess
from itertools import dropwhile
import re

def CallBWA(str_f_fwdfq, str_f_revfq, str_refgem, str_of_psam, 
            bln_long=True, int_thrds=1):
    str_f_revfq = str_f_revfq.strip()
    str_bwa_l = r'bwa mem -t %d ${BWADB}/%s.fa %s %s > %s'% (int_thrds, str_refgem, str_f_fwdfq, str_f_revfq, str_of_psam)

    func_bwa_aln = lambda str_fq, str_ref=str_refgem, int_td=int_thrds: r'bwa aln -t %d ${BWADB}/%s.fa %s'% (int_td, str_ref, str_fq)
    
    str_bwa_s = ''
    if str_f_revfq == '':
        str_bwa_s = r'bwa samse ${BWADB}/%s.fa <(%s) %s > %s'% (str_refgem, func_bwa_aln(str_f_fwdfq), str_f_fwdfq, str_of_psam)
    else:
        str_bwa_s = r'bwa sampe ${BWADB}/%s.fa <(%s) <(%s) %s %s > %s'% (str_refgem, func_bwa_aln(str_f_fwdfq), func_bwa_aln(str_f_revfq), str_f_fwdfq, str_f_revfq, str_of_psam)

    if bln_long == True:
        int_status = subprocess.call(str_bwa_l, shell=True, executable='/bin/bash')
    else:
        int_status = subprocess.call(str_bwa_s, shell=True, executable='/bin/bash')
    return int_status

def __TestSam(lst_rc, re_pat, bln_pair):
    idx_flag = 1
    idx_quality = 4
    idx_cigar = 5
    idx_pair = 6

    if int(lst_rc[idx_quality]) < 10:
        return False

    if re_pat.search(lst_rc[idx_cigar])is not None:
        return False

    if bln_pair == True:
        if int(lst_rc[idx_flag]) not in (83, 163, 99, 147):
            return False
        if lst_rc[idx_pair] != '=':
            return False
    else:
        if int(lst_rc[idx_flag]) not in (0, 16):
            return False
    return True
    

def FilterSam(str_f_psam, str_of_sam, bln_pair=True):
    f_psam = open(str_f_psam, 'r')
    g_f_psam = (line.strip().split('\t') for line in dropwhile(lambda line: line.startswith('@'), f_psam) if line.strip() != '')
    f_sam = open(str_of_sam, 'w')
    
    idx_idx = 0
    re_cigar = re.compile('H')
    g_psam = (lst_rc for lst_rc in g_f_psam if __TestSam(lst_rc, re_cigar, bln_pair) == True)

    str_prevrc = ''
    if bln_pair == True:
        lst_fwd = []
        lst_rev = []
        bln_rc = False
        for lst_rc in g_psam:
            if bln_rc == False:
                str_prevrc = lst_rc[idx_idx]
                lst_fwd = lst_rc
                bln_rc = True
            else:
                if lst_rc[idx_idx] != str_prevrc:
                    str_prevrc = lst_rc[idx_idx]
                    lst_fwd = lst_rc
                else:
                    lst_rev = lst_rc
                    bln_rc = False
                    f_sam.write('\t'.join(lst_fwd) + '\n')
                    f_sam.write('\t'.join(lst_rev) + '\n')
    else:
        for lst_rc in g_psam:
            if lst_rc[idx_idx] == str_prevrc:
                continue
            else:
                str_prevrc = lst_rc[idx_idx]
                f_sam.write('\t'.join(lst_rc) + '\n')
    f_sam.close()
                
                
            
    
        
    

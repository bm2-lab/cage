
from lxml import etree
from collections import namedtuple

Fs = namedtuple('Fs', ['reg', 'cord', 'nt'])
Cord = namedtuple('Cord', ['cord', 'tick'])
Fsc = namedtuple('Fsc', ['ups', 'spa', 'pam', 'dws'])

def __ParseFsRep(str_f_fesrep):
    etree_fes = etree.parse(str_f_fesrep)
    enode_root = etree_fes.getroot()
    enode_fes = enode_root[0]
    lst_fs = [ed.tag for ed in enode_fes]
    int_ups = int(enode_fes.get('ups'))
    int_dws = int(enode_fes.get('dws'))
    return (lst_fs, int_ups, int_dws)

def __ParseFs(lst_fs):
    dict_nt = dict(A=1, C=2, G=3, T=4)
    lst_ups = []
    lst_spa = []
    lst_pam = []
    lst_dws = []
    gn_lfs = (lfs.split('_') for lfs in lst_fs)
    gn_fs = (Fs(reg=fs[0], cord=int(fs[1]), nt=fs[2]) for fs in gn_lfs)
    for fs in gn_fs:
        if fs.reg == 'ups':
            lst_ups.append(fs)
        elif fs.reg == 'spa':
            lst_spa.append(fs)
        elif fs.reg == 'pam':
            lst_pam.append(fs)
        elif fs.reg == 'dws':
            lst_dws.append(fs)
    cord_last = 0
    lst_ups_cord = []
    lst_ups_tick = []
    i = 0
    for fs in lst_ups:
        if fs.cord != cord_last:
            i += 1
            lst_ups_tick.append('%d/%d'% (fs.cord, i))
            cord_last = fs.cord
        lst_ups_cord.append('(%s-%d-%d.center) {%s}'% (fs.reg, dict_nt[fs.nt], i, fs.nt))
        
    cord_last = 0
    lst_spa_cord = []
    lst_spa_tick = []
    i = 0
    for fs in lst_spa:
        if fs.cord != cord_last:
            i += 1
            lst_spa_tick.append('%d/%d'% (fs.cord, i))
            cord_last = fs.cord
        lst_spa_cord.append('(%s-%d-%d.center) {%s}'% (fs.reg, dict_nt[fs.nt], i, fs.nt))

    cord_last = 0
    lst_pam_cord = []
    lst_pam_tick = []
    i = 0
    for fs in lst_pam:
        if fs.cord != cord_last:
            i += 1
            lst_pam_tick.append('%d/%d'% (fs.cord, i))
            cord_last = fs.cord
        lst_pam_cord.append('(%s-%d-%d.center) {%s}'% (fs.reg, dict_nt[fs.nt], i, fs.nt))

    cord_last = 0
    lst_dws_cord = []
    lst_dws_tick = []
    i = 0
    for fs in lst_dws:
        if fs.cord != cord_last:
            i += 1
            lst_dws_tick.append('%d/%d'% (fs.cord, i))
            cord_last = fs.cord
        lst_dws_cord.append('(%s-%d-%d.center) {%s}'% (fs.reg, dict_nt[fs.nt], i, fs.nt))
        
    return Fsc(ups=Cord(cord=lst_ups_cord, tick=lst_ups_tick),
               spa=Cord(cord=lst_spa_cord, tick=lst_spa_tick),
               pam=Cord(cord=lst_pam_cord, tick=lst_pam_tick),
               dws=Cord(cord=lst_dws_cord, tick=lst_dws_tick))
            
def VisualizeFeature(str_f_fesrep, str_of_lax, flt_ampfct=1.00):
    lst_fs, int_ups, int_dws = __ParseFsRep(str_f_fesrep)
    fsc = __ParseFs(lst_fs)
    func_mat = lambda x: r'\rowzero' if x==0 else (r'\rowone' if x==1 else r'\row{%d}'%x)
    lst_lax = []
    str_prem = r'''
\documentclass[tikz, border=10pt]{standalone}

\usepackage{tikz}
\usetikzlibrary{arrows}
\usetikzlibrary{calc}
\usetikzlibrary{matrix}
\usetikzlibrary{math}
\usepackage{etoolbox}

\begin{document}

\newcommand{\row}[1]
{
  \tikzmath{int \x; \x = #1 - 1;}
  \foreach \i in {1,...,4}
  {
  \foreach \i in {1,...,\x} 
  {
    \gappto\mat{\expandonce{
    A \&
    }}
  }
  \gappto\mat{A \\}
  }
}

\newcommand{\rowone}
{
  \foreach \i in {1,...,4}
  {
  \gappto\mat{A \\}
  }
}

\newcommand{\rowzero}
{
  \foreach \i in {1,...,4}
  {
  \gappto\mat{\\}
  }
}
'''
    lst_lax.append(str_prem)
    
    str_pd = r'''
%% Parameter Definition 
\tikzmath
{
  coordinate \begpt;
  \begpt = (0,0);
  let \lbsep = 1ex;
  let \regsep = 5ex;
  let \axsep = 2em;
  let \fct = %f;
}
'''% flt_ampfct
    lst_lax.append(str_pd)

    str_preups = r'''
\begin{tikzpicture}
  \tikzstyle{sty_nt}=[fill=white, font=\sffamily\footnotesize]
  \tikzstyle{sty_pmx}=[nodes={rectangle, fill=black, inner ysep=1pt, inner
    xsep=0, font=\ttfamily\footnotesize}]
  \tikzstyle{sty_mx}=[matrix of nodes, sty_pmx, row sep=1ex, column sep=1ex,
  ampersand replacement=\&, inner sep=0, anchor=north west]
  \tikzstyle{sty_lb}=[align=right, font=\ttfamily, inner sep=0,
  anchor=east, rotate=90]
'''
    lst_lax.append(str_preups)
    
    str_ups_cord = '\n'.join([r'\node[sty_nt] at %s;'%c for c in fsc.ups.cord])
    str_ups_tick = ','.join(fsc.ups.tick)
    
    str_ups = r'''
  \begin{scope}  
  \let\mat\empty
  %s
  \matrix (ups) [sty_mx, nodes={fill=blue, text=blue}] at (\begpt)
  {
    \mat
  };
  %s
  \foreach \x/\y in {%s}
  {  
    \node[sty_lb] at ($(ups-4-\y.south)-(0,\lbsep)$) {\footnotesize{-\x}};
  }
  \end{scope}
'''% (func_mat(len(fsc.ups.tick)), str_ups_cord, str_ups_tick)
    lst_lax.append(str_ups)

    str_spa_cord = '\n'.join([r'\node[sty_nt] at %s;'%c for c in fsc.spa.cord])
    str_spa_tick = ','.join(fsc.spa.tick)
    str_spa = r'''
  \tikzmath
  {
    coordinate \cspa;
    \cspa = (ups.north east);
    \cspa = (\cspa) + (\regsep,0);
  }
  \begin{scope}
  \let\mat\empty
  %s
  \matrix (spa) [sty_mx, nodes={fill=red, text=red}] at (\cspa)
  {
    \mat
  };
  %s
  \foreach \x/\y in {%s}
  {  
    \node[sty_lb] at ($(spa-4-\y.south)-(0,\lbsep)$) {\footnotesize{\x}};
  }
  \end{scope}
'''% (func_mat(len(fsc.spa.tick)), str_spa_cord, str_spa_tick)
    lst_lax.append(str_spa)

    str_pam_cord = '\n'.join([r'\node[sty_nt] at %s;'%c for c in fsc.pam.cord])
    str_pam_tick = ','.join(fsc.pam.tick)
    str_pam = r'''
  \tikzmath
  {
    coordinate \cpam;
    \cpam = (spa.north east);
    \cpam = (\cpam) + (\regsep,0);
  }
  \begin{scope}
  \let\mat\empty
  %s
  \matrix (pam) [sty_mx, nodes={fill=orange, text=orange}] at (\cpam)
  {
    \mat
  };
  %s
  \foreach \x/\y in {%s}
  {
    
    \node[sty_lb] at ($(pam-4-\y.south)-(0,\lbsep)$) {\footnotesize{\x}};
  }
  \end{scope}
'''% (func_mat(len(fsc.pam.tick)), str_pam_cord, str_pam_tick)
    lst_lax.append(str_pam)

    str_dws_cord = '\n'.join([r'\node[sty_nt] at %s;'%c for c in fsc.dws.cord])
    str_dws_tick = ','.join(fsc.dws.tick)
    str_dws = r'''
  \tikzmath
  {
    coordinate \cdws;
    \cdws = (pam.north east);
    \cdws = (\cdws) + (\regsep,0);
  }
  \begin{scope}
  \let\mat\empty
  %s
  \matrix (dws) [sty_mx, nodes={fill=cyan, text=cyan}] at (\cdws)
  {
    \mat
  };
  %s
  \foreach \x/\y in {%s}
  {
    
    \node[sty_lb] at ($(dws-4-\y.south)-(0,\lbsep)$) {\footnotesize{\x}};
  }
  \end{scope}
'''% (func_mat(len(fsc.dws.tick)), str_dws_cord, str_dws_tick)
    lst_lax.append(str_dws)

    str_ups_line = '' if len(fsc.ups.tick)==0 else r'''
  \draw[dline, blue] (\bups) -- (ups.north west);
  \draw[dline, blue] (\bspa) -- (ups.north east);
'''
    str_spa_line = '' if len(fsc.spa.tick)==0 else r'''
  \draw[dline, red] (\bspa) -- (spa.north west);
  \draw[dline, red] (\bpam) -- (spa.north east);
'''
    str_pam_line = '' if len(fsc.pam.tick)==0 else r'''
  \draw[dline, orange] (\bpam) -- (pam.north west);
  \draw[dline, orange] (\bdws) -- (pam.north east);
'''
    str_dws_line = '' if len(fsc.dws.tick)==0 else r'''
  \draw[dline, cyan] (\bdws) -- (dws.north west);
  \draw[dline, cyan] (\edws) -- (dws.north east);
'''
    str_suf = r'''
%% Axis

  \tikzmath
  {
    coordinate \cinit, \bups, \bspa, \bpam, \bdws, \edws;
    \cinit = (ups.north west) + (0,\axsep);
    coordinate \tmpa, \tmpb, \tmpm;
    \tmpa = (ups.north west);
    \tmpb = (dws.north east);
    \tmpm = 1/2*(\tmpa) + 1/2*(\tmpb);
    %%let \tmp = \tmpbx - \tmpax;
    let \tmp = 420pt;
    \bups = (\cinit);
    \bspa = (\bups) + \fct*%d/90*(\tmp,0);
    \bpam = (\bspa) + \fct*20/90*(\tmp,0);
    \bdws = (\bpam) + \fct*3/90*(\tmp,0);
    \edws = (\bdws) + \fct*%d/90*(\tmp,0);
    coordinate \cmid;
    \cmid = 1/2*(\bups) + 1/2*(\edws);
    let \dtmp = \tmpmx - \cmidx;
    \bups = (\bups) + (\dtmp,0);
    \bspa = (\bspa) + (\dtmp,0);
    \bpam = (\bpam) + (\dtmp,0);
    \bdws = (\bdws) + (\dtmp,0);
    \edws = (\edws) + (\dtmp,0);
  }
  \tikzstyle{line}=[ very thick]
  \tikzstyle{axlb1}=[midway, above=3ex, sloped]
  \tikzstyle{axlb2}=[midway, above=0.5ex, sloped]
  \draw[line,blue] (\bups) -- (\bspa)
  node[axlb1] {Upstream}  node[axlb2] {%dnt};
  \draw[line,red] (\bspa) -- (\bpam)
  node[axlb1] {Spacer} node[axlb2] {20nt};
  \draw[line,orange] (\bpam) -- (\bdws)
  node[axlb1] {PAM} node[axlb2] {3nt};
  \draw[line,cyan] (\bdws) -- (\edws)
  node[axlb1] {Downstream} node[axlb2] {%dnt};

  \tikzstyle{dline}=[dotted,thick]
  %s
  %s
  %s
  %s
\end{tikzpicture}
\end{document}
'''% (int_ups, int_dws, int_ups, int_dws, str_ups_line, str_spa_line, str_pam_line, str_dws_line)
    lst_lax.append(str_suf)
    f_lax = open(str_of_lax, 'w')
    f_lax.write('\n'.join(lst_lax))
    f_lax.close()

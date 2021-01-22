#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu Jan 30 17:50:43 2020

@author: jackaraz
contact: jackaraz@gmail.com
"""
from CutFlowReader import *
from FoM import FoM
import os


class CutFlowTable:
    def __init__(self, *args,**kwargs):
        """
        Transforms MadAnalysis 5 CutFlows into LaTeX table.

        Parameters
        ----------
        *args : list of SR Collection
            This list contains SR collections i.e. background and signal. It can 
            have multiple collections but all collections has to have same cutflow.
        **kwargs : 
            ref_sample : INT
                The index of the reference sample in the SR collection.
            sample_names : LIST
                Names of the samples.
            notes : STR
                Notes to be written in the caption. Default ''
            SR_list : LIST
                List of the SRs to be written. Default all in the ref. input.

        Returns
        -------
        None.
        """
        samples = [x for x in args if type(x) == Collection]
        sample_names = kwargs.get('sample_names',[])
        if len(sample_names) == len(samples):
            self.sample_names = sample_names
        else:
            self.sample_names = ['Sample '+str(x) for x in range(len(samples))]
        self.SR_list = kwargs.get('SR_list',[])
        self.notes   = kwargs.get('notes','')
        ref_sample = kwargs.get('ref_sample',0)
        self.ref_name = self.sample_names[ref_sample]
        self.ref_sample = samples[ref_sample]
        samples.remove(self.ref_sample)
        self.sample_names.remove(self.ref_name)
        self.samples = samples


    def _sorter(self,x):
        return self.ref_sample[x].get_final_cut().Nentries


    def write_comparison_table(self,*args,**kwargs):
        """
        Writes sample comparison table.

        Parameters
        ----------
        *args : FILE
            Optional, if there is a file input, tables will be written in the 
            file otherwise all will be printed on the screen.
        **kwargs : 
            only_alive : BOOLEAN (default True)
                only write the SRs which has more than zero yield for reference
                collection.
            make : BOOL
                Write the Makefile -> (default, True)
            raw  : BOOL optional
                Generate table with raw number of entries. Default False.
            event_style : STR optional
                Decimal style of the events, default '{:.1f}'
            eff_style : STR optional
                Decimal style of the efficiencies, default '{:.3f}'
            ratio_style : STR optional
                Decimal style of the ref/input ratio, default '{:.1f}'
            mcunc : BOOL
                Monte Carlo uncertainty of the cut efficiency. Default False.
            finalMCunc : BOOL
                Write Monte Carlo uncertainty for the last cut. Default False.

        Returns
        -------
        LaTeX tables of signal regions.
        """
        if self.SR_list == []:
            SR_list = self.ref_sample.keys()
            if kwargs.get('only_alive',True): 
                SR_list = [x for x in SR_list if self.ref_sample[x].isAlive()]
            SR_list.sort(key=self._sorter, reverse=True)
        else:
            SR_list = self.SR_list

        # Generate table with number of entries
        raw = kwargs.get('raw',False)
        # Get table style
        event_style = kwargs.get('event_style','{:.1f}')
        if raw: event_style = '{:.0f}'
        eff_style   = kwargs.get('eff_style','{:.3f}')
        ratio_style = kwargs.get('ratio_style','{:.1f}')
        MCunc       = kwargs.get('mcunc',False)
        finalMCunc  = kwargs.get('finalMCunc',False)

        TeX=None
        if any([x for x in args if isinstance(x,file)]):
            TeX = [x for x in args if isinstance(x,file)][0]
            TeX.write(r'\documentclass[12pt]{article}'+'\n'+\
                       r'\usepackage{pdflscape,slashed}'+'\n'+\
                       r'\begin{document}'+'\n'+\
                       r'\begin{landscape}'+'\n\n\n\n'+\
                       '%%%%%% \\delta := |Ref. smp - smp_i| / ref_smp\n\n\n')
            for line in self.notes.split('\n'):
                TeX.write('%%%% '+line+'\n')
            if MCunc:
                TeX.write('\n%%%% MC Unc = Nevt * sqrt((1-eff)/NMC)\n')
            TeX.write('\n\n\n\n')

        for SR in SR_list:
            txt = '\n\n%% '+SR+'\n\n'
            txt+='\\begin{table}[h]\n'
            txt+='  \\begin{center}\n'
            txt+='    \\renewcommand{\\arraystretch}{1.}\n'
            # txt+='    \\setlength\\tabcolsep{2pt}\n'
            n_rows = len(self.samples)
            txt+='    \\begin{tabular}{l||cc|'+'|'.join(['ccc']*(n_rows))+'}\n'
            txt+='      & '

            # Write header of the table
            txt += '\\multicolumn{2}{c|}{'+self.ref_name+'} '
            for smp in self.sample_names:
                txt += '& \\multicolumn{3}{c'+(self.sample_names.index(smp) != len(self.sample_names)-1)*'|'+'}{'+smp+'} '
                # if not self.sample_names.index(smp) == len(self.sample_names)-1:
                #     txt += '&'
                # else:
            txt += '\\\ \\hline\\hline\n'
            txt +='      & '+(not raw)*'Events'+(raw)*'Entries'+' & $\\varepsilon$'
            for smp in self.sample_names:
                txt += ' & '+(not raw)*'Events'+(raw)*'Entries'+' & $\\varepsilon$ & $\\delta$ [\%]'
                # if not self.sample_names.index(smp) == len(self.sample_names)-1:
                #     txt += ' & '
                # else:
            txt += '\\\ \\hline\n'
            # write cutflow
            for cutID, cut in self.ref_sample[SR].items():
                name = cut.Name
                if '$' not in name:
                    name = name.replace('_',' ')
                txt += '      '+name.ljust(40,' ') + '& '
                if cutID == 0:
                    tmp = '{}'+' & - '
                    if raw:
                        txt += tmp.format(scientific_LaTeX(cut.Nentries,sty=event_style))
                    else:
                        txt += tmp.format(scientific_LaTeX(cut.Nevents,sty=event_style))
                else:
                    tmp = '{}'+(MCunc and cut.Nentries>0)*(' $ \pm $ '+event_style)+' & '+eff_style
                    if raw:
                        txt += tmp.format(scientific_LaTeX(cut.Nentries,sty=event_style),cut.raw_rel_eff)
                    else:
                        if not (MCunc and cut.Nentries>0):
                            txt += tmp.format(scientific_LaTeX(cut.Nevents,sty=event_style),cut.rel_eff)
                        else:
                            txt += tmp.format(scientific_LaTeX(cut.Nevents,sty=event_style),cut.MCunc,cut.rel_eff)

                for sample in self.samples:
                    smp = sample[SR]
                    if cutID == 0:
                        tmp = ' & {} & - & - '
                        if raw:
                            txt += tmp.format(scientific_LaTeX(smp[cutID].Nentries,sty=event_style))
                        else:
                            txt += tmp.format(scientific_LaTeX(smp[cutID].Nevents,sty=event_style))
                    elif cutID > 0 and cut.rel_eff == 0:
                        tmp = ' & {}'+(MCunc and smp[cutID].Nentries>0)*(' $ \pm $ '+event_style)+\
                              ' & '+eff_style+' & - '
                        if raw:
                            txt += tmp.format(scientific_LaTeX(smp[cutID].Nentries,sty=event_style),smp[cutID].raw_rel_eff)
                        else:
                            if not (MCunc and smp[cutID].Nentries>0):
                                txt += tmp.format(scientific_LaTeX(smp[cutID].Nevents,sty=event_style),smp[cutID].rel_eff)
                            else:
                                txt += tmp.format(scientific_LaTeX(smp[cutID].Nevents,sty=event_style),smp[cutID].MCunc,smp[cutID].rel_eff)
                    else:
                        tmp = ' & {}'+(MCunc and smp[cutID].Nentries>0)*(' $ \pm $ '+event_style)+\
                              ' & '+eff_style+' & '+ratio_style+' '
                        if raw:
                            rel_eff =abs(1-(smp[cutID].raw_rel_eff/cut.raw_rel_eff))
                            txt  += tmp.format(scientific_LaTeX(smp[cutID].Nentries,sty=event_style),smp[cutID].raw_rel_eff,rel_eff*100.)
                        else:
                            rel_eff =abs(1-(smp[cutID].rel_eff/cut.rel_eff))
                            if not (MCunc and smp[cutID].Nentries>0):
                                txt  += tmp.format(scientific_LaTeX(smp[cutID].Nevents,sty=event_style),smp[cutID].rel_eff,rel_eff*100.)
                            else:
                                txt += tmp.format(scientific_LaTeX(smp[cutID].Nevents,sty=event_style),smp[cutID].MCunc,smp[cutID].rel_eff,rel_eff*100.)
                    # if smp != self.samples[-1][SR]:
                    #     txt += ' & '  
                    # else:
                txt += r'\\'
                txt += '\n'
            
            if finalMCunc:
                tmp = '$ '+event_style+' \\pm '+event_style+' $'
                finalMCunc = [tmp.format(smp.Nevents,smp.MCunc) for smp in [self.ref_sample[SR].get_final_cut()]+\
                                                                           [sample[SR].get_final_cut() for sample in self.samples]]
            else:
                finalMCunc = ''
            entries = [(x.Nentries,r' ($\Delta_{MC}'+r'={:.2f}\%$)'.format(100.*x.MCunc/max(x.Nevents,1e-10)))\
                       for x in [self.ref_sample[SR].get_final_cut()]+[sample[SR].get_final_cut() for sample in self.samples]]
            txt+='    \\end{tabular}\n'
            txt+='    \\caption{'+SR.replace('_',' ')+\
            (any([x[0]<100 for x in entries]))*(' (This region might need more event $\\to$ MC event count = '+\
                                             ', '.join([(x[0]<1e99)*(str(x[0])+x[1])+(x[0]==1e99)*' - ' for x in entries])+') ')+\
                 (self.notes != '')*self.notes+(finalMCunc!='')*('   ($N \\pm \\Delta_{\\rm MC} = $ '+', '.join(finalMCunc)+')')+'}\n'
            txt+='  \\end{center}\n'
            txt+='\\end{table}\n'
            if TeX != None:
                TeX.write(txt)
            else:
                print(txt)
        if TeX != None:
            TeX.write('\n\n\n\n'+r'\end{landscape}'+'\n'+r'\end{document}'+'\n')
            if kwargs.get('make',True):
                self.WriteMake(TeX,make=kwargs.get('make',True))


    def write_signal_comparison_table(self,*args,**kwargs):
        """
        Writes Signal vs Bkg comparison table.

        Parameters
        ----------
        *args : FILE
            Optional, if there is a file input, tables will be written in the 
            file otherwise all will be printed on the screen.
        **kwargs : 
            sys : FLOAT ]0,1]
                Systematic uncertainty, default 20%
            only_alive : BOOLEAN (default True)
                only write the SRs which has more than zero yield for reference
                collection.
            sys_sig : BOOL
                Calculate S/sqrt(B+(B*sys)^2) -> (default False)
            ZA : BOOL
                Calculate Assimov significance -> (default False)
            make : BOOL
                Write the Makefile -> (default, True)

        Returns
        -------
        Signal over Background comparison table.

        """
        sys = kwargs.get('sys',0.2)
        SR_list = self.ref_sample.keys()
        if kwargs.get('only_alive',True): 
            SR_list = [x for x in SR_list if self.ref_sample[x].isAlive()]
        SR_list.sort(key=self._sorter, reverse=True)
        file = None
        if len(args) > 0:
            file = args[0]
            file.write(r'\documentclass[12pt]{article}'+'\n'+\
                       r'\usepackage{pdflscape,slashed}'+'\n'+\
                       r'\begin{document}'+'\n'+\
                       r'\begin{landscape}'+'\n\n\n\n')
            if kwargs.get('ZA',False):
                file.write(r'%%%    Z_A=\sqrt{ 2\left('+'\n')
                file.write(r'%%%    (S+B)\ln\left[\frac{(S+B)(S+\sigma^2_B)}{B^2+(S+B)\sigma^2_B}\right] -'+'\n')
                file.write(r'%%%    \frac{B^2}{\sigma^2_B}\ln\left[1+\frac{\sigma^2_BS}{B(B+\sigma^2_B)}\right]'+'\n')
                file.write(r'%%%    \right)}'+'\n\n\n\n\n\n')
        for SR in SR_list:
            txt = '\n\n%% '+SR+'\n\n'
            txt+='\\begin{table}[h]\n'
            txt+='  \\begin{center}\n'
            txt+='  \\renewcommand{\\arraystretch}{1.}\n'
            # txt+='  \\setlength\\tabcolsep{2pt}\n'
            n_rows = len(self.samples)
            txt+='    \\begin{tabular}{l||cc|'+'|'.join(['cc']*(n_rows))+'}\n'
            txt+='      & '

            # Write header of the table
            txt += '\\multicolumn{2}{c|}{'+self.ref_name+'} &'
            for smp in self.sample_names:
                txt += '\\multicolumn{2}{c'+(self.sample_names.index(smp) != len(self.sample_names)-1)*'|'+'}{'+smp+'} '
                if not self.sample_names.index(smp) == len(self.sample_names)-1:
                    txt += '&'
                else:
                    txt += '\\\ \\hline\\hline\n'
            txt +='      & Events & $\\varepsilon$ &'
            for smp in self.sample_names:
                txt += 'Events & $\\varepsilon$ '
                if not self.sample_names.index(smp) == len(self.sample_names)-1:
                    txt += ' & '
                else:
                    txt += '\\\ \\hline\n'
            # write cutflow
            for cutID, cut in self.ref_sample[SR].items():
                name = cut.Name
                if '$' not in name:
                    name = name.replace('_',' ')
                txt += '      '+name.ljust(40,' ') + '& '
                if cutID == 0:
                    txt += '{:.1f} & - &'.format(cut.Nevents)
                else:
                    txt += '{:.1f} & {:.3f} &'.format(cut.Nevents,cut.rel_eff)
                
                for sample in self.samples:
                    smp = sample[SR]
                    if cutID == 0:
                        txt += '{:.1f} & - '.format(smp[cutID].Nevents)
                    elif cutID > 0 and cut.rel_eff == 0:
                        txt += '{:.1f} & {:.3f} '.format(smp[cutID].Nevents,smp[cutID].rel_eff)
                    else:
                        txt += '{:.1f} & {:.3f} '.format(smp[cutID].Nevents,smp[cutID].rel_eff)
                    if smp != self.samples[-1][SR]:
                        txt += ' & '  
                    else:
                        txt += r'\\'

                if cut == self.ref_sample[SR].get_final_cut():
                    txt += r'\hline\hline'
                    txt += '\n     \\multicolumn{3}{c}{$S/B$} &'
                    for sample in self.samples:
                        smp = sample[SR]
                        fom = FoM(smp[cutID].Nevents,cut.Nevents,sys=sys)
                        txt += '\\multicolumn{2}{c}{'+'{:.3f}\\%'.format(100.*fom.S_B)+'}'
                        if smp != self.samples[-1][SR]:
                            txt += ' & ' 
                        else:
                            txt += r'\\'
 
                    txt += '\n     \\multicolumn{3}{c}{$S/S+B$} &'
                    for sample in self.samples:
                        smp = sample[SR]
                        fom = FoM(smp[cutID].Nevents,cut.Nevents,sys=sys)
                        txt += '\\multicolumn{2}{c}{'+'{:.3f}\\%'.format(100.*fom.S_SB)+'}'
                        if smp != self.samples[-1][SR]:
                            txt += ' & ' 
                        else:
                            txt += r'\\'

                               
                    txt += '\n     \\multicolumn{3}{c}{$S/\sqrt{B}$}  &'
                    for sample in self.samples:
                        smp = sample[SR]
                        fom = FoM(smp[cutID].Nevents,cut.Nevents,sys=sys)
                        txt += '\\multicolumn{2}{c}{'+\
                                '{:.3f}'.format(fom.sig)+'}'
                        if smp != self.samples[-1][SR]:
                            txt += ' & ' 
                        else:
                            txt += r'\\'

                    if kwargs.get('sig_sys',False):
                        txt += '\n     \\multicolumn{3}{c}{$S/\sqrt{B+(B\Delta_{sys})^2}$}  &'
                        for sample in self.samples:
                            smp = sample[SR]
                            fom = FoM(smp[cutID].Nevents,cut.Nevents,sys=sys)
                            txt += '\\multicolumn{2}{c}{'+\
                                    '{:.3f}'.format(fom.sig_sys)+'}'
                            if smp != self.samples[-1][SR]:
                                txt += ' & ' 
                            else:
                                txt += r'\\'

                    if kwargs.get('ZA',False):
                        txt += '\n     \\multicolumn{3}{c}{$Z_A$} &'
                        for sample in self.samples:
                            smp = sample[SR]
                            fom = FoM(smp[cutID].Nevents,cut.Nevents,sys=sys)
                            txt += '\\multicolumn{2}{c}{'+\
                                    '${:.3f} \\pm {:.3f} $'.format(fom.ZA,fom.ZA_err)+'}'
                            if smp != self.samples[-1][SR]:
                                txt += ' & ' 
                            else:
                                txt += r'\\'
                    
                txt += '\n'

            txt+='    \\end{tabular}\n'
            txt+='    \\caption{'+SR.replace('_',' ')+\
            (cut.Nentries<100)*'(This SR needs more event:: MC event count = {:.0f})'.format(cut.Nentries)+'}\n' 
            txt+='  \\end{center}\n'
            txt+='\\end{table}\n'
            if file != None:
                file.write(txt)
            else:
                print(txt)
        if file != None:
            file.write('\n\n\n\n'+r'\end{landscape}'+'\n'+r'\end{document}'+'\n')
            if kwargs.get('make',True):
                self.WriteMake(file,make=kwargs.get('make',True))


    def WriteMake(self,file,make=True):
        """
        Writes make file for given tex file.

        Parameters
        ----------
        file : FILE
            TeX file to write the Makefile for.
        make : BOOL, optional
            Compile or not. The default is True.

        Raises
        ------
        ValueError
            Can not find the file.

        Returns
        -------
        None.

        """
        if not file.name.endswith('.tex'):
            raise ValueError('Input does not have .tex extention.')
        if os.path.isfile(file.name):
            make = open('Makefile','w')
            make.write('all:\n')
            make.write('\tpdflatex '+file.name[:-4]+'\n'+\
                       '\tpdflatex '+file.name[:-4]+'\n'+\
                       '\trm -f *.aux *.log *.out *.toc *.blg *.dvi *.t1 *.1 *.mp *spl\n'+\
                       'clean:\n'+\
                       '\trm -f *.aux *.log *.out *.toc *.blg *.dvi *.t1 *.1 *.mp *spl *.lol *Notes.bib\n')
            if make:
                try:
                    file.close()
                    os.system('make')
                except:
                    print 'Compilation failed.'
        else:
            raise ValueError('Can not find '+file.name)


def scientific_LaTeX(val,sty='{:.1f}'):
    if val >= 1e5:
        tmp = '{:.1e}'.format(val)
        tmp = [float(x) for x in tmp.split('e+')]
        tmp = r'${:.1f} \times 10^'.format(tmp[0]) + '{' + '{:.0f}'.format(tmp[1])+'}$'
    elif val < 1e-3 and val > 0.:
        tmp = '{:.1e}'.format(val)
        tmp = [float(x) for x in tmp.split('e-')]
        tmp = r'${:.1f} \times 10^'.format(tmp[0]) + '{-' + '{:.0f}'.format(tmp[1])+'}$'
    else:
        tmp = sty.format(val)
    return tmp

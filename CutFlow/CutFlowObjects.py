#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Fri Jan 31 10:37:07 2020

@author  : jackaraz
@contact : Jack Y. Araz <jackaraz@gmail.com>
"""
from math import sqrt

class Cut(object):
    def __init__(self,Name='NaN',Nentries=-1,sumw=None,sumw2=-1, precut=None,cut_0=None, xsec=1.,Nevents=None):
        """

        Parameters
        ----------
        Name : STR, optional
            Name of the cut. The default is -1.
        Nentries : FLOAT, optional
            Number of entries. The default is -1.
        sumw : FLOAT, optional
            Sum of weights. The default is -1.
        sumw2 : FLOAT, optional
            Square of sum of weights. The default is -1.
        precut : Cut, optional
            Previous cut. The default is None.
        cut_0 : Cut, optional
            Initial cut. The default is None.
        xsec : FLOAT, optional
            Cross section. The default is 1..

        Returns
        -------
        None.

        """
        self.Name        = str(Name)
        self.Nentries    = Nentries
        self.sumw        = sumw
        self.sumw2       = sumw2

        if Nevents == None:
            if cut_0 == None:
                self.eff         = 1.
                self.raw_eff     = 1.
            else:
                if cut_0.sumw > 0 and cut_0.Nentries > 0:
                    self.eff         = round(sumw/cut_0.sumw,8)
                    self.raw_eff     = round(float(Nentries)/float(cut_0.Nentries),8)
                else:
                    self.eff         = 1.
                    self.raw_eff     = 1.
            if precut == None:
                self.rel_eff     = 1.
                self.raw_rel_eff = 1.
            else:
                if precut.sumw > 0 and precut.Nentries>0:
                    self.rel_eff     = round(sumw/precut.sumw,8)
                    self.raw_rel_eff = round(float(Nentries)/float(precut.Nentries),8)
                else:
                    self.rel_eff     = 1.
                    self.raw_rel_eff = 1.
            self.nevt        = round(self.eff*xsec,8)
            self.Nevents     = round(self.eff*xsec,8)
        else:
            self.nevt        = round(Nevents, 8)
            self.Nevents     = self.nevt
            if cut_0 == None:
                self.eff         = 1.
                self.raw_eff     = 1.
            else:
                if cut_0.sumw > 0:
                    self.eff         = round(sumw/cut_0.sumw,8)
                elif cut_0.sumw == None:
                    self.eff         = round(Nevents/cut_0.Nevents,8)
                else:
                    self.eff         = 1.
                if cut_0.Nentries > 0:
                    self.raw_eff     = round(float(Nentries)/float(cut_0.Nentries),8)
                else:
                    self.raw_eff     = 1.
            if precut == None:
                self.rel_eff     = 1.
                self.raw_rel_eff = 1.
            else:
                if precut.sumw > 0:
                    self.rel_eff     = round(sumw/precut.sumw,8)
                elif precut.sumw == None:
                    self.rel_eff      = round(Nevents/precut.Nevents,8)
                else:
                    self.rel_eff     = 1.
                if precut.Nentries>0:
                    self.raw_rel_eff = round(float(Nentries)/float(precut.Nentries),8)
                else:
                    self.raw_rel_eff = 1.
            if self.Nevents < 0.:
                self.nevt        = 0.
                self.Nevents     = 0.
        if self.eff <= 1:
            self.MCunc = self.Nevents*sqrt(self.eff*(1.-self.eff)/float(self.Nentries)) if self.Nentries>0 else 0.
        else:
            self.MCunc = -1


    def __mul__(self,lumi):
        # in ifb
        self.Nevents *= 1000.*lumi
        self.nevt     = self.Nevents
        if self.Nentries > 0:
            self.MCunc = self.Nevents*sqrt((1.-self.eff)/float(self.Nentries))
        else:
            self.MCunc = 0.
        return self
    
    def set_xsec(self,xsec):
        self.nevt    = round(self.eff*xsec, 8)
        self.Nevents = self.nevt
        if self.Nentries > 0:
            self.MCunc = self.Nevents*sqrt((1.-self.eff)/float(self.Nentries))
        else:
            self.MCunc = 0.
        return self
    
    def __str__(self):
        if self.eff < 1:
            return  '   '+self.Name+'\n'+\
                    '      Nentries: {:.0f}\n      Nevents : {:.3f} ± {:.3f}(ΔMC)\n'.format(self.Nentries,self.Nevents,self.MCunc)+\
                    '      Cut Eff : {:.5f}\n      Rel Eff : {:.5f}'.format(self.eff,self.rel_eff)
        else:
            return  '   '+self.Name+'\n'+\
                    '      Nentries: {:.0f}\n      Nevents : {:.3f}\n'.format(self.Nentries,self.Nevents)+\
                    '      Cut Eff : {:.5f}\n      Rel Eff : {:.5f}'.format(self.eff,self.rel_eff)



class SignalRegion(object):
    def __init__(self,name):
        self.name = name
        self.cutlist = []

    def __getitem__(self,cut_num):
        return self.cutlist[cut_num]

    @classmethod
    def __type__(self):
        return __name__

    def __len__(self):
        return len(self.cutlist)

    def items(self):
        return [(i,self.cutlist[i]) for i in range(len(self.cutlist))]

    def add_cut(self,cut):
        self.cutlist.append(cut)
    
    def get_names(self):
        return [x.Name for x in self.cutlist]

    def get_name(self,n):
        return self.cutlist[n].Name

    def get_cut(self,n):
        return self.cutlist[n]

    def get_final_cut(self):
        return self.cutlist[len(self)-1]

    def isAlive(self):
        return self.get_final_cut().Nentries > 0

    def __mul__(self,lumi):
        self.cutlist = [cut * lumi for cut in self.cutlist]
        return self

    def set_xsec(self,xsec):
        self.cutlist = [cut.set_xsec(xsec) for cut in self.cutlist]
        return self

    def regiondata(self):
        return {self.name : {'Nf' : self.get_final_cut().sumw,
                             'N0' : self.get_cut(0).sumw}}

    def __str__(self):
        return '\n'.join(['   '+str(i)+'. '+str(self.cutlist[i]) for i in range(len(self.cutlist))])
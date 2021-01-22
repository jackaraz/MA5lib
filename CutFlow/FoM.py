#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Fri Jan 31 11:44:54 2020

@author  : jackaraz
@contact : Jack Y. Araz <jackaraz@gmail.com>
"""

from numpy import sqrt, log, power

class FoM:
    def __init__(self,nsignal,nbkg,sys=0.):
        self.nsignal = nsignal
        self.nbkg    = nbkg
        self.sys     = sys
        if nbkg == 0.:
            self.ZA      = -1
            self.ZA_err  = -1
            self.sig_sys = -1
            self.sig     = -1
            self.S_B     = -1
            self.S_SB    = -1
            self.S_sqSB  = -1
        else:
            if sys > 0.:
                self.ZA      = self.asimovZ()
                self.ZA_err  = self.asimovError()
                self.sig_sys = self.significance()
            self.sig  = nsignal/sqrt(nbkg)
            self.S_B  = nsignal/nbkg
            self.S_SB = nsignal/(nbkg+nsignal)
            self.S_sqSB = nsignal/sqrt(nbkg+nsignal)
    
    def asimovZ(self):
        """
            arXiv:1007.1727
        """
        try:
            varb = self.nbkg*self.sys*self.nbkg*self.sys
            tot  = self.nsignal + self.nbkg
            asimovsig = sqrt(2*(tot*log((tot*(varb+self.nbkg))/((self.nbkg*self.nbkg)+\
                        tot*varb))-(self.nbkg*self.nbkg/varb)*log(1+(varb*self.nsignal)/(self.nbkg*(self.nbkg+varb)))))
        except:
            return 0.0
        return asimovsig
    
    def asimovError(self):
        sig = self.sys
        es  = sqrt(self.nsignal)
        eb  = sqrt(self.nbkg)
        s   = self.nsignal
        b   = self.nbkg
        try:
            err = power(-(eb*eb)/( 1.0/(sig*sig)*log( b/( b+(b*b)*(sig*sig))*(sig*sig)*s+1.0)-\
                           ( b+s)*log(( b+s)*( b+(b*b)*(sig*sig))/( (b*b)+( b+s)*(b*b)*(sig*sig))))*\
                    power( 1.0/( b/( b+(b*b)*(sig*sig))*(sig*sig)*s+1.0)/(sig*sig)*( 1.0/( b+(b*b)*(sig*sig))*\
                                (sig*sig)*s-b/power( b+(b*b)*(sig*sig),2.0)*(sig*sig)*( 2.0*b*(sig*sig)+1.0)*s)-\
                    ( ( b+s)*( 2.0*b*(sig*sig)+1.0)/( (b*b)+( b+s)*(b*b)*(sig*sig))+( b+(b*b)*(sig*sig))/( (b*b)+\
                       ( b+s)*(b*b)*(sig*sig))-( b+s)*( 2.0*( b+s)*b*(sig*sig)+2.0*b+(b*b)*(sig*sig))*( b+\
                       (b*b)*(sig*sig))/power( (b*b)+( b+s)*(b*b)*(sig*sig),2.0))/( b+(b*b)*(sig*sig))*( \
                    (b*b)+( b+s)*(b*b)*(sig*sig))-log(( b+s)*( b+(b*b)*(sig*sig))/( (b*b)+( b+s)*(b*b)*(sig*sig))),\
                    2.0)/2.0-1.0/( 1.0/(sig*sig)*log( b/( b+(b*b)*(sig*sig))*(sig*sig)*s+1.0)-\
                                 ( b+s)*log(( b+s)*( b+(b*b)*(sig*sig))/( (b*b)+( b+s)*(b*b)*(sig*sig))))*\
                    power( log(( b+s)*( b+(b*b)*(sig*sig))/( (b*b)+( b+s)*(b*b)*(sig*sig)))+1.0/( 
                            b+(b*b)*(sig*sig))*( ( b+(b*b)*(sig*sig))/( (b*b)+( b+s)*(b*b)*(sig*sig))-\
                              ( b+s)*(b*b)*( b+(b*b)*(sig*sig))*(sig*sig)/power( (b*b)+( b+s)*(b*b)*(sig*sig),
                               2.0))*( (b*b)+( b+s)*(b*b)*(sig*sig))-1.0/( b/( b+(b*b)*(sig*sig))*(sig*sig)*s+\
                               1.0)*b/( b+(b*b)*(sig*sig)),2.0)*(es*es)/2.0,(1.0/2.0))
        except:
            return 0.
        return err


    def significance(self):
        varb    = self.nbkg*self.sys*self.nbkg*self.sys
        stopsig = self.nsignal / sqrt(self.nbkg+varb)
        return stopsig

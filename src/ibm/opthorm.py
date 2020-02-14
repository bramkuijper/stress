#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Feb 14 22:49:08 2020

@author: bram
"""


class OptHorm:
    
    """
    Calculates optimal hormone levels in an environment
    without any autocorrelationa
    """
    
    def __init__(self, ad, ap, mu_bg, sP2NP, sNP2P, p_att):
        """
        Initialize the optimal hormone level class with
        parameters for the calculation

        Parameters
        ----------
        ad : float
            power of hormone damage cost function.
        ap : float
            power of hormone-dependent predation function. 
        mu_bg : float
            background mortality rate.
        sP2NP : float
            switch probability from predator (i.e., stress) to non predator. 
        sNP2P : float
            switch probability from non predator to predator.
        p_att : float
            attack probability of predator.
        
        Returns
        -------
        None.

        """
        
        self.ad = ad
        self.ap = ap
        self.mu_bg = mu_bg
        self.pr_P = sNP2P / (sP2NP + sNP2P)
        self.p_att = p_att


    def mort(self, h):
        """
        Calculates mortality probability, which depends on hormone level h

        Parameters
        ----------
        h : float
            hormone level.

        Returns
        -------
        a float reflecting the average mortality probability

        """
        return(self.mu_bg +\
                (1.0 - self.mu_bg) * self.pr_P * self.p_att *\
                    (1 - h**self.ap))

    def repr(self, h):
        """
        Reproductive success, dependent on hormone level h

        Parameters
        ----------
        h : float
            hormone level.

        Returns
        -------
        fecundity (float).

        """
        return(1.0 - h**self.ad)

    def lrs_neg(self, h):
        """
        Lifetime reproductive success, dependent on hormone level h
        multiplied by -1. We take this negative as python only has 
        minimization, rather than maximization

        Parameters
        ----------
        h : float
            hormone level.

        Returns
        -------
        Life time reproductive success (float) multiplied by -1

        """
        
        # reproductive success is (1 - mort) * repr
        # lifespan = 1/mort(h)
        return(-1*(1.0 - self.mort(h)) * self.repr(h) / self.mort(h))

    def print_lrs_vs_hormone(self):
        """
        Prints lifetime reproductive success versus hormone level
        to see whether we are indeed minimizing correctly (i.e., at a
        valley of negative lifetime RS
        
        Returns
        -------
        None.

        """

        hval = list(np.linspace(0,1,50))

        the_str = "h;lrs;\n"

        for hval_i in hval:
            the_lrs = -self.lrs_neg(hval_i)
            the_str += f"{hval_i};{the_lrs}\n"

        print(the_str)

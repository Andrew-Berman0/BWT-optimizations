# bidirectional_fm.py

from bwtConst import bwtConst
from collections import defaultdict

class BiFMIndex:
    def __init__(self, bwtConst):
        """
        Bidirectional FM-index using the BwtConst class.
        Stores BWT and auxiliary data for forward/backward extension.
        """
        # Use bwtConst to construct BWT, suffix array, C array, and OCC table
        self.bwt_struct = bwtConst
        self.text = self.bwt_struct.text
        self.bwt = self.bwt_struct.bwt
        self.C = self.bwt_struct.C
        self.OCC = self.bwt_struct.OCC
        self.alphabet = sorted(set(self.bwt))

    # -----------------------------
    # LF-mapping for backward extension
    def lf_backward(self, l, r, c):
        """
        Extend pattern backward with character c.
        l, r = current interval in BWT
        """
        l_new = self.C[c] + (self.OCC[c][l-1] if l > 0 else 0)
        r_new = self.C[c] + self.OCC[c][r] - 1
        return l_new, r_new

    # LF-mapping for forward extension
    def lf_forward(self, l, r, c):
        """
        Extend pattern forward with character c.
        Uses the concept of inverse BWT interval.
        """
        # For simplicity, using full scan (can be optimized with wavelet trees)
        positions = [i for i in range(len(self.bwt)) if self.bwt[i] == c]
        # Find positions that map into current interval [l,r]
        l_new = next((i for i in positions if l <= i <= r), None)
        r_new = next((i for i in reversed(positions) if l <= i <= r), None)
        if l_new is None or r_new is None:
            return None, None
        return l_new, r_new

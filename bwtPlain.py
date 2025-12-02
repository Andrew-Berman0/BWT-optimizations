# fm_index_plain.py

from bwtConst import bwtConst
from collections import defaultdict

class FMIndexPlain:
    def __init__(self, bwtConst):
        """
        Plain FM-index using only C array and full OCC table.
        """
        # Build BWT and auxiliary structures using bwtConst
        self.bwt_struct = bwtConst
        self.bwt = self.bwt_struct.bwt
        self.C = self.bwt_struct.C
        self.OCC = self.bwt_struct.OCC
        self.alphabet = sorted(set(self.bwt))

    def backward_search(self, pattern):
        """
        Backward search for a pattern using C array and OCC table.
        Returns True if pattern exists in the text.
        """
        l, r = 0, len(self.bwt) - 1
        for c in reversed(pattern):
            if c not in self.C:
                return False
            l = self.C[c] + (self.OCC[c][l-1] if l > 0 else 0)
            r = self.C[c] + self.OCC[c][r] - 1
            if l > r:
                return False
        return True

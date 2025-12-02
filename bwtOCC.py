# occ_checkpoint.py

from bwtConst import bwtConst
from collections import defaultdict

class OCCCompressed:
    def __init__(self, bwtConst, checkpoint_interval=4):
        """
        OCC table with checkpoints to save space.
        Uses BWT from bwtConst.
        :param text: original text to build BWT from
        :param checkpoint_interval: distance between checkpoints
        """
        # Use bwtConst to generate BWT
        self.bwt_struct = bwtConst
        self.bwt = self.bwt_struct.bwt
        self.checkpoint_interval = checkpoint_interval
        self.chars = sorted(set(self.bwt))
        self.checkpoints = defaultdict(list)  # checkpoints[c] = list of counts at multiples of interval

        # Initialize counts and save checkpoints
        counts = {c: 0 for c in self.chars}
        for i, ch in enumerate(self.bwt):
            counts[ch] += 1
            if (i + 1) % checkpoint_interval == 0:
                for c in self.chars:
                    self.checkpoints[c].append(counts[c])

    def occ(self, c, idx):
        """
        Return number of occurrences of character c in bwt[0..idx]
        """
        if c not in self.chars:
            return 0
        
        # Find nearest checkpoint
        cp_idx = idx // self.checkpoint_interval
        count = self.checkpoints[c][cp_idx - 1] if cp_idx > 0 else 0

        # Count occurrences between checkpoint and idx
        start = cp_idx * self.checkpoint_interval
        end = idx + 1
        count += sum(1 for i in range(start, end) if self.bwt[i] == c)
        return count

from array import array

class WaveletTreeNode:
    def __init__(self, alphabet, bwt):
        self.left = None
        self.right = None
        self.alphabet = alphabet
        
        # Optimization 1: Use 'array' for 1-byte integers (vs 28-byte Python ints)
        # 'b' = signed char (1 byte)
        self.bitmap = array('b') 
        
        # Optimization 2: Sparse Checkpoints (store rank only every K bits)
        self.checkpoint_interval = 32
        self.checkpoints = array('I') # 'I' = unsigned int (4 bytes)

        if len(alphabet) == 0 or len(bwt) == 0:
            return

        if len(alphabet) == 1:
            return

        # Split alphabet
        mid_idx = len(alphabet) // 2
        self.left_alphabet = alphabet[:mid_idx]
        self.right_alphabet = alphabet[mid_idx:]
        left_set = set(self.left_alphabet)

        # Build bitmap and Checkpoints in one pass
        running_ones = 0
        left_bwt = []
        right_bwt = []

        for i, c in enumerate(bwt):
            # 1. Record Checkpoint every 32 bits
            if i % self.checkpoint_interval == 0:
                self.checkpoints.append(running_ones)

            # 2. Determine bit
            if c in left_set:
                bit = 0
                left_bwt.append(c)
            else:
                bit = 1
                right_bwt.append(c)
                running_ones += 1
            
            self.bitmap.append(bit)
        
        # Capture final count for range boundary checks
        if len(bwt) % self.checkpoint_interval == 0:
            self.checkpoints.append(running_ones)

        # Recursively create children
        if left_bwt:
            self.left = WaveletTreeNode(self.left_alphabet, left_bwt)
        if right_bwt:
            self.right = WaveletTreeNode(self.right_alphabet, right_bwt)

    def rank(self, c, idx):
        """
        Rank using Sparse Checkpoints.
        Complexity: O(1) [Lookup + max 32 iterations]
        """
        if idx < 0:
            return 0
            
        if len(self.alphabet) == 1:
            return idx + 1 if self.alphabet[0] == c else 0

        # 1. Jump to nearest checkpoint
        chunk_idx = (idx + 1) // self.checkpoint_interval
        ones_count = self.checkpoints[chunk_idx]

        # 2. Scan the remaining bits (max 31 iterations)
        start_scan = chunk_idx * self.checkpoint_interval
        
        # This slice is tiny, so the loop is very fast
        # Note: We scan up to idx (inclusive)
        limit = idx + 1
        for k in range(start_scan, limit):
             if self.bitmap[k] == 1:
                 ones_count += 1

        # 3. Calculate Zeros
        zeros_count = (idx + 1) - ones_count

        if c in self.left_alphabet:
            return self.left.rank(c, zeros_count - 1) if self.left else zeros_count
        else:
            return self.right.rank(c, ones_count - 1) if self.right else ones_count

class WaveletTree:
    def __init__(self, bwt_const_obj):
        self.bwt_struct = bwt_const_obj
        self.bwt = self.bwt_struct.bwt
        self.alphabet = sorted(list(self.bwt_struct.C.keys()))
        self.root = WaveletTreeNode(self.alphabet, list(self.bwt))

    def rank(self, c, idx):
        return self.root.rank(c, idx)
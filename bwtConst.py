# bwtConst.py

class bwtConst:
    def __init__(self, text):
        # Ensure text is terminated (standard for BWT)
        if not text.endswith('$'):
            text += '$'
        
        self.text = text
        
        # 1. Build Suffix Array (Optimized)
        self.sa = self._build_suffix_array_prefix_doubling(text)
        
        # 2. Build BWT from SA
        # BWT[i] = Text[SA[i] - 1]
        self.bwt = ''.join([text[i - 1] if i > 0 else '$' for i in self.sa])
        
        # 3. Build C table and OCC table (Optimized)
        self.C = self._build_C_table()
        self.OCC = self._build_OCC_table()

    def _build_suffix_array_prefix_doubling(self, text):
        """
        Constructs Suffix Array using Prefix Doubling (Manber-Myers).
        Complexity: O(N log^2 N)
        Avoids string slicing entirely.
        """
        n = len(text)
        
        # Initial ranking based on first character
        # rank[i] is the rank of the suffix starting at i
        rank = [ord(c) for c in text]
        
        # Tuple to store (rank of first half, rank of second half, original index)
        # We need n tuples
        # k is the length of the prefix we are currently sorting by
        k = 1
        
        while k < n:
            # Create a list of tuples: (rank[i], rank[i+k], i)
            # If i+k is out of bounds, use -1 (smaller than any char)
            tuples = []
            for i in range(n):
                first = rank[i]
                second = rank[i + k] if i + k < n else -1
                tuples.append((first, second, i))
            
            # Sort based on the tuple values (first, then second)
            # Python's sort is stable and efficient for tuples of ints
            tuples.sort()
            
            # Recompute ranks based on the sorted order
            new_rank = [0] * n
            sa = [0] * n
            
            # Assign first rank
            sa[0] = tuples[0][2]
            new_rank[tuples[0][2]] = 0
            
            for i in range(1, n):
                sa[i] = tuples[i][2]
                prev_tuple = tuples[i-1]
                curr_tuple = tuples[i]
                
                # If tuples are same, they get same rank; otherwise increment rank
                if prev_tuple[0] == curr_tuple[0] and prev_tuple[1] == curr_tuple[1]:
                    new_rank[sa[i]] = new_rank[sa[i-1]]
                else:
                    new_rank[sa[i]] = new_rank[sa[i-1]] + 1
            
            rank = new_rank
            k *= 2
            
            # Optimization: If all ranks are distinct, we are done
            if rank[sa[n-1]] == n - 1:
                break
                
        return sa

    def _build_C_table(self):
        """Builds the C table: count of characters lexicographically smaller than c"""
        counts = {}
        for c in self.bwt:
            counts[c] = counts.get(c, 0) + 1
        
        sorted_chars = sorted(counts.keys())
        C = {}
        total = 0
        for c in sorted_chars:
            C[c] = total
            total += counts[c]
        return C

    def _build_OCC_table(self):
        """
        Builds the full 2D OCC table.
        Note: For very large genomes, this consumes O(Sigma * N) RAM.
        If memory is tight, this logic moves to the compressed classes.
        """
        # Get unique characters from C table to ensure consistent order
        chars = sorted(self.C.keys())
        
        # Initialize dictionary of lists
        # occ[char] = [count_at_0, count_at_1, ...]
        occ = {c: [0] * (len(self.bwt) + 1) for c in chars}
        
        # Running totals
        current_counts = {c: 0 for c in chars}
        
        # Single pass O(N)
        for i, char in enumerate(self.bwt):
            # Update running count for the current char
            if char in current_counts:
                current_counts[char] += 1
            
            # Copy current state to the table
            # (Note: This is still memory heavy, but faster than multiple passes)
            for c in chars:
                occ[c][i] = current_counts[c]
                
        # Fill the last index (for range queries that go up to len(bwt))
        for c in chars:
            occ[c][len(self.bwt)] = current_counts[c]
            
        return occ
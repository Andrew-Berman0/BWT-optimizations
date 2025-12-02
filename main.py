from bwtConst import bwtConst
from bwtWaveletTree import WaveletTree
from bwtOCC import OCCCompressed
from bwtBD import BiFMIndex
from bwtPlain import FMIndexPlain

import time as t
import random as r
import tracemalloc  # <--- NEW LIBRARY

# -----------------------------
# Helper for Measuring Memory & Time
# -----------------------------
def measure_construction(class_constructor, *args):
    """
    Runs a constructor, measures execution time, 
    and captures the net memory increase (Object Size).
    """
    # 1. Start Memory Tracking
    tracemalloc.start()
    
    # 2. Snapshot time
    start_time = t.perf_counter()
    
    # 3. Build the object
    obj = class_constructor(*args)
    
    # 4. Snapshot time
    end_time = t.perf_counter()
    
    # 5. Get Memory Usage (current, peak)
    # We care about 'current' (how much memory the object is holding right now)
    current_mem, peak_mem = tracemalloc.get_traced_memory()
    
    # 6. Stop Tracking
    tracemalloc.stop()
    
    time_taken = end_time - start_time
    return obj, time_taken, current_mem

# -----------------------------
# Search Functions (Unchanged)
# -----------------------------
def random_pattern(length):
    return ''.join(r.choice('ACGT') for _ in range(length))

def backward_search_bwt(pattern, bwt_obj):
    l, r = 0, len(bwt_obj.bwt) - 1
    for c in reversed(pattern):
        if c not in bwt_obj.C: return False
        # Note: This relies on the full OCC table being in bwt_obj
        l = bwt_obj.C[c] + (bwt_obj.OCC[c][l-1] if l > 0 else 0) 
        r = bwt_obj.C[c] + bwt_obj.OCC[c][r] - 1
        if l > r: return False
    return True

def backward_search_wavelet(pattern, wt_obj):
    l, r = 0, len(wt_obj.bwt) - 1
    C = wt_obj.bwt_struct.C
    for c in reversed(pattern):
        if c not in C: return False
        l = C[c] + wt_obj.rank(c, l - 1)
        r = C[c] + wt_obj.rank(c, r) - 1
        if l > r: return False
    return True

def backward_search_occ(pattern, occ_obj):
    l, r = 0, len(occ_obj.bwt) - 1
    for c in reversed(pattern):
        l = occ_obj.bwt_struct.C[c] + (occ_obj.occ(c, l-1) if l > 0 else 0)
        r = occ_obj.bwt_struct.C[c] + occ_obj.occ(c, r) - 1
        if l > r: return False
    return True

def backward_search_bifm(pattern, bifm_obj):
    l, r = 0, len(bifm_obj.bwt) - 1
    for c in reversed(pattern):
        l, r = bifm_obj.lf_backward(l, r, c)
        if l > r: return False
    return True

def backward_search_plain(pattern, plain_obj):
    # This relies on the Plain object having the backward_search method
    return plain_obj.backward_search(pattern)

# -----------------------------
# Main Execution
# -----------------------------
if __name__ == "__main__":
    genome_sizes = [10000, 100000, 1000000] 
    pattern_length = 10
    num_patterns = 100

    # Table Header
    print(f"{'Structure':<15} | {'Size':<8} | {'Build Time':<10} | {'Memory (KB)':<12} | {'Search Time':<10}")
    print("-" * 75)

    for size in genome_sizes:
        # Generate random genome string
        s = ''.join(r.choice('ACGT') for _ in range(size))

        # --- CORRECTLY MEASURE THE PLAIN FM INDEX COST (bwtConst) ---
        # 1. Measure the BWT Construction Time and Memory
        tracemalloc.start()
        start_build_plain = t.perf_counter()

        # This allocation creates the huge OCC table (the Plain FM Index data)
        bwt_obj = bwtConst(s) 

        end_build_plain = t.perf_counter()
        plain_mem_usage, _ = tracemalloc.get_traced_memory() # Current memory usage is the key
        tracemalloc.stop()

        plain_build_time = end_build_plain - start_build_plain
        # -------------------------------------------------------------

        # List of strategies to test
        # Note: We now treat bwt_obj as the base Plain FM Index object
        strategies = [
            # 1. PLAIN INDEX (Base object is the index)
            ("FMIndexPlain", lambda: bwt_obj, (None,), backward_search_bwt, plain_build_time, plain_mem_usage),
            # 2. OPTIMIZED INDICES (Must be built on top of bwt_obj)
            ("OCCCompressed", OCCCompressed, (bwt_obj, 16), backward_search_occ, 0, 0),
            ("WaveletTree", WaveletTree, (bwt_obj,), backward_search_wavelet, 0, 0),
            ("BiFMIndex", BiFMIndex, (bwt_obj,), backward_search_bifm, 0, 0),
        ]

        for name, cls, args, search_func, pre_time, pre_mem in strategies:
            
            # A. Measure Construction (Time + Memory)
            if name == "FMIndexPlain":
                # Use the pre-measured values for the base object
                obj = bwt_obj
                build_time = pre_time
                mem_usage = pre_mem
            else:
                # Measure the new memory and time added by the compressed index
                obj, build_time, mem_usage = measure_construction(cls, *args)
                # The memory reported here is just the COMPRESSED structure's new size
                
            # B. Measure Search Speed
            matches_found = 0
            patterns = [random_pattern(pattern_length) for _ in range(num_patterns)]
            
            start_search = t.perf_counter()
            for pat in patterns:
                # NOTE: For FMIndexPlain, we must use backward_search_bwt 
                # because the base bwt_obj does not have the 'backward_search' method
                if name == "FMIndexPlain" and search_func == backward_search_plain:
                    # Reroute to the function that uses the bwtConst's OCC table
                    if backward_search_bwt(pat, obj):
                        matches_found += 1
                else:
                    if search_func(pat, obj):
                        matches_found += 1
                        
            end_search = t.perf_counter()
            search_time = end_search - start_search

            # Convert Bytes to KB
            mem_kb = mem_usage / 1024

            print(f"{name:<15} | {size:<8} | {build_time:.5f}s   | {mem_kb:<12.2f} | {search_time:.5f}s")
        
        print("-" * 75)
from bwtConst import bwtConst
from bwtWaveletTree import WaveletTree
from bwtOCC import OCCCompressed
from bwtBD import BiFMIndex
from bwtPlain import FMIndexPlain

import time as t
import random as r
import tracemalloc
import matplotlib.pyplot as plt
import numpy as np
import collections

# -----------------------------
# Helper for Measuring Memory and Time (Unchanged)
# -----------------------------
def measure_construction(class_constructor, *args):
    tracemalloc.start()
    start_time = t.perf_counter()
    obj = class_constructor(*args)
    end_time = t.perf_counter()
    current_mem, peak_mem = tracemalloc.get_traced_memory()
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
    return plain_obj.backward_search(pattern)

# -----------------------------
# Plotting Function (FIXED)
# -----------------------------
def plot_results(results_data, filename="bwt_optimization_benchmarks.png", num_patterns=100):
    """Generates and saves a three-panel plot of the benchmark results."""
    
    # Organize data by structure name
    data_by_structure = collections.defaultdict(lambda: collections.defaultdict(list))
    for entry in results_data:
        name = entry['Structure']
        size = entry['Size']
        data_by_structure[name]['Size'].append(size)
        data_by_structure[name]['BuildTime'].append(entry['BuildTime'])
        data_by_structure[name]['MemoryKB'].append(entry['MemoryKB'])
        data_by_structure[name]['SearchTime'].append(entry['SearchTime'])

    # Setup the figure
    fig, axes = plt.subplots(3, 1, figsize=(10, 15))
    fig.suptitle('BWT Index Performance & Resource Usage', fontsize=16)

    # --- Plot 1: Build Time ---
    ax = axes[0]
    for name, data in data_by_structure.items():
        ax.plot(data['Size'], data['BuildTime'], marker='o', label=name)
    ax.set_title('Construction Time vs. Genome Size')
    ax.set_xlabel('Genome Size (N)', fontsize=10)
    ax.set_ylabel('Build Time (seconds)', fontsize=10)
    ax.legend(fontsize=8)
    ax.grid(True, linestyle='--')
    ax.ticklabel_format(axis='x', style='plain')

    # --- Plot 2: Memory Usage ---
    ax = axes[1]
    for name, data in data_by_structure.items():
        # Mask zero memory values for log scale
        memory_data = np.array(data['MemoryKB'])
        memory_data[memory_data <= 0] = 0.01 
        ax.plot(data['Size'], memory_data, marker='o', label=name)
        
    ax.set_title('Memory Usage vs. Genome Size (Log Scale)')
    ax.set_xlabel('Genome Size (N)', fontsize=10)
    ax.set_ylabel('Peak Memory (Kilobytes)', fontsize=10)
    ax.set_yscale('log') # Use log scale to clearly show difference between Plain and Compressed
    ax.legend(fontsize=8)
    ax.grid(True, linestyle='--')
    ax.ticklabel_format(axis='x', style='plain')

    # --- Plot 3: Search Time (FIXED Y-AXIS LABEL) ---
    ax = axes[2]
    for name, data in data_by_structure.items():
        ax.plot(data['Size'], data['SearchTime'], marker='o', label=name)
    ax.set_title('Average Search Time vs. Genome Size')
    ax.set_xlabel('Genome Size (N)', fontsize=10)
    ax.set_ylabel(f'Search Time (s) for {num_patterns} Patterns', fontsize=10) # <-- FIX HERE
    ax.legend(fontsize=8)
    ax.grid(True, linestyle='--')
    ax.ticklabel_format(axis='x', style='plain')
    
    plt.tight_layout(rect=[0, 0.03, 1, 0.97])
    plt.savefig(filename)
    print(f"\n[SUCCESS] Benchmark results saved to {filename}")

# -----------------------------
# Main Execution (FIXED PLOT CALL)
# -----------------------------
if __name__ == "__main__":
    genome_sizes = [10000, 100000, 500000] # Increased 500k back for better visualization
    pattern_length = 10
    num_patterns = 100
    
    all_results = [] # Stores data for plotting

    # Table Header
    print(f"{'Structure':<15} | {'Size':<8} | {'Build Time':<10} | {'Memory (KB)':<12} | {'Search Time':<10}")
    print("-" * 75)

    for size in genome_sizes:
        # Generate random genome string
        s = ''.join(r.choice('ACGT') for _ in range(size))

        # --- MEASURE THE PLAIN FM INDEX COST (bwtConst) ---
        tracemalloc.start()
        start_build_plain = t.perf_counter()
        bwt_obj = bwtConst(s) 
        end_build_plain = t.perf_counter()
        plain_mem_usage, _ = tracemalloc.get_traced_memory() 
        tracemalloc.stop()

        plain_build_time = end_build_plain - start_build_plain
        # -------------------------------------------------------------

        # List of strategies to test
        strategies = [
            ("FMIndexPlain", lambda: bwt_obj, (None,), backward_search_bwt, plain_build_time, plain_mem_usage),
            ("OCCCompressed", OCCCompressed, (bwt_obj, 16), backward_search_occ, 0, 0),
            ("WaveletTree", WaveletTree, (bwt_obj,), backward_search_wavelet, 0, 0),
            ("BiFMIndex", BiFMIndex, (bwt_obj,), backward_search_bifm, 0, 0),
        ]

        for name, cls, args, search_func, pre_time, pre_mem in strategies:
            
            result_entry = {'Structure': name, 'Size': size}
            
            # A. Measure Construction (Time + Memory)
            if name == "FMIndexPlain":
                obj = bwt_obj
                build_time = pre_time
                mem_usage = pre_mem
            else:
                obj, build_time, mem_usage = measure_construction(cls, *args)
                
            result_entry['BuildTime'] = build_time
            result_entry['MemoryKB'] = mem_usage / 1024

            # B. Measure Search Speed
            matches_found = 0
            patterns = [random_pattern(pattern_length) for _ in range(num_patterns)]
            
            start_search = t.perf_counter()
            for pat in patterns:
                # Handle search function routing
                if name == "FMIndexPlain":
                    if backward_search_bwt(pat, obj):
                        matches_found += 1
                else:
                    if search_func(pat, obj):
                        matches_found += 1
                        
            end_search = t.perf_counter()
            search_time = end_search - start_search
            result_entry['SearchTime'] = search_time

            all_results.append(result_entry)
            
            # Print to console
            print(f"{name:<15} | {size:<8} | {build_time:.5f}s   | {mem_usage / 1024:<12.2f} | {search_time:.5f}s")
        
        print("-" * 75)
        
    # Generate the graphical output (FIXED CALL)
    if all_results:
        plot_results(all_results, num_patterns=num_patterns)
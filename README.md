1. Background
Modern sequence aligners, such as BWA, Bowtie, and SOAP2, utilize the Burrows–Wheeler Transform (BWT) and its FM-index to perform substring searches within genomes. The FM-index enables pattern matching in time proportional to the query length rather than the genome size, which is important for large-scale alignment tasks. However, when genomes become larger and data grows, even FM-index–based searches can become bottlenecks. Optimizing the data structures and access patterns used in BWT-based searching can significantly reduce alignment times and memory usage. This project aims to find performance optimizations beyond the base FM-index implementation, benchmark them, and visualize their effects on alignment speed across genomes of different sizes.
2. Project Objectives
1.	Implement a baseline BWT +FM-index to perform exact substring searches on DNA sequences.
2.	Implement selected optimizations to improve search efficiency, such as:
o	Wavelet tree–based rank/select operations
o	Compressed occurrence tables (OCC arrays with checkpointing)
o	Bidirectional BWT for flexible pattern extension
3.	Benchmark and visualize performance across multiple genome sizes (e.g., 100k, 1M, and 10M base pairs).
4.	Analyze trade-offs between speed, memory, and index construction time.
3. Methods
1. Baseline Implementation
•	Build a basic FM-index implementation using the BWT of a DNA string.
•	Implement backward search to count occurrences of a query substring.
2. Optimizations
At least two optimizations will be added incrementally:
•	Optimization A:
Implement checkpointed OCC tables (store character occurrence counts at intervals and compute intermediate counts using bitwise operations).
•	Optimization B:
Implement a bidirectional FM-index to support extensions in both directions without rebuilding the index.
(MAYBE) explore a wavelet tree–based rank/select structure for compressed access.
3. Benchmarking
•	Measure runtime for substring searches across increasing genome sizes.
•	Measure memory usage of each variant.
•	Perform multiple runs and record average runtime.
4. Visualization
•	Plot runtime vs. genome size for:
o	Baseline FM-index
o	FM-index + Optimization A
o	FM-index + Optimization A + Optimization B
•	Use matplotlib to generate graphs.

5. Expected Outcomes
•	A working FM-index search implementation.
•	At least two optimized variants with measurable performance improvements.
•	Graphs demonstrating speed and scalability improvements.
•	A short analysis report discussing results and potential further optimizations.
6. Timeline
		
Week 1	Literature review on BWT, FM-index, and common optimizations. Implement baseline FM-index and backward search.
Week 2	Implement Optimization A (checkpointed OCC tables) and test for correctness.
Week 3	Implement Optimization B (bidirectional FM-index) and expand benchmarks.
Week 4	Benchmark all versions using different genome lengths. Record runtimes and memory usage.
Week 5	Visualize results and prepare final analysis/report.

8. References
•	Li, H., & Durbin, R. (2009). Fast and accurate short read alignment with Burrows–Wheeler transform. Bioinformatics, 25(14), 1754–1760.
•	Ferragina, P., & Manzini, G. (2000). Opportunistic data structures with applications. Proceedings of FOCS.
•	Langmead, B., & Salzberg, S. L. (2012). Fast gapped-read alignment with Bowtie 2. Nature Methods, 9(4), 357–359.



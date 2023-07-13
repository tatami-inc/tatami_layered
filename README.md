# tatami helpers for creating layered matrices

Layered matrices are a space optimization of sparse matrices containing small positive counts,
where we store rows in different "layers" depending on whether their maximum count is large enough to fit into an unsigned 8-bit, 16-bit or 32-bit integer.
This reduces the memory usage compared to naively storing counts for all rows in the largest integer size across the entire matrix.
It is intended to be used with gene expression data where different genes (rows) can vary widely in their expression.

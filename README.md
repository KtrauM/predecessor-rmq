# Predecessor and Range Minimum Query Data Structures
- Implementation of Predecessor and Range Minimum Queries.
- Predecessor queries are implemented using Elias Fano coding using succinct bit vector data structure with rank/select query support.
- For RMQ three different data structures with O(N^2), O(NlogN) and O(N) time and space complexities are used. 
- To achieve linear complexity in RMQ, input array is divided into subblocks, and encoded using the corresponding cartesian tree of the block.

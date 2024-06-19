# wglib
Library for tunnelled Wheeler Graph (WG).

File src/lib.rs contains class for tunnelled WG with structures supporting FM-index. The tunnelled version is created using heuristics described in work "Efficient construction of compressed index for large text collections" written by Klára Sládečková, available in the file named Masters_thesis.pdf. The file comprises the class, and implemented heuristics as well as functions for reconstruction of the tunnelled WG and various tests showing the functionality of the functions as well as their usage.

The package DATA contains some of the data used in the tests. Other data are not provided, because of their size.   

For full functionality, please download simple_sds library using the last link.

## useful links
https://github.com/andynet/pfg
https://github.com/andynet/pfp_wg
https://github.com/jltsiren/simple-sds


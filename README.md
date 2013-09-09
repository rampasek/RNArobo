RNArobo
=======

An RNA structural motif searching tool. RNArobo can search sequence databases in FASTA format for a motif defined by a "descriptor", which can specify primary and secondary structure constraints.

Compile RNArobo:
----------------
* `tar -zxf rnarobo-2.1.0.tar.gz`
* `cd rnarobo-2.1.0/`
* `make`


Run RNArobo by the command:
---------------------------
`./rnarobo [options] <descriptor-file> <sequence-file>`  
e.g.: `./rnarobo -cu motif.des db.fa > occurrences.txt`


Available RNArobo Options:
--------------------------
`-c`&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;search both strands of the database  
`-u`&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;report only non-overlapping occurrences  
`-f`&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;print output in plain FASTA format  
`-s`&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;print output in FASTA format with element separators  
`--nratio FLOAT`  &nbsp;set max allowed ratio of “N”s in reported occurrences to their length;
must be within `<0,1>`

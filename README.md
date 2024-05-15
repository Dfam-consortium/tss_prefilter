
# Tensor Slide Sketch - Prefilter For Alignment (tss_prefilter)

## Description

This is an experimental Rust implementation of an alignment prefilter.  
The prefilter is based on the Tensor Slide Sketch (TSS) algorithm [https://www.biorxiv.org/content/10.1101/2020.11.13.381814v5.abstract]
and nearest neighbor search using the FAISS library [https://github.com/facebookresearch/faiss]. Used in conjunction 
with a more expensive alignment step, the prefilter reduces the number of sequences that need to be aligned. 

*TSS Method*
Fast Alignment-Free Similarity Estimation By Tensor Sketching
Amir Joudaki, Gunnar Rätsch, André Kahles
bioRxiv 2020.11.13.381814; doi: https://doi.org/10.1101/2020.11.13.381814

*TSS Align Method*
Joudaki, A., Meterez, A., Mustafa, H., Koerkamp, R.G., Kahles, A. and Rätsch, G., 2023. 
Aligning distant sequences to graphs using long seed sketches. Genome Research, 33(7), pp.1208-1217.


## Dependencies

FAIS 1.8.0  : https://github.com/facebookresearch/faiss/archive/refs/tags/v1.8.0.tar.gz

```
% tar zxvf v1.8.0.tar.gz
% cd faiss-1.8.0
% /usr/local/cmake/bin/cmake -DFAISS_ENABLE_GPU=OFF -DFAISS_ENABLE_C_API=ON -DBUILD_SHARED_LIBS=ON -B build .
% make -C build -j faiss
% make -C build install
```
### Does not appear to install C library files by default

```
% cd build/c_api
% cp lib* /usr/local/lib
```

## Build

```
% cargo build --release
```

## Example

### Build a small index on a subset of the Dfam TE consensus sequences

```
% head -n 1000 data/Dfam36-human.fa > sequences.fasta
% /bin/time target/release/tss_prefilter build
```

### Search for the nearest neighbors to query sequence(s)

First create a file called 'query.fasta' containing one or more query sequences in FASTA format.
Here we use mir.fa from the data directory.
Then run the search command.

```
% cp data/mir.fa query.fasta
% target/release/tss_prefilter search
```

### Build a larger index using a chunk of Dfam36 representing all of the human TEs (plus some more)

```
% rm faiss.index
% cp data/Dfam36-human.fa > sequences.fasta
% /bin/time target/release/tss_prefilter --build
```

### Search for the nearest neighbors to SVA
This TE family is known to contain a bit of Alu, HERVK, and LTR5 elements

```
% cp data/sva_a.fa query.fasta
% target/release/tss_prefilter --search  
```

### Align candidates and lookup TE names for identifiers
This is just an example script and will not work out-of-the-box.  It is meant to
demonstrate how to use the output of the search command to align the candidates.
It also performs a lookup of the TE names for the identifiers.  NOTE: It can be
made to work by changing the path to the 'align.pl' script if you have
installed a RepeatMasker/RepeatModeler distribution on your system.

```
/bin/time target/release/tss_prefilter --search | scripts/alignQueryToSet.pl
```

## Implementation Details and Work in Progress

 * This currently uses a 4 character alphabet.  All IUB codes are converted to 'A's.
 * The sketch parameters and k-mer window size were chosen to be compatible with the
   TSS Align method. 
 * Took 3hr 15min (single threaded) to build an index with 44,294 TE sequences (71MB).  The
   index is 684MB in size.
      - With 10 nearest neighbors for each plus-strand K-mer in SVA (stride=1), and 
        10 nearest neighbors for each minus-strand K-mer in SVA (stride=1), the search
        took < 1s and produced 10,847 candidates.
      - Aligning (threads=4) SVA to the 10,847 candidates took 12s.
      - Aligning (threads=4) SVA to the full 44,294 TE set took 39s.


Robert Hubley 5/2024


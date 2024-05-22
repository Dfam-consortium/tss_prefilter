
# Tensor Slide Sketch - Prefilter For Alignment (tss_prefilter)

## Description

This is an experimental Rust implementation of an alignment prefilter. 
The prefilter is based on the Tensor Slide Sketch (TSS) algorithm 
and nearest neighbor search using the FAISS library 
[https://github.com/facebookresearch/faiss]. Used in conjunction with 
a more expensive alignment step, the prefilter reduces the number 
of sequences that need to be aligned. 

*TSS Method*  
Fast Alignment-Free Similarity Estimation By Tensor Sketching
Amir Joudaki, Gunnar Rätsch, André Kahles
bioRxiv 2020.11.13.381814; doi: https://doi.org/10.1101/2020.11.13.381814

*TSS Align Method*  
Joudaki, A., Meterez, A., Mustafa, H., Koerkamp, R.G., Kahles, A. and Rätsch, G., 2023. 
Aligning distant sequences to graphs using long seed sketches. Genome Research, 33(7), pp.1208-1217.
https://genome.cshlp.org/content/33/7/1208.full.pdf


## Dependencies
You will need the FAISS dynamic libraries installed to build this project.
FAIS 1.8.0  : https://github.com/facebookresearch/faiss/archive/refs/tags/v1.8.0.tar.gz

```
% tar zxvf v1.8.0.tar.gz
% cd faiss-1.8.0
% /usr/local/cmake/bin/cmake -DFAISS_ENABLE_GPU=OFF -DFAISS_ENABLE_C_API=ON -DBUILD_SHARED_LIBS=ON -B build .
% make -C build -j faiss
% make -C build install
```

Does not appear to install C library files by default
```
% cd build/c_api
% cp lib* /usr/local/lib
```

## Build

```
% cargo build --release
```

## Example

#### Build a small index on a subset of the Dfam TE consensus sequences

```
% head -n 1000 data/Dfam36-human.fa > sequences.fasta
% /bin/time target/release/tss_prefilter build
```

#### Search for the nearest neighbors to query sequence(s)

First create a file called 'query.fasta' containing one or more query sequences in FASTA format.
Here we use mir.fa from the data directory.
Then run the search command.

```
% cp data/mir.fa query.fasta
% target/release/tss_prefilter search
```

#### Build a larger index using a chunk of Dfam36 representing all of the human TEs (plus some more)

```
% rm faiss.index
% cp data/Dfam36-human.fa > sequences.fasta
% /bin/time target/release/tss_prefilter --build
```

#### Search for the nearest neighbors to SVA
This TE family is known to contain a bit of Alu, HERVK, and LTR5 elements

```
% cp data/sva_a.fa query.fasta
% target/release/tss_prefilter --search  
```

#### Align candidates and lookup TE names for identifiers
This is just an example script and will not work out-of-the-box.  It is meant to
demonstrate how to use the output of the search command to align the candidates.
It also performs a lookup of the TE names for the identifiers.  NOTE: It can be
made to work by changing the path to the 'align.pl' script if you have
installed a RepeatMasker/RepeatModeler distribution on your system.

```
/bin/time target/release/tss_prefilter --search | scripts/alignQueryToSet.pl
```

## Stats

 UPDATE: The following stats are out of date.  The current implementation is much faster at building
         the index.  Turns out that FAISS.add_with_ids() is parallelized, so calling this function
         with a large batch of vectors is much faster than adding them one at a time. 

 * Small Index ( 44,294 TE sequences, 71 mbases ):
    * Index Build: 3hr 15min
    * Index Size: 684MB
    * With 10 nearest neighbors for each plus-strand K-mer in SVA (stride=1), and 
      10 nearest neighbors for each minus-strand K-mer in SVA (stride=1), the search
      took < 1s and produced 10,847 candidates.
    * Aligning (threads=4) SVA to the 10,847 candidates took 12s.
    * Aligning (threads=4) SVA to the full 44,294 TE set took 39s.
    * At this rate it would take 10 days to build a monolithic index of all Dfam_3.8 3.5M TE sequences with
      an estimated size of 54GB for a monolithic index buitl on a single cpu.

 * Dfam_3.6 ( 732,993 TE sequences, 968 mbases )
    * Index Build: 65hr 15min ( monolithic index built on single cpu )
    * Index Size: 8.9GB
    * Searching SVA took 8.3 seconds ( TODO: need to break out index search vs result [eg. sorting and duplicate removal] processing time ).  This produced 19,874 candidates.
    * Aligning (threads=4) SVA to the 19,874 candidates took 24s (692 alignments)
    * Aligning (threads=4) SVA to the full 732,993 TE set took 9 min 51s (3112 alignments)
    * NOTE: Expected results of Alu/HERVK were not returned in this search as there are many results that are equally distant in the larger index. Choice of nearest neighbor parameter will impact this but also increase the number of candidates requiring alignment.  



## Implementation Details and Work in Progress

 * This currently uses a 4 character alphabet.  All IUB codes are converted to 'A's.
 * The sketch parameters and k-mer window size were chosen to be compatible with the
   TSS Align method. 
 * Unlike tss_align we do not store a unique id per k-mer, rather we store the ID of the family from 
   which the k-mer originated.  This simplification allows for quickly referencing the family, but 
   does not allow for easily locating the position of the indexed k-mer within the family.
 * Explore building index in parallel on independent segments of the database and join
   them in the end to create a monolithic index *or* search them independently and join
   the results. E.g. https://github.com/facebookresearch/faiss/wiki/Indexing-1T-vectors#building-the-index
 * The current implementation collapses all k-mer matches to the set of distinct families identified by the
   k-mers.  The positional information of the kmers nor the relationship between them are utilized to further
   cull the set of candidate families.  That means that a single k-mer (80bp) match between query/index is
   sufficient to return a family.  We should experiment with increasing the nearest neighbor parameter and 
   further filtering the return set to increase sensitivity while reducing the number of candidates.
 * How costly is it to delete or incrementally add new families to the index?



Robert Hubley 5/2024


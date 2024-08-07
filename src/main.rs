use std::collections::HashMap;
use std::fs::File;
use std::io::{BufReader, BufRead, Result};
use std::time::Instant;

use clap::Args;
use clap::Command;
use clap::FromArgMatches;

extern crate faiss;
use faiss::{Index, index_factory, IdMap};
use faiss::index::io::write_index;
use faiss::index::io::read_index;
use faiss::index::Idx;

pub mod tensor_sketch;  
use tensor_sketch::TensorSketch;



fn parse_fasta(filename: &str) -> Result<Vec<Vec<u8>>> {
    let file = File::open(filename)?;
    let reader = BufReader::new(file);

    let mut sequences = Vec::new();
    let mut current_sequence = Vec::new();

    for line in reader.lines() {
        let line = line?;
        if line.starts_with('>') {
            // New sequence started, store previous sequence if any
            if !current_sequence.is_empty() {
                sequences.push(current_sequence.clone());
                current_sequence.clear();
            }
        } else {
            // Process sequence data
            let converted_sequence: Vec<u8> = line
                .chars()
                .map(|c| match c {
                    'A' => 0,
                    'C' => 1,
                    'G' => 2,
                    'T' => 3,
                    'a' => 0,
                    'c' => 1,
                    'g' => 2,
                    't' => 3,
                    _ => 0, // For all other IUB codes
                })
                .collect();
            current_sequence.extend(converted_sequence);
        }
    }

    // Store the last sequence
    if !current_sequence.is_empty() {
        sequences.push(current_sequence);
    }

    Ok(sequences)
}

fn reverse_complement(dna: &[u8]) -> Vec<u8> {
    // Create a mapping for complement bases
    let complement = [3, 2, 1, 0]; // A <-> T, C <-> G

    // Convert each base to its complement
    let complemented: Vec<u8> = dna.iter().map(|&base| complement[base as usize]).collect();

    // Reverse the complemented vector
    let mut reversed_complement = complemented;
    reversed_complement.reverse();

    reversed_complement
}

#[derive(Debug)]
struct Neighbor {
  min_distance: f32,
  label_count: u32,
}


#[derive(Args, Debug)]                
#[command(author, version, about, long_about = None)] // Read from `Cargo.toml`
struct CLI {                    
    /// Build the index 'faiss.index' from the current hard-coded input filename 'sequences.fasta'
    #[arg(short, long)]
    build: bool,
                
    /// Search the index 'faiss.index' using a hard coded search query filename 'query.fasta'
    #[arg(short, long)]               
    search: bool,           
   
    /// Experiment #1: Compute the minimum distance between any 1-bp shifted window (k_param) of A to any non-overlapping window of B
    #[arg(short, long)]               
    exp1: bool,           
}

fn main() {

    // Using builder interface to support a custom help template
    let cli = Command::new("tss_prefilter").help_template(
"
{before-help}{name} {version}
{author-with-newline}{about-with-newline}
tss_prefilter is a tool for indexing/searching a DNA sequence collection for alignment candidates
...
The index is comprised of Tensor Slide Sketches stored in Facebook AI Similarity Search system (FAISS).

{usage-heading} {usage}

  Examples:

    # Build the index
    ./tss_prefilter --build

    # Search the index
    ./tss_prefilter --search

    # Experiment #1 (expects a.fa and b.fa to already exist )
    #   - Compute the minimum distance between any 1-bp shifted window (k_param) of A to any non-overlapping window of B
    ./tss_prefilter --exp1 

{all-args}{after-help}
",
    );
    // NOTE: This is quite convoluted...I am not sure if there is
    //       a simpler way to mix the Build/Derive APIs
    let cli = CLI::augment_args(cli);
    let args = CLI::from_arg_matches(&cli.get_matches())
        .map_err(|err| err.exit())
        .unwrap();

    // Parameters for sketching (TODO: Expose and set defaults)
    let alphabet_size = 4; // NOTE: all IUB codes are mapped to 'A'

    // 
    // Metagraph Parameters 0.629 Sp. Corr @ len 126:
    //
    //let k_param = 80;
    //let sketch_dim = 14;
    //let s_param = 8;   // AKA stride length
    //let subsequence_len = 6; // AKA tuple length
    //let w_param = 16;  // AKA windows size

    // Set #2: 80,15,7,1,20 0.892 Sp. Corr @ len 135
    //let k_param = 80;
    //let sketch_dim = 15;
    //let s_param = 7;
    //let subsequence_len = 1;
    //let w_param = 20;

    let k_param = 80;
    let sketch_dim = 4;
    let s_param = 4;
    let subsequence_len = 1;
    let w_param = 13;

    // #2-a 80,7,4,1,13  0.894381 Sp. Corr @ len 119
    // #2-b 80,9,4,1,14  0.89689 Sp. Corr @ len 153
    // Set #3 80,7,4,1,13 0.889 Sp. Corr @ len 119
    // Set #4 80,14,4,1,13  0.901 Sp. Corr @ len 238!

    let seed = 0; // Seed for the random number generator
 

    // Create a TensorSketch instance
    let mut tensor: TensorSketch<u8> = TensorSketch::new(alphabet_size, sketch_dim, subsequence_len, seed);

    if subsequence_len == 1 && sketch_dim == 4 {
      println!("WARNING: Using hardcoded translation hashes for 1-tuple sketching.");
      let h = vec![ vec![0, 1, 2, 3] ];
      let s = vec![ vec![true, false, false, false] ];
      tensor.set_hashes_for_testing(h, s);
    }

    // Experiment #1
    // Take two sequences A and B, and compute the minimum distance between any 1-bp shifted
    // window (k_param) of A to any non-overlapping window of B.
    if args.exp1 {
        let mut b_sketches = Vec::new();
        println!("Reading sequences from b.fa and computing non-overlapping sketches");
        match parse_fasta("b.fa") {
            Ok(sequences) => {
               b_sketches = tensor.compute_slide_sketch_2d(&sequences[0], k_param, k_param, w_param, s_param);
            }
            Err(e) => eprintln!("Error parsing FASTA file: {:?}", e),
        }
        println!("Reading sequences from a.fa and computing non-overlapping sketches");
        match parse_fasta("a.fa") {
            Ok(sequences) => {
               let a_sketches = tensor.compute_slide_sketch_2d(&sequences[0], k_param, 1, w_param, s_param);
               for i in 0..a_sketches.len() {
                   let mut min_distance = std::f64::MAX;
                   let mut min_distance_index = 0;
                   for j in 0..b_sketches.len() {
                       let distance = TensorSketch::<u8>::l2_dist(&a_sketches[i], &b_sketches[j]);
                       if distance < min_distance {
                           min_distance = distance;
                           min_distance_index = j;
                       }
                   }
                   //println!("Minimum distance between A[{}] and B[{}]: {:?}", i, min_distance_index, min_distance);
                   println!("{},{},{}", i, (min_distance_index*80), min_distance);
                }
            }
            Err(e) => eprintln!("Error parsing FASTA file: {:?}", e),
        }
    }

    if args.build {
        // Create a new index
        println!("Initializing the index...");
        let mut t = Instant::now();
        let num_sketches_in_kmer = ((k_param - w_param + 1) as f64 / s_param as f64).ceil() as u32; 
        let index = index_factory(sketch_dim as u32 * num_sketches_in_kmer, "HNSW32", faiss::MetricType::L2).unwrap();
        let sketch_vec_size = num_sketches_in_kmer as usize * sketch_dim as usize;
        let mut index = IdMap::new(index).unwrap();
        println!("  Index initialized in {:?}", t.elapsed());

        // Read sequences from the FASTA file
        println!("Reading sequences from the FASTA file...");
        match parse_fasta("sequences.fasta") {
            Ok(sequences) => {
                let mut batch = Vec::new();
                let mut ids = Vec::new();
                for (i, sequence) in sequences.iter().enumerate() {

                    // For debugging/timing
                    //println!("Sequence {}: {:?}", i + 1, sequence);
                    //for j in 0..sequence.len() {
                    //    if sequence[j] == 0 {
                    //        print!("A");
                    //    } else if sequence[j] == 1 {
                    //        print!("C");
                    //    } else if sequence[j] == 2 {
                    //        print!("G");
                    //    } else if sequence[j] == 3 {
                    //        print!("T");
                    //    } else {
                    //        println!("OOOPS....this sequence has a bad base! {:?}", sequence[j]);
                    //    }
                    // }
                    //println!();
                    //print!("Sketching Sequence {} ", i + 1);
                    //t = Instant::now();

                    // This should be < window size (80) however it failed with sequences < 100
                    // TODO: figure out why
                    if sequence.len() < 100 {
                            continue;
                    }
                    let sketches = tensor.compute_slide_sketch_1d(&sequence, k_param, k_param, w_param, s_param);

                    //println!("Sketches of the sequence: {:?}", sketches);
    
                    let kmer_count = (sequence.len() as f64 / 80.0).floor() as usize;

                    // For debugging/timing TODO remove
                    if sketches.len() / sketch_vec_size != kmer_count {
                        println!("OOOPS....this sketch is not the right length! {:?}", sketches.len());
                    }

                    batch.extend(sketches);
                    for _j in 0..kmer_count {
                        ids.push(Idx::new(i as u64));
                    }
                    //println!("in {:?}", t.elapsed());

                    // Add vectors to the index once we reach an adequate batch size
                    if i > 0 && i % 100000 == 0 {
                        println!("Adding a batch of sketches ( sequences {} .. {} ) to the index...", i - 100000, i);
                        t = Instant::now();
                        // This form of index does not need training.  This probably just did
                        // nothing anyway, but for efficiency I removed the call altogether.
                        //index.train(&sketches).unwrap();
                        index.add_with_ids(&batch, &ids).unwrap();
                        batch.clear();
                        ids.clear();
                        println!("Added sketches to the index in {:?}", t.elapsed());
                    }
                }
                if batch.len() > 0  {
                    println!("Adding the last batch of sketches to the index...");
                    // This form of index does not need training.  This probably just did
                    // nothing anyway, but for efficiency I removed the call altogether.
                    //index.train(&sketches).unwrap();
                    index.add_with_ids(&batch, &ids).unwrap();
                    println!("Added sketches to the index in {:?}", t.elapsed());
                }
            }
            Err(e) => eprintln!("Error parsing FASTA file: {:?}", e),
        }
        write_index(&index, "index.faiss").unwrap();
        println!("Index created successfully!");
    }

    if args.search {
        let mut t = Instant::now();
        let mut index = read_index("index.faiss").unwrap();
        println!("# Read index into memory in {:?}", t.elapsed());

        let nearest_neighbors = 20;
        t = Instant::now();
        match parse_fasta("query.fasta") {
            Ok(sequences) => {
                for (i, sequence) in sequences.iter().enumerate() {
                    println!("Sketching Query {}", i + 1);
                    
                    let sketches = tensor.compute_slide_sketch_1d(&sequence, k_param, 1, w_param, s_param);

                    //println!("Sketches of the query: {:?}", sketches.len());
                    //println!("Sketches of the query: {:?}", sketches);
                    
                    let result = index.search(&sketches, nearest_neighbors).unwrap();

                    println!("total forward results: {:?}", result.distances.len());
                    //println!("forward nodes: {:?}", result.labels);
                    println!("forward distances: {:?}", result.distances);

                    let rc_sequence = reverse_complement(&sequence);
                    let rc_sketches = tensor.compute_slide_sketch_1d(&rc_sequence, k_param, 1, w_param, s_param);
                    let rc_result = index.search(&rc_sketches, nearest_neighbors).unwrap();

                    println!("total reverse results: {:?}", rc_result.distances.len());
                    //println!("reverse nodes: {:?}", rc_result.labels);
                    println!("reverse distances: {:?}", rc_result.distances);
    
                    let mut all_results = result;
                    all_results.distances.extend(rc_result.distances);
                    all_results.labels.extend(rc_result.labels);
    
                    // Print them out
                    let mut label_counts: HashMap<usize, Neighbor> = HashMap::new();
                    for (distance, label) in all_results.distances.iter().zip(all_results.labels.iter()) {
                        let label_val = match label.get() {
                                                Some(value) => value as usize,
                                                None => {
                                                    continue;
                                                }
                                            };
                        let label_count = label_counts.entry(label_val).or_insert(Neighbor{min_distance: *distance, label_count: 1});
                        if distance < &label_count.min_distance {
                            label_count.min_distance = *distance;
                        }
                        label_count.label_count += 1;
                    }
    
                    // Collect the key-value pairs into a vector of tuples
                    let mut sorted_distances: Vec<_> = label_counts.iter().collect();
    
                    // Sort the vector by distances
                    sorted_distances.sort_by(|(_, neighbor1), (_, neighbor2)| {
                        neighbor1.min_distance.partial_cmp(&neighbor2.min_distance).unwrap()
                    });
    
                    // Iterate over the sorted vector
                    for (key, neighbor) in sorted_distances {
                        println!("Key: {}, Neighbor: {:?}", key, neighbor);
                    }
    
                }
            }
            Err(e) => eprintln!("Error parsing FASTA file: {:?}", e),
        }
        println!("# Total search time {:?}", t.elapsed());
    }
}


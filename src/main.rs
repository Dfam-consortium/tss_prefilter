use std::collections::HashMap;
use std::fs::File;
use std::io::{BufReader, BufRead, Result};

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
    let sketch_dim = 14;
    let subsequence_len = 6;
    let k_param = 80;
    let w_param = 16;
    let s_param = 8;
    let seed = 0; // Seed for the random number generator

    // Create a TensorSketch instance
    let tensor: TensorSketch<u8> = TensorSketch::new(alphabet_size, sketch_dim, subsequence_len, seed);

    if args.build {
        // Create a new index
        let num_sketches_in_kmer = ((k_param - w_param + 1) as f64 / s_param as f64).ceil() as u32; 
        let index = index_factory(sketch_dim as u32 * num_sketches_in_kmer, "HNSW32", faiss::MetricType::L2).unwrap();
        let mut index = IdMap::new(index).unwrap();

        match parse_fasta("sequences.fasta") {
            Ok(sequences) => {
                for (i, sequence) in sequences.iter().enumerate() {

                    //println!("Sequence {}: {:?}", i + 1, sequence);
                    println!("Sketching Sequence {}", i + 1);
                    if sequence.len() < 100 {
                            continue;
                    }
                    let sketches = tensor.compute_slide_sketch_1d(&sequence, 80, 80, 16, 8);
                    //println!("Sketches of the sequence: {:?}", sketches);
    
                    let kmer_count = (sequence.len() as f64 / 80.0).floor() as usize;
    
                    let mut ids = Vec::with_capacity(kmer_count);
                    for _j in 0..kmer_count {
                        ids.push(Idx::new(i as u64));
                    }
                
                    // Add vectors to the index
                    index.train(&sketches).unwrap();
                    index.add_with_ids(&sketches, &ids).unwrap();
                }
            }
            Err(e) => eprintln!("Error parsing FASTA file: {:?}", e),
        }
        write_index(&index, "index.faiss").unwrap();
        println!("Index created successfully!");
    }

    if args.search {
        let mut index = read_index("index.faiss").unwrap();

        let nearest_neighbors = 10;
        match parse_fasta("query.fasta") {
            Ok(sequences) => {
                for (i, sequence) in sequences.iter().enumerate() {
                    println!("Sketching Query {}", i + 1);
                    if sequence.len() < 100 {
                        continue;
                    }
                    let sketches = tensor.compute_slide_sketch_1d(&sequence, 80, 1, 16, 8);
                    let result = index.search(&sketches, nearest_neighbors).unwrap();
    
                    let rc_sequence = reverse_complement(&sequence);
                    let rc_sketches = tensor.compute_slide_sketch_1d(&rc_sequence, 80, 1, 16, 8);
                    let rc_result = index.search(&rc_sketches, nearest_neighbors).unwrap();
    
                    let mut all_results = result;
                    all_results.distances.extend(rc_result.distances);
                    all_results.labels.extend(rc_result.labels);
    
                    // Print them out
                    let mut label_counts: HashMap<usize, Neighbor> = HashMap::new();
                    for (distance, label) in all_results.distances.iter().zip(all_results.labels.iter()) {
                        let label_count = label_counts.entry(label.get().unwrap() as usize).or_insert(Neighbor{min_distance: *distance, label_count: 1});
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
    }
}




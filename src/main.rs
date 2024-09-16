use std::cmp::Ordering;
use std::collections::HashMap;
use std::collections::BinaryHeap;
use std::fs::File;
use std::io::{BufReader, BufRead, Result};
use std::time::Instant;
use rand::Rng;

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


fn edit_distance<T: PartialEq>(s1: &[T], s2: &[T]) -> usize {
    let m = s1.len();
    let n = s2.len();

    if m == 0 {
        return n;
    }
    if n == 0 {
        return m;
    }

    let mut costs = vec![0; n + 1];

    for k in 0..=n {
        costs[k] = k;
    }

    for (i, it1) in s1.iter().enumerate() {
        costs[0] = i + 1;
        let mut corner = i;

        for (j, it2) in s2.iter().enumerate() {
            let upper = costs[j + 1];
            if it1 == it2 {
                costs[j + 1] = corner;
            } else {
                let t = if upper < corner { upper } else { corner };
                costs[j + 1] = if costs[j] < t { costs[j] } else { t } + 1;
            }

            corner = upper;
        }
    }
    costs[n]
}

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
   
    /// Experiment #2: Compute the minimum distance between any 1-bp shifted window (k_param) of query.fasta to any non-overlapping window of sequences.fasta
    #[arg(short, long)]               
    exp2: bool,           
    
    /// Experiment #3: Compute the average magnitude of the sketches over a sequence
    #[arg(short, long)]               
    exp3: bool,           
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

    # Experiment #2 (expects query.fasta and sequences.fasta to already exist )
    #   - Compute the minimum distance between any 1-bp shifted window (k_param) of query to any non-overlapping window of sequences
    ./tss_prefilter --exp2
    
    # Experiment #3: Compute the average magnitude of the sketches over a sequence
    ./tss_prefilter --exp3

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
    let k_param = 80;
    let sketch_dim = 14;
    let s_param = 8;   // AKA stride length
    let subsequence_len = 6; // AKA tuple length
    let w_param = 16;  // AKA windows size

    // Set #1.....just lowering tuple size
    //let k_param = 80;
    //let sketch_dim = 14;
    //let s_param = 8;   // AKA stride length
    //let subsequence_len = 3; // AKA tuple length
    //let w_param = 16;  // AKA windows size


    // Set #2: 80,15,7,1,20 0.892 Sp. Corr @ len 135
    //let k_param = 80;
    //let sketch_dim = 15;
    //let s_param = 7;
    //let subsequence_len = 1;
    //let w_param = 20;

    // my preferred set
    //let k_param = 80;
    //let sketch_dim = 4;
    //let s_param = 4;
    //let subsequence_len = 1;
    //let w_param = 13;

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

    //
    // Experiment #3: 
    //   Given query.fasta, sketch all subsequences up to a given number of times (replicates),
    //   and calculate the average magnitude of sketches for each subsequence.
    //
    if args.exp3 {
        let replicates = 1000;
        let mut rng = rand::thread_rng();
        println!("Reading sequences from query.fa and computing average sketch magnitudes"); 
        match parse_fasta("query.fasta") {
            Ok(sequences) => {
                let mut mags: Vec<f64> = vec![0.0; sequences[0].len()];
                for i in 0..replicates {
                    let seed: u64 = rng.gen();
                    // Create a TensorSketch instance
                    //println!("Replicate {}: seed={:?}", i, seed);
                    tensor = TensorSketch::new(alphabet_size, sketch_dim, subsequence_len, seed);
                    let q_sketches = tensor.compute_slide_sketch_2d(&sequences[0], k_param, 1, w_param, s_param);
                    for (i, q_sketch) in q_sketches.iter().enumerate() {
                        let v_magnitude = TensorSketch::<u8>::vector_magnitude(&q_sketch);
                        mags[i] += v_magnitude;
                    }
                }   
                for i in 0..sequences[0].len() {
                    println!("avgmag,{},{:?}", i, mags[i] / replicates as f64);
                }
            }
            Err(e) => eprintln!("Error parsing FASTA file: {:?}", e),
        }
    }

    // Experiment #2
    //   Identify the N closest neighbors to a query sequence (query.fasta) subsequences given a 
    //   small subject sequence (sequences.fasta, single sequence only) using brute force ( not FAISS ).
    if args.exp2 {

        #[derive(Debug)]
        struct DistanceInfo {
            distance: f64,
            //l2_dist: f64,
            //dot_prod: f64,
            //cos_sim: f64,
            ed: usize,
            seq_index: usize,
            pos_index: usize,
        }

        //const metric: &str = "dot_product";
        //const metric: &str = "cosine_similarity";
        //const metric: &str = "l2_distance";
        //const metric: &str = "normalized_dot_product";
        const metric: &str = "edit_distance";

        // Sort order for distance metrics
        impl Ord for DistanceInfo {
            fn cmp(&self, other: &Self) -> Ordering {
                if metric == "dot_product" || metric == "cosine_similarity" || metric == "normalized_dot_product" ||
                   metric == "edit_distance" { 
                    // Find maximum metric ( closest matching )
                    self.distance.partial_cmp(&other.distance).unwrap()
                } else {
                    // Find minimum metric ( closest matching )
                    other.distance.partial_cmp(&self.distance).unwrap()
                }
            }
        }

        impl PartialOrd for DistanceInfo {
            fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
                Some(self.cmp(other))
            }
        }

        impl PartialEq for DistanceInfo {
            fn eq(&self, other: &Self) -> bool {
                self.distance == other.distance
            }
        }

        impl Eq for DistanceInfo {}
        let neighbors = 20;
        let mut q_sketches = Vec::new();
        println!("Running experiment #3: Metric={}", metric);
        println!("Reading query.fasta and computing sketches for each subsequence...");
        let mut q_sequence = Vec::new();
        match parse_fasta("query.fasta") {
            Ok(sequences) => {
               q_sketches = tensor.compute_slide_sketch_2d(&sequences[0], k_param, 1, w_param, s_param);
               q_sequence = sequences[0].clone();
            }
            Err(e) => eprintln!("Error parsing FASTA file: {:?}", e),
        }

        let mut closest_distances: Vec<BinaryHeap<DistanceInfo>> = Vec::with_capacity(q_sketches.len());
        for _ in 0..q_sketches.len() {
            closest_distances.push(BinaryHeap::new());
        }

        //for (i, q_sketch) in q_sketches.iter().enumerate() {
        //    let v_magnitude = TensorSketch::<u8>::vector_magnitude(&q_sketch);
        //    println!("vmag,{},{:?}", i, v_magnitude);
        //}

        match parse_fasta("sequences.fasta") {
            Ok(sequences) => {
                // Lazy....this should only process the first sequence if multiple are present
                for (s, sequence) in sequences.iter().enumerate() {
                    println!("Sketching sequence {}", s + 1);
                    let mut s_sketches = tensor.compute_slide_sketch_2d(&sequence, k_param, 1, w_param, s_param);
                    for (i, mut q_sketch) in q_sketches.iter_mut().enumerate() {
                        if metric == "normalized_dot_product" {
                            TensorSketch::<u8>::normalize_vector(&mut q_sketch);
                        }
                        for (j, mut s_sketch) in s_sketches.iter_mut().enumerate() {
                            let qseq = &q_sequence[i..i + k_param];
                            let sseq = &sequence[j..j + k_param];

                            let ed = edit_distance(qseq, sseq);

                            let mut distance = 0.0;
                            if metric == "l2_distance" {
                                distance = TensorSketch::<u8>::l2_dist(&mut q_sketch, &mut s_sketch);
                            } else if metric == "dot_product" {
                                distance = TensorSketch::<u8>::dot_product(&mut q_sketch, &mut s_sketch);
                            } else if metric == "cosine_similarity" {
                                distance = TensorSketch::<u8>::cosine_similarity(&mut q_sketch, &mut s_sketch);
                            } else if metric == "normalized_dot_product" {
                                TensorSketch::<u8>::normalize_vector(&mut s_sketch);
                                distance = TensorSketch::<u8>::dot_product(&mut q_sketch, &mut s_sketch);
                            } else if metric == "edit_distance" {
                                distance = (80 as f64 - ed as f64)/80 as f64;
                            } else {
                                panic!("Unknown metric: {:?}", metric);
                            }
                            //let distance = (80 as f64 - ed as f64)/80 as f64;
                            //println!("Distance between q_sketch pos {} and s_sketch-{} pos {}: {:?}", i, s, /j, distance);
                            //println!("Distance between q_sketch pos {} and s_sketch-{} pos {}: ed={:?}, l2={:?}, dot={:?}, cos={:?}", i, s, j, ed, l2_dist, dot_prod, cos_sim);

                            //let dist_info = DistanceInfo {
                            //    distance: distance, seq_index: s, pos_index: j };
                            let dist_info = DistanceInfo {
                                //distance: distance, l2_dist: l2_dist, dot_prod: dot_prod, 
                                //          cos_sim: cos_sim, ed: ed, seq_index: s, pos_index: j;
                                distance: distance, ed: ed, seq_index: s, pos_index: j 
                              };
                            let heap = &mut closest_distances[i];

                            if  heap.len() < neighbors {
                                heap.push(dist_info);
                            }else if metric == "dot_product" || metric == "cosine_similarity" || metric == "normalized_dot_product" ||
                                     metric == "edit_distance" {
                                // Maximum distances 
                                if heap.peek().unwrap().distance < distance {
                                    heap.pop();
                                    heap.push(dist_info);
                                }
                            }else if heap.peek().unwrap().distance > distance {
                                // Minimium distances
                                heap.pop();
                                heap.push(dist_info);
                            }
                        }
                   }
                }
            }
            Err(e) => eprintln!("Error parsing FASTA file: {:?}", e),
        }

        for (i, heap) in closest_distances.iter().enumerate() {
            println!("{},{},{:?},{:?}", metric, i, heap.peek().unwrap().distance,((80 as f64 - heap.peek().unwrap().ed as f64)/80 as f64));
            //println!("Closest distances for q_sketch {}: {:?}", i, heap);
        }
    }


    // For FAISS testing
    let use_cosine_similarity = 1;
 
    if args.build {
        let mut normalize: bool = false;
        if ( use_cosine_similarity == 1 ) {
          println!("Using cosine similarity");
          normalize = true;
        }else {
            println!("Using L2 distance");
        }

        // Create a new index
        println!("Initializing the index...");
        let mut t = Instant::now();
        let num_sketches_in_kmer = ((k_param - w_param + 1) as f64 / s_param as f64).ceil() as u32; 
        let index;
        if ( use_cosine_similarity == 1 ) {
          index = index_factory(sketch_dim as u32 * num_sketches_in_kmer, "HNSW32", faiss::MetricType::InnerProduct).unwrap();
        }else {
          index = index_factory(sketch_dim as u32 * num_sketches_in_kmer, "HNSW32", faiss::MetricType::L2).unwrap();
        }
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
                    let sketches = tensor.compute_slide_sketch_1d(&sequence, k_param, k_param, w_param, s_param, normalize);

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

        let mut normalize: bool = false;
        if ( use_cosine_similarity == 1 ) {
          println!("Using cosine similarity");
          normalize = true;
        }else {
            println!("Using L2 distance");
        }

        let nearest_neighbors = 40;
        t = Instant::now();
        match parse_fasta("query.fasta") {
            Ok(sequences) => {
                for (i, sequence) in sequences.iter().enumerate() {
                    println!("Sketching Query {}", i + 1);
                    
                    let sketches = tensor.compute_slide_sketch_1d(&sequence, k_param, 1, w_param, s_param, normalize);

                    //println!("Sketches of the query: {:?}", sketches.len());
                    //println!("Sketches of the query: {:?}", sketches);
                    
                    let result = index.search(&sketches, nearest_neighbors).unwrap();

                    println!("total forward results: {:?}", result.distances.len());
                    //println!("forward nodes: {:?}", result.labels);
                    println!("forward distances: {:?}", result.distances);

                    let rc_sequence = reverse_complement(&sequence);
                    let rc_sketches = tensor.compute_slide_sketch_1d(&rc_sequence, k_param, 1, w_param, s_param, normalize);
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


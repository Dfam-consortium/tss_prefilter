use rand::distributions::{Distribution, Uniform};
use rand::{rngs::StdRng, SeedableRng};
use std::f64;

pub struct TensorSketch<SeqType> {
    alphabet_size: SeqType,
    sketch_dim: usize,
    subsequence_len: usize,
    rng: StdRng,
    hashes: Vec<Vec<SeqType>>,
    signs: Vec<Vec<bool>>,
}

impl<SeqType> TensorSketch<SeqType>
where
    SeqType: Copy + From<u8> + Into<usize> + std::marker::Sized,
{
    pub fn new(
        alphabet_size: SeqType,
        sketch_dim: usize,
        subsequence_len: usize,
        seed: u64,
    ) -> Self {
        let mut tensor = TensorSketch {
            alphabet_size,
            sketch_dim,
            subsequence_len,
            rng: StdRng::seed_from_u64(seed),
            hashes: vec![vec![SeqType::from(0); alphabet_size.into()]; subsequence_len],
            signs: vec![vec![false; alphabet_size.into()]; subsequence_len],
        };
        tensor.init();
        tensor
    }

    fn init(&mut self) {
        let hash_dist = Uniform::from(0..self.sketch_dim);
        let bool_dist = Uniform::from(0..2);

        for h in 0..self.subsequence_len {
            for c in 0..self.alphabet_size.into() {
                self.hashes[h][c] = SeqType::from(hash_dist.sample(&mut self.rng) as u8);
                self.signs[h][c] = bool_dist.sample(&mut self.rng) != 0;
            }
        }
    }

    pub fn compute_sketch(&self, seq: &[SeqType]) -> Vec<f64> {
        let mut tp: Vec<Vec<f64>> = vec![vec![0.0; self.sketch_dim]; self.subsequence_len + 1];
        let mut tm: Vec<Vec<f64>> = vec![vec![0.0; self.sketch_dim]; self.subsequence_len + 1];

        tp[0][0] = 1.0;                                                                

        for (i, &c) in seq.iter().enumerate() {
            let c = c.into();
            if c >= self.alphabet_size.into() {
                continue;                     
            }        

            for p in (1..=std::cmp::min(i + 1, self.subsequence_len)).rev() {
                let z = p as f64 / (i as f64 + 1.0);
                let r = self.hashes[p - 1][c].into();
                let s = self.signs[p - 1][c];
                let len = self.sketch_dim;
                if s {
                    for i in 0..len {
                        //println!(" - a[{}] = {}, b[{}] = {}, z = {}, shift = {}, len = {}\n", i, tp[p][i], i, tp[p-1][(len + i - r) % len], z, r, len );
                        tp[p][i] = z * tp[p-1][(len + i - r) % len] + (1.0 - z) * tp[p][i];
                    }
                    for i in 0..len {
                        tm[p][i] = z * tm[p-1][(len + i - r) % len] + (1.0 - z) * tm[p][i];
                    }
                    //println!("S: i: {}, p: {}, c: {}, r: {}, s: {}, z: {}, tp: {:?}, tm: {:?}", i, p, c, r, s, z, tp[p], tm[p]);
                } else {
                    for i in 0..len {
                        tp[p][i] = z * tm[p-1][(len + i - r) % len] + (1.0 - z) * tp[p][i];
                    }
                    for i in 0..len {
                        tm[p][i] = z * tp[p-1][(len + i - r) % len] + (1.0 - z) * tm[p][i];
                    }
                    //println!("!S: i: {}, p: {}, c: {}, r: {}, s: {}, z: {}, tp: {:?}, tm: {:?}", i, p, c, r, s, z, tp[p], tm[p]);
                }
            }                  
        }                                                                                  

        let mut sketch = vec![0.0; self.sketch_dim];
        for m in 0..self.sketch_dim {      
            sketch[m] = tp[self.subsequence_len][m] - tm[self.subsequence_len][m];
            //println!("m: {}, tp[{}][]: {}, tm: {}, sketch: {}", m, self.subsequence_len, tp[self.subsequence_len][m], tm[self.subsequence_len][m], sketch[m]);
        }                       
                                                 
        sketch
    }

            
    pub fn compute_slide_sketch_2d(
        &self,
        seq: &[SeqType],
        k_size: usize,
        k_stride: usize,
        t_size: usize,
        t_stride: usize,
    ) -> Vec<Vec<f64>> {

        let sketches_per_kmer = (((k_size-t_size) as f64 / t_stride as f64).floor() as usize) + 1;
        let kmer_count = (((seq.len()-k_size) as f64 / k_stride as f64).floor() as usize) + 1;
        let mut tensors = vec![vec![0.0; self.sketch_dim * sketches_per_kmer]; kmer_count];
 
        for i in (0..=seq.len()-k_size).step_by(k_stride) {
            let outer_subseq = &seq[i..i + k_size];
            let window_index = i / k_stride;
            
            println!();
            print!("KMER[{}]: ", window_index);
            for k in 0..outer_subseq.len() {
                if outer_subseq[k].into() == 0  {
                  print!("A");
                } else if outer_subseq[k].into() == 1 {
                  print!("C");
                } else if outer_subseq[k].into() == 2 {
                  print!("G");
                } else if outer_subseq[k].into() == 3 {
                  print!("T");
                }
            } 
            println!();

            for j in (0..=k_size - t_size).step_by(t_stride) {
                let inner_subseq = &outer_subseq[j..j + t_size];

                print!("T        ");
                for k in 0..j {
                  print!(" ");
                }
                for k in 0..inner_subseq.len() {
                    if inner_subseq[k].into() == 0  {
                      print!("A");
                    } else if inner_subseq[k].into() == 1 {
                      print!("C");
                    } else if inner_subseq[k].into() == 2 {
                      print!("G");
                    } else if inner_subseq[k].into() == 3 {
                      print!("T");
                    }
                }
                println!();

                let tensor = self.compute_sketch(inner_subseq);

                // Determine the position to store the tensor in the 2D vector
                let tensor_index = (j / t_stride) * self.sketch_dim;
      
                // Store the tensor in the appropriate position
                for (idx, &value) in tensor.iter().enumerate() {
                    tensors[window_index][tensor_index + idx] = value;
                }
            }

            let mut tid = 0;
            for j in (0..=k_size - t_size).step_by(t_stride) {
                //print!("S: ");
                for k in 0..self.sketch_dim {
                    print!("{:+.4}, ", tensors[window_index][tid]);
                    tid = tid + 1;
                }
                println!();
            }
            
        }
        tensors
    }

    pub fn compute_slide_sketch_1d(
        &self,
        seq: &[SeqType],
        k_size: usize,
        k_stride: usize,
        t_size: usize,
        t_stride: usize,
    ) -> Vec<f32> {

        let sketches_per_kmer = (((k_size-t_size) as f64 / t_stride as f64).floor() as usize) + 1;
        let kmer_count = (((seq.len()-k_size) as f64 / k_stride as f64).floor() as usize) + 1;

        let mut tensors = vec![0.0; self.sketch_dim * sketches_per_kmer * kmer_count];
 
        let mut tensor_index = 0;
        for i in (0..=seq.len()-k_size).step_by(k_stride) {
            let outer_subseq = &seq[i..i + k_size];

            //println!();
            //print!("KMER[{}]: ", i / k_stride);
            //for k in 0..outer_subseq.len() {
            //    if outer_subseq[k].into() == 0  {
            //      print!("A");
            //    } else if outer_subseq[k].into() == 1 {
            //      print!("C");
            //    } else if outer_subseq[k].into() == 2 {
            //      print!("G");
            //    } else if outer_subseq[k].into() == 3 {
            //      print!("T");
            //    }
           // }
            //println!();


            for j in (0..=k_size - t_size).step_by(t_stride) {
                let inner_subseq = &outer_subseq[j..j + t_size];

                //print!("Subsequence: ");
                //for k in 0..inner_subseq.len() {
                //    if inner_subseq[k].into() == 0  {
                //      print!("A");
                //    } else if inner_subseq[k].into() == 1 {
                //      print!("C");
                //    } else if inner_subseq[k].into() == 2 {
                //      print!("G");
                //    } else if inner_subseq[k].into() == 3 {
                //      print!("T");
                //    }
               // }
               // println!();

                let tensor = self.compute_sketch(inner_subseq);

                //println!("Tensor: {:?}", tensor);

                // Store the tensor in the appropriate position
                for (idx, &value) in tensor.iter().enumerate() {
                    tensors[tensor_index + idx] = value as f32;
                }
                tensor_index += self.sketch_dim;
            }
        }
        tensors
    }

    pub fn set_hashes_for_testing(&mut self, h: Vec<Vec<SeqType>>, s: Vec<Vec<bool>>) {
        self.hashes = h;
        self.signs = s;
    }

    pub fn l2_dist(a: &[f64], b: &[f64]) -> f64 {
        a.iter().zip(b.iter())
            .map(|(x, y)| (x - y) * (x - y))
            .sum::<f64>()
            .sqrt()
    }

}



#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    // Compare implementation with output of graph_align
    fn test_compute_sketch() {
        let alphabet_size = 4;
        let sketch_dim = 14;
        let subsequence_len = 6;
        let seed = 0;

        let mut tensor: TensorSketch<u8> = TensorSketch::new(alphabet_size, sketch_dim, subsequence_len, seed);

        // Fixed pairwise hashes (generated randomly by graph_align)
        let h = vec![
            vec![7, 10, 8, 7],
            vec![5, 9, 6, 12],
            vec![13, 5, 11, 7],
            vec![7, 12, 0, 1],
            vec![0, 11, 10, 12],
            vec![13, 11, 6, 10],
        ];
        let s = vec![
            vec![true, true, true, true],
            vec![true, false, false, false],
            vec![false, false, true, false],
            vec![false, true, false, true],
            vec![false, true, false, true],
            vec![false, true, true, true],
        ];

        tensor.set_hashes_for_testing(h, s);

        // Example for which we have output from graph_align
        let sequence = vec![1, 3, 3, 0, 3, 1, 0, 0, 3, 0, 1, 3, 2, 2, 1, 1];
        let sketch = tensor.compute_sketch(&sequence);

        // graph_align output
        let expected_sketch = vec![
            0.0161088911088911,
           -0.026598401598401596,
            0.014110889110889105,
           -0.007742257742257744,
           -0.02122877122877123,
            0.026598401598401596,
           -0.016608391608391615,
           -0.009990009990009985,
            0.021728271728271728,
           -0.029220779220779213,
            0.008991008991008995,
            0.019480519480519477,
           -0.04832667332667335,
            0.022977022977022983
        ];

        for (a, b) in sketch.iter().zip(expected_sketch.iter()) {
            assert!((a - b).abs() < f64::EPSILON);
        }
    }

    #[test]
    // Compare implementation with simple example that we can validate by hand
    fn test_compute_paper_sketch() {
        let alphabet_size = 4;
        let sketch_dim = 3;
        let subsequence_len = 2;
        let seed = 0;

        let mut tensor: TensorSketch<u8> = TensorSketch::new(alphabet_size, sketch_dim, subsequence_len, seed);

        // Pairwise hash functions
        let h = vec![
            vec![2, 1, 0, 2],
            vec![0, 0, 1, 2],
        ];
        let s = vec![
            vec![true, true, false, true],
            vec![true, false, false, false],
        ];

        tensor.set_hashes_for_testing(h, s);

        // Alphabet encoding: ACGT = 0, 1, 2, 3
        // Example string: GACGC
        let sequence = vec![2, 0, 1, 2, 1];
        let sketch = tensor.compute_sketch(&sequence);

        let expected_sketch = vec![
            0.1,
            0.0,
            -0.3,
        ];

        //for (a, b) in sketch.iter().zip(expected_sketch.iter()) {
        //    print!("a: {}, b: {}\n", a, b);
        //}

        for (a, b) in sketch.iter().zip(expected_sketch.iter()) {
            assert!((a - b).abs() < f64::EPSILON);
        }

    }

}

pub mod hyperloglog;
pub mod minhash;
pub mod nodegraph;

use std::any::Any;

use dyn_clonable::*;

use crate::encodings::HashFunctions;
use crate::encodings::{aa_to_dayhoff, aa_to_hp, revcomp, to_aa, VALID};
use crate::Error;
use crate::HashIntoType;

//https://github.com/dtolnay/typetag
#[typetag::serde(tag = "type")]
#[clonable]
pub trait Sketch: Clone {
    fn size(&self) -> usize;
    fn to_vec(&self) -> Vec<u64>;
    fn ksize(&self) -> usize;
    fn check_compatible(&self, other: &Box<dyn Sketch>) -> Result<(), Error>;
    fn seed(&self) -> u64;

    fn hash_function(&self) -> HashFunctions;

    fn add_hash(&mut self, hash: HashIntoType);

    fn add_sequence(&mut self, seq: &[u8], force: bool) -> Result<(), Error> {
        let ready_hashes = SeqToHashes::new(
            seq,
            self.ksize(),
            force,
            false,
            self.hash_function(),
            self.seed(),
        );

        for hash_value in ready_hashes {
            match hash_value {
                Ok(0) => continue,
                Ok(x) => self.add_hash(x),
                Err(err) => return Err(err),
            }
        }

        // Should be always ok
        Ok(())
    }

    fn add_protein(&mut self, seq: &[u8]) -> Result<(), Error> {
        let ready_hashes = SeqToHashes::new(
            seq,
            self.ksize(),
            false,
            true,
            self.hash_function(),
            self.seed(),
        );

        for hash_value in ready_hashes {
            match hash_value {
                Ok(0) => continue,
                Ok(x) => self.add_hash(x),
                Err(err) => return Err(err),
            }
        }

        // Should be always ok
        Ok(())
    }

    fn similarity(&self, other: &Box<dyn Sketch>) -> Result<f64, Error>;
    fn containment(&self, other: &Box<dyn Sketch>) -> Result<f64, Error>;
    fn md5sum(&self) -> String;

    fn as_any(&self) -> &dyn Any;
}

impl std::fmt::Debug for Box<dyn Sketch> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        //write!(f, "{{ x: {}, y: {} }}", self.x, self.y)
        write!(f, "Sketch")
    }
}

impl PartialEq for Box<dyn Sketch> {
    fn eq(&self, other: &Box<dyn Sketch>) -> bool {
        unimplemented!()
    }
}

// Iterator for converting sequence to hashes
pub struct SeqToHashes {
    sequence: Vec<u8>,
    kmer_index: usize,
    k_size: usize,
    max_index: usize,
    force: bool,
    is_protein: bool,
    hash_function: HashFunctions,
    seed: u64,
    hashes_buffer: Vec<u64>,

    dna_configured: bool,
    dna_rc: Vec<u8>,
    dna_ksize: usize,
    dna_len: usize,
    dna_last_position_check: usize,

    prot_configured: bool,
    aa_seq: Vec<u8>,
}

impl SeqToHashes {
    pub fn new(
        seq: &[u8],
        k_size: usize,
        force: bool,
        is_protein: bool,
        hash_function: HashFunctions,
        seed: u64,
    ) -> SeqToHashes {
        let mut ksize: usize = k_size;

        // Divide the kmer size by 3 if protein
        if is_protein || !hash_function.dna() {
            ksize = k_size / 3;
        }

        // By setting _max_index to 0, the iterator will return None and exit
        let _max_index: usize;
        if seq.len() >= ksize {
            _max_index = seq.len() - ksize + 1;
        } else {
            _max_index = 0;
        }

        SeqToHashes {
            // Here we convert the sequence to upper case
            sequence: seq.to_ascii_uppercase(),
            k_size: ksize,
            kmer_index: 0,
            max_index: _max_index as usize,
            force,
            is_protein,
            hash_function,
            seed,
            hashes_buffer: Vec::with_capacity(1000),
            dna_configured: false,
            dna_rc: Vec::with_capacity(1000),
            dna_ksize: 0,
            dna_len: 0,
            dna_last_position_check: 0,
            prot_configured: false,
            aa_seq: Vec::new(),
        }
    }
}

impl Iterator for SeqToHashes {
    type Item = Result<u64, Error>;

    fn next(&mut self) -> Option<Self::Item> {
        // TODO: Remove the hashes buffer
        // Priority for flushing the hashes buffer

        if (self.kmer_index < self.max_index) || !self.hashes_buffer.is_empty() {
            // Processing DNA or Translated DNA
            if !self.is_protein {
                // Setting the parameters only in the first iteration
                if !self.dna_configured {
                    self.dna_ksize = self.k_size as usize;
                    self.dna_len = self.sequence.len();
                    if self.dna_len < self.dna_ksize
                        || (!self.hash_function.dna() && self.dna_len < self.k_size * 3)
                    {
                        return None;
                    }
                    // pre-calculate the reverse complement for the full sequence...
                    self.dna_rc = revcomp(&self.sequence);
                    self.dna_configured = true;
                }

                // Processing DNA
                if self.hash_function.dna() {
                    let kmer = &self.sequence[self.kmer_index..self.kmer_index + self.dna_ksize];

                    for j in std::cmp::max(self.kmer_index, self.dna_last_position_check)
                        ..self.kmer_index + self.dna_ksize
                    {
                        if !VALID[self.sequence[j] as usize] {
                            if !self.force {
                                return Some(Err(Error::InvalidDNA {
                                    message: String::from_utf8(kmer.to_vec()).unwrap(),
                                }));
                            } else {
                                self.kmer_index += 1;
                                // Move the iterator to the next step
                                return Some(Ok(0));
                            }
                        }
                        self.dna_last_position_check += 1;
                    }

                    // ... and then while moving the k-mer window forward for the sequence
                    // we move another window backwards for the RC.
                    //   For a ksize = 3, and a sequence AGTCGT (len = 6):
                    //                   +-+---------+---------------+-------+
                    //   seq      RC     |i|i + ksize|len - ksize - i|len - i|
                    //  AGTCGT   ACGACT  +-+---------+---------------+-------+
                    //  +->         +->  |0|    2    |       3       |   6   |
                    //   +->       +->   |1|    3    |       2       |   5   |
                    //    +->     +->    |2|    4    |       1       |   4   |
                    //     +->   +->     |3|    5    |       0       |   3   |
                    //                   +-+---------+---------------+-------+
                    // (leaving this table here because I had to draw to
                    //  get the indices correctly)

                    let krc = &self.dna_rc[self.dna_len - self.dna_ksize - self.kmer_index
                        ..self.dna_len - self.kmer_index];
                    let hash = crate::_hash_murmur(std::cmp::min(kmer, krc), self.seed);
                    self.kmer_index += 1;
                    Some(Ok(hash))
                } else if self.hashes_buffer.is_empty() {
                    // Processing protein by translating DNA
                    // TODO: make it a real iterator not a buffer

                    // Three frames
                    for i in 0..3 {
                        let substr: Vec<u8> = self
                            .sequence
                            .iter()
                            .cloned()
                            .skip(i)
                            .take(self.sequence.len() - i)
                            .collect();

                        let aa = to_aa(
                            &substr,
                            self.hash_function.dayhoff(),
                            self.hash_function.hp(),
                        )
                        .unwrap();

                        aa.windows(self.k_size as usize).for_each(|n| {
                            let hash = crate::_hash_murmur(n, self.seed);
                            self.hashes_buffer.push(hash);
                        });

                        let rc_substr: Vec<u8> = self
                            .dna_rc
                            .iter()
                            .cloned()
                            .skip(i)
                            .take(self.dna_rc.len() - i)
                            .collect();
                        let aa_rc = to_aa(
                            &rc_substr,
                            self.hash_function.dayhoff(),
                            self.hash_function.hp(),
                        )
                        .unwrap();

                        aa_rc.windows(self.k_size as usize).for_each(|n| {
                            let hash = crate::_hash_murmur(n, self.seed);
                            self.hashes_buffer.push(hash);
                        });
                    }
                    self.kmer_index = self.max_index;
                    Some(Ok(self.hashes_buffer.remove(0)))
                } else {
                    let first_element: u64 = self.hashes_buffer.remove(0);
                    Some(Ok(first_element))
                }
            } else {
                // Processing protein
                // The kmer size is already divided by 3

                if self.hash_function.protein() {
                    let aa_kmer = &self.sequence[self.kmer_index..self.kmer_index + self.k_size];
                    let hash = crate::_hash_murmur(aa_kmer, self.seed);
                    self.kmer_index += 1;
                    Some(Ok(hash))
                } else {
                    if !self.prot_configured {
                        self.aa_seq = match self.hash_function {
                            HashFunctions::murmur64_dayhoff => {
                                self.sequence.iter().cloned().map(aa_to_dayhoff).collect()
                            }
                            HashFunctions::murmur64_hp => {
                                self.sequence.iter().cloned().map(aa_to_hp).collect()
                            }
                            invalid => {
                                return Some(Err(Error::InvalidHashFunction {
                                    function: format!("{}", invalid),
                                }));
                            }
                        };
                    }

                    let aa_kmer = &self.aa_seq[self.kmer_index..self.kmer_index + self.k_size];
                    let hash = crate::_hash_murmur(aa_kmer, self.seed);
                    self.kmer_index += 1;
                    Some(Ok(hash))
                }
            }
        } else {
            // End the iterator
            None
        }
    }
}

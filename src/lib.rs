#![allow(unused_variables, dead_code)]
/* Class for tunnelled BWT
 * implemented: 
 *  - construction of BWT
 *  - tunneling BWT using the heurisitcs
 *  - reconstructing the TBWT
 *  - reading TBWT from files
 *  - tests for employed functions
 */

use bio::data_structures::rank_select::RankSelect;
use bv::BitsMut;
use bv::BitVec;
use::std::cmp;
use std::collections::{HashMap, VecDeque};
use std::fs::File;
use std::io::{self, Read};
use simple_sds::bit_vector::BitVector;
use simple_sds::raw_vector::{RawVector, AccessRaw, PushRaw};
use simple_sds::ops::{BitVec as bitvec, Rank, Select};




pub struct TBWT {
    bwt: Vec<u8>,  // for access to bwt
    o: BitVector, // outgoing nodes
    i: BitVector, // incoming nodes

    o_raw: RawVector, //for tunneling
    i_raw: RawVector, //for tunneling

    l: Vec<RankSelect>, // for ranks on bwt used in lf mapping
    f: Vec<usize>, // f[c] is the number of letters smaller than c, or the first position of letter
                   // c in array f

    char_map: Vec<u8>, //mapping original characters to characters 0, 1, 2, ...
}


impl TBWT {

    /* Function reads a bitvector from file
     * Result is raw version of the bitvector and Bitvector from the raw vesrion
     */
    fn bitvector_from_file (filename: &str, size: usize) -> io::Result<(BitVector, RawVector)> {
        let mut file = File::open(filename).expect("Unable to open file ");
        let mut vector = Vec::new();
        file.read_to_end(&mut vector).expect("Unable to read file");

        let mut bv_raw = RawVector::with_len(size+1, false);

        let mut index = 0;
        for byte in vector {
            for i in (0..8).rev() {
                if byte & (1 << i) != 0 && index < size {
                    bv_raw.set_bit(index, true);
                }
                index += 1;
            }
        }

        bv_raw.set_bit(size, true);

        let bv: BitVector = BitVector::from(bv_raw.clone());

        Ok((bv, bv_raw))
    }

    /* Function reads a tunneled Wheeler graph from three files,
     * containing tunneled l vector, din vector and dout vector
     * Results instance of TBWT
     */
    pub fn from_binary(lfile_name: &str, ofile_name: &str, ifile_name: &str) -> io::Result<TBWT> {
        let mut lfile = File::open(lfile_name).expect("Unable to open file bwt");
        let mut bwt = Vec::new();
        lfile.read_to_end(&mut bwt).expect("Unable to read file bwt");

        let (din, in_raw) = Self::bitvector_from_file(ifile_name, bwt.len()).unwrap();
        let (dout, out_raw) = Self::bitvector_from_file(ofile_name, bwt.len()).unwrap();

        Ok(TBWT::new(&bwt, dout, din, out_raw, in_raw))
    }

    /* Function  construcitng instance of class TBWT from a simple bwt
     * sigma is the size of the alphabet
     * used only for data from alphabet {0, 1, , ..., sigma-1}
     */
    pub fn from_bwt(bwt: &[u8], sigma: usize) -> TBWT {

        // construct l for ranks an selcets on the bwt
        let mut wm: Vec<RankSelect> = Vec::new();
        for i in 0..sigma {
            let mut bitvec: BitVec<u8> = BitVec::new_fill(false, bwt.len() as u64);
            for j in 0..bwt.len() {
                if bwt[j] == i as u8 {
                    bitvec.set_bit(j as u64, true);
                }
            }
            wm.push(RankSelect::new(bitvec, 1));
        }

        // count the letter occurances to construct the array f
        let mut char_occ: Vec<usize> = vec![0; sigma];
        for c in bwt {
            char_occ[*c as usize] += 1;
        }
        let mut f: Vec<usize> = vec![0; sigma];
        for i in 1..sigma {
            f[i] = f[i - 1] + char_occ[i - 1];
        }

        let o_raw = RawVector::with_len(bwt.len()+1, true);
        let mut o = BitVector::from(o_raw.clone());
        o.enable_rank(); o.enable_select();

        let i_raw = RawVector::with_len(bwt.len()+1, true);
        let mut i = BitVector::from(i_raw.clone());
        i.enable_rank(); i.enable_select();

        TBWT {
            bwt: bwt.to_vec(),
            l: wm,
            i,
            o,
            o_raw,
            i_raw,
            f,
            char_map: (0..sigma as u8).collect(),
        }
    }

    /* Function  construcitng instance of class TBWT from a simple bwt
     * used for general data, enriched with the function mapping the characters from any alphabet to alphabet {0, 1, 2, ..., k}
    */
    pub fn from_bwt2(l: &[u8]) -> TBWT {
        //find size of alphabet sigma
        let mut alphabet = HashMap::new();
        for c in l {
            alphabet.entry(c).or_insert(1);
        }
        let sigma: usize = alphabet.keys().len();

        // map the real characters in alphabet to characters {0, 1, ..., sigma-1} and remember the
        // mapping for later reconstruction
        let mut char_map = Vec::new();
        for (c, i) in alphabet {
            char_map.push(*c);
        }
        char_map.sort();

        //get necessary info from bwt
        let mut bwt: Vec<u8> = vec![0; l.len()];
        for i in 0..l.len() {
            bwt[i] = char_map.iter().position(|&x| x == l[i]).unwrap() as u8;
        }

        // construct l for ranks an selcets on the bwt
        let mut wm: Vec<RankSelect> = Vec::new();
        for i in 0..sigma {
            let mut bitvec: BitVec<u8> = BitVec::new_fill(false, bwt.len() as u64);
            for j in 0..bwt.len() {
                if bwt[j] == i as u8 {
                    bitvec.set_bit(j as u64, true);
                }
            }
            wm.push(RankSelect::new(bitvec, 1));
        }

        // count the letter occurances to construct the array f
        let mut char_occ: Vec<usize> = vec![0; sigma];
        for c in bwt.clone() {
            char_occ[c as usize] += 1;
        }
        let mut f: Vec<usize> = vec![0; sigma];
        for i in 1..sigma {
            f[i] = f[i - 1] + char_occ[i - 1];
        }

        let o_raw = RawVector::with_len(bwt.len()+1, true);
        let mut o = BitVector::from(o_raw.clone());
        o.enable_rank(); o.enable_select();

        let i_raw = RawVector::with_len(bwt.len()+1, true);
        let mut i = BitVector::from(i_raw.clone());
        i.enable_rank(); i.enable_select();

        TBWT {
            bwt,
            l: wm,
            i,
            o,
            i_raw,
            o_raw,
            f,
            char_map,
        }
    }

    /* Function to create a new TBWT instant from a tunneled bwt
     * input - tunneled bwt l with corresponding bitvectors din i an dout o an their raw version
     */
    pub fn new(l: &[u8], mut o: BitVector, mut i: BitVector, o_raw: RawVector, i_raw: RawVector) -> TBWT {
        let mut wm: Vec<RankSelect> = Vec::new();

        //find size of alphabet sigma
        let mut alphabet = HashMap::new();
        for c in l {
            alphabet.entry(c).or_insert(1);
        }
        let sigma: usize = alphabet.keys().len();

        // map the real characters in alphabet to characters {0, 1, ..., sigma-1} and remember the
        // mapping for later reconstruction
        let mut char_map = Vec::new();
        for (c, i) in alphabet {
            char_map.push(*c);
        }
        char_map.sort();

        //get necessary info from bwt
        let mut bwt: Vec<u8> = vec![0; l.len()];
        for i in 0..l.len() {
            bwt[i] = char_map.iter().position(|&x| x == l[i]).unwrap() as u8;
        }

        // construct l for ranks an selcets on the bwt
        for i in 0..sigma {
            let mut bitvec: BitVec<u8> = BitVec::new_fill(false, bwt.len() as u64);
            for j in 0..bwt.len() {
                if bwt[j] == i as u8 {
                    bitvec.set_bit(j as u64, true);
                }
            }
            wm.push(RankSelect::new(bitvec, 1));
        }

        // count the letter occurances to construct the array f
        let mut char_occ: Vec<usize> = vec![0; sigma];
        for c in bwt.clone() {
            char_occ[c as usize] += 1;
        }
        let mut f: Vec<usize> = vec![0; sigma];
        for i in 1..sigma {
            f[i] = f[i - 1] + char_occ[i - 1];
        }

        let k: usize = o.len().pow(2) as usize / 32;

        o.enable_select();
        i.enable_rank();
        i.enable_select();

        TBWT {
            bwt,
            l: wm,
            o,
            i,
            o_raw,
            i_raw,
            f,
            char_map,
        }
    }

    /* Function to determine if any tunnel ends at the position of the node
     */
    fn tunnel_ends(&self, node: usize) -> bool { 
        // find where the bitrun for the node starts
        let l_pos = self.o.select(node).unwrap();

        // if second position is 0, it means there is > 2 edges outgoing and tunnel ends 
        return !self.o.get(l_pos + 1);
    }

    /* Fucntion to determine if any tunnel starts at the position of the node
     */ 
    fn tunnel_starts(&self, node: usize) -> bool { 
        // same logic as tunnel ends, only on incoming bitvector
        let f_pos = self.i.select(node).unwrap();

        // if the second position is 0, it means there starts a tunnel
        return !self.i.get(f_pos + 1);
    }

    /* Function to navigate through tunneled bwt
     * results the previous node of the input node in the graph, or previous character in the initial string
     * stack serves for remembering offsets of tunnels
     */
    fn backward_step(&self, node: usize, stack: &mut Vec<usize>) -> usize {
        // backwardstep in Baier does not use node, but edge

        let node_o_bits_start = self.o.select(node).unwrap() as usize;

        let mut offset = 0;
        if self.tunnel_ends(node) && !stack.is_empty() { offset = stack.pop().unwrap(); }

        let l_pos = node_o_bits_start + offset;
        let letter = self.bwt[l_pos] as usize;

        // this rank is inclusive, we want exclusive, thus -1;
        let letter_rank = self.l[letter].rank(l_pos as u64).unwrap() - 1;
        let f_pos = self.f[letter] as u64 + letter_rank;

        // this rank is inclusive, we want exclusive, thus -1;
        let new_node = self.i.rank(f_pos as usize +1) -1;

        if self.tunnel_starts(new_node) {
            let new_node_i_bits_start = self.i.select(new_node).unwrap();
            let new_offset = f_pos as usize - new_node_i_bits_start;
            stack.push(new_offset);
        }

        return new_node as usize;
    }

    /* Function to determine the letter at position of the node in f
     */
    fn incoming_letter(&self, node: usize) -> u8 {
        let f_pos =  self.i.select(node).unwrap();
        for i in (0..self.f.len()).rev() {
            if self.f[i] <= f_pos {
                return i as u8;
            }
        }
        return 0;
    }

    /* Function that reconstruts the original string from the TBWT
     * uses function backward_step an incomming_letter
     * returns the reconstructed original string
     */
    pub fn reconstruct(&self) -> Vec<u8> {
        let mut node = 0;
        let mut offsets = vec![0];

        let mut result = Vec::with_capacity(self.bwt.len());

        loop {
            result.push(self.char_map[self.incoming_letter(node) as usize]);
            node = self.backward_step(node, &mut offsets);
            if node == 0 { break; }
        }

        return result;
    }

    /* Function that checks if the reconstructed string from the TBWT is equal to the expected
     * string
     * using functions backward_step an incomming_letter
     */
    fn reconstruct_with_expected(&self, expected: Vec<u8>) {
        let mut node = 0;
        let mut offsets = vec![0];

        let mut counter = expected.len()-1;

        loop {
            let x = self.char_map[self.incoming_letter(node) as usize];
            assert_eq!(expected[counter], x);
 
            if counter == 0 {break;}

            node = self.backward_step(node, &mut offsets);
            if node == 0 { break; }

            counter = counter-1;
        }
    }

    /* LF-mapping
     * returns the position in f of the letter at pos in l
     */
    fn lf(&self, pos: usize) -> usize {
        let c = self.bwt[pos] as usize;
        let rank = self.l[c].rank(pos as u64).unwrap() as usize - 1;
        let less = self.f[c];
        return less + rank;
    }

    /* inverse LF-mapping
     * returns the position in l of the letter at pos in f
     */
    fn inverse_lf (&self, pos:usize) -> usize {
        let mut c = 0;
        for i in (0..self.f.len()).rev() {
            if self.f[i] as u64 <= pos as u64 {
                c = i;
                break;
            }
        }

        return self.l[c].select((pos-self.f[c]) as u64 +1).unwrap() as usize;
    }

    /* Function to update l and f array of the freshly tunneled TBWT 
     * according to the tunneled bwt, i an o bitvectors 
     */
    fn update_tunneled(&mut self, bwt: Vec<u8>, o: BitVector, i: BitVector) {
        let mut wm: Vec<RankSelect> = Vec::new();

        //find size of alphabet
        let sigma = self.f.len();

        // update l
        for i in 0..sigma {
            let mut bitvec: BitVec<u8> = BitVec::new_fill(false, bwt.len() as u64);
            for j in 0..bwt.len() {
                if bwt[j] == i as u8 {
                    bitvec.set_bit(j as u64, true);
                }
            }
            wm.push(RankSelect::new(bitvec, 1));
        }

        // update f
        let mut char_occ: Vec<usize> = vec![0; sigma];
        for c in bwt.clone() {
            char_occ[c as usize] += 1;
        }
        let mut f: Vec<usize> = vec![0; sigma];
        for i in 1..sigma {
            f[i] = f[i - 1] + char_occ[i - 1];
        }

        self.l = wm;
        self.f = f;
    }

    /* Function to mark the tunnel of the block
     * argument - block as a triplet of width, start and end position in bwt
     * marks the tunnel in raw version of i and o array, as only they can be changed
     */
    pub fn mark_tunnel(&mut self, block: (usize, usize, usize)) {
        let w = block.0; // block width
        let i = block.1; // starting bwt position
        let j = block.2; // ending bwt position

        for k in i + 1..=j {
            let mut pos = k;
            for _ in 0..w - 1 {
                self.i_raw.set_bit(pos, false);
                pos = self.lf(pos);
                self.o_raw.set_bit(pos, false);
            }
        }
    }

    /* Function to tunnel self according to the marked tunneles in arrays i and o
     */
    pub fn tunnel(&mut self) {
        // L = 5$1422334    -> L:	5 $ 1 4 2 3   4
        // O = 1110101111   -> O:	1 1 1 0 1 1   1 1
        // I = 1111101011   -> I:	1 1 1   1 1 0 1 1
        // F = $12233445    -> F:   $ 1 2   3 4 4 5
        //                     C:  [0,1,2,  3,4,  6]

        let mut last = Vec::with_capacity(self.bwt.len());
        let mut outgoing = RawVector::new();
        let mut incoming = RawVector::new();

        let n = self.i.len() - 1;
        for i in 0..n {
            let ibit = self.i_raw.bit(i);
            let obit = self.o_raw.bit(i);
            match (ibit, obit) {
                (true, true) => {
                    // if 1 1 write 1 1 and letter  - no tunnel
                    incoming.push_bit(true);
                    outgoing.push_bit(true);
                    last.push(self.bwt[i as usize]);
                }
                (true, false) => {
                    // if 1 0 write - 0 and letter  - beginning of a tunnel
                    outgoing.push_bit(false);
                    last.push(self.bwt[i as usize]);
                }
                (false, true) => {
                    // if 0 1 write 0 -   - end of a tunnel
                    incoming.push_bit(false);
                } 
                (false, false) => {} // if 0 0 skip   -  in a tunnel
            }
        }
        incoming.push_bit(true);
        outgoing.push_bit(true);

        // update bwt, i and o vectors
        self.bwt = last;
        self.i = BitVector::from(incoming);
        self.i.enable_rank(); self.i.enable_select();
        self.o = BitVector::from(outgoing);
        self.o.enable_rank(); self.o.enable_select();
 
        self.update_tunneled(self.bwt.clone(), self.o.clone(), self.i.clone());
    }

    /* Function to find maximal blocks
     * return the set of maximal blocks as triplets (width, start position, end position)
     */
    fn find_maximal_blocks(&self) -> Vec<(usize, usize, usize)>{
        // for fast computation, memorize the inverse lf in an arra
        let mut inverse_lf = vec![0; self.bwt.len()];
        for i in 0..self.bwt.len() {
            inverse_lf[self.lf(i)] = i;
        }

        // compute the special lcs array - longest common suffix array
        let mut lcs = vec![0; self.bwt.len()];

        lcs[1] = 0;
        let mut j = inverse_lf[1];
        let mut l = 0;

        for i in 2..self.bwt.len() {
            if j > 0 && self.bwt[j] == self.bwt[j-1] {
                l += 1;
            } else {
                l = 0;
            }
            if j > 0 && self.bwt[inverse_lf[j]] == self.bwt[inverse_lf[j-1]] {
                lcs[j] = l+1;
            } else {
                lcs[j] = 0;
            }

            j = inverse_lf[j];
        }

        // Find all maximal blocks
        // Store beginnings of potential blocks as (start position, width) on stack
        let mut stack: VecDeque<(usize, usize)> = VecDeque::new();
        stack.push_back((1, 0));

        let mut all_blocks: Vec<(usize, usize, usize)> = Vec::new();

        for i in 1..lcs.len() {
            let mut start = stack.back().cloned().unwrap();

            // End of the block of width of size start.1
            while start.1 > lcs[i] {
                stack.pop_back();

                //if big enough and right-maximal, report the block
                if i- start.0 > 1 { // the block is high enough
                    let mut max = false;
                    for j in start.0..i {
                        if self.bwt[inverse_lf[inverse_lf[j]]] != self.bwt[inverse_lf[inverse_lf[start.0]]] {
                            max = true;
                            break;
                        }
                    }

                    //check if the block expanded one column to the right is also a block or not,
                    //in order to ensure the right maximality
                    //the block is right maximal if characters to the right do not
                    //form a consecutive part in the bwt with characters two times right being
                    //equally labeled (then there woul be a bigger block)
                    if !(!max 
                        && inverse_lf[i-1] > inverse_lf[start.0]
                        && inverse_lf[i-1] - inverse_lf[start.0] == i-1-start.0) {

                        if start.1 > 1 {
                            all_blocks.push((start.1, start.0, i-1));
                        }
                    }
                }

                // Assurance of the height-maximality of blocks
                if lcs[i] >= 1 && (stack.len() <= 1 || stack.back().unwrap().1 < lcs[i]) {
                    stack.push_back((start.0, lcs[i]));
                }
                start = stack.back().cloned().unwrap();
            }

            // Store a possible start of a block
            if start.1 < lcs[i] {
                stack.push_back((i - 1, lcs[i]));
            }
        }

        // check for blocks ending at the end of the bwt - assurance of reporting all blocks
        let bwt_size = self.bwt.len();
        while let Some(start) = stack.pop_back() {
            if bwt_size - start.0 > 1 {
                let mut max = false;
                for j in (start.0+1)..bwt_size {
                    if self.bwt[inverse_lf[inverse_lf[j]]] != self.bwt[inverse_lf[inverse_lf[start.0]]]{
                        max = true;
                        break;
                    }
                }

                //assurance of the right maximality
                if !(!max
                    && inverse_lf[bwt_size-1] > inverse_lf[start.0]
                    && inverse_lf[bwt_size-1] - inverse_lf[start.0] == bwt_size-1-start.0) {
                        if start.1 > 1 {
                            all_blocks.push((start.1, start.0, bwt_size-1));
                        }
                }
            }

        }

        return all_blocks;
    }


    /* Function to merge a set of overlapping intervals into a minimal set of 
     * non overlapping intervals that cover the original set
     */
    fn merge_intervals(&self, mut intervals: Vec<(usize, usize)>) -> Vec<(usize, usize)> {
        intervals.sort_by(|a, b| if a.0 == b.0 {
                (b.1).cmp(&a.1)
            } else {
                (a.0).cmp(&b.0)
        });

        let mut start = intervals[0].0;
        let mut end = intervals[0].1;
        let mut result = Vec::new();

        for i in intervals {
            if i.0 <= end {
                continue;
            }

            result.push((start, end));
            start = i.0;
            end = i.1;
        }

        result.push((start, end));

        return result;
    }

    /* Function that finds inner blocks of width one and outputs them
     * first the intervals in l an f occupied by the blocks are computed
     * next they are merged
     * finally, their overlaps are computed and returned as blocks of width one
     */
    fn one_column_overlappings(&self, mut blocks: Vec<(usize, usize, usize)>) -> Vec<(usize, usize, usize)> {
        let mut f_blocks = Vec::new(); //occupied intervals in f
        let mut l_blocks = Vec::new(); //occupied intervals in l

        // remember occupied interval of each block
        for i in 0..blocks.len() {
            l_blocks.push((blocks[i].1, blocks[i].2));
            let mut start = blocks[i].1;
            for _ in 0..blocks[i].0-1 {
                start = self.lf(start);
            }
            f_blocks.push((start, start +blocks[i].2-blocks[i].1));
        }

        //merge the intervals
        let f_interval = self.merge_intervals(f_blocks);
        let l_interval = self.merge_intervals(l_blocks);

        //compute the overlaps and save them as blocks of width one
        let mut fi = 0;
        let mut li = 0;
        while li < l_interval.len() {
            if fi >= f_interval.len() {
                fi -= 1;
                if l_interval[li].0 > f_interval[fi].1 || li == l_interval.len()-1{
                    break;
                }
                li += 1;
                continue;
            }

            if f_interval[fi].0 > l_interval[li].1 {
                li += 1;
                continue;
            }
            if l_interval[li].0 > f_interval[fi].1 {
                fi += 1;
                continue;
            }

            blocks.push((1, cmp::min(l_interval[li].0, f_interval[fi].0), l_interval[li].1/*cmp::max(l_interval[li].1, f_interval[fi].1)*/));
            li += 1;
        }

        return blocks;
    }

    /* Function to resolve right-aligned collision via vertical division
     * input: set of blocks
     * output: set of blocks free of right-aligned (and corner) collisions
     */
    fn vertical_division(&self, mut blocks: Vec<(usize, usize, usize)>) -> Vec<(usize, usize, usize)> {
        // sort the blocks
        blocks.sort_by(|a, b| if a.1 == b.1 && a.2 == b.2 {
                    (a.0).cmp(&b.0)
                } else if a.1 == b.1 {
                    (b.2).cmp(&a.2)
                } else {
                    (a.1).cmp(&b.1)
                });

        let mut cycle = true; //variable defining whether there are still right-aligned collisions
                              //present in the set

        while cycle {
            cycle = false;
            let mut stack: VecDeque<(usize, usize, usize)> = VecDeque::new(); // path to the root
                                                                              // of the block tree
            let mut tmp = Vec::new(); //for newly created blocks by the division
            let mut i = 0;
 
            while i < blocks.len() {
                // remove false blocks in the path to the root
                while !stack.is_empty() && stack.back().cloned().unwrap().2 < blocks[i].1 {
                    stack.pop_back();
                }

                if !stack.is_empty() {
                    cycle = true;

                    let parent = stack.back().cloned().unwrap();

                    // if the parent is wider, merge the block and the parent
                    if parent.0 >= blocks[i].0 {
                        if blocks[i].0 != parent.0 || parent.1 < blocks[i].1 || parent.2 > blocks[i].2 {
                            tmp.push((blocks[i].0, parent.1, cmp::max(parent.2, blocks[i].2)));
                        } 
                    } else if parent.0+1 < blocks[i].0 { //create new block that is divided by the
                                                         //with of the parent
                        let mut start = blocks[i].1;
                        for _ in 0..(parent.0) {
                            start = self.lf(start);
                        }
                        tmp.push((blocks[i].0-parent.0, start, start+blocks[i].2-blocks[i].1));
                    }
                    stack.push_back(blocks[i]);

                    blocks.remove(i);
                    i -= 1;
                } else {
                    stack.push_back(blocks[i]);
                }

                i += 1;
            }

            tmp.sort_by(|a, b| if a.1 == b.1 && a.2 == b.2 {
                    (a.0).cmp(&b.0)
                } else if a.1 == b.1 {
                    (b.2).cmp(&a.2)
                } else {
                    (a.1).cmp(&b.1)
                });

            // merge the leftover set of blocks with the sorted set of newly created blocks
            if cycle {
                let mut ind1 = 0;
                let mut ind2 = 0;
                let mut tmp2 = Vec::new();
                for _ in 0..(tmp.len() + blocks.len()) {
                    if ind2 >= tmp.len() || (ind1 < blocks.len() && 
                                             (blocks[ind1].1 < tmp[ind2].1 || 
                                              (blocks[ind1].1 == tmp[ind2].1 && 
                                               (blocks[ind1].2 > tmp[ind2].2 || (blocks[ind1].0 == tmp[ind2].2 && blocks[ind1].0 < tmp[ind2].0))))) {
                        tmp2.push(blocks[ind1]);
                        ind1 += 1;
                    } else {
                        tmp2.push(tmp[ind2]);
                        ind2 += 1;
                    }
                }
                blocks = tmp2;
            }
        }

        return blocks;
    }

    /* Function to remove all self-colliding blocks in the set
     * and blocks of width less than 2
     */
    fn remove_selfcolliding_blocks (&self, all_blocks : Vec<(usize, usize, usize)>) -> Vec<(usize, usize, usize)> {
        //save only non-selfcolliding blocks of with at least two
        let mut blocks: Vec<(usize, usize, usize)> = Vec::new();

        for block in all_blocks {
            if block.0 < 2 {
                continue;
            }
            let mut run_starts = Vec::new();
            let mut start = block.1;

            for i in 0..block.0 {
                run_starts.push(start);
                start = self.lf(start);
            }

            run_starts.sort();

            let mut non_selfcolliding = true;

            // if the first row positoin are closer in the f than is the height, the block is
            // self-colliding
            for i in 1..block.0 {
                if run_starts[i]-run_starts[i-1] <= (block.2-block.1) {
                    non_selfcolliding = false;
                    break;
                }
            }

            if non_selfcolliding {
                blocks.push(block);
            }
        }

        return blocks;
    }

    /* Function to solve the left-aligne collisions via shortening
     * input: set of blocks
     * output: set of blocks free of left-aligned collisions
     */
    fn shortening_blocks(&self, blocks: Vec<(usize, usize, usize)>) -> Vec<(usize, usize, usize)> {

        //align blocks to the other side of the bwt matrix and sort them
        let mut obverse_blocks : Vec<(usize, usize, usize, usize)> = Vec::new();
        let mut ind = 0;
        for block in blocks.clone() {
            let mut start = block.1;
            for i in 0..block.0-1 {
                start = self.lf(start);
            }
            obverse_blocks.push((block.0, start, start + (block.2-block.1), ind));
            ind += 1;
        }
 
        obverse_blocks.sort_by(|a, b| if a.1 == b.1 && a.2 == b.2 {
                    (a.0).cmp(&b.0)
                } else if a.1 == b.1 {
                    (b.2).cmp(&a.2)
                } else {
                    (a.1).cmp(&b.1)
                });

        let mut changed = true;
        while changed {
            let mut blocks_subtract = vec![0; obverse_blocks.len()]; //vector for remembering the
                                                                     //with of the shortening of
                                                                     //each block
            let mut stack: VecDeque<usize> = VecDeque::new(); // remember the path to a leaf
            let mut i = 0;

            while i < obverse_blocks.len() {
                stack.push_back(i);

                let mut j = i+1;
                let mut depth = 0;

                // push blocks until a leaf is reached
                while j < obverse_blocks.len() && obverse_blocks[j].1 <= obverse_blocks[j-1].2 {
                    stack.push_back(j);
                    j += 1;
                }

                // update the size of the shortening of all block on the path, according to the
                // depth
                while !stack.is_empty() && (j >= obverse_blocks.len() || obverse_blocks[j].1 > obverse_blocks[stack.back().cloned().unwrap()].2) {
                    let index = stack.pop_back().unwrap();
                    blocks_subtract[index] = cmp::max(depth, blocks_subtract[index]);
                    depth += 1;
                }
                for k in (0..stack.len()).rev() {
                    blocks_subtract[stack[k]] = cmp::max(depth, blocks_subtract[stack[k]]);
                    depth += 1;
                }

                i += j-i;
            }

            let mut shorten_blocks = Vec::new();
            let mut initial_blocks = Vec::new();

            // shorten the blocks
            i = 0;
            changed = false;
            for block in obverse_blocks.clone() {
                if blocks_subtract[i] > 0 {
                    if block.0 > blocks_subtract[i] + 1 {
                        changed = true;

                        let mut start = block.1;
                        for _ in 0..blocks_subtract[i] {
                            start = self.inverse_lf(start);
                        }
                        shorten_blocks.push((block.0-blocks_subtract[i], start, start + block.2-block.1, block.3));
                    }
                } else {
                    initial_blocks.push(block);
                }
                i += 1;
            }

            //merge the leftover set of blocks with sorted set of the newly created blocks
            if changed {
                shorten_blocks.sort_by(|a, b| if a.1 == b.1 && a.2 == b.2 {
                    (a.0).cmp(&b.0)
                } else if a.1 == b.1 {
                    (b.2).cmp(&a.2)
                } else {
                    (a.1).cmp(&b.1)
                });

                let mut tmp = Vec::new();
                let mut sb_index = 0;
                let mut ib_index = 0;
                for _ in 0..initial_blocks.len() + shorten_blocks.len() {
                    if sb_index >= shorten_blocks.len() || 
                        (ib_index < initial_blocks.len() && 
                         (initial_blocks[ib_index].1 < shorten_blocks[sb_index].1 ||
                          (initial_blocks[ib_index].1 == shorten_blocks[sb_index].1 && initial_blocks[ib_index].2 > shorten_blocks[sb_index].2))) {
                        tmp.push(initial_blocks[ib_index]);
                        ib_index += 1;
                    } else {
                        tmp.push(shorten_blocks[sb_index]);
                        sb_index += 1;
                    }
                }
                obverse_blocks = tmp;
            } else {
                obverse_blocks = initial_blocks;
                break;
            }
        }

        //remap the blocks to the l column
        let mut result = Vec::new();
        for block in obverse_blocks {
            result.push((block.0, blocks[block.3].1, blocks[block.3].2));
        }

        return result;
    }

    /* Function to find a set of blocks with no critical collisions
     * right-aligned collisions are solved with vertical division
     * left-aligned collisions are solved with shortening method
     */
    fn heuristic(&self) -> Vec<(usize, usize, usize)> {
        let mut blocks = self.find_maximal_blocks();

        if blocks.len() <= 1 {
            blocks = self.remove_selfcolliding_blocks(blocks);
            return blocks;
        }

        blocks = self.one_column_overlappings(blocks);

        blocks = self.vertical_division(blocks);
 
        blocks = self.remove_selfcolliding_blocks(blocks);

        blocks = self.shortening_blocks(blocks);

        return blocks;
    }

    /* Function to tunnel self using the heuristic for block choice
     */
    pub fn heuristic_tunnel(&mut self) {
        let blocks = self.heuristic();

        if blocks.len() == 0 {
            return;
        }

        for i in 0..blocks.len() {
            if blocks[i].0 > 1 {  
                self.mark_tunnel((blocks[i].0, blocks[i].1, blocks[i].2));
            }
        }

        self.tunnel();
    }

    /* Print the TBWT
     */
    pub fn print(&self) {
        for c in &self.bwt {
            print!("{}", *c);
        }
        println!();
        for i in 0..self.o.len() {
            print!("{}", self.o.get(i) as u8);
        }
        println!();
        for i in 0..self.i.len() {
            print!("{}", self.i.get(i) as u8);
        }
        println!();
    }
}


//   TESTS   ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#[cfg(test)]
mod tests {
    use super::*;
    use bv::bit_vec;
    use bio::data_structures::bwt::{bwt as bwtransform, self};
    use bio::data_structures::suffix_array::suffix_array;
    use proptest::prelude::*;
    use std::process::Command;
    use std::time::Instant;
    use std::io::Write;

    #[test]
    #[rustfmt::skip]
    fn test_simple_case() {
        // input
        let seq = [1, 2, 3, 4, 2, 3, 4, 5, 0];
        let sa  = suffix_array(&seq);
        let bwt = bwtransform(&seq, &sa);
        let wg  = TBWT::from_bwt(&bwt, 6);

        // expected output
        let expected_nodes   = [0, 8, 7, 5, 3, 6, 4, 2, 1]; // 0, lf(0), lf(lf(0))...
        let expected_offsets = [0, 0, 0, 0, 0, 0, 0, 0, 0];
        let expected_letters = seq.into_iter().rev().collect::<Vec<u8>>();

        // test code
        wg.print();

        let mut node = 0;
        let mut offsets = vec![0];

        for i in 0..expected_nodes.len() {
            let c = wg.incoming_letter(node);

            assert_eq!(expected_letters[i], c);
            assert_eq!(  expected_nodes[i], node);
            assert_eq!(expected_offsets[i], *offsets.last().unwrap());

            node = wg.backward_step(node, &mut offsets);
        }
    }

    #[test]
    #[rustfmt::skip]
    fn test_tunneled_case() {
        // input
        let test_seq = [5, 0, 1, 4, 2, 3, 4];
        let bit_o = vec![true, true, true, false, true, true, true, true];
        let bit_i = vec![true, true, true, true, true, false, true, true];
 
        let mut raw_o = RawVector::with_len(bit_o.len(), true);
        let mut raw_i = RawVector::with_len(bit_i.len(), true);

        for i in 0..bit_o.len() {
            if !bit_o[i] {  raw_o.set_bit(i, false);  }
        }
        for i in 0..bit_i.len() {
            if !bit_i[i] {  raw_i.set_bit(i, false);  }
        }

        let test_o: BitVector = BitVector::from(raw_o.clone());
        let test_i: BitVector = BitVector::from(raw_i.clone());                                           

        let wg: TBWT = TBWT::new(&test_seq, test_o, test_i, raw_o, raw_i);

        // expected output
        let expected_nodes   = [0, 5, 4, 3, 2, 4, 3, 2, 1];
        let expected_offsets = [0, 0, 1, 1, 1, 0, 0, 0, 0];
        let expected_letters = [0, 5, 4, 3, 2, 4, 3, 2, 1];
 
        // test code
        wg.print();

        let mut node = 0;
        let mut offsets = vec![0];

        for i in 0..expected_nodes.len() {
            let c = wg.incoming_letter(node);
            //println!("{} {} {}", c, node, offsets.last().unwrap());

            assert_eq!(expected_letters[i], c);
            assert_eq!(  expected_nodes[i], node);
            assert_eq!(expected_offsets[i], *offsets.last().unwrap());

            node = wg.backward_step(node, &mut offsets);
        }
    }

    #[test]
    fn test_reconstruction1() {
        let test_seq = [6, 0, 1, 2, 3, 4, 5];
        let bit_o = vec![true, true, true, true, true, true, true, true];
        let bit_i = vec![true, true, true, true, true, true, true, true];
  
        let mut raw_o = RawVector::with_len(bit_o.len(), true);
        let mut raw_i = RawVector::with_len(bit_i.len(), true);

        for i in 0..bit_o.len() {
            if !bit_o[i] {  raw_o.set_bit(i, false);  }
        }
        for i in 0..bit_i.len() {
            if !bit_i[i] {  raw_i.set_bit(i, false);  }
        }

        let test_o: BitVector = BitVector::from(raw_o.clone());
        let test_i: BitVector = BitVector::from(raw_i.clone());                                           

        let wg: TBWT = TBWT::new(&test_seq, test_o, test_i, raw_o, raw_i);
 
        let expected_letters = [0, 6, 5, 4, 3, 2, 1];

        let letters = wg.reconstruct();
        for i in 0..expected_letters.len() { assert_eq!(letters[i], expected_letters[i]); }
    }

    #[test]
    fn test_reconstruction2() {
        let test_seq = [5, 0, 2, 2, 3, 1, 1, 5, 1, 2, 2, 3, 3, 4, 4];
        let bit_o = vec![true, true, false, true, true, true, true, true, false, true, true, true, true, true, true, true];
        let bit_i = vec![true, true, true, true, true, false, true, true, true, true, false, true, true, true, true, true];

        let mut raw_o = RawVector::with_len(bit_o.len(), true);
        let mut raw_i = RawVector::with_len(bit_i.len(), true);

        for i in 0..bit_o.len() {
            if !bit_o[i] {  raw_o.set_bit(i, false);  }
        }
        for i in 0..bit_i.len() {
            if !bit_i[i] {  raw_i.set_bit(i, false);  }
        }

        let test_o: BitVector = BitVector::from(raw_o.clone());
        let test_i: BitVector = BitVector::from(raw_i.clone());                                           

        let wg: TBWT = TBWT::new(&test_seq, test_o, test_i, raw_o, raw_i);

        let expected_letters = [0, 5, 4, 3, 2, 5, 4, 3, 2, 1, 3, 2, 1, 2, 1, 2, 1];

        let letters = wg.reconstruct();
        for i in 0..expected_letters.len() { assert_eq!(letters[i], expected_letters[i]); }
    }

    #[test]
    fn test_reconstruction3() {
        let test_seq = [6, 0, 1, 9, 5, 2, 3, 3, 3, 9, 4, 5, 5, 5, 9, 6, 6, 7, 8];
        let bit_o = vec![true, true, true, true, true, true, true, false, false, false, true, true, true, true, true, false, false, true, true, true];
        let bit_i = vec![true, true, true, true, true, true, true, true, false, false, false, true, true, true, true, true, true, false, false, true];

        let mut raw_o = RawVector::with_len(bit_o.len(), true);
        let mut raw_i = RawVector::with_len(bit_i.len(), true);

        for i in 0..bit_o.len() {
            if !bit_o[i] {  raw_o.set_bit(i, false);  }
        }
        for i in 0..bit_i.len() {
            if !bit_i[i] {  raw_i.set_bit(i, false);  }
        }

        let test_o: BitVector = BitVector::from(raw_o.clone());
        let test_i: BitVector = BitVector::from(raw_i.clone());                                           

        let wg: TBWT = TBWT::new(&test_seq, test_o, test_i, raw_o, raw_i);

        let expected_letters = [0, 6, 5, 4, 3, 5, 4, 3, 9, 8, 7, 9, 8, 7, 6, 5, 4, 9, 8, 7, 6, 5, 4, 3, 2, 1];

        let letters = wg.reconstruct();
        for i in 0..expected_letters.len() { assert_eq!(letters[i], expected_letters[i]); }
    }

    #[test]
    fn test_reconstruction4() {
        let test_seq = [6, 6, 0, 6, 1, 6, 6, 6, 1, 2, 3, 3, 3, 3, 3, 3, 4, 4, 5, 5, 5, 5, 5, 4];
        let bit_o = vec![true, true, true, true, false, false, false, false, false, true, true, true, false, false, false, false, true, true, true, true, true, true, true, true, true];
        let bit_i = vec![true, true, true, true, true, false, false, false, false, false, true, true, true, true, false, false, false, false, true, true, true, true, true, true, true];

        let mut raw_o = RawVector::with_len(bit_o.len(), true);
        let mut raw_i = RawVector::with_len(bit_i.len(), true);

        for i in 0..bit_o.len() {
            if !bit_o[i] {  raw_o.set_bit(i, false);  }
        }
        for i in 0..bit_i.len() {
            if !bit_i[i] {  raw_i.set_bit(i, false);  }
        }

        let test_o: BitVector = BitVector::from(raw_o.clone());
        let test_i: BitVector = BitVector::from(raw_i.clone());                                           

        let wg: TBWT = TBWT::new(&test_seq, test_o, test_i, raw_o, raw_i);

        let expected_letters = [0, 6, 5, 4, 3, 2, 1, 6, 5, 4, 3, 2, 6, 5, 4, 3, 2, 6, 4, 4, 3, 2, 6, 5, 4, 3, 2, 6, 5, 4, 3, 2, 1];

        let letters = wg.reconstruct();
        for i in 0..expected_letters.len() { assert_eq!(letters[i], expected_letters[i]); }
    }

    #[test]
    fn test_reconstruction5() {
        let test_seq = [1, 2, 3, 4, 5, 1, 2, 2, 0];
        let bit_o = vec![true, true, true, true, true, false, false, true, false, true];
        let bit_i = vec![true, true, false, true, false, false, true, true, true, true];
 
        let mut raw_o = RawVector::with_len(bit_o.len(), true);
        let mut raw_i = RawVector::with_len(bit_i.len(), true);

        for i in 0..bit_o.len() {
            if !bit_o[i] {  raw_o.set_bit(i, false);  }
        }
        for i in 0..bit_i.len() {
            if !bit_i[i] {  raw_i.set_bit(i, false);  }
        }

        let test_o: BitVector = BitVector::from(raw_o.clone());
        let test_i: BitVector = BitVector::from(raw_i.clone());                                           

        let wg: TBWT = TBWT::new(&test_seq, test_o, test_i, raw_o, raw_i);

        let expected_letters = [0, 1, 2, 3, 4, 5, 2, 3, 4, 2, 3, 4, 1, 2, 3, 4, 5];

        let letters = wg.reconstruct();
        for i in 0..expected_letters.len() { assert_eq!(letters[i], expected_letters[i]); }
    }

    #[test]
    fn test_reconstruction_character_mapping() {
        let test_seq = [12, 13, 12, 13, 16, 16, 13, 11, 11, 11, 11, 14, 12, 11, 12, 16, 16, 16, 16, 16, 14, 14, 14, 10, 13, 15, 14, 13];
        let bit_o = vec![true, true, true, true, true, true, true, true, true, true, false, true, true, true, true, true, true, true, true, true, true, true, true, true, true, true, true, true, true];
        let bit_i = vec![true, true, true, true, true, true, true, true, true, true, true, true, true, true, false, true, true, true, true, true, true, true, true, true, true, true, true, true, true];

        let mut raw_o = RawVector::with_len(bit_o.len(), true);
        let mut raw_i = RawVector::with_len(bit_i.len(), true);

        for i in 0..bit_o.len() {
            if !bit_o[i] {  raw_o.set_bit(i, false);  }
        }
        for i in 0..bit_i.len() {
            if !bit_i[i] {  raw_i.set_bit(i, false);  }
        }

        let test_o: BitVector = BitVector::from(raw_o.clone());
        let test_i: BitVector = BitVector::from(raw_i.clone());                                           

        let wg: TBWT = TBWT::new(&test_seq, test_o, test_i, raw_o, raw_i);

        let expected_letters = [10, 12, 13, 11, 16, 14, 16, 14, 16, 13, 12, 11, 16, 14, 16, 15, 14, 16, 13, 12, 11, 13, 12, 11, 12, 11, 13, 14, 16];

        let letters = wg.reconstruct();
        for i in 0..expected_letters.len() { assert_eq!(letters[i], expected_letters[i]); }
    }

    #[test]
    fn test_reconstruciton_from_tunnel_hierarchy() {
        let test_seq = [2, 3, 6, 4, 0, 2, 6, 1, 5];
        let bit_o : Vec<bool> = vec![true, true, true, true, false, true, false, true, true, true];
        let bit_i : Vec<bool> = vec![true, true, true, false, true, true, true, true, false, true];

        let mut raw_o = RawVector::with_len(bit_o.len(), true);
        let mut raw_i = RawVector::with_len(bit_i.len(), true);

        for i in 0..bit_o.len() {
            if !bit_o[i] {  raw_o.set_bit(i, false);  }
        }
        for i in 0..bit_i.len() {
            if !bit_i[i] {  raw_i.set_bit(i, false);  }
        }

        let test_o: BitVector = BitVector::from(raw_o.clone());
        let test_i: BitVector = BitVector::from(raw_i.clone());                                           

        let wg: TBWT = TBWT::new(&test_seq, test_o, test_i, raw_o, raw_i);

        let expected_letters = [0, 2, 6, 5, 1, 3, 4, 2, 6, 5, 1, 3, 4, 6, 5, 1, 3];

        let letters = wg.reconstruct();
        for i in 0..expected_letters.len() { assert_eq!(letters[i], expected_letters[i]); }
    }

    #[test]
    fn showcase_how_select_is_used() {
        let bits = bit_vec![true, true, true, false, true, true, true, true];
        let rs = RankSelect::new(bits.clone(), 1);

        assert_eq!(rs.select(1), Some(0)); // first true is at position 0
        assert_eq!(rs.select(2), Some(1));
        assert_eq!(rs.select(4), Some(4)); // fourth true is as position 4 

    }

    proptest! {
        fn proptest_tunnel_binaryfile_reconstruct(s in "[a-z][a-z][a-z][a-z][a-z][a-z]+0") {
            let seq = s.into_bytes();
            let mut file = File::create("input").expect("Unable to create file");
            file.write_all(&seq.clone()).expect("unable to write to the file");

            let args = vec!["/home/klara/Documents/DIPLOMOVKA/wglib/input", "/home/klara/Documents/DIPLOMOVKA/wglib/output"];

            let output = Command::new("./../../BAKALARIS/BWT-Tunneling/seqana/tfm_index_construct.x")
                .arg("/home/klara/Documents/DIPLOMOVKA/wglib/input")
                .arg("/home/klara/Documents/DIPLOMOVKA/wglib/output")
                .output()
                .expect("failed to execute process");

            // Print the output of the command
            /*println!("Status: {}", output.status);
            if !output.stdout.is_empty() {
                println!("Output:\n{}", String::from_utf8_lossy(&output.stdout));
            }
            if !output.stderr.is_empty() {
                eprintln!("Error:\n{}", String::from_utf8_lossy(&output.stderr));
            }*/

            let bwt_file = "/home/klara/Documents/DIPLOMOVKA/wglib/bwt";
            let din_file = "/home/klara/Documents/DIPLOMOVKA/wglib/bwt.din";
            let dout_file = "/home/klara/Documents/DIPLOMOVKA/wglib/bwt.dout";

            match TBWT::from_binary(bwt_file, dout_file, din_file) {
                Ok(wg) => {
                    let letters = wg.reconstruct().iter().rev().copied().collect::<Vec<u8>>();
                    for i in 1..seq.len() { assert_eq!(letters[i], seq[i]); }
                }
                Err(e) => eprint!("Error read in file: {}", e),
            }
        }
    }

    #[test]
    fn test_non_selfcolliding_blocks() {
        let seq = [3, 1, 5, 6, 4, 3, 1, 5, 6, 2, 4, 3, 1, 5, 6, 2, 0];
        let sa  = suffix_array(&seq);
        let bwt = bwtransform(&seq, &sa);
        let wg  = TBWT::from_bwt2(&bwt);

        let mut blocks = wg.find_maximal_blocks();
        blocks = wg.remove_selfcolliding_blocks(blocks);

        assert_eq!(blocks, [(6, 4, 5), (4, 14, 16)]);
    }

    #[test]
    fn test_find_nonselfcolliding_blocks_on_example_with_selfcolliding_blocks() {
        let seq = [2, 1, 3, 4, 2, 1, 3, 4, 2, 1, 3, 5, 0];
        let sa  = suffix_array(&seq);
        let bwt = bwtransform(&seq, &sa);
        let wg  = TBWT::from_bwt2(&bwt);

        let mut blocks = wg.find_maximal_blocks();
        blocks = wg.remove_selfcolliding_blocks(blocks);

        assert_eq!(blocks, [(3, 7, 9)]);
    }

    #[test]
    fn test_heuristic_simple_example_with_one_tunnel() {
        let seq = [1, 2, 3, 4, 2, 3, 4, 5, 0];
        let sa = suffix_array(&seq);
        let bwt = bwtransform(&seq, &sa);
        let mut wg = TBWT::from_bwt2(&bwt);

        wg.heuristic_tunnel();

        let expected_l = [5, 0, 1, 4, 2, 3, 4];
        let expected_o = [true, true, true, false, true, true, true, true];
        let expected_i = [true, true, true, true, true, false, true, true];

        for i in 0..wg.i.len() {
            assert_eq!(wg.i.get(i), expected_i[i]);
        }

        for i in 0..wg.o.len() {
            assert_eq!(wg.o.get(i), expected_o[i]);
        }

        assert_eq!(wg.bwt, expected_l);

        let result = wg.reconstruct().iter().rev().copied().collect::<Vec<u8>>();
        assert_eq!(result, seq);
    }

    #[test]
    fn test_heuristic_simple_colliding_example() {
        let seq = [5, 4, 3, 2, 1, 4, 3, 2, 4, 3, 2, 5, 4, 3, 2, 1, 0];
        let sa = suffix_array(&seq);
        let bwt = bwtransform(&seq, &sa);
        let mut wg = TBWT::from_bwt2(&bwt);

        //wg.heuristic_tunnel();
        wg.mark_tunnel((5, 1, 2));
        wg.mark_tunnel((3, 3, 6));

        wg.tunnel();

        let expected_l = [1, 2, 3, 4, 5, 1, 2, 2, 0];
        let expected_o = [true, true, true, true, true, false, false, true, false, true];
        let expected_i = [true, true, false, true, false, false, true, true, true, true];

        for i in 0..wg.i.len() {
            assert_eq!(wg.i.get(i), expected_i[i]);
        }

        for i in 0..wg.o.len() {
            assert_eq!(wg.o.get(i), expected_o[i]);
        }

        assert_eq!(wg.bwt, expected_l);

        let result = wg.reconstruct().iter().rev().copied().collect::<Vec<u8>>();
        assert_eq!(result, seq);
}

    proptest! {
        #[test]
        fn proptest_heuristic_and_reconstruct(s in "[a-d][a-d][a-d][a-d][a-d][a-d]+0") {
            let seq = s.clone().into_bytes();
            let sa = suffix_array(&seq);
            let bwt = bwtransform(&seq, &sa);
            let mut wg = TBWT::from_bwt2(&bwt);
 
            wg.heuristic_tunnel();

            let expected_letters = seq.clone().into_iter().rev().collect::<Vec<u8>>();
            wg.reconstruct_with_expected(seq);
        }
    }


    #[test]
    fn test_heuristics() {
        let mut file = File::open("data/example.txt").expect("Unable to open file ");
        //let mut file = File::open("data/protein.fasta").expect("Unable to open file ");
        //let mut file = File::open("data/zinc_fingers.fa").expect("Unable to open file ");
        //let mut file = File::open("data/bacteriophage.fasta").expect("Unable to open file ");
        //let mut file = File::open("data/S-cereale.fasta").expect("Unable to open file ");
        //let mut file = File::open("data/huYchr.fasta").expect("Unable to open file ");
        //let mut file = File::open("data/repetitive.txt").expect("Unable to open file ");
        //let mut file = File::open("data/HIV.txt").expect("Unable to open file");

        //let mut file = File::open("data/corona_virus.fasta").expect("Unable to open file");
        //let mut file = File::open("data/prion_protein.fasta").expect("Unable to open file");
        //let mut file = File::open("data/drosophila_protein.fasta").expect("Unable to open file"); 
        //let mut file = File::open("data/ecoli_plasmid.fasta").expect("Unable to open file"); 
        //let mut file = File::open("data/HIV.fasta").expect("Unable to open file"); 
        //let mut file = File::open("data/Caplha.fasta").expect("Unable to open file");
        //let mut file = File::open("data/G_suppressor.fasta").expect("Unable to open file");
        //let mut file = File::open("data/penicilium.fasta").expect("Unable to open file");
        //let mut file = File::open("data/zinc_fingers.fasta").expect("Unable to open file");
        //let mut file = File::open("data/stachybotris_elegans.fasta").expect("Unable to open file");
        //let mut file = File::open("data/nematocida.fasta").expect("Unable to open file");

        let mut seq = Vec::new();
        file.read_to_end(&mut seq).expect("Unable to read file");
        seq.push(0);

        let sa = suffix_array(&seq);
        let bwt = bwtransform(&seq, &sa);
        let mut wg = TBWT::from_bwt2(&bwt);

        let start = Instant::now();
        wg.heuristic_tunnel();
        println!("{:?}", start.elapsed());

        println!("{:?}", seq.len());
        println!("{:?}", wg.bwt.len());

        wg.reconstruct_with_expected(seq.clone());
    }

    #[test]
    fn bitvector_usage() {
        let data :Vec<bool> = vec![false, true, true, true, false];
        let mut bv: BitVector = data.iter().cloned().collect();
        assert_eq!(bv.len(), data.len(), "Invalid bitvector length");

        for i in 0..bv.len() {
            assert_eq!(bv.get(i), data[i], "Invalid bit {}", i);
        }

        bv.enable_select();
        println!("{:?}", bv.select(1).unwrap());
        bv.enable_rank();
        println!("{:?}", bv.rank(4));
    }
}


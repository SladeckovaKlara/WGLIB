use bio::data_structures::suffix_array::suffix_array;
use bio::data_structures::bwt::bwt as bwtransform;
use std::collections::HashMap;
use bio::data_structures::wavelet_matrix::WaveletMatrix;
use bv::BitVec;
use bv::BitsMut;
use bv::bit_vec;

use bio::data_structures::rank_select::RankSelect;

fn bit_rank(vec: BitVec, idx: u64) -> u64 {
    let mut rk: u64 = 0;
    for i in 0..idx {
        if vec.get(i) {
            rk += 1;
        }
    }
    rk
}


struct TBWT {
    bwt: Vec<u8>,           // for access to bwt
    o: RankSelect,          // outgoing nodes
    i: RankSelect,          // incoming nodes

    l: Vec<RankSelect>,     // for ranks on bwt used in lf mapping
    f: Vec<usize>           // f[c] is the number of letters smaller than c
}

impl TBWT {
    fn new(l: &[u8], o: BitVec<u8>, i: BitVec<u8>) -> TBWT {
        let mut wm: Vec<RankSelect> = Vec::new();

        //find size of alphabet
        let mut alphabet = HashMap::new();
        for c in l {
            alphabet.entry(c).or_insert(1);
        }
        let sigma: usize = alphabet.keys().len();

        //get necessary info from bwt
        let bwt = l;
        for i in 0..sigma {
            let mut bitvec: BitVec<u8> = BitVec::new_fill(false, bwt.len() as u64);
            for j in 0..bwt.len() {
                if bwt[j] == i as u8 { bitvec.set_bit(j as u64, true); }
            }
            wm.push(RankSelect::new(bitvec, 1));
        }

        let mut char_occ: Vec<usize> = vec![0; sigma];
        for c in bwt { char_occ[*c as usize] += 1; }
        let mut f: Vec<usize> = vec![0; sigma];
        for i in 1..sigma { f[i] = f[i-1] + char_occ[i-1]; }

        let k : usize = o.len().pow(2) as usize / 32;
        let o_wt = RankSelect::new(o,k);
        let i_wt = RankSelect::new(i,k);
        TBWT {bwt: bwt.to_vec() , l: wm, o: o_wt, i: i_wt, f }
    }


    fn from_bwt(bwt: &[u8], sigma: usize) -> TBWT {
        let mut wm: Vec<RankSelect> = Vec::new();
        for i in 0..sigma {
            let mut bitvec: BitVec<u8> = BitVec::new_fill(false, bwt.len() as u64);
            for j in 0..bwt.len() {
                if bwt[j] == i as u8 { bitvec.set_bit(j as u64, true); }
            }
            wm.push(RankSelect::new(bitvec, 1));
        }

        let mut char_occ: Vec<usize> = vec![0; sigma];
        for c in bwt { char_occ[*c as usize] += 1; }
        let mut f: Vec<usize> = vec![0; sigma];
        for i in 1..sigma { f[i] = f[i-1] + char_occ[i-1]; }

        let k : usize = (bwt.len()+1).pow(2) as usize / 32;
        let o_wt = RankSelect::new(BitVec::new_fill(true, (bwt.len()+1) as u64),k);
        let i_wt = RankSelect::new(BitVec::new_fill(true, (bwt.len()+1) as u64),k);

        TBWT { 
            bwt: bwt.to_vec(),
            l: wm,
            o: o_wt,
            i: i_wt,
            f
        }
    }

    fn backward_step(&self, node: usize, offset: usize) -> (usize, usize) {
        let mut letter = 0;
        for i in 0..self.l.len() {
            if self.l[i].get(node as u64) {
                letter = i;
                break
            }
        }
        let mut idx = self.f[letter] as u64 + self.l[letter].rank(node as u64).unwrap()-1;
        let node_rank = self.i.rank(idx).unwrap();
        let mut e = offset;

        if self.i.bits().len() > idx + 1 {
            if self.i.get(idx) == false || self.i.get(idx + 1) == false {
                e = (idx - self.i.select(node_rank).unwrap()) as usize;
            }
        }

        idx = self.o.select(node_rank).unwrap();
        if self.o.bits().len() > idx + 1 {
            if self.o.get(idx + 1) == false {
                idx += e as u64;
                e = 0;
            }
        }
        
        return (idx as usize, e);
    } 

    fn incoming_letter(&self, node:usize, offset: usize) -> u8 {
        let idx = (self.o.select(node as u64 + 1)).unwrap() % self.l[0].bits().len();
        for i in 0..self.l.len() {
            if self.l[i].get(idx)  {
                return i as u8
            }
        }
        return 0;
    }

    fn mark_tunnel(&mut self, block: (usize, usize, usize)) {}
    fn tunnel(&mut self) {}
    fn print(&self) {
        println!("L:");
        print!("[");
        for i in 0..self.l[0].bits().len() {
            if i != 0 {
                print!(", ");
            }
            for j in 0..self.l.len() {
                if self.l[j].get(i) {
                    print!("{}", j);
                    break
                }
            }
        }
        println!("]");
        println!("O:");
        print!("[");
        for i in 0..self.o.bits().len() {
            if i != 0 {
                print!(", ");
            }
            print!("{}", self.o.get(i));
        }
        println!("]");
        println!("I:");
        print!("[");
        for i in 0..self.i.bits().len() {
            if i != 0 {
                print!(", ");
            }
            print!("{}", self.i.get(i));
        }
        println!("]");
        println!("F:");
        println!("{:?}", self.f);
    }
}

fn main() {
    //let seq = b"AACCAGCGGATCTGTTAACCAGCGGATCTGTTAACCAGCGGATCTGTTA$";
    // let seq = b"AAAAAA$";

    //let sa = suffix_array(seq);
    //let bwt = bwtransform(seq, &sa);
    
    
    let seq2: [u8;50] = [1,1,2,2,1,3,2,3,3,1,4,2,4,3,4,4,1,1,2,2,1,3,2,3,3,1,4,2,4,3,4,4,1,1,2,2,1,3,2,3,3,1,4,2,4,3,4,4,1,0];
    let sa2 = suffix_array(&seq2);
    let bwt2 = bwtransform(&seq2, &sa2);

    //let mut wg: TBWT = TBWT::from_bwt(&bwt2,6);

    let test_seq : [u8;9] = [3,2,2,2,1,1,1,0,2];
    let test_o : BitVec<u8> = bit_vec![true,true,false,false,true,false,false,false,true,true];
    let test_i : BitVec<u8> = bit_vec![true,true,false,false,true,false,false,false,true,true];

    let mut test_wg: TBWT = TBWT::new(&test_seq, test_o, test_i);


    /*let blocks = [(16, 1, 3)];
    for b in blocks {
        wg.mark_tunnel(b);
    }
    wg.tunnel();*/
    

}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn simple_case() {
        // input
        let seq: [u8;9] = [1,2,3,4,2,3,4,5,0];
        let sa = suffix_array(&seq);
        let bwt = bwtransform(&seq,&sa);
        let mut wg: TBWT = TBWT::from_bwt(&bwt,6);
        // expected output
        let result_edges = [8,7,5,3,6,4,2,1,0];
        let result_offs = [0,0,0,0,0,0,0,0,0];

        // what is done
        wg.print();
        println!("{}", wg.backward_step(4,0).0);

        let mut edge = 0;
        let mut off = 0;
        println!("{}", wg.incoming_letter(edge,off));

        let mut i = 0;

        while true {
            println!("{} {}", edge, off);
            (edge,off) = wg.backward_step(edge,off);
            assert_eq!(result_edges[i],edge);
            assert_eq!(result_offs[i],off);
            i+=1;
            if edge == 0 && off == 0 {
                break;
            }
        }
        println!();
    }

    #[test]
    fn tunneled_case() {
        // input
        let test_seq : [u8;7] = [5, 0, 1, 4, 2, 3, 4];
        let test_o : BitVec<u8> = bit_vec![true, true, true, false, true, true, true];
        let test_i : BitVec<u8> = bit_vec![true, true, true, true, true, false, true];
        let mut wg: TBWT = TBWT::new(&test_seq, test_o, test_i);
        // expected output
        let result_edges = [6,5,4,3,5,4,2,1,0];
        let result_offs = [0,1,1,0,0,0,0,0,0];

        // what is done

        let mut edge = 0;
        let mut off = 0;
        println!("{}", wg.incoming_letter(edge,off));

        let mut i = 0;

        while true {
            println!("{} {}", edge, off);
            (edge,off) = wg.backward_step(edge,off);
            assert_eq!(result_edges[i],edge);
            assert_eq!(result_offs[i],off);
            i+=1;
            if edge == 0 && off == 0 {
                break;
            }
        }
        println!();
    }
}

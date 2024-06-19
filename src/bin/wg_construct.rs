use bio::data_structures::bwt::bwt as bwtransform;
use bio::data_structures::suffix_array::{lcp as lcp_array, suffix_array};
use std::collections::{HashMap, HashSet};
use std::str;
use wglib::TBWT;

fn lcs_array(text: &[u8], sa: &Vec<usize>) -> Vec<isize> {
    let mut lcs = Vec::new();
    lcs.push(-1);

    let n = sa.len();
    for row in 1..n {
        let i = sa[row - 1];
        let j = sa[row];

        let mut w = 0;
        let mut i_is_dollar = text[(n + i - 1 - w) % n] == b'$';
        let mut j_is_dollar = text[(n + j - 1 - w) % n] == b'$';
        let mut i_equals_j = text[(n + i - 1 - w) % n] == text[(n + j - 1 - w) % n];
        while !i_is_dollar && !j_is_dollar && i_equals_j {
            w += 1;
            i_is_dollar = text[(n + i - 1 - w) % n] == b'$';
            j_is_dollar = text[(n + j - 1 - w) % n] == b'$';
            i_equals_j = text[(n + i - 1 - w) % n] == text[(n + j - 1 - w) % n];
        }
        lcs.push(w as isize);
    }
    lcs.push(-1);
    return lcs;
}

fn psv_array(arr: &[isize]) -> Vec<isize> {
    let mut res: Vec<isize> = Vec::new();
    res.push(-1);
    let mut stack = Vec::new();
    stack.push((arr[0], 0));

    for i in 1..arr.len() - 1 {
        let (mut v, mut pos) = stack.pop().unwrap();
        while v >= arr[i] {
            (v, pos) = stack.pop().unwrap();
        }
        res.push(pos);
        stack.push((v, pos));
        stack.push((arr[i], i as isize));
    }
    res.push(-1);

    return res;
}

fn nsv_array(arr: &Vec<isize>) -> Vec<isize> {
    let mut lcs_rev = arr.clone();
    lcs_rev.reverse();
    let mut res = psv_array(&lcs_rev);
    res.reverse();
    let n = res.len();
    res.iter().map(|x| (n - 1) as isize - x).collect::<Vec<_>>()
}

fn get_blocks(lcs: &Vec<isize>) -> Vec<(usize, usize, usize)> {
    let psv = psv_array(lcs);
    let nsv = nsv_array(lcs);
    let n = lcs.len() - 1;
    let mut blocks = Vec::new();

    for i in 1..n {
        if lcs[i] > 0 {
            blocks.push((lcs[i] as usize, psv[i] as usize, (nsv[i] - 1) as usize));
        }
    }
    return blocks;
}

fn translate_blocks(blocks: &[(usize, usize, usize)], sa: &[usize]) -> Vec<(usize, usize, usize)> {
    let mut isa = vec![0; sa.len()];
    for i in 0..sa.len() {
        isa[sa[i]] = i;
    }

    let mut res = Vec::new();
    for block in blocks {
        let b = (
            block.0,
            isa[sa[block.1] - block.0],
            isa[sa[block.2] - block.0],
        );
        res.push(b);
    }
    return res;
}

fn translate_back(blocks: &[(usize, usize, usize)], sa: &[usize]) -> Vec<(usize, usize, usize)> {
    let mut isa = vec![0; sa.len()];
    for i in 0..sa.len() {
        isa[sa[i]] = i;
    }

    let mut res = Vec::new();
    for block in blocks {
        let b = (
            block.0,
            isa[sa[block.1] + block.0],
            isa[sa[block.2] + block.0],
        );
        res.push(b);
    }
    return res;
}

fn get_maximal_blocks(blocks: &[(usize, usize, usize)]) -> Vec<(usize, usize, usize)> {
    let mut map: HashMap<(usize, usize), usize> = HashMap::new();
    for block in blocks {
        let l = map.get(&(block.1, block.2));
        match l {
            Some(l) => {
                map.insert((block.1, block.2), usize::max(*l, block.0));
            }
            None => {
                map.insert((block.1, block.2), block.0);
            }
        }
    }
    let mut res = Vec::new();
    for (k, v) in map.iter() {
        res.push((*v, k.0, k.1));
    }
    return res;
}

fn is_self_colliding(block: (usize, usize, usize), lf: &[usize]) -> bool {
    let mut set = HashSet::new();
    let (l, mut i, mut j) = block;
    for _ in 0..l {
        if set.contains(&i) {
            return true;
        } else {
            set.insert(i);
        }
        if set.contains(&j) {
            return true;
        } else {
            set.insert(j);
        }
        i = lf[i];
        j = lf[j];
    }
    return false;
}

fn compute_lf(sa: &[usize]) -> Vec<usize> {
    let mut isa = vec![0; sa.len()];
    for i in 0..sa.len() {
        isa[sa[i]] = i;
    }

    let mut lf = Vec::new();
    for i in 0..sa.len() {
        let x;
        if sa[i] != 0 {
            x = isa[sa[i] - 1];
        } else {
            x = 0;
        }
        lf.push(x);
    }
    return lf;
}

fn get_alphabet2num(seq: &[u8]) -> HashMap<u8, u8> {
    let mut set = HashSet::new();
    for i in 0..seq.len() {
        set.insert(seq[i]);
    }
    let mut vec: Vec<u8> = set.into_iter().collect();
    vec.sort();
    let mut map = HashMap::new();
    for i in 0..vec.len() {
        map.insert(vec[i], i as u8);
    }
    return map;
}

fn normalize_dna(seq: &[u8], map: &HashMap<u8, u8>) -> Vec<u8> {
    let mut res = Vec::with_capacity(seq.len());
    for i in 0..seq.len() {
        match map.get(&seq[i]) {
            Some(v) => {
                res.push(*v);
            }
            None => {
                panic!("Unexpected letter {}", seq[i]);
            }
        }
    }
    return res;
}

fn print_table(
    seq: &[u8],
    blocks: &Vec<(usize, usize, usize)>,
    sa: &[usize],
    lcp: &[isize],
    lcs: &[isize],
) {
    let n = seq.len();
    let mut tmp_blocks = blocks.clone();
    tmp_blocks.reverse();
    let mut b = tmp_blocks.pop().unwrap();
    for i in 0..n {
        print!("{:2} {:2} {:2} ", i, sa[i], lcp[i]);
        let s = str::from_utf8(&seq[sa[i]..n]).unwrap().to_owned()
            + str::from_utf8(&seq[0..sa[i]]).unwrap();
        print!("{} {:2}", s, lcs[i]);
        if lcs[i] == b.0 as isize {
            print!(" {:?}", b);
            if let Some(top) = tmp_blocks.pop() {
                b = top;
            };
        }
        println!();
    }
}

fn filter_blocks(blocks: &[(usize, usize, usize)], sa: &[usize]) -> Vec<(usize, usize, usize)> {
    let lf = compute_lf(&sa);
    let mut filtered_blocks = Vec::new();
    for b in blocks {
        if is_self_colliding(*b, &lf) {
            println!("{:?} is self-colliding.", b);
            continue;
        }
        if b.0 < 2 {
            println!("{:?} is too short.", b);
            continue;
        }
        filtered_blocks.push(*b);
    }

    return filtered_blocks;
}

fn main() {
    let seq = b"AACCAGCGGATCTGTTAACCAGCGGATCTGTTAACCAGCGGATCTGTTA$";
    
    let map = get_alphabet2num(seq);
    println!("{:?}", map);
    let norm_seq = normalize_dna(seq, &map);
    println!("{:?}", norm_seq);

    let sa = suffix_array(seq);
    let lcp = lcp_array(seq, &sa).decompress();
    let bwt = bwtransform(&norm_seq, &sa);
    let lcs = lcs_array(seq, &sa);
    let blocks = get_blocks(&lcs);
    print_table(seq, &blocks, &sa, &lcp, &lcs);

    println!("Original blocks: {:?}", blocks);
    let blocks = translate_blocks(&blocks, &sa);
    println!("Translated to text order: {:?}", blocks);
    let blocks = get_maximal_blocks(&blocks);
    println!("Maximal blocks {:?}", blocks);
    let blocks = translate_back(&blocks, &sa);
    println!("Translated to bwt order {:?}", blocks);
    let blocks = filter_blocks(&blocks, &sa);
    println!("Filtered blocks: {:?}", blocks);
    println!();


    // tunelling starts
    let mut wg: TBWT = TBWT::from_bwt(&bwt, map.len());
    
    let seq_reconstructed = wg.reconstruct(); // reconstructed sequence is reversed?
    println!("{:?}", seq_reconstructed);

    for b in blocks {
        wg.mark_tunnel(b);
    }
    wg.tunnel();
    wg.print();

    // let seq3 = wg.reconstruct();
    // println!("{:?}", seq3);

}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn tunnelling_works() {
        let blocks = [(2, 6, 7)];
        let bwt = b"E$ADBBCCD";
        let map = get_alphabet2num(bwt);
        println!("{:?}", map);
        let norm_bwt = normalize_dna(bwt, &map);
        println!("{:?}", norm_bwt);

        let mut wg: TBWT = TBWT::from_bwt(&norm_bwt, map.len());
        for b in blocks {
            wg.print();
            wg.mark_tunnel(b);
        }
        wg.print();
        wg.tunnel();
        wg.print();
    }
}




struct WheelerGraph {
    adj_list: Vec<Vec<(usize, u8)>>
}

impl WheelerGraph {
    fn new(n: usize) -> WheelerGraph {
        WheelerGraph {
            adj_list: vec![Vec::new(); n]
        }
    }

    fn is_valid(&self) -> bool {
        // 0 in-degree at the beginning
        let mut in_deg = vec![0; self.adj_list.len()];
        for (_u, list) in self.adj_list.iter().enumerate() {
            for (v, _l) in list {
                in_deg[*v] += 1;
            }
        }
        for i in 1..self.adj_list.len() {
            if in_deg[i] == 0 && in_deg[i-1] > 0 {
                println!("{:?}", in_deg);
                return false;
            }
        }

        // a < a' => v < v'
        // a = a' && u < u' => v =< v'
        for (u1, list1) in self.adj_list.iter().enumerate() {
            for (v1, l1) in list1 {
                for (u2, list2) in self.adj_list.iter().enumerate() {
                    for (v2, l2) in list2 {
                        if l1 < l2 && v1 >= v2 {
                            println!("({} {} {}) and ({} {} {}) violate a < a' => v < v'", u1, v1, l1, u2, v2, l2);
                            return false;
                        }
                        if l1 == l2 && u1 < u2 && v1 > v2 {
                            println!("({} {} {}) and ({} {} {}) violate a = a' AND u < u' => v =< v'", u1, v1, l1, u2, v2, l2);
                            return false;
                        }
                    }
                }
            }
        }

        return true;
    }

    fn add_edge(&mut self, u: usize, v: usize, l: u8) {
        self.adj_list[u].push((v, l));
    }

    // fn edges(&self) -> Iter<(usize, usize, u8)> {}

    fn _print(&self) {
        println!("digraph G {{");
        for (u, list) in self.adj_list.iter().enumerate() {
            for (v, l) in list {
                println!("    {} -> {} [label=\"{}\"]", u, v, *l as char);
            }
        }
        println!("}}");
    }
}

fn main() {
    let n = 10;
    let mut wg = WheelerGraph::new(n);
    let edge_list = [
        (3, 6, b'G'), 
        (9, 3, b'A'),
        (8, 9, b'T'),
        (1, 8, b'T'),
        (4, 1, b'A'),
        (2, 4, b'C'),
        (7, 2, b'A'),
        (0, 7, b'T'),
        (6, 0, b'$'),
        (5, 1, b'A'),
        (2, 5, b'G'),
    ];
    for (u, v, l) in edge_list {
        wg.add_edge(u, v, l);
    }

    println!("{:?}", edge_list);

    // wg._print();

    if wg.is_valid() {
        println!("Wheeler graph is valid.");
    }

    let edge_list = [
        (3, 6, b'G'), 
        (9, 3, b'A'),
        (8, 9, b'T'),
        (1, 8, b'T'),
        (4, 1, b'A'),
        (2, 4, b'C'),
        (7, 2, b'A'),
        (0, 7, b'T'),
        (6, 0, b'$'),
        (5, 1, b'A'),
        (2, 5, b'T'),
    ];
    let mut wg = WheelerGraph::new(n);
    for (u, v, l) in edge_list { wg.add_edge(u, v, l); }
    wg._print();
    if wg.is_valid() { println!("Wheeler graph is valid."); }
}
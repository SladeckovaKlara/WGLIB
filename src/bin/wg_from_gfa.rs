use gfa::parser::GFAParser;
use gfa::optfields::*;

fn main() {
    let parser: GFAParser<Vec<u8>, Vec<OptField>> = GFAParser::new();
    let gfa = parser.parse_file("data/lil.gfa").expect("Error parsing gfa");
    for seg in gfa.segments {
        println!("{:?}", seg);
    }
    for link in gfa.links {
        println!("{:?}", link);
    }
}

use bio::io::fasta;
use crate::helper;
use std::collections;

pub fn all_records (fasta_file_path: &str)
    -> Result<collections::HashMap<String,String>, helper::PublicError>
{
    let mut r = collections::HashMap::<String,String>::new ();
    let mut records = fasta::Reader::from_file (fasta_file_path)?.records ();

    while let Some (Ok(record)) = records.next () {
        r.insert (record.id ().to_string (), String::from_utf8 (record.seq ().to_vec ()).unwrap ());
    }
    Ok (r)
}

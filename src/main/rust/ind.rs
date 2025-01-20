
use crate::file;
use crate::helper;
use itertools::Itertools;
use log::debug;
use std::collections;
use std::io;
use std::iter;
use std::str;

/// Implementation is inspired by https://doi.org/10.1101/082214
pub fn revcomp (seq: String)
    -> String
{
    seq.bytes ()
        .rev ()
        .map (|c| { let xb = [if c & 2 != 0 { c ^ 4 } else { c ^ 21 }];str::from_utf8 (&xb).unwrap ().to_owned () })
        .collect::<String> ()
}

pub fn ni_to_seq (nucleotides: &Vec<String>, ni: &Vec<usize>)
    -> String
{
    ni.iter ().map (|&x| { nucleotides.get (x).unwrap () }).cloned ().collect::<String> ()
}

pub fn seq_to_ni (nucleotides: &Vec<String>, seq: &str)
    -> Vec<usize>
{
    seq.bytes ().map (|x| { let xb = [x];let xs = str::from_utf8 (&xb).unwrap ();nucleotides.iter ().position (|y| { y == xs }).unwrap () }).collect::<Vec<_>> ()
}

pub fn create_indices (nucleotides: &Vec<String>, index_length: u64)
    -> Vec<Vec<usize>>
{
    let multi_prod = (0..index_length).map (|_| 0..nucleotides.len ())
        .multi_cartesian_product()
        .collect::<Vec<_>> ();

    /*
    for p in &multi_prod
    {
        debug! ("p: {:?} {}", p, ni_to_seq (&nucleotides, &p));
    }
    */
    multi_prod
}

pub fn pick_indices (adapter_three_prime: String, nucleotides: &Vec<String>, skip_first: bool, multi_prod: &Vec<Vec<usize>>)
    -> Result<Vec<Vec<usize>>, helper::PublicError>
{

    let mut not_liked = (0..nucleotides.len ()).map (|i| vec![i,i,i]).collect::<Vec<_>> ();
    not_liked.push (seq_to_ni (&nucleotides, "ACA"));
    not_liked.push (seq_to_ni (&nucleotides, "CAC"));
    not_liked.push (seq_to_ni (&nucleotides, "GTG"));
    not_liked.push (seq_to_ni (&nucleotides, "TGT"));

    for p in &not_liked
    {
        debug! ("not_liked: {:?}", p);
    }

    let gg = seq_to_ni (&nucleotides, "GG");

    let reasonable = multi_prod.iter ()
        .filter (|x| {
            !x.windows (3).any (|y| not_liked.contains (&y.to_vec ()))
        })
        .filter (|&x| {
            (0..nucleotides.len ()).map (|y| {
                    x.iter ().filter (|&z| z == &y ).count ()
                })
                .filter (|&y| {
                    y == 0usize
                })
                .count () < 2
        })
        .filter (|&x| {
            let mm = iter::zip (gg.iter (), x.iter ())
                            .filter (|(y,z)| y == z)
                            .count ();
            /*
            if mm == 2
            {
                debug! ("when testing {} found gg in pick mm: {}", ni_to_seq (&nucleotides, &x), mm);
            }
            */
            mm != 2
        })
        .collect::<Vec<_>> ();

    if reasonable.is_empty ()
    {
        Err (helper::PublicError::ApplicationError (String::from ("No reasonable indices to pick from")))
    }
    else
    {
        let mut picks: Vec<Vec<usize>> = Vec::with_capacity (reasonable.len ());
        let mut skips: Vec<Vec<usize>> = Vec::with_capacity (reasonable.len ());
        //let mut oi = 0;

        debug! ("picking from {}", reasonable.len ());

        if skip_first
        {
            let adapter_nucleotide = revcomp (adapter_three_prime).chars ().take (1).collect::<String> ();
            debug! ("adapter_nucleotide: '{}'", adapter_nucleotide);

            for p in &reasonable
            {
                //oi += 1;
                //debug! ("pr: {:?} {}", p, ni_to_seq (&nucleotides, &p));
                let search_pick = picks.iter ()
                    .rev ()
                    .find (|pick| {
                        let mm = iter::zip (pick.iter (), p.iter ())
                            .filter (|(x,y)| x != y)
                            .count ();
                        /*
                        if mm <= 2
                        {
                            debug! ("when searching for {} found hit {} in pick mm: {}", ni_to_seq (&nucleotides, &p), ni_to_seq (&nucleotides, &pick), mm);
                        }
                        */
                        mm <= 2
                    });

                let search_skip = skips.iter ()
                    .find (|pick| {
                        let mm = iter::zip (pick.iter (), p.iter ())
                            .filter (|(x,y)| x != y)
                            .count ();
                        /*
                        if mm <= 1
                        {
                            debug! ("when searching for {} found hit {} in skip mm: {}", ni_to_seq (&nucleotides, &p), ni_to_seq (&nucleotides, &pick), mm);
                        }
                        */
                        mm <= 1
                    });

                if search_pick.is_none () && search_skip.is_none ()
                {
                    let psf = p.iter ()
                        .skip (1)
                        .copied ()
                        .chain (seq_to_ni (&nucleotides, &adapter_nucleotide))
                        .collect::<Vec<_>> ();

                    let search_psf = skips.iter ()
                        .chain (picks.iter ().rev ())
                        .find (|pick| {
                            let mm = iter::zip (pick.iter (), psf.iter ())
                                .filter (|(x,y)| x != y)
                                .count ();
                            /*
                            if mm <= 1
                            {
                                debug! ("when searching for {} found hit {} in psf mm: {}", ni_to_seq (&nucleotides, &psf), ni_to_seq (&nucleotides, &pick), mm);
                            }
                            */
                            mm <= 1
                        });
                    if search_psf.is_none ()
                    {
                        skips.push (psf);
                        picks.push ((*p).to_vec ());
                    }
                }
            }
        }
        else
        {
            for p in &reasonable
            {
                //oi += 1;
                //debug! ("pr: {:?} {}", p, ni_to_seq (&nucleotides, &p));
                let search_pick = picks.iter ()
                    .rev ()
                    .find (|pick| {
                        let mm = iter::zip (pick.iter (), p.iter ())
                            .filter (|(x,y)| x != y)
                            .count ();
                        /*
                        if mm <= 2
                        {
                            debug! ("when searching for {} found hit {} in pick mm: {}", ni_to_seq (&nucleotides, &p), ni_to_seq (&nucleotides, &pick), mm);
                        }
                        */
                        mm <= 2
                    });

                if search_pick.is_none ()
                {
                    picks.push ((*p).to_vec ());
                }
            }
        }
        Ok (picks)
    }
}


pub fn pick_indices_florian (adapter_three_prime: String, nucleotides: &Vec<String>, skip_first: bool, multi_prod: &Vec<Vec<usize>>)
    -> Result<Vec<Vec<usize>>, helper::PublicError>
{

    let mut not_liked = (0..nucleotides.len ()).map (|i| vec![i,i,i]).collect::<Vec<_>> ();
    not_liked.push (seq_to_ni (&nucleotides, "ACA"));
    not_liked.push (seq_to_ni (&nucleotides, "CAC"));
    not_liked.push (seq_to_ni (&nucleotides, "GTG"));
    not_liked.push (seq_to_ni (&nucleotides, "TGT"));

    for p in &not_liked
    {
        debug! ("not_liked: {:?}", p);
    }

    let reasonable = multi_prod.iter ()
        .filter (|x| {
            !x.windows (3).any (|y| not_liked.contains (&y.to_vec ()))
        })
        .filter (|&x| {
            (0..nucleotides.len ()).map (|y| {
                    x.iter ().filter (|&z| z == &y ).count ()
                })
                .filter (|&y| {
                    y == 0usize
                })
                .count () < 2
        })
        .collect::<Vec<_>> ();

/*
    // test case (only first four should be picked)
    let reasonable = vec![
        seq_to_ni (&nucleotides, "AACCATGC"),
        seq_to_ni (&nucleotides, "AACCGCCA"),
        seq_to_ni (&nucleotides, "AACGCCAT"),
        seq_to_ni (&nucleotides, "AATAAGGA"),
        seq_to_ni (&nucleotides, "AATTAAGG"),
        seq_to_ni (&nucleotides, "ACCAGGCA"),
    ];
    */

    if reasonable.is_empty ()
    {
        Err (helper::PublicError::ApplicationError (String::from ("No reasonable indices to pick from")))
    }
    else
    {
        let mut picks: Vec<Vec<usize>> = Vec::with_capacity (reasonable.len ());
        let mut skips: Vec<Vec<usize>> = Vec::with_capacity (reasonable.len ());
        //let mut oi = 0;

        debug! ("picking from {}", reasonable.len ());

        if skip_first
        {
            let adapter_nucleotide = revcomp (adapter_three_prime).chars ().take (1).collect::<String> ();
            debug! ("adapter_nucleotide: '{}'", adapter_nucleotide);

            for p in &reasonable
            {
                //oi += 1;
                //debug! ("pr: {:?} {}", p, ni_to_seq (&nucleotides, &p));
                let search_pick = picks.iter ()
                    .rev ()
                    .find (|pick| {
                        let mm = iter::zip (pick.iter (), p.iter ())
                            .filter (|(x,y)| x != y)
                            .count ();
                        /*
                        if mm <= 2
                        {
                            debug! ("when searching for {} found hit {} in pick mm: {}", ni_to_seq (&nucleotides, &p), ni_to_seq (&nucleotides, &pick), mm); 
                        }
                        */
                        mm <= 2
                    });

                let search_skip = skips.iter ()
                    .find (|pick| {
                        let mm = iter::zip (pick.iter (), p.iter ())
                            .filter (|(x,y)| x != y)
                            .count ();
                        /*
                        if mm <= 1
                        {
                            debug! ("when searching for {} found hit {} in skip mm: {}", ni_to_seq (&nucleotides, &p), ni_to_seq (&nucleotides, &pick), mm); 
                        }
                        */
                        mm <= 1
                    });

                if search_pick.is_none () && search_skip.is_none ()
                {
                    let psf = p.iter ()
                        .skip (1)
                        .copied ()
                        .chain (seq_to_ni (&nucleotides, &adapter_nucleotide))
                        .collect::<Vec<_>> ();

                    let search_psf = skips.iter ()
                        .chain (picks.iter ().rev ())
                        .find (|pick| {
                            let mm = iter::zip (pick.iter (), psf.iter ())
                                .filter (|(x,y)| x != y)
                                .count ();
                            /*
                            if mm <= 1
                            {
                                debug! ("when searching for {} found hit {} in psf mm: {}", ni_to_seq (&nucleotides, &psf), ni_to_seq (&nucleotides, &pick), mm); 
                            }
                            */
                            mm <= 1
                        });
                    if search_psf.is_none ()
                    {
                        skips.push (psf);
                        picks.push ((*p).to_vec ());
                    }
                }
            }
        }
        else
        {
            for p in &reasonable
            {
                //oi += 1;
                //debug! ("pr: {:?} {}", p, ni_to_seq (&nucleotides, &p));
                let search_pick = picks.iter ()
                    .rev ()
                    .find (|pick| {
                        let mm = iter::zip (pick.iter (), p.iter ())
                            .filter (|(x,y)| x != y)
                            .count ();
                        /*
                        if mm <= 2
                        {
                            debug! ("when searching for {} found hit {} in pick mm: {}", ni_to_seq (&nucleotides, &p), ni_to_seq (&nucleotides, &pick), mm); 
                        }
                        */
                        mm <= 2
                    });

                if search_pick.is_none ()
                {
                    picks.push ((*p).to_vec ());
                }
            }
        }
        Ok (picks)
    }
}

pub fn output_indices (adapter_five_prime: String, adapter_three_prime: String, nucleotides: &Vec<String>, index_length: u64, picks: &Vec<Vec<usize>>)
{
    for (i,p) in picks.iter ().enumerate ()
    {
        let rc = revcomp (ni_to_seq (&nucleotides, p)).to_ascii_lowercase ();
        println! ("{}\tindex_{}nt_{}\t{}{}{}", ni_to_seq (&nucleotides, p), index_length, i+1, adapter_three_prime, rc, adapter_five_prime);
    }
}

pub fn filter_indices (nucleotides: Vec<String>, index_file_paths: Vec<String>, cutoff: usize)
    -> Result<(), helper::PublicError>
{
        let mut filter = collections::HashMap::<Vec<usize>, collections::HashSet<(String,String)>>::new ();
        for ifp in index_file_paths
        {
            for (k,v) in file::fasta::all_records (&ifp)?
            {
                filter.entry (seq_to_ni (&nucleotides, &v)).or_insert (collections::HashSet::<(String,String)>::new ()).insert ( (k, ifp.clone ()) );
            }
        }
        //debug! ("filter: {:?}", filter);

        let lines = io::stdin ().lines ();
        for line in lines
        {
            //println!("got a line: {}", line.unwrap ());
            let ldata: Vec<String> = line.unwrap ().trim_end_matches ("\n").split ("\t").map (String::from).collect ();
            //debug! ("ldata: {:?}", ldata);
            if let Some (s) = ldata.first ()
            {
                let sk = seq_to_ni (&nucleotides, &s);
                if let Some ( (m, hit) ) = filter.iter ().find (|(k,_)| {
                    //iter::zip (sk.iter (), k.iter ()).all (|(x,y)| x == y)
                    let mm = iter::zip (sk.iter (), k.iter ())
                        .filter (|(x,y)| x != y)
                        .count ();
                    mm <= cutoff
                })
                {
                    eprintln! ("discarding '{}'. Matched: {} from {}", s, ni_to_seq (&nucleotides, m), hit.iter ().map (|(sn,sf)| { format! ("{}:{}", sn, sf) }).collect::<String> ());
                }
                else
                {
                    println! ("{}", ldata.join ("\t"));
                }
            }
        }
        Ok (())
}



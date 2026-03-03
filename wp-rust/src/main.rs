use std::{env,ptr,thread};
use std::collections::VecDeque;
use std::sync::mpsc;

mod bindings;
use crate::bindings::StaticData;
use crate::bindings::DynamicData;
use crate::bindings::ps_static_init;
use crate::bindings::ps_static_clear;
use crate::bindings::ps_dynamic_init;
use crate::bindings::ps_dynamic_clear;
use crate::bindings::ps_dynamic_split;
use crate::bindings::next_pol;

// Wrapper types for raw pointers. These are required to convince Rust
// to allow these to be passed between threads.
#[derive(Copy, Clone)]
struct StaticPtr{ ptr: *const StaticData }
unsafe impl Send for StaticPtr {}

struct DynamicPtr{ ptr: *mut DynamicData }
unsafe impl Send for DynamicPtr {}

// Record polynomials.
fn record(res: Vec<i64>) {
    for j in res { print!("{j} "); } 
    println!();
}

fn main() {
    // Read command line arguments and count threads.
    let args: Vec<String> = env::args().collect();
    let d0 = &args[1].parse::<i64>().unwrap();
    let q = &args[2].parse::<i64>().unwrap();
    let lead = &args[3].parse::<i64>().unwrap();
    let max_threads = 200; // Recommended value is 2n^2 where n = # of available cores
    eprintln!("Computing Weil polynomials with d = {d0}, q = {q}, lead = {lead} (threads: {max_threads})");

    let max_steps = 1000; // Maximum steps in a single call to next_pol
    let d = d0/2;
    let d_size = (d0+1) as usize;
    let d32 = d as i32;
    let mut ans_count = 0;

    // Initialize static data used by the C code, a deque (double-ended queue) of work packets,
    // and a deque of reserve packets.
    let st_data;
    let mut work: VecDeque<DynamicPtr> = VecDeque::new();
    unsafe {
        use flint3_sys::*;
        let mut temp_lead: fmpz = 0;
        let mut temp_q: fmpz = 0;
        let temp_array;
        fmpz_init(&mut temp_lead);
        fmpz_set_si(&mut temp_lead, *lead);
        fmpz_init(&mut temp_q);
        fmpz_set_si(&mut temp_q, *q);
        temp_array = _fmpz_vec_init(d+1);
        fmpz_zero(temp_array.add(d as usize));
        for i in 1..=(d as usize) { fmpz_one(temp_array.add(i)); }
        st_data = StaticPtr{ptr: ps_static_init(d32, &mut temp_q, &mut temp_lead, temp_array, -1, 0) };
        fmpz_clear(&mut temp_lead);
        fmpz_clear(&mut temp_q);
        _fmpz_vec_zero(temp_array, d);
        fmpz_set(temp_array.add(d as usize), &temp_lead);
        work.push_back(DynamicPtr{ptr: ps_dynamic_init(d32, temp_array) });
        _fmpz_vec_clear(temp_array, d+1);
    }

    // Construct deque of empty packets.
    let mut reserve: VecDeque<DynamicPtr> = VecDeque::new();
    for _ in 1..max_threads {
        reserve.push_back(DynamicPtr{ptr: unsafe { ps_dynamic_init(d32, ptr::null_mut()) }});
    }

    // Construct deque of Senders to dispatch work to threads.
    let mut dispatch: VecDeque<mpsc::Sender<DynamicPtr>> = VecDeque::new();

    // Construct channel to return answers.
    let (tx_answers, rx_answers) = mpsc::channel::<Vec<i64>>();

    // Construct channel to return packets.
    let (tx_data, rx_data) = mpsc::channel::<DynamicPtr>();

    // Run the outer loop while there is still outstanding work.
    while reserve.len() < max_threads {
        println!("{ans_count}");
        // Spawn threads with queued packets.
        for data in work.drain(..) {
            let tx_answers_clone = tx_answers.clone();
            let tx_data_clone = tx_data.clone();
            let (tx_dispatch, rx_dispatch) = mpsc::channel::<DynamicPtr>();
            dispatch.push_back(tx_dispatch);
            thread::spawn(move || {
                let _ = st_data.clone(); // Makes st_data.ptr available in the spawned thread
                let (mut flag, mut sympol);
                loop {
                    unsafe { flag = next_pol(st_data.ptr, data.ptr, max_steps); }
                    if flag == 2 { // Return a polynomial.
                        let mut ans: Vec<i64> = Vec::with_capacity(d_size);
                        sympol = unsafe { (*data.ptr).sympol };
                        for j in 0..d_size { ans.push( unsafe { *sympol.add(j) }); }
                        tx_answers_clone.send(ans).unwrap();
                    }

                    // Check the termination conditions.
                    let x = rx_dispatch.try_recv();
                    if flag == 0 || x.is_ok() {
                        let data2 = if x.is_ok() { x.unwrap() } else { rx_dispatch.recv().unwrap() };
                        unsafe { ps_dynamic_split(st_data.ptr, data.ptr, data2.ptr); }
                        tx_data_clone.send(data).unwrap();
                        tx_data_clone.send(data2).unwrap();
                        break;
                    }
                }
            });
        }

        // Collect answers. Assumes all fmpz's represent slong's.
        for res in rx_answers.try_iter() {
            record(res);
            ans_count += 1;
        }

        // Collect split and terminated processes.
        for data in rx_data.try_iter() {
            let p = data.ptr;
            if !p.is_null() { // NULL is possible here due to dummy splits
                let deq = if unsafe { (*p).flag } == 0 { &mut reserve } else { &mut work };
                deq.push_back(data);
            }
        }

        // Distribute reserve processes.
        while let Some(tx) = dispatch.pop_front() {
            if let Some(data) = reserve.pop_front() { tx.send(data).unwrap(); }
            else { // Create one dummy split to prevent deadlock.
                tx.send(DynamicPtr{ ptr: ptr::null_mut() }).unwrap();
                break;
            }
        }
    }

    // Unblock the answer channel, then collect remaining answers.
    drop(tx_answers);
    for res in rx_answers {
        record(res);
        ans_count += 1;
    }
    eprintln!("Number of polynomials found: {ans_count}");

    // Release allocated memory.
    unsafe {
       ps_static_clear(st_data.ptr);
       for data in reserve.drain(..) { ps_dynamic_clear(data.ptr); }
    }
}

use std::{env,ptr,thread};
use std::alloc::{alloc_zeroed, dealloc, Layout};
use std::collections::VecDeque;
use std::sync::mpsc;

mod bindings;
use crate::bindings::*;

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

    let max_steps = 1000; // Maximum steps in a single call to ps_next_pol
    let d = d0/2;
    let d_size = (d0+1) as usize;
    let d32 = d as i32;
    let mut ans_count = 0;

    // Initialize static data used by the C code.
    let st_data = unsafe { StaticPtr{ptr: ps_static_init(d32, q as *const i64, lead as *const i64, ptr::null_mut(), -1, 0)} };
    
    // Construct deques for loaded work packets, empty packets, and Senders to dispatch work to threads.
    let mut work: VecDeque<DynamicPtr> = VecDeque::with_capacity(max_threads);
    let mut reserve: VecDeque<DynamicPtr> = VecDeque::with_capacity(max_threads);
    let mut dispatch: VecDeque<mpsc::Sender<DynamicPtr>> = VecDeque::with_capacity(max_threads);

    let layout = Layout::array::<i64>((d+1) as usize).unwrap();
    unsafe{
        let ptr = alloc_zeroed(layout) as *mut i64;
        *ptr.add(d as usize) = *lead;
        work.push_back(DynamicPtr{ptr: ps_dynamic_init(d32, ptr) });
        dealloc(ptr as *mut u8, layout);
    }
    for _ in 1..max_threads { reserve.push_back(DynamicPtr{ptr: unsafe { ps_dynamic_init(d32, ptr::null_mut()) }}); }

    // Construct channel to return answers.
    let (tx_answers, rx_answers) = mpsc::channel::<Vec<i64>>();

    // Construct channel to return packets.
    let (tx_data, rx_data) = mpsc::channel::<DynamicPtr>();

    // Run the outer loop while there is still outstanding work.
    while reserve.len() < max_threads {
        let i = max_threads - reserve.len();
        eprintln!("Active threads: {i}; answer count: {ans_count}");
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
                    unsafe { flag = ps_next_pol(st_data.ptr, data.ptr, max_steps); }
                    if flag == 2 { // Return a polynomial.
                        let mut ans: Vec<i64> = Vec::with_capacity(d_size);
                        sympol = unsafe { (*data.ptr).sympol };
                        for j in 0..d_size { ans.push( unsafe { *sympol.add(j) }); }
                        tx_answers_clone.send(ans).unwrap();
                    }

                    // Check for termination triggers.
                    let x = rx_dispatch.try_recv();
                    if flag == 0 || x.is_ok() {
                        let data2 = if x.is_ok() { x.unwrap() } else { rx_dispatch.recv().unwrap() };
                        unsafe {
                          ps_dynamic_split(st_data.ptr, data.ptr, data2.ptr);
                          ps_cleanup(1); // Clear FLINT cache for this thread.
                        }
                        tx_data_clone.send(data).unwrap();
                        tx_data_clone.send(data2).unwrap();
                        break;
                    }
                }
            });
        }

        // Collect answers. Assumes all fmpz's represent slong's.
        let ans_count0 = ans_count;
        for res in rx_answers.try_iter() {
            record(res);
            ans_count += 1;
        }

        // Collect split and terminated processes.
        for data in rx_data.try_iter() {
            let p = data.ptr;
            if !p.is_null() { // null is possible here due to dummy splits.
                let deq = if unsafe { (*p).flag } == 0 { &mut reserve } else { &mut work };
                deq.push_back(data);
            }
        }

        // Distribute reserve processes.
        while let Some(tx) = dispatch.pop_front() {
            if let Some(data) = reserve.pop_front() { tx.send(data).unwrap(); }
            else { // If the answer queue was empty, create a dummy split to prevent deadlock.
                if ans_count > ans_count0 { dispatch.push_front(tx); }
                else { tx.send(DynamicPtr{ ptr: ptr::null_mut() }).unwrap(); }
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
       ps_cleanup(0);
    }
}

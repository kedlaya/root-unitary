use std::{env, ptr, thread};
use std::alloc::{alloc_zeroed, dealloc, Layout};
use std::sync::mpsc::{channel, Sender};
use cty::{c_int,c_long};

mod bindings;
use crate::bindings::*;

// Wrapper types for raw pointers. These are required to convince Rust
// to allow these to be passed between threads.
#[derive(Copy, Clone)]
struct StaticPtr{ ptr: *const StaticData }
unsafe impl Send for StaticPtr {}

struct DynamicPtr{ ptr: *mut DynamicData, flag: i32 }
unsafe impl Send for DynamicPtr {}

// Record polynomials taken from an iterator.
// Assumes all fmpz's returned are single limbs.
fn record(res: Vec<i64>) {
    for j in res { print!("{j} "); }
    println!();
}

fn main() {
    // Read command line arguments and count threads.
    let args: Vec<String> = env::args().collect();
    let d0 = &args[1].parse::<c_int>().unwrap();
    let q = &args[2].parse::<c_long>().unwrap();
    let lead = &args[3].parse::<c_long>().unwrap();
    let max_threads = 200; // Recommended value is n^2 where n = # of available cores
    eprintln!("Computing Weil polynomials with d = {d0}, q = {q}, lead = {lead} (threads: {max_threads})");

    let max_steps = 10000; // Maximum (weighted) steps in a single call to ps_next_pol
    let mut steps = 1;
    let d = d0/2;
    let d_size = (d0+1) as usize;
    let d32 = d as c_int;
    let answer_quota = 100000;
    let mut ans_count = 0;

    // Initialize static data used by the C code.
    let null = ptr::null_mut();
    let st_data = StaticPtr{ptr: unsafe {ps_static_init(d32, q as *const c_long, lead as *const c_long, null, -1, 0) }};
    
    // Construct deques for loaded work packets, empty packets, and Senders to dispatch work to threads.
    let mut work: Vec<DynamicPtr> = Vec::with_capacity(max_threads);
    let mut reserve: Vec<DynamicPtr> = Vec::with_capacity(max_threads);
    let mut dispatch: Vec<Sender<Option<DynamicPtr>>> = Vec::with_capacity(max_threads);

    let layout = Layout::array::<c_long>((d+1) as usize).unwrap();
    unsafe{
        let ptr = alloc_zeroed(layout) as *mut c_long;
        *ptr.add(d as usize) = *lead;
        work.push(DynamicPtr{ptr: ps_dynamic_init(d32, ptr), flag: 1});
        dealloc(ptr as *mut u8, layout);
    }
    for _ in 1..max_threads { reserve.push(DynamicPtr{ptr: unsafe { ps_dynamic_init(d32, null) }, flag: 0}); }

    // Construct channel to return answers.
    let (tx_answers, rx_answers) = channel::<Vec<c_long>>();

    // Construct channel to return packets.
    let (tx_data, rx_data) = channel::<DynamicPtr>();

    // Run the outer loop while there is still outstanding work.
    while reserve.len() < max_threads {
        eprintln!("Active workers: {}; answer count: {}", max_threads - reserve.len(), ans_count);

        // Spawn threads with queued packets.
        for mut data in work.drain(..) {
            let tx_answers_clone = tx_answers.clone();
            let tx_data_clone = tx_data.clone();
            let (tx_dispatch, rx_dispatch) = channel::<Option<DynamicPtr>>();
            dispatch.push(tx_dispatch);
            thread::spawn(move || {
                let _ = st_data.clone(); // Makes st_data.ptr available in the spawned thread
                let sympol = unsafe { (*data.ptr).sympol };
                let x = loop {
                    // Step this process forward. If we find a polynomial, return it via a channel.
                    data.flag = unsafe { ps_next_pol(st_data.ptr, data.ptr, steps) };
                    if data.flag == 2 {
                        let iter = (0..d_size).map(|x| unsafe { *sympol.add(x) });
                        tx_answers_clone.send(Vec::from_iter(iter)).unwrap();
                    }
                    // Check exit/split conditions.
                    if data.flag == 0 { break rx_dispatch.recv().unwrap(); }
                    else if let Ok(x) = rx_dispatch.try_recv() { break x; }
                };
                // Split if necessary, clean up, and return packets.
                if let Some(mut data2) = x {
                    data2.flag = unsafe { ps_dynamic_split(st_data.ptr, data.ptr, data2.ptr) };
                    tx_data_clone.send(data2).unwrap();
                }
                tx_data_clone.send(data).unwrap();
                unsafe { ps_cleanup(1); } // Clear FLINT cache for this thread.
            });
        }

        // Collect split and terminated processes.
        for data in rx_data.try_iter() {
            if data.flag == 0 { reserve.push(data); }
            else { work.push(data); }
        }

        // Distribute reserve processes to trigger splitting.
        // If reserve gets exhausted, we force one dummy split to avert deadlock.
        while let Some(tx) = dispatch.pop() {
            if steps < max_steps { steps += 1; } // Recalibrate steps after runup
            let data = reserve.pop(); // Do not unwrap! Pass as an option
            let exhausted = data.is_none(); // Must check this before data is moved
            tx.send(data).unwrap();
            if exhausted { break; }
        }

        // Collect and count answers.
        // We cap the number of answers collected at once to keep this loop circulating.
        ans_count += rx_answers.try_iter().take(answer_quota).map(|x| record(x)).count();
    }

    // Unblock the answer channel, then collect remaining answers.
    drop(tx_answers);
    ans_count += rx_answers.iter().map(|x| record(x)).count();
    eprintln!("Found {ans_count} Weil polynomials with d = {d0}, q = {q}, lead = {lead}");

    // Release allocated memory.
    unsafe {
       ps_static_clear(st_data.ptr);
       for data in reserve { ps_dynamic_clear(data.ptr); }
       ps_cleanup(0);
    }
}

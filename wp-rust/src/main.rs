use std::{env,ptr,thread,time};
use std::collections::VecDeque;
use std::mem::zeroed;
use std::sync::mpsc;
use flint3_sys::*;

mod bindings;
use crate::bindings::StaticData;
use crate::bindings::DynamicData;
use crate::bindings::ps_static_init;
use crate::bindings::ps_static_clear;
use crate::bindings::ps_dynamic_init;
use crate::bindings::ps_dynamic_clear;
use crate::bindings::ps_dynamic_split;
use crate::bindings::next_pol;

#[derive(Copy, Clone)]
struct StaticPtr{ ptr: *const StaticData }
unsafe impl Send for StaticPtr {}
unsafe impl Sync for StaticPtr {}

struct DynamicPtr { ptr: *mut DynamicData }
unsafe impl Send for DynamicPtr {}
unsafe impl Sync for DynamicPtr {}

fn main() {
    // Read command line arguments and count threads.
    let args: Vec<String> = env::args().collect();
    let d0 = &args[1].parse::<i64>().unwrap();
    let q = &args[2].parse::<i64>().unwrap();
    let lead = &args[3].parse::<i64>().unwrap();
    let max_threads = 1000;
    eprintln!("Computing Weil polynomials with d = {d0}, q = {q}, lead = {lead} (threads: {max_threads})");

    // Initialize static data used by the C code.
    let mut steps = 1;
    let max_steps = 10000;
    let d = d0/2;
    let d_size = (d0+1) as usize;
    let mut ans_count = 0;

    let mut temp_lead: fmpz = unsafe { zeroed() };
    let mut temp_q: fmpz = unsafe { zeroed() };
    let temp_array: *mut fmpz = unsafe { _fmpz_vec_init(d+1) };

    unsafe {
        fmpz_init(&mut temp_lead);
        fmpz_set_si(&mut temp_lead, *lead);
        fmpz_init(&mut temp_q);
        fmpz_set_si(&mut temp_q, *q);
        fmpz_zero(temp_array.add(d as usize));
        for i in 1..(d as usize)+1 { fmpz_one(temp_array.add(i)); }
    }

    let st_data = StaticPtr{ptr: unsafe { ps_static_init(d as i32, &mut temp_q, &mut temp_lead, temp_array, -1, 0) }};

    unsafe {
        fmpz_clear(&mut temp_lead);
        fmpz_clear(&mut temp_q);
        _fmpz_vec_zero(temp_array, d);
        fmpz_set(temp_array.add(d as usize), &temp_lead);
    }

    // Construct deque (double-ended queue) of work packets, initially with only one term.
    let mut work: VecDeque<DynamicPtr> = VecDeque::new();
    work.push_back(DynamicPtr{ptr: unsafe { ps_dynamic_init(d as i32, temp_array) }});
    unsafe { _fmpz_vec_clear(temp_array, d+1); }

    // Construct deque of empty packets.
    let mut reserve: VecDeque<DynamicPtr> = VecDeque::new();
    for _ in 1..max_threads {
        reserve.push_back(DynamicPtr{ptr: unsafe { ps_dynamic_init(d as i32, ptr::null_mut()) }});
    }

    // Construct deque of Senders to dispatch work to threads.
    let mut dispatch: VecDeque<mpsc::Sender<DynamicPtr>> = VecDeque::new();

    // Construct channel for threads to return answers.
    let (tx_answers, rx_answers) = mpsc::channel::<Vec<i64>>();

    // Construct channel for threads to return packets.
    let (tx_data, rx_data) = mpsc::channel::<DynamicPtr>();

    // Run the outer loop while there is still outstanding work.
    while reserve.len() < max_threads {

        // Spawn threads with queued packets.
        for data in work.drain(..) {
            if steps < max_steps { steps *= 2; } // The first few threads should poll more frequently.
            let tx_answers_clone = tx_answers.clone();
            let tx_data_clone = tx_data.clone();
            let (tx_dispatch, rx_dispatch) = mpsc::channel::<DynamicPtr>();
            dispatch.push_back(tx_dispatch);
            thread::spawn(move || {
                let _ = st_data.clone(); // Makes st_data available in the spawned thread
                let mut flag;
                let pause = time::Duration::from_millis(10);
                let mut sympol;
                loop {
                    // Do some work. If we find a polynomial, send it back.
                    // If there is no more work, pause to allow other threads to catch up.
                    unsafe { flag = next_pol(st_data.ptr, data.ptr, steps); }
                    if flag == 0 { thread::sleep(pause); }
                    else if flag == 2 {
                        let mut ans: Vec<i64> = Vec::with_capacity(d_size);
                        sympol = unsafe { (*data.ptr).sympol };
                        for j in 0..d_size { ans.push( unsafe { *sympol.add(j) }); }
                        tx_answers_clone.send(ans).unwrap();
                    }

                    // When we get a dispatch, split work and exit the loop.
                    if let Ok(data2) = rx_dispatch.try_recv() {
                        unsafe { ps_dynamic_split(st_data.ptr, data.ptr, data2.ptr); }
                        tx_data_clone.send(data).unwrap();
                        tx_data_clone.send(data2).unwrap();
                        break;
                    }
               }
            });
        }

        // Collect answers.
        for res in rx_answers.try_iter() {
            for j in res { print!("{j} "); } println!();
            ans_count += 1;
        }

        // Collect split and terminated processes.
        for data in rx_data.try_iter() {
            let flag = unsafe { (*data.ptr).flag };
            if flag == 0 { reserve.push_back(data); }
            else { work.push_back(data); }
        }

        // Create a dummy split to avoid blocking.
        if reserve.len() == 0 {
           let tx = dispatch.pop_front().unwrap(); 
           let data = DynamicPtr{ ptr: ptr::null_mut() };
           tx.send(data).unwrap();
        }

        // Distribute reserve processes.
        while dispatch.len() > 0 && reserve.len() > 0 {
           let tx = dispatch.pop_front().unwrap(); // Can only be used *once*
           let data = reserve.pop_front().unwrap();
           tx.send(data).unwrap();
        }
    }

    // Collect remaining answers.
    drop(tx_answers);
    for res in rx_answers.iter() {
        for j in res { print!("{j} "); } println!();
        ans_count += 1;
    }
    eprintln!("Number of polynomials found: {ans_count}");

    // Release allocated memory.
    unsafe {
       ps_static_clear(st_data.ptr);
       for data in reserve.drain(..) { ps_dynamic_clear(data.ptr); }
    }
}

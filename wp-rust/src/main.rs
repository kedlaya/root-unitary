use std::collections::VecDeque;
use std::env;
use std::mem::zeroed;
use std::ptr;
use std::sync::mpsc;
use std::thread;
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
    // Read command line arguments
    let args: Vec<String> = env::args().collect();
    let d0 = &args[1].parse::<i64>().unwrap();
    let q = &args[2].parse::<i64>().unwrap();
    let lead = &args[3].parse::<i64>().unwrap();

    // Initialize stuff
    let max_threads = 1000;
    let mut steps = 1;
    let max_steps = 10000;
    let d = d0/2;
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
        for i in 1..d+1 { fmpz_one(temp_array.add(i as usize)); }
    }
    
    let st_data = StaticPtr{ptr: unsafe { ps_static_init(d as i32, &mut temp_q, &mut temp_lead, temp_array, -1, 0) }};
    println!("d = {d0}, q = {q}, lead = {lead}");
    
    unsafe {
        fmpz_clear(&mut temp_lead);
        fmpz_clear(&mut temp_q);
        _fmpz_vec_zero(temp_array, d);
        fmpz_set(temp_array.add(d as usize), &temp_lead);
    }
    
    let mut dispatch: VecDeque<mpsc::Sender<DynamicPtr>> = VecDeque::new();
    let mut work: VecDeque<DynamicPtr> = VecDeque::new();
    let mut reserve: VecDeque<DynamicPtr> = VecDeque::new();
    work.push_back(DynamicPtr{ptr: unsafe { ps_dynamic_init(d as i32, temp_array) }});
    unsafe { _fmpz_vec_clear(temp_array, d+1); }

    for _ in 2..max_threads+1 { 
        reserve.push_back(DynamicPtr{ptr: unsafe { ps_dynamic_init(d as i32, ptr::null_mut()) }});
    }    

    let (tx_answers, rx_answers) = mpsc::channel::<i32>();
    let (tx_data, rx_data) = mpsc::channel::<DynamicPtr>();

    // Run the outer loop while there is still outstanding work.
    while reserve.len() < max_threads {
        let i = reserve.len();
        if steps < max_steps && max_threads - reserve.len() > (steps as usize) { steps *= 2; }
        print!("{i}, {steps}, {ans_count}\n");
        
        // Launch queued processes.
        for data in work.drain(..) {
            let tx_answers_clone = tx_answers.clone();
            let tx_data_clone = tx_data.clone();
            let (tx_dispatch, rx_dispatch) = mpsc::channel::<DynamicPtr>();
            dispatch.push_back(tx_dispatch); // Can only be used *once*
            thread::spawn(move || {
                let st_data_local = st_data.clone();
                let mut flag;
                loop {
                    // Do some work.
                    unsafe { flag = next_pol(st_data_local.ptr, data.ptr, steps); }
    
                    // Report an answer if we have one.
                    if flag == 2 as i32 {
                       tx_answers_clone.send(1).unwrap();
                    }
                    
                    let x = rx_dispatch.try_recv();
                    if let Ok(data2) = x { // Split work and terminate.
                        unsafe { ps_dynamic_split(st_data_local.ptr, data.ptr, data2.ptr); }
                        tx_data_clone.send(data).unwrap();
                        tx_data_clone.send(data2).unwrap();
                        break;
                    }
               }
            });
        }
        
        // Collect answers.
        for _res in rx_answers.try_iter() {
            ans_count = ans_count + 1;
        }
        
        // Collect split and terminated processes.
        for data in rx_data.try_iter() {
            let flag = unsafe { (*(data.ptr)).flag };
            if flag == (0 as i32) { reserve.push_back(data); }
            else { work.push_back(data); }
        }

        // Create a dummy split to avoid blocking.
        if reserve.len() == 0 {
           let tx = dispatch.pop_front().unwrap(); 
           let data = DynamicPtr{ptr: ptr::null_mut() };
           tx.send(data).unwrap();
        }
        
        // Distribute reserve processes;
        while dispatch.len() > 0 && reserve.len() > 0 {
           let tx = dispatch.pop_front().unwrap(); // Can only be used *once*
           let data = reserve.pop_front().unwrap();
           tx.send(data).unwrap();
        }

    }
    
    // Collect remaining answers.
    drop(tx_answers);
    for _res in rx_answers.iter() {
        ans_count = ans_count + 1;
    }
    
    // Release allocated memory.
    unsafe {
       ps_static_clear(st_data.ptr);
       for data in reserve.drain(..) { ps_dynamic_clear(data.ptr); }
    }
}

use cty::c_int;
use cty::c_long;
use flint3_sys::*; // Needed for linking

#[repr(C)]
pub struct StaticData {
    pub d: c_int,
    pub force_squarefree: c_int,
    pub node_limit: c_long,
    pub q: *const i64,
    pub modlist: *const i64,
    pub binom_mat: *const i64,
    pub sum_mats: *const i64,
    pub eval_pm2_mats: *const i64
}

#[repr(C)]
pub struct DynamicData {
    pub d: cty::c_int,
    pub n: cty::c_int,
    pub ascend: cty::c_int,
    pub flag: cty::c_int,
    pub node_count: cty::c_long,
    pub pol: *mut i64,
    pub sympol: *mut i64,
    pub upper: *mut i64,
    pub power_sums_num: *mut i64,
    pub hankel_dets: *mut i64,
    pub w: *mut i64,
    pub wlen: cty::c_long
}

unsafe extern "C" {
    pub fn ps_static_init(
        d: c_int,
        q: *const i64,
        lead: *const i64,
        modlist: *const i64,
        node_limit: c_long,
        force_squarefree: c_int,
    ) -> *const StaticData;
}

unsafe extern "C" {
    pub fn ps_dynamic_init(
        d: c_int,
        coefflist: *mut i64,
    ) ->  *mut DynamicData;
}

unsafe extern "C" {
    pub fn ps_static_clear(st_data: *const StaticData);
}

unsafe extern "C" {
    pub fn ps_dynamic_clear(dy_data: *mut DynamicData);
}

unsafe extern "C" {
    pub fn ps_cleanup(n: c_int);
}

unsafe extern "C" {
    pub fn ps_dynamic_split(
        st_data: *const StaticData,
        dy_data: *mut DynamicData,
        dy_data2: *mut DynamicData,
    ) -> c_int;
}

unsafe extern "C" {
    pub fn ps_next_pol(
        st_data: *const StaticData,
        dy_data: *mut DynamicData,
        max_steps: c_int,
    ) -> c_int;
}



use cty::c_int;
use cty::c_long;

#[allow(unused)]
use flint3_sys::*; // Needed for linking

#[repr(C)]
pub struct StaticData {
    pub d: c_int,
    pub force_squarefree: c_int,
    pub node_limit: c_long,
    pub q: *const c_long,
    pub q_sqrt: *const c_long,
    pub c0: *const c_long,
    pub c1: *const c_long,
    pub q_is_1: c_int,
    pub q_is_square: c_int,
    pub lead_is_1: c_int,    
    pub num_constraints: c_int,
    pub constraint_lens: *const c_int,
    pub modlist: *const c_long,
    pub binom_mat: *const c_long,
    pub sum_mats: *const c_long,
    pub eval_pm2_mats: *const c_long,
    pub ranges: *const c_long,
    pub lead_pows: *const c_long,
    pub constraints: *const c_long,
    pub pol_to_sym: fmpz_mat_t
}

#[repr(C)]
pub struct DynamicData {
    pub d: cty::c_int,
    pub n: cty::c_int,
    pub ascend: cty::c_int,
    pub flag: cty::c_int,
    pub node_count: cty::c_long,
    pub pol: *mut c_long,
    pub sympol: *mut c_long,
    pub upper: *mut c_long,
    pub power_sums_num: *mut c_long,
    pub hankel_dets: *mut c_long,
    pub w: *mut c_long,
    pub wlen: cty::c_long
}

unsafe extern "C" {
    pub fn ps_static_init(
        d: c_int,
        q: *const c_long,
        lead: *const c_long,
        modlist: *const c_long,
        num_constraints: c_int,
        constraints: *const c_long,
        node_limit: c_long,
        force_squarefree: c_int,
    ) -> *const StaticData;
}

unsafe extern "C" {
    pub fn ps_dynamic_init(
        d: c_int,
        coefflist: *mut c_long,
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
    pub fn ascend_step_forward(
        st_data: *const StaticData,
        dy_data: *mut DynamicData,
        n: c_int,
    ) -> c_int;
}

unsafe extern "C" {
    pub fn reduce_range_from_rolle(
        st_data: *const StaticData,
        dy_data: *mut DynamicData,
        n: c_int,
    ) -> c_int;
}

unsafe extern "C" {
    pub fn set_range_from_power_sums(
        st_data: *const StaticData,
        dy_data: *mut DynamicData,
        n: c_int,
    ) -> c_int;
}

unsafe extern "C" {
    pub fn reciprocal_transform(
        st_data: *const StaticData,
        dy_data: *mut DynamicData,
    );
}

unsafe extern "C" {
    pub fn ps_dynamic_split(
        st_data: *const StaticData,
        dy_data: *mut DynamicData,
        dy_data2: *mut DynamicData,
    ) -> c_int;
}


use cty::c_int;
use cty::c_long;
use flint3_sys::*;

#[repr(C)]
pub struct StaticData {
    pub d: c_int,
    pub force_squarefree: c_int,
    pub node_limit: c_long,
    pub q: fmpz_t,
    pub modlist: *mut fmpz,
    pub binom_mat: *mut fmpz,
    pub sum_mats: *mut fmpz,
    pub eval_pm2_mats: *mut fmpz
}

#[repr(C)]
pub struct DynamicData {
    pub d: cty::c_int,
    pub n: cty::c_int,
    pub ascend: cty::c_int,
    pub flag: cty::c_int,
    pub node_count: cty::c_long,
    pub pol: *mut fmpz,
    pub sympol: *mut fmpz,
    pub upper: *mut fmpz,
    pub power_sums_num: *mut fmpz,
    pub hankel_dets: *mut fmpz,
    pub w: *mut fmpz,
    pub wlen: cty::c_long
}

unsafe extern "C" {
    pub fn ps_static_init(
        d: c_int,
        q: *mut fmpz,
        lead: *mut fmpz,
        modlist: *mut fmpz,
        node_limit: c_long,
        force_squarefree: c_int,
    ) -> *const StaticData;
}

unsafe extern "C" {
    pub fn ps_dynamic_init(
        d: c_int,
        coefflist: *mut fmpz,
    ) ->  *mut DynamicData;
}

unsafe extern "C" {
    pub fn ps_static_clear(st_data: *const StaticData);
}

unsafe extern "C" {
    pub fn ps_dynamic_clear(dy_data: *mut DynamicData);
}

unsafe extern "C" {
    pub fn ps_dynamic_split(
        st_data: *const StaticData,
        dy_data: *mut DynamicData,
        dy_data2: *mut DynamicData,
    ) -> c_int;
}

unsafe extern "C" {
    pub fn next_pol(
        st_data: *const StaticData,
        dy_data: *mut DynamicData,
        max_steps: c_int,
    ) -> c_int;
}



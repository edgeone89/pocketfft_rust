/*
 * This file is Rust rewrite of pocketfft.c from https://gitlab.mpcdf.mpg.de/mtr/pocketfft
 * Licensed under a 3-clause BSD style license - see LICENSE
 */

use core::ffi::c_void;
use libc::free;
use libc::malloc;
use libc::memcpy;
use std::mem::forget;
use std::mem::size_of;
use std::ptr::null;
use std::ptr::null_mut;
use std::slice::from_raw_parts;
use std::slice::from_raw_parts_mut;

fn my_sincosm1pi(a: f64, res: &mut [f64]) {
    let s = a * a;
    /* Approximate cos(pi*x)-1 for x in [-0.25,0.25] */
    let r: f64 = -1.0369917389758117e-4;
    let r1: f64 = r.mul_add(s, 1.9294935641298806e-3);

    let r: f64 = r1.mul_add(s, -2.5806887942825395e-2);

    let r1: f64 = r.mul_add(s, 2.3533063028328211e-1);

    let r: f64 = r1.mul_add(s, -1.3352627688538006e+0);

    let r1: f64 = r.mul_add(s, 4.0587121264167623e+0);

    let r: f64 = r1.mul_add(s, -4.9348022005446790e+0);

    let c = r * s;

    /* Approximate sin(pi*x) for x in [-0.25,0.25] */
    let r: f64 = 4.6151442520157035e-4;
    let r1: f64 = r.mul_add(s, -7.3700183130883555e-3);

    let r: f64 = r1.mul_add(s, 8.2145868949323936e-2);

    let r1: f64 = r.mul_add(s, -5.9926452893214921e-1);

    let r: f64 = r1.mul_add(s, 2.5501640398732688e+0);

    let r1: f64 = r.mul_add(s, -5.1677127800499516e+0);

    let r: f64 = r1 * s * a;
    let s1 = a.mul_add(3.1415926535897931e+0, r);
    res[0] = c;
    res[1] = s1;
}

fn calc_first_octant(den: usize, res: &mut [f64]) {
    let n = (den + 4) >> 3;
    if n != 0 && n != 1 {
        res[0] = 1.0;
        res[1] = 0.0;
        let l1 = (n as f64).sqrt() as usize;

        for i in 1..l1 {
            let a = (2.0 * i as f64) / (den as f64);
            my_sincosm1pi(a, &mut res[2 * i..]);
        }

        let mut start = l1;

        while start < n {
            let mut cs: [f64; 2] = [0.0; 2];
            let a = (2.0 * (start as f64)) / (den as f64);
            my_sincosm1pi(a, &mut cs);
            res[2 * start] = cs[0] + 1.0;
            res[2 * start + 1] = cs[1];
            let mut end = l1;
            if start + end > n {
                end = n - start;
            }
            for k in 1..end {
                let mut csx: [f64; 2] = [0.0; 2];
                csx[0] = res[2 * k];
                csx[1] = res[2 * k + 1];
                res[2 * (start + k)] = ((cs[0] * csx[0] - cs[1] * csx[1] + cs[0]) + csx[0]) + 1.0;
                res[2 * (start + k) + 1] = (cs[0] * csx[1] + cs[1] * csx[0]) + cs[1] + csx[1];
            }
            start += l1;
        }
        for i in 1..l1 {
            res[2 * i] += 1.0;
        }
    }
}

fn calc_first_quadrant(n: usize, res: &mut [f64]) {
    //let mut p = &mut res[n..];
    //calc_first_octant(n<<1, &mut p);
    calc_first_octant(n << 1, &mut res[n..]);
    let ndone = (n + 2) >> 2;
    let mut i: usize = 0;
    let mut idx1: usize = 0;
    let mut idx2: usize = 2 * ndone - 2;
    while (i + 1) < ndone {
        //res[idx1]   = p[2*i];
        res[idx1] = res[n + 2 * i];
        //res[idx1+1] = p[2*i+1];
        res[idx1 + 1] = res[n + 2 * i + 1];
        //res[idx2]   = p[2*i+3];
        res[idx2] = res[n + 2 * i + 3];
        //res[idx2+1] = p[2*i+2];
        res[idx2 + 1] = res[n + 2 * i + 2];
        i += 2;
        idx1 += 2;
        idx2 -= 2;
    }
    if i != ndone {
        //res[idx1  ] = p[2*i];
        res[idx1] = res[n + 2 * i];
        //res[idx1+1] = p[2*i+1];
        res[idx1 + 1] = res[n + 2 * i + 1];
    }
}

fn calc_first_half(n: usize, res: &mut [f64]) {
    let ndone = ((n + 1) >> 1) as usize;
    //double * p = res+n-1;
    //calc_first_octant(n<<2, p);
    calc_first_octant(n << 2, &mut res[(n - 1)..]);
    let mut i4 = 0;
    let mut iN = n;
    let mut i = 0;
    while i4 <= iN - i4
    /* octant 0 */
    {
        //res[2*i] = p[2*i4];
        res[2 * i] = res[n - 1 + 2 * i4];
        //res[2*i+1] = p[2*i4+1];
        res[2 * i + 1] = res[n - 1 + 2 * i4 + 1];
        i += 1;
        i4 += 4;
    }
    while i4 - iN <= 0
    /* octant 1 */
    {
        let xm = iN - i4;
        //res[2*i] = p[2*xm+1];
        res[2 * i] = res[n - 1 + 2 * xm + 1];
        //res[2*i+1] = p[2*xm];
        res[2 * i + 1] = res[n - 1 + 2 * xm];
        i += 1;
        i4 += 4;
    }
    while i4 <= 3 * iN - i4
    /* octant 2*/
    {
        let xm = i4 - iN;
        //res[2*i] = -p[2*xm+1];
        res[2 * i] = -1.0 * res[n - 1 + 2 * xm + 1];
        //res[2*i+1] = p[2*xm];
        res[2 * i + 1] = res[n - 1 + 2 * xm];
        i += 1;
        i4 += 4;
    }
    while i < ndone
    /* octant 3 */
    {
        let xm = 2 * iN - i4;
        //res[2*i] = -p[2*xm];
        res[2 * i] = -1.0 * res[n - 1 + 2 * xm];
        //res[2*i+1] = p[2*xm+1];
        res[2 * i + 1] = res[n - 1 + 2 * xm + 1];
        i += 1;
        i4 += 4;
    }
}

fn fill_first_quadrant(n: usize, res: &mut [f64]) {
    let hsqt2 = 0.707106781186547524400844362104849;
    let quart = n >> 2;
    if (n & 7) == 0 {
        res[quart + 1] = hsqt2;
        res[quart] = hsqt2
    }
    let mut i = 2;
    let mut j = 2 * quart - 2;
    while i < quart {
        res[j] = res[i + 1];
        res[j + 1] = res[i];
        i += 2;
        j -= 2;
    }
}

fn fill_first_half(n: usize, res: &mut [f64]) {
    let half = n >> 1;
    if (n & 3) == 0 {
        let mut i = 0;
        while i < half {
            res[i + half] = -res[i + 1];
            res[i + half + 1] = res[i];
            i += 2;
        }
    } else {
        let mut i = 2;
        let mut j = 2 * half - 2;
        while i < half {
            res[j] = -res[i];
            res[j + 1] = res[i + 1];
            i += 2;
            j -= 2;
        }
    }
}

fn fill_second_half(n: usize, res: &mut [f64]) {
    if (n & 1) == 0 {
        for i in 0..n {
            res[i + n] = -res[i];
        }
    } else {
        let mut i = 2;
        let mut j = 2 * n - 2;
        while i < n {
            res[j] = res[i];
            res[j + 1] = -res[i + 1];
            i += 2;
            j -= 2;
        }
    }
}

fn sincos_2pibyn_half(n: usize, res: &mut [f64]) {
    if (n & 3) == 0 {
        calc_first_octant(n, res);
        fill_first_quadrant(n, res);
        fill_first_half(n, res);
    } else if (n & 1) == 0 {
        calc_first_quadrant(n, res);
        fill_first_half(n, res);
    } else {
        calc_first_half(n, res);
    }
}

fn sincos_2pibyn(n: usize, res: &mut [f64]) {
    sincos_2pibyn_half(n, res);
    fill_second_half(n, res);
}

fn largest_prime_factor(n: usize) -> usize {
    let mut n_temp = n;
    let mut res: usize = 1;
    let mut tmp = n_temp >> 1;
    while (tmp << 1) == n_temp {
        res = 2;
        n_temp = tmp;
        tmp = n_temp >> 1;
    }

    let mut limit = ((n_temp as f64) + 0.01).sqrt() as usize;
    let mut x = 3;
    while x <= limit {
        tmp = n_temp / x;
        while tmp * x == n_temp {
            res = x;
            n_temp = tmp;
            limit = ((n_temp as f64) + 0.01).sqrt() as usize;
        }
        x += 2;
    }
    if n_temp > 1 {
        res = n_temp;
    }

    return res;
}

fn cost_guess(n: usize) -> f64 {
    let lfp: f64 = 1.1; // penalty for non-hardcoded larger factors
    let mut n_temp = n;
    let ni = n;
    let mut result: f64 = 0.0;
    let mut tmp = n_temp >> 1;
    while (tmp << 1) == n_temp {
        result += 2.0;
        n_temp = tmp;
        tmp = n_temp >> 1;
    }

    let mut limit = ((n_temp as f64) + 0.01).sqrt() as usize;
    let mut x: usize = 3;
    while x <= limit {
        tmp = n_temp / x;
        while (tmp * x) == n_temp {
            if x <= 5 {
                result += x as f64;
            } else {
                result += lfp * (x as f64);
            }
            // penalize larger prime factors
            n_temp = tmp;
            limit = ((n_temp as f64) + 0.01).sqrt() as usize;
            tmp = n_temp / x;
        }
        x += 2;
    }
    if n_temp > 1 {
        if n_temp <= 5 {
            result += (n_temp as f64);
        } else {
            result += lfp * (n_temp as f64);
        }
    }

    return result * (ni as f64);
}

/* returns the smallest composite of 2, 3, 5, 7 and 11 which is >= n */
fn good_size(n: usize) -> usize {
    if n <= 6 {
        return n;
    }

    let mut bestfac: usize = 2 * n;
    let mut f2: usize = 1;
    while f2 < bestfac {
        let mut f23 = f2;
        while f23 < bestfac {
            let mut f235 = f23;
            while f235 < bestfac {
                let mut f2357 = f235;
                while f2357 < bestfac {
                    let mut f235711 = f2357;
                    while f235711 < bestfac {
                        if f235711 >= n {
                            bestfac = f235711;
                        }
                        f235711 *= 11;
                    }
                    f2357 *= 7;
                }
                f235 *= 5;
            }
            f23 *= 3;
        }
        f2 *= 2;
    }
    return bestfac;
}

/*#[repr(C)]
struct cmplx {
    r: f64,
    i: f64,
}*/

const NFCT: usize = 25;

#[repr(C)]
struct cfftp_fctdata {
    fct: usize,
    tw: Vec<f64>,
    tws: Vec<f64>,
}

#[repr(C)]
struct cfftp_plan_i {
    length: usize,
    nfct: usize,
    mem: Vec<f64>,
    fct: Vec<cfftp_fctdata>,
}

pub type cfftp_plan = *mut cfftp_plan_i;

fn PMC(a: &mut [f64], b: &mut [f64], c: &[f64], d: &[f64]) {
    a[0] = c[0] + d[0];
    a[1] = c[1] + d[1];
    b[0] = c[0] - d[0];
    b[1] = c[1] - d[1];
}

fn ADDC(a: &mut [f64], b: &[f64], c: &[f64]) {
    a[0] = b[0] + c[0];
    a[1] = b[1] + c[1];
}

fn SCALEC(a: &mut [f64], b: f64) {
    a[0] *= b;
    a[1] *= b;
}

fn ROT90(a: &mut [f64]) {
    let mut tmp_: f64 = a[0];
    a[0] = -a[1];
    a[1] = tmp_;
}

fn ROTM90(a: &mut [f64]) {
    let mut tmp_: f64 = -a[0];
    a[0] = a[1];
    a[1] = tmp_;
}

fn A_EQ_B_MUL_C(a: &mut [f64], b: &[f64], c: &[f64]) {
    a[0] = b[0] * c[0] - b[1] * c[1];
    a[1] = b[0] * c[1] + b[1] * c[0];
}

fn A_EQ_CB_MUL_C(a: &mut [f64], b: &[f64], c: &[f64]) {
    a[0] = b[0] * c[0] + b[1] * c[1];
    a[1] = b[0] * c[1] - b[1] * c[0];
}

fn PMSIGNC(a: &mut [f64], b: &mut [f64], c: &[f64], d: &[f64], sign: f64) {
    a[0] = c[0] + sign * d[0];
    a[1] = c[1] + sign * d[1];
    b[0] = c[0] - sign * d[0];
    b[1] = c[1] - sign * d[1];
}

fn MULPMSIGNC(a: &mut [f64], b: &[f64], c: &[f64], sign: f64) {
    a[0] = b[0] * c[0] - sign * b[1] * c[1];
    a[1] = b[0] * c[1] + sign * b[1] * c[0];
}

fn MULPMSIGNCEQ(a: &mut [f64], b: &mut [f64], sign: f64) {
    let xtmp = a[0];
    a[0] = b[0] * a[0] - sign * b[1] * a[1];
    a[1] = b[0] * a[1] + sign * b[1] * xtmp;
}

fn pass2b(ido: usize, l1: usize, cc: &[f64], ch: &mut [f64], wa: &[f64]) {
    let cdim: usize = 2;

    if ido == 1 {
        for k in 0..l1 {
            //PMC (CH(0,k,0),CH(0,k,1),CC(0,0,k),CC(0,1,k))
            ch[(0) + ido * ((k) + l1 * (0))] =
                cc[(0) + ido * ((0) + cdim * (k))] + cc[(0) + ido * ((1) + cdim * (k))];
            ch[(0) + ido * ((k) + l1 * (0)) + 1] =
                cc[(0) + ido * ((0) + cdim * (k)) + 1] + cc[(0) + ido * ((1) + cdim * (k)) + 1];
            ch[(0) + ido * ((k) + l1 * (1))] =
                cc[(0) + ido * ((0) + cdim * (k))] - cc[(0) + ido * ((1) + cdim * (k))];
            ch[(0) + ido * ((k) + l1 * (1)) + 1] =
                cc[(0) + ido * ((0) + cdim * (k)) + 1] - cc[(0) + ido * ((1) + cdim * (k)) + 1];
        }
    } else {
        for k in 0..l1 {
            // CH(a,b,c) ch[(a)+ido*((b)+l1*(c))]
            // CC(a,b,c) cc[(a)+ido*((b)+cdim*(c))]
            //PMC (CH(0,k,0),CH(0,k,1),CC(0,0,k),CC(0,1,k))
            ch[(0) + ido * ((k) + l1 * (0))] =
                cc[(0) + ido * ((0) + cdim * (k))] + cc[(0) + ido * ((1) + cdim * (k))];
            ch[(0) + ido * ((k) + l1 * (0)) + 1] =
                cc[(0) + ido * ((0) + cdim * (k)) + 1] + cc[(0) + ido * ((1) + cdim * (k)) + 1];

            ch[(0) + ido * ((k) + l1 * (1))] =
                cc[(0) + ido * ((0) + cdim * (k))] - cc[(0) + ido * ((1) + cdim * (k))];
            ch[(0) + ido * ((k) + l1 * (1)) + 1] =
                cc[(0) + ido * ((0) + cdim * (k)) + 1] - cc[(0) + ido * ((1) + cdim * (k)) + 1];
            for i in 1..ido {
                let mut t: [f64; 2] = [0.0; 2]; //cmplx { r: 0.0, i: 0.0 };
                                                //PMC (CH(i,k,0),t,CC(i,0,k),CC(i,1,k))
                ch[(i) + ido * ((k) + l1 * (0))] =
                    cc[(i) + ido * ((0) + cdim * (k))] + cc[(i) + ido * ((1) + cdim * (k))];
                ch[(i) + ido * ((k) + l1 * (0)) + 1] =
                    cc[(i) + ido * ((0) + cdim * (k)) + 1] + cc[(i) + ido * ((1) + cdim * (k)) + 1];
                t[0] = cc[(i) + ido * ((0) + cdim * (k))] - cc[(i) + ido * ((1) + cdim * (k))];
                t[1] =
                    cc[(i) + ido * ((0) + cdim * (k)) + 1] - cc[(i) + ido * ((1) + cdim * (k)) + 1];

                //A_EQ_B_MUL_C(a,b,c) { a.r=b.r*c.r-b.i*c.i; a.i=b.r*c.i+b.i*c.r; }
                //A_EQ_B_MUL_C (CH(i,k,1),WA(0,i),t)
                // WA(x,i) wa[(i)-1+(x)*(ido-1)]
                ch[(i) + ido * ((k) + l1 * (1))] =
                    wa[(i) - 1 + (0) * (ido - 1)] * t[0] - wa[(i) - 1 + (0) * (ido - 1) + 1] * t[1];
                ch[(i) + ido * ((k) + l1 * (1)) + 1] =
                    wa[(i) - 1 + (0) * (ido - 1)] * t[1] + wa[(i) - 1 + (0) * (ido - 1) + 1] * t[0];
            }
        }
    }
}

fn pass2f(ido: usize, l1: usize, cc: &[f64], ch: &mut [f64], wa: &[f64]) {
    let cdim: usize = 2;

    if ido == 1 {
        for k in 0..l1 {
            //PMC (CH(0,k,0),CH(0,k,1),CC(0,0,k),CC(0,1,k))
            // ch[(a)+ido*((b)+l1*(c))]
            // cc[(a)+ido*((b)+cdim*(c))]
            ch[(0) + ido * ((k) + l1 * (0))] =
                cc[(0) + ido * ((0) + cdim * (k))] + cc[(0) + ido * ((1) + cdim * (k))];
            ch[(0) + ido * ((k) + l1 * (0)) + 1] =
                cc[(0) + ido * ((0) + cdim * (k)) + 1] + cc[(0) + ido * ((1) + cdim * (k)) + 1];
            ch[(0) + ido * ((k) + l1 * (1))] =
                cc[(0) + ido * ((0) + cdim * (k))] - cc[(0) + ido * ((1) + cdim * (k))];
            ch[(0) + ido * ((k) + l1 * (1)) + 1] =
                cc[(0) + ido * ((0) + cdim * (k)) + 1] - cc[(0) + ido * ((1) + cdim * (k)) + 1];
        }
    } else {
        for k in 0..l1 {
            //PMC (CH(0,k,0),CH(0,k,1),CC(0,0,k),CC(0,1,k))
            //{ a.r=c.r+d.r; a.i=c.i+d.i; b.r=c.r-d.r; b.i=c.i-d.i; }
            ch[(0) + ido * ((k) + l1 * (0))] =
                ch[(0) + ido * ((0) + l1 * (k))] + ch[(0) + ido * ((1) + l1 * (k))];
            ch[(0) + ido * ((k) + l1 * (0)) + 1] =
                ch[(0) + ido * ((0) + l1 * (k)) + 1] + ch[(0) + ido * ((1) + l1 * (k)) + 1];
            ch[(0) + ido * ((k) + l1 * (1))] =
                cc[(0) + ido * ((0) + cdim * (k))] - cc[(0) + ido * ((1) + cdim * (k))];
            ch[(0) + ido * ((k) + l1 * (1)) + 1] =
                cc[(0) + ido * ((0) + cdim * (k)) + 1] - cc[(0) + ido * ((1) + cdim * (k)) + 1];
            for i in 1..ido {
                let mut t: [f64; 2] = [0.0; 2]; //cmplx { r: 0.0, i: 0.0 };
                                                //PMC (CH(i,k,0),t,CC(i,0,k),CC(i,1,k))
                ch[(i) + ido * ((k) + l1 * (0))] =
                    cc[(i) + ido * ((0) + cdim * (k))] + cc[(i) + ido * ((1) + cdim * (k))];
                ch[(i) + ido * ((k) + l1 * (0)) + 1] =
                    cc[(i) + ido * ((0) + cdim * (k)) + 1] + cc[(i) + ido * ((1) + cdim * (k)) + 1];
                t[0] = cc[(i) + ido * ((0) + cdim * (k))] - cc[(i) + ido * ((1) + cdim * (k))];
                t[1] =
                    cc[(i) + ido * ((0) + cdim * (k)) + 1] - cc[(i) + ido * ((1) + cdim * (k)) + 1];

                //A_EQ_CB_MUL_C (CH(i,k,1),WA(0,i),t)
                //A_EQ_B_MUL_C(a,b,c) { a.r=b.r*c.r-b.i*c.i; a.i=b.r*c.i+b.i*c.r; }
                ch[(i) + ido * ((k) + l1 * (1))] =
                    wa[(i) - 1 + (0) * (ido - 1)] * t[0] - wa[(i) - 1 + (0) * (ido - 1) + 1] * t[1];
                ch[(i) + ido * ((k) + l1 * (1)) + 1] =
                    wa[(i) - 1 + (0) * (ido - 1)] * t[1] + wa[(i) - 1 + (0) * (ido - 1) + 1] * t[0];
            }
        }
    }
}

fn pass3b(ido: usize, l1: usize, cc: &[f64], ch: &mut [f64], wa: &[f64]) {
    let cdim: usize = 3;
    let tw1r: f64 = -0.5;
    let tw1i: f64 = 0.86602540378443864676;

    if ido == 1 {
        for k in 0..l1 {
            //PREP3(0)
            let mut t0: [f64; 2] = [0.0; 2]; //&cc[(0) + ido * ((0) + cdim * (k))];
            t0[0] = cc[(0) + ido * ((0) + cdim * (k))];
            t0[1] = cc[(0) + ido * ((0) + cdim * (k)) + 1];
            let mut t1: [f64; 2] = [0.0; 2]; //cmplx { r: 0.0, i: 0.0 };
            let mut t2: [f64; 2] = [0.0; 2]; //cmplx { r: 0.0, i: 0.0 };
                                             //PMC (t1,t2,CC(idx,1,k),CC(idx,2,k));
            t1[0] = cc[(0) + ido * ((1) + cdim * (k))] + cc[(0) + ido * ((2) + cdim * (k))];
            t1[1] = cc[(0) + ido * ((1) + cdim * (k)) + 1] + cc[(0) + ido * ((2) + cdim * (k)) + 1];
            t2[0] = cc[(0) + ido * ((1) + cdim * (k))] - cc[(0) + ido * ((2) + cdim * (k))];
            t2[1] = cc[(0) + ido * ((1) + cdim * (k)) + 1] - cc[(0) + ido * ((2) + cdim * (k)) + 1];
            ch[(0) + ido * ((k) + l1 * (0))] = t0[0] + t1[0];
            ch[(0) + ido * ((k) + l1 * (0)) + 1] = t0[1] + t1[1];

            //PARTSTEP3a(1,2,tw1r,tw1i)
            //PARTSTEP3a(u1,u2,twr,twi)

            let mut ca: [f64; 2] = [0.0; 2]; //cmplx { r: 0.0, i: 0.0 };
            let mut cb: [f64; 2] = [0.0; 2]; //cmplx { r: 0.0, i: 0.0 };
            ca[0] = t0[0] + tw1r * t1[0];
            ca[1] = t0[1] + tw1r * t1[1];
            cb[1] = tw1i * t2[0];
            cb[0] = -1.0 * (tw1i * t2[1]);

            ch[(0) + ido * ((k) + l1 * (1))] = ca[0] + cb[0];
            ch[(0) + ido * ((k) + l1 * (1)) + 1] = ca[1] + cb[1];
            ch[(0) + ido * ((k) + l1 * (2))] = ca[0] - cb[0];
            ch[(0) + ido * ((k) + l1 * (2)) + 1] = ca[1] - cb[1];
        }
    } else {
        for k in 0..l1 {
            {
                //PREP3(0)
                let mut t0: [f64; 2] = [0.0; 2];
                t0[0] = cc[(0) + ido * ((0) + cdim * (k))];
                t0[1] = cc[(0) + ido * ((0) + cdim * (k)) + 1];
                let mut t1: [f64; 2] = [0.0; 2]; //cmplx { r: 0.0, i: 0.0 };
                let mut t2: [f64; 2] = [0.0; 2]; //cmplx { r: 0.0, i: 0.0 };

                //PMC (t1,t2,CC(idx,1,k),CC(idx,2,k))
                t1[0] = cc[(0) + ido * ((1) + cdim * (k))] + cc[(0) + ido * ((2) + cdim * (k))];
                t1[1] =
                    cc[(0) + ido * ((1) + cdim * (k)) + 1] + cc[(0) + ido * ((2) + cdim * (k)) + 1];
                t2[0] = cc[(0) + ido * ((1) + cdim * (k))] - cc[(0) + ido * ((2) + cdim * (k))];
                t2[1] =
                    cc[(0) + ido * ((1) + cdim * (k)) + 1] - cc[(0) + ido * ((2) + cdim * (k)) + 1];
                //CH(idx,k,0).r=t0.r+t1.r;
                ch[(0) + ido * ((k) + l1 * (0))] = t0[0] + t1[0];

                //CH(idx,k,0).i=t0.i+t1.i;
                ch[(0) + ido * ((k) + l1 * (0)) + 1] = t0[1] + t1[1];

                //PARTSTEP3a(1,2,tw1r,tw1i)
                let mut ca: [f64; 2] = [0.0; 2]; //cmplx { r: 0.0, i: 0.0 };
                let mut cb: [f64; 2] = [0.0; 2]; //cmplx { r: 0.0, i: 0.0 };
                ca[0] = t0[0] + tw1r * t1[0];
                ca[1] = t0[1] + tw1r * t1[1];
                cb[1] = tw1i * t2[0];
                cb[0] = -1.0 * (tw1i * t2[1]);
                ch[(0) + ido * ((k) + l1 * (1))] = ca[0] + cb[0];
                ch[(0) + ido * ((k) + l1 * (1)) + 1] = ca[1] + cb[1];
                ch[(0) + ido * ((k) + l1 * (2))] = ca[0] - cb[0];
                ch[(0) + ido * ((k) + l1 * (2)) + 1] = ca[1] - cb[1];
            }
            for i in 1..ido {
                //PREP3(i)
                let mut t0: [f64; 2] = [0.0; 2];
                t0[0] = cc[(i) + ido * ((0) + cdim * (k))];
                t0[1] = cc[(i) + ido * ((0) + cdim * (k)) + 1];
                let mut t1: [f64; 2] = [0.0; 2]; //cmplx { r: 0.0, i: 0.0 };
                let mut t2: [f64; 2] = [0.0; 2]; //cmplx { r: 0.0, i: 0.0 };
                t1[0] = cc[(i) + ido * ((1) + cdim * (k))] + cc[(i) + ido * ((2) + cdim * (k))];
                t1[1] =
                    cc[(i) + ido * ((1) + cdim * (k)) + 1] + cc[(i) + ido * ((2) + cdim * (k)) + 1];
                t2[0] = cc[(i) + ido * ((1) + cdim * (k))] - cc[(i) + ido * ((2) + cdim * (k))];
                t2[1] =
                    cc[(i) + ido * ((1) + cdim * (k)) + 1] - cc[(i) + ido * ((2) + cdim * (k)) + 1];
                ch[(i) + ido * ((k) + l1 * (0))] = t0[0] + t1[0];
                ch[(i) + ido * ((k) + l1 * (0)) + 1] = t0[1] + t1[1];

                //PARTSTEP3b(1,2,tw1r,tw1i)
                let mut ca: [f64; 2] = [0.0; 2];
                let mut cb: [f64; 2] = [0.0; 2];
                let mut da: [f64; 2] = [0.0; 2];
                let mut db: [f64; 2] = [0.0; 2];
                ca[0] = t0[0] + tw1r * t1[0];
                ca[1] = t0[1] + tw1r * t1[1];
                cb[1] = tw1i * t2[0];
                cb[0] = -(tw1i * t2[1]);
                da[0] = ca[0] + cb[0];
                da[1] = ca[1] + cb[1];
                db[0] = ca[0] - cb[0];
                db[1] = ca[1] - cb[1];

                ch[(i) + ido * ((k) + l1 * (1))] = wa[(i) - 1 + (1 - 1) * (ido - 1)] * da[0]
                    - wa[(i) - 1 + (1 - 1) * (ido - 1) + 1] * da[1];
                ch[(i) + ido * ((k) + l1 * (1)) + 1] = wa[(i) - 1 + (1 - 1) * (ido - 1)] * da[1]
                    + wa[(i) - 1 + (1 - 1) * (ido - 1) + 1] * da[0];

                ch[(i) + ido * ((k) + l1 * (2))] = wa[(i) - 1 + (2 - 1) * (ido - 1)] * db[0]
                    - wa[(i) - 1 + (2 - 1) * (ido - 1) + 1] * db[1];
                ch[(i) + ido * ((k) + l1 * (2)) + 1] = wa[(i) - 1 + (2 - 1) * (ido - 1)] * db[1]
                    + wa[(i) - 1 + (2 - 1) * (ido - 1) + 1] * db[0];
            }
        }
    }
}

fn pass3f(ido: usize, l1: usize, cc: &[f64], ch: &mut [f64], wa: &[f64]) {
    let cdim: usize = 3;
    let tw1r: f64 = -0.5;
    let tw1i: f64 = -0.86602540378443864676;

    if ido == 1 {
        for k in 0..l1 {
            let mut t0: [f64; 2] = [0.0; 2];
            t0[0] = cc[(0) + ido * ((0) + cdim * (k))];
            t0[1] = cc[(0) + ido * ((0) + cdim * (k)) + 1];
            let mut t1: [f64; 2] = [0.0; 2];
            let mut t2: [f64; 2] = [0.0; 2];
            {
                t1[0] = cc[(0) + ido * ((1) + cdim * (k))] + cc[(0) + ido * ((2) + cdim * (k))];
                t1[1] =
                    cc[(0) + ido * ((1) + cdim * (k)) + 1] + cc[(0) + ido * ((2) + cdim * (k)) + 1];
                t2[0] = cc[(0) + ido * ((1) + cdim * (k))] - cc[(0) + ido * ((2) + cdim * (k))];
                t2[1] =
                    cc[(0) + ido * ((1) + cdim * (k)) + 1] - cc[(0) + ido * ((2) + cdim * (k)) + 1];
            }
            ch[(0) + ido * ((k) + l1 * (0))] = t0[0] + t1[0];
            ch[(0) + ido * ((k) + l1 * (0)) + 1] = t0[1] + t1[1];
            {
                let mut ca: [f64; 2] = [0.0; 2];
                let mut cb: [f64; 2] = [0.0; 2];
                ca[0] = t0[0] + tw1r * t1[0];
                ca[1] = t0[1] + tw1r * t1[1];
                cb[1] = tw1i * t2[0];
                cb[0] = -(tw1i * t2[1]);
                {
                    ch[(0) + ido * ((k) + l1 * (1))] = ca[0] + cb[0];
                    ch[(0) + ido * ((k) + l1 * (1)) + 1] = ca[1] + cb[1];
                    ch[(0) + ido * ((k) + l1 * (2))] = ca[0] - cb[0];
                    ch[(0) + ido * ((k) + l1 * (2)) + 1] = ca[1] - cb[1];
                }
            }
        }
    } else {
        for k in 0..l1 {
            {
                let mut t0: [f64; 2] = [0.0; 2];
                t0[0] = cc[(0) + ido * ((0) + cdim * (k))];
                t0[1] = cc[(0) + ido * ((0) + cdim * (k)) + 1];
                let mut t1: [f64; 2] = [0.0; 2];
                let mut t2: [f64; 2] = [0.0; 2];
                {
                    t1[0] = cc[(0) + ido * ((1) + cdim * (k))] + cc[(0) + ido * ((2) + cdim * (k))];
                    t1[1] = cc[(0) + ido * ((1) + cdim * (k)) + 1]
                        + cc[(0) + ido * ((2) + cdim * (k)) + 1];
                    t2[0] = cc[(0) + ido * ((1) + cdim * (k))] - cc[(0) + ido * ((2) + cdim * (k))];
                    t2[1] = cc[(0) + ido * ((1) + cdim * (k)) + 1]
                        - cc[(0) + ido * ((2) + cdim * (k)) + 1];
                }
                ch[(0) + ido * ((k) + l1 * (0))] = t0[0] + t1[0];
                ch[(0) + ido * ((k) + l1 * (0)) + 1] = t0[1] + t1[1];
                {
                    let mut ca: [f64; 2] = [0.0; 2];
                    let mut cb: [f64; 2] = [0.0; 2];
                    ca[0] = t0[0] + tw1r * t1[0];
                    ca[1] = t0[1] + tw1r * t1[1];
                    cb[1] = tw1i * t2[0];
                    cb[0] = -(tw1i * t2[1]);
                    {
                        ch[(0) + ido * ((k) + l1 * (1))] = ca[0] + cb[0];
                        ch[(0) + ido * ((k) + l1 * (1)) + 1] = ca[1] + cb[1];
                        ch[(0) + ido * ((k) + l1 * (2))] = ca[0] - cb[0];
                        ch[(0) + ido * ((k) + l1 * (2)) + 1] = ca[1] - cb[1];
                    }
                }
            }
            for i in 1..ido {
                let mut t0: [f64; 2] = [0.0; 2];
                t0[0] = cc[(i) + ido * ((0) + cdim * (k))];
                t0[1] = cc[(i) + ido * ((0) + cdim * (k)) + 1];
                let mut t1: [f64; 2] = [0.0; 2];
                let mut t2: [f64; 2] = [0.0; 2];
                {
                    t1[0] = cc[(i) + ido * ((1) + cdim * (k))] + cc[(i) + ido * ((2) + cdim * (k))];
                    t1[1] = cc[(i) + ido * ((1) + cdim * (k)) + 1]
                        + cc[(i) + ido * ((2) + cdim * (k)) + 1];
                    t2[0] =
                        cc[(i) + ido * ((1) + cdim * (k))] - cc[(i) + ido * ((2) + cdim * (k)) + 1];
                    t2[1] = cc[(i) + ido * ((1) + cdim * (k)) + 1]
                        - cc[(i) + ido * ((2) + cdim * (k)) + 1];
                }
                ch[(i) + ido * ((k) + l1 * (0))] = t0[0] + t1[0];
                ch[(i) + ido * ((k) + l1 * (0)) + 1] = t0[1] + t1[1];
                {
                    let mut ca: [f64; 2] = [0.0; 2];
                    let mut cb: [f64; 2] = [0.0; 2];
                    let mut da: [f64; 2] = [0.0; 2];
                    let mut db: [f64; 2] = [0.0; 2];
                    ca[0] = t0[0] + tw1r * t1[0];
                    ca[1] = t0[1] + tw1r * t1[1];
                    cb[1] = tw1i * t2[0];
                    cb[0] = -(tw1i * t2[1]);
                    {
                        da[0] = ca[0] + cb[0];
                        da[1] = ca[1] + cb[1];
                        db[0] = ca[0] - cb[0];
                        db[1] = ca[1] - cb[1];
                    }
                    {
                        ch[(i) + ido * ((k) + l1 * (1))] = wa[(i) - 1 + (1 - 1) * (ido - 1)]
                            * da[0]
                            + wa[(i) - 1 + (1 - 1) * (ido - 1) + 1] * da[1];
                        ch[(i) + ido * ((k) + l1 * (1)) + 1] = wa[(i) - 1 + (1 - 1) * (ido - 1)]
                            * da[1]
                            - wa[(i) - 1 + (1 - 1) * (ido - 1) + 1] * da[0];
                    }
                    {
                        ch[(i) + ido * ((k) + l1 * (2))] = wa[(i) - 1 + (2 - 1) * (ido - 1)]
                            * db[0]
                            + wa[(i) - 1 + (2 - 1) * (ido - 1) + 1] * db[1];
                        ch[(i) + ido * ((k) + l1 * (2)) + 1] = wa[(i) - 1 + (2 - 1) * (ido - 1)]
                            * db[1]
                            - wa[(i) - 1 + (2 - 1) * (ido - 1) + 1] * db[0];
                    }
                }
            }
        }
    }
}
fn pass4b(ido: usize, l1: usize, cc: &mut [f64], ch: &mut [f64], wa: &[f64]) {
    let cdim: usize = 4;

    if ido == 1 {
        for k in 0..l1 {
            let mut t1: [f64; 2] = [0.0; 2];
            let mut t2: [f64; 2] = [0.0; 2];
            let mut t3: [f64; 2] = [0.0; 2];
            let mut t4: [f64; 2] = [0.0; 2];
            {
                t2[0] = cc[(0) + ido * ((0) + cdim * (k))] + cc[(0) + ido * ((2) + cdim * (k))];
                t2[1] =
                    cc[(0) + ido * ((0) + cdim * (k)) + 1] + cc[(0) + ido * ((2) + cdim * (k)) + 1];
                t1[0] = cc[(0) + ido * ((0) + cdim * (k))] - cc[(0) + ido * ((2) + cdim * (k))];
                t1[1] =
                    cc[(0) + ido * ((0) + cdim * (k)) + 1] - cc[(0) + ido * ((2) + cdim * (k)) + 1];
            }
            {
                t3[0] = cc[(0) + ido * ((1) + cdim * (k))] + cc[(0) + ido * ((3) + cdim * (k))];
                t3[1] =
                    cc[(0) + ido * ((1) + cdim * (k)) + 1] + cc[(0) + ido * ((3) + cdim * (k)) + 1];
                t4[0] = cc[(0) + ido * ((1) + cdim * (k))] - cc[(0) + ido * ((3) + cdim * (k))];
                t4[1] =
                    cc[(0) + ido * ((1) + cdim * (k)) + 1] - cc[(0) + ido * ((3) + cdim * (k)) + 1];
            }
            {
                let tmp_ = t4[0];
                t4[0] = -t4[1];
                t4[1] = tmp_;
            }
            {
                ch[(0) + ido * ((k) + l1 * (0))] = t2[0] + t3[0];
                ch[(0) + ido * ((k) + l1 * (0)) + 1] = t2[1] + t3[1];
                ch[(0) + ido * ((k) + l1 * (2))] = t2[0] - t3[0];
                ch[(0) + ido * ((k) + l1 * (2)) + 1] = t2[1] - t3[1];
            }
            {
                ch[(0) + ido * ((k) + l1 * (1))] = t1[0] + t4[0];
                ch[(0) + ido * ((k) + l1 * (1)) + 1] = t1[1] + t4[1];
                ch[(0) + ido * ((k) + l1 * (3))] = t1[0] - t4[0];
                ch[(0) + ido * ((k) + l1 * (3)) + 1] = t1[1] - t4[1];
            }
        }
    } else {
        for k in 0..l1 {
            {
                let mut t1: [f64; 2] = [0.0; 2];
                let mut t2: [f64; 2] = [0.0; 2];
                let mut t3: [f64; 2] = [0.0; 2];
                let mut t4: [f64; 2] = [0.0; 2];
                {
                    t2[0] = cc[(0) + ido * ((0) + cdim * (k))] + cc[(0) + ido * ((2) + cdim * (k))];
                    t2[1] = cc[(0) + ido * ((0) + cdim * (k)) + 1]
                        + cc[(0) + ido * ((2) + cdim * (k)) + 1];
                    t1[0] = cc[(0) + ido * ((0) + cdim * (k))] - cc[(0) + ido * ((2) + cdim * (k))];
                    t1[1] = cc[(0) + ido * ((0) + cdim * (k)) + 1]
                        - cc[(0) + ido * ((2) + cdim * (k)) + 1];
                }
                {
                    t3[0] = cc[(0) + ido * ((1) + cdim * (k))] + cc[(0) + ido * ((3) + cdim * (k))];
                    t3[1] = cc[(0) + ido * ((1) + cdim * (k)) + 1]
                        + cc[(0) + ido * ((3) + cdim * (k)) + 1];
                    t4[0] = cc[(0) + ido * ((1) + cdim * (k))] - cc[(0) + ido * ((3) + cdim * (k))];
                    t4[1] = cc[(0) + ido * ((1) + cdim * (k)) + 1]
                        - cc[(0) + ido * ((3) + cdim * (k)) + 1];
                }
                {
                    let tmp_ = t4[0];
                    t4[0] = -t4[1];
                    t4[1] = tmp_;
                }
                {
                    ch[(0) + ido * ((k) + l1 * (0))] = t2[0] + t3[0];
                    ch[(0) + ido * ((k) + l1 * (0)) + 1] = t2[1] + t3[1];
                    ch[(0) + ido * ((k) + l1 * (2))] = t2[0] - t3[0];
                    ch[(0) + ido * ((k) + l1 * (2)) + 1] = t2[1] - t3[1];
                }
                {
                    ch[(0) + ido * ((k) + l1 * (1))] = t1[0] + t4[0];
                    ch[(0) + ido * ((k) + l1 * (1)) + 1] = t1[1] + t4[1];
                    ch[(0) + ido * ((k) + l1 * (3))] = t1[0] - t4[0];
                    ch[(0) + ido * ((k) + l1 * (3)) + 1] = t1[1] - t4[1];
                }
            }
            for i in 1..ido {
                let mut c2: [f64; 2] = [0.0; 2];
                let mut c3: [f64; 2] = [0.0; 2];
                let mut c4: [f64; 2] = [0.0; 2];
                let mut t1: [f64; 2] = [0.0; 2];
                let mut t2: [f64; 2] = [0.0; 2];
                let mut t3: [f64; 2] = [0.0; 2];
                let mut t4: [f64; 2] = [0.0; 2];
                let mut cc0: [f64; 2] = [0.0; 2];
                cc0[0] = cc[(i) + ido * ((0) + cdim * (k))];
                cc0[1] = cc[(i) + ido * ((0) + cdim * (k)) + 1];
                let mut cc1: [f64; 2] = [0.0; 2];
                cc1[0] = cc[(i) + ido * ((1) + cdim * (k))];
                cc1[1] = cc[(i) + ido * ((1) + cdim * (k)) + 1];
                let mut cc2: [f64; 2] = [0.0; 2];
                cc2[0] = cc[(i) + ido * ((2) + cdim * (k))];
                cc2[1] = cc[(i) + ido * ((2) + cdim * (k)) + 1];
                let mut cc3: [f64; 2] = [0.0; 2];
                cc3[0] = cc[(i) + ido * ((3) + cdim * (k))];
                cc3[1] = cc[(i) + ido * ((3) + cdim * (k)) + 1];
                {
                    t2[0] = cc0[0] + cc2[0];
                    t2[1] = cc0[1] + cc2[1];
                    t1[0] = cc0[0] - cc2[0];
                    t1[1] = cc0[1] - cc2[1];
                }
                {
                    t3[0] = cc1[0] + cc3[0];
                    t3[1] = cc1[1] + cc3[1];
                    t4[0] = cc1[0] - cc3[0];
                    t4[1] = cc1[1] - cc3[1];
                }
                {
                    let tmp_ = t4[0];
                    t4[0] = -t4[1];
                    t4[1] = tmp_;
                }
                let mut wa0: [f64; 2] = [0.0; 2];
                wa0[0] = wa[(i) - 1 + (0) * (ido - 1)];
                wa0[1] = wa[(i) - 1 + (0) * (ido - 1) + 1];
                let mut wa1: [f64; 2] = [0.0; 2];
                wa1[0] = wa[(i) - 1 + (1) * (ido - 1)];
                wa1[1] = wa[(i) - 1 + (1) * (ido - 1) + 1];
                let mut wa2: [f64; 2] = [0.0; 2];
                wa2[0] = wa[(i) - 1 + (2) * (ido - 1)];
                wa2[1] = wa[(i) - 1 + (2) * (ido - 1) + 1];
                {
                    ch[(i) + ido * ((k) + l1 * (0))] = t2[0] + t3[0];
                    ch[(i) + ido * ((k) + l1 * (0)) + 1] = t2[1] + t3[1];
                    c3[0] = t2[0] - t3[0];
                    c3[1] = t2[1] - t3[1];
                }
                {
                    c2[0] = t1[0] + t4[0];
                    c2[1] = t1[1] + t4[1];
                    c4[0] = t1[0] - t4[0];
                    c4[1] = t1[1] - t4[1];
                }
                {
                    ch[(i) + ido * ((k) + l1 * (1))] = wa0[0] * c2[0] - wa0[1] * c2[1];
                    ch[(i) + ido * ((k) + l1 * (1)) + 1] = wa0[0] * c2[1] + wa0[1] * c2[0];
                }
                {
                    ch[(i) + ido * ((k) + l1 * (2))] = wa1[0] * c3[0] - wa1[1] * c3[1];
                    ch[(i) + ido * ((k) + l1 * (2)) + 1] = wa1[0] * c3[1] + wa1[1] * c3[0];
                }
                {
                    ch[(i) + ido * ((k) + l1 * (3))] = wa2[0] * c4[0] - wa2[1] * c4[1];
                    ch[(i) + ido * ((k) + l1 * (3)) + 1] = wa2[0] * c4[1] + wa2[1] * c4[0];
                }
            }
        }
    }
}

fn pass4f(ido: usize, l1: usize, cc: &[f64], ch: &mut [f64], wa: &[f64]) {
    let cdim: usize = 4;

    if ido == 1 {
        for k in 0..l1 {
            let mut t1: [f64; 2] = [0.0; 2];
            let mut t2: [f64; 2] = [0.0; 2];
            let mut t3: [f64; 2] = [0.0; 2];
            let mut t4: [f64; 2] = [0.0; 2];
            {
                t2[0] = cc[(0) + ido * ((0) + cdim * (k))] + cc[(0) + ido * ((2) + cdim * (k))];
                t2[1] =
                    cc[(0) + ido * ((0) + cdim * (k)) + 1] + cc[(0) + ido * ((2) + cdim * (k)) + 1];
                t1[0] = cc[(0) + ido * ((0) + cdim * (k))] - cc[(0) + ido * ((2) + cdim * (k))];
                t1[1] =
                    cc[(0) + ido * ((0) + cdim * (k)) + 1] - cc[(0) + ido * ((2) + cdim * (k)) + 1];
            }
            {
                t3[0] = cc[(0) + ido * ((1) + cdim * (k))] + cc[(0) + ido * ((3) + cdim * (k))];
                t3[1] =
                    cc[(0) + ido * ((1) + cdim * (k)) + 1] + cc[(0) + ido * ((3) + cdim * (k)) + 1];
                t4[0] = cc[(0) + ido * ((1) + cdim * (k))] - cc[(0) + ido * ((3) + cdim * (k))];
                t4[1] =
                    cc[(0) + ido * ((1) + cdim * (k)) + 1] - cc[(0) + ido * ((3) + cdim * (k)) + 1];
            }
            {
                let tmp_ = -t4[0];
                t4[0] = t4[1];
                t4[1] = tmp_;
            }
            {
                ch[(0) + ido * ((k) + l1 * (0))] = t2[0] + t3[0];
                ch[(0) + ido * ((k) + l1 * (0)) + 1] = t2[1] + t3[1];
                ch[(0) + ido * ((k) + l1 * (2))] = t2[0] - t3[0];
                ch[(0) + ido * ((k) + l1 * (2)) + 1] = t2[1] - t3[1];
            }
            {
                ch[(0) + ido * ((k) + l1 * (1))] = t1[0] + t4[0];
                ch[(0) + ido * ((k) + l1 * (1)) + 1] = t1[1] + t4[1];
                ch[(0) + ido * ((k) + l1 * (3))] = t1[0] - t4[0];
                ch[(0) + ido * ((k) + l1 * (3)) + 1] = t1[1] - t4[1];
            }
        }
    } else {
        for k in 0..l1 {
            {
                let mut t1: [f64; 2] = [0.0; 2];
                let mut t2: [f64; 2] = [0.0; 2];
                let mut t3: [f64; 2] = [0.0; 2];
                let mut t4: [f64; 2] = [0.0; 2];
                {
                    t2[0] = cc[(0) + ido * ((0) + cdim * (k))] + cc[(0) + ido * ((2) + cdim * (k))];
                    t2[1] = cc[(0) + ido * ((0) + cdim * (k)) + 1]
                        + cc[(0) + ido * ((2) + cdim * (k)) + 1];
                    t1[0] = cc[(0) + ido * ((0) + cdim * (k))] - cc[(0) + ido * ((2) + cdim * (k))];
                    t1[1] = cc[(0) + ido * ((0) + cdim * (k)) + 1]
                        - cc[(0) + ido * ((2) + cdim * (k)) + 1];
                }
                {
                    t3[0] = cc[(0) + ido * ((1) + cdim * (k))] + cc[(0) + ido * ((3) + cdim * (k))];
                    t3[1] = cc[(0) + ido * ((1) + cdim * (k)) + 1]
                        + cc[(0) + ido * ((3) + cdim * (k)) + 1];
                    t4[0] = cc[(0) + ido * ((1) + cdim * (k))] - cc[(0) + ido * ((3) + cdim * (k))];
                    t4[1] = cc[(0) + ido * ((1) + cdim * (k)) + 1]
                        - cc[(0) + ido * ((3) + cdim * (k)) + 1];
                }
                {
                    let tmp_ = -t4[0];
                    t4[0] = t4[1];
                    t4[1] = tmp_;
                }
                {
                    ch[(0) + ido * ((k) + l1 * (0))] = t2[0] + t3[0];
                    ch[(0) + ido * ((k) + l1 * (0)) + 1] = t2[1] + t3[1];
                    ch[(0) + ido * ((k) + l1 * (2))] = t2[0] - t3[0];
                    ch[(0) + ido * ((k) + l1 * (2)) + 1] = t2[1] - t3[1];
                }
                {
                    ch[(0) + ido * ((k) + l1 * (1))] = t1[0] + t4[0];
                    ch[(0) + ido * ((k) + l1 * (1)) + 1] = t1[1] + t4[1];
                    ch[(0) + ido * ((k) + l1 * (3))] = t1[0] - t4[0];
                    ch[(0) + ido * ((k) + l1 * (3)) + 1] = t1[1] - t4[1];
                }
            }
            for i in 1..ido {
                let mut c2: [f64; 2] = [0.0; 2];
                let mut c3: [f64; 2] = [0.0; 2];
                let mut c4: [f64; 2] = [0.0; 2];
                let mut t1: [f64; 2] = [0.0; 2];
                let mut t2: [f64; 2] = [0.0; 2];
                let mut t3: [f64; 2] = [0.0; 2];
                let mut t4: [f64; 2] = [0.0; 2];
                let mut cc0: [f64; 2] = [0.0; 2];
                cc0[0] = cc[(i) + ido * ((0) + cdim * (k))];
                cc0[1] = cc[(i) + ido * ((0) + cdim * (k)) + 1];
                let mut cc1: [f64; 2] = [0.0; 2];
                cc1[0] = cc[(i) + ido * ((1) + cdim * (k))];
                cc1[1] = cc[(i) + ido * ((1) + cdim * (k)) + 1];
                let mut cc2: [f64; 2] = [0.0; 2];
                cc2[0] = cc[(i) + ido * ((2) + cdim * (k))];
                cc2[1] = cc[(i) + ido * ((2) + cdim * (k)) + 1];
                let mut cc3: [f64; 2] = [0.0; 2];
                cc3[0] = cc[(i) + ido * ((3) + cdim * (k))];
                cc3[1] = cc[(i) + ido * ((3) + cdim * (k)) + 1];
                {
                    t2[0] = cc0[0] + cc2[0];
                    t2[1] = cc0[1] + cc2[1];
                    t1[0] = cc0[0] - cc2[0];
                    t1[1] = cc0[1] - cc2[1];
                }
                {
                    t3[0] = cc1[0] + cc3[0];
                    t3[1] = cc1[1] + cc3[1];
                    t4[0] = cc1[0] - cc3[0];
                    t4[1] = cc1[1] - cc3[1];
                }
                {
                    let tmp_ = -t4[0];
                    t4[0] = t4[1];
                    t4[1] = tmp_;
                }
                let mut wa0: [f64; 2] = [0.0; 2];
                wa0[0] = wa[(i) - 1 + (0) * (ido - 1)];
                wa0[1] = wa[(i) - 1 + (0) * (ido - 1) + 1];
                let mut wa1: [f64; 2] = [0.0; 2];
                wa1[0] = wa[(i) - 1 + (1) * (ido - 1)];
                wa1[1] = wa[(i) - 1 + (1) * (ido - 1) + 1];
                let mut wa2: [f64; 2] = [0.0; 2];
                wa2[0] = wa[(i) - 1 + (2) * (ido - 1)];
                wa2[1] = wa[(i) - 1 + (2) * (ido - 1) + 1];
                {
                    ch[(i) + ido * ((k) + l1 * (0))] = t2[0] + t3[0];
                    ch[(i) + ido * ((k) + l1 * (0)) + 1] = t2[1] + t3[1];
                    c3[0] = t2[0] - t3[0];
                    c3[1] = t2[1] - t3[1];
                }
                {
                    c2[0] = t1[0] + t4[0];
                    c2[1] = t1[1] + t4[1];
                    c4[0] = t1[0] - t4[0];
                    c4[1] = t1[1] - t4[1];
                }
                {
                    ch[(i) + ido * ((k) + l1 * (1))] = wa0[0] * c2[0] + wa0[1] * c2[1];
                    ch[(i) + ido * ((k) + l1 * (1)) + 1] = wa0[0] * c2[1] - wa0[1] * c2[0];
                }
                {
                    ch[(i) + ido * ((k) + l1 * (2))] = wa1[0] * c3[0] + wa1[1] * c3[1];
                    ch[(i) + ido * ((k) + l1 * (2)) + 1] = wa1[0] * c3[1] - wa1[1] * c3[0];
                }
                {
                    ch[(i) + ido * ((k) + l1 * (3))] = wa2[0] * c4[0] + wa2[1] * c4[1];
                    ch[(i) + ido * ((k) + l1 * (3)) + 1] = wa2[0] * c4[1] - wa2[1] * c4[0];
                }
            }
        }
    }
}

fn pass5b(ido: usize, l1: usize, cc: &[f64], ch: &mut [f64], wa: &[f64]) {
    let cdim: usize = 5;
    let tw1r: f64 = 0.3090169943749474241;
    let tw1i: f64 = 0.95105651629515357212;
    let tw2r: f64 = -0.8090169943749474241;
    let tw2i: f64 = 0.58778525229247312917;

    if ido == 1 {
        for k in 0..l1 {
            let mut t0: [f64; 2] = [0.0; 2];
            t0[0] = cc[(0) + ido * ((0) + cdim * (k))];
            t0[1] = cc[(0) + ido * ((0) + cdim * (k)) + 1];
            let mut t1: [f64; 2] = [0.0; 2];
            let mut t2: [f64; 2] = [0.0; 2];
            let mut t3: [f64; 2] = [0.0; 2];
            let mut t4: [f64; 2] = [0.0; 2];
            {
                t1[0] = cc[(0) + ido * ((1) + cdim * (k))] + cc[(0) + ido * ((4) + cdim * (k))];
                t1[1] =
                    cc[(0) + ido * ((1) + cdim * (k)) + 1] + cc[(0) + ido * ((4) + cdim * (k)) + 1];
                t4[0] = cc[(0) + ido * ((1) + cdim * (k))] - cc[(0) + ido * ((4) + cdim * (k))];
                t4[1] =
                    cc[(0) + ido * ((1) + cdim * (k)) + 1] - cc[(0) + ido * ((4) + cdim * (k)) + 1];
            }
            {
                t2[0] = cc[(0) + ido * ((2) + cdim * (k))] + cc[(0) + ido * ((3) + cdim * (k))];
                t2[1] =
                    cc[(0) + ido * ((2) + cdim * (k)) + 1] + cc[(0) + ido * ((3) + cdim * (k)) + 1];
                t3[0] = cc[(0) + ido * ((2) + cdim * (k))] - cc[(0) + ido * ((3) + cdim * (k))];
                t3[1] =
                    cc[(0) + ido * ((2) + cdim * (k)) + 1] - cc[(0) + ido * ((3) + cdim * (k)) + 1];
            }
            ch[(0) + ido * ((k) + l1 * (0))] = t0[0] + t1[0] + t2[0];
            ch[(0) + ido * ((k) + l1 * (0)) + 1] = t0[1] + t1[1] + t2[1];
            {
                let mut ca: [f64; 2] = [0.0; 2];
                let mut cb: [f64; 2] = [0.0; 2];
                ca[0] = t0[0] + tw1r * t1[0] + tw2r * t2[0];
                ca[1] = t0[1] + tw1r * t1[1] + tw2r * t2[1];
                cb[1] = tw1i * t4[0] + tw2i * t3[0];
                cb[0] = -(tw1i * t4[1] + tw2i * t3[1]);
                {
                    ch[(0) + ido * ((k) + l1 * (1))] = ca[0] + cb[0];
                    ch[(0) + ido * ((k) + l1 * (1)) + 1] = ca[1] + cb[1];
                    ch[(0) + ido * ((k) + l1 * (4))] = ca[0] - cb[0];
                    ch[(0) + ido * ((k) + l1 * (4)) + 1] = ca[1] - cb[1];
                }
            }
            {
                let mut ca: [f64; 2] = [0.0; 2];
                let mut cb: [f64; 2] = [0.0; 2];
                ca[0] = t0[0] + tw2r * t1[0] + tw1r * t2[0];
                ca[1] = t0[1] + tw2r * t1[1] + tw1r * t2[1];
                cb[1] = tw2i * t4[0] - tw1i * t3[0];
                cb[0] = -(tw2i * t4[1] - tw1i * t3[1]);
                {
                    ch[(0) + ido * ((k) + l1 * (2))] = ca[0] + cb[0];
                    ch[(0) + ido * ((k) + l1 * (2)) + 1] = ca[1] + cb[1];
                    ch[(0) + ido * ((k) + l1 * (3))] = ca[0] - cb[0];
                    ch[(0) + ido * ((k) + l1 * (3)) + 1] = ca[1] - cb[1];
                }
            }
        }
    } else {
        for k in 0..l1 {
            {
                let mut t0: [f64; 2] = [0.0; 2];
                t0[0] = cc[(0) + ido * ((0) + cdim * (k))];
                t0[1] = cc[(0) + ido * ((0) + cdim * (k)) + 1];
                let mut t1: [f64; 2] = [0.0; 2];
                let mut t2: [f64; 2] = [0.0; 2];
                let mut t3: [f64; 2] = [0.0; 2];
                let mut t4: [f64; 2] = [0.0; 2];
                {
                    t1[0] = cc[(0) + ido * ((1) + cdim * (k))] + cc[(0) + ido * ((4) + cdim * (k))];
                    t1[1] = cc[(0) + ido * ((1) + cdim * (k)) + 1]
                        + cc[(0) + ido * ((4) + cdim * (k)) + 1];
                    t4[0] = cc[(0) + ido * ((1) + cdim * (k))] - cc[(0) + ido * ((4) + cdim * (k))];
                    t4[1] = cc[(0) + ido * ((1) + cdim * (k)) + 1]
                        - cc[(0) + ido * ((4) + cdim * (k)) + 1];
                }
                {
                    t2[0] = cc[(0) + ido * ((2) + cdim * (k))] + cc[(0) + ido * ((3) + cdim * (k))];
                    t2[1] = cc[(0) + ido * ((2) + cdim * (k)) + 1]
                        + cc[(0) + ido * ((3) + cdim * (k)) + 1];
                    t3[0] = cc[(0) + ido * ((2) + cdim * (k))] - cc[(0) + ido * ((3) + cdim * (k))];
                    t3[1] = cc[(0) + ido * ((2) + cdim * (k)) + 1]
                        - cc[(0) + ido * ((3) + cdim * (k)) + 1];
                }
                ch[(0) + ido * ((k) + l1 * (0))] = t0[0] + t1[0] + t2[0];
                ch[(0) + ido * ((k) + l1 * (0)) + 1] = t0[1] + t1[1] + t2[1];
                {
                    let mut ca: [f64; 2] = [0.0; 2];
                    let mut cb: [f64; 2] = [0.0; 2];
                    ca[0] = t0[0] + tw1r * t1[0] + tw2r * t2[0];
                    ca[1] = t0[1] + tw1r * t1[1] + tw2r * t2[1];
                    cb[1] = tw1i * t4[0] + tw2i * t3[0];
                    cb[0] = -(tw1i * t4[1] + tw2i * t3[1]);
                    {
                        ch[(0) + ido * ((k) + l1 * (1))] = ca[0] + cb[0];
                        ch[(0) + ido * ((k) + l1 * (1)) + 1] = ca[1] + cb[1];
                        ch[(0) + ido * ((k) + l1 * (4))] = ca[0] - cb[0];
                        ch[(0) + ido * ((k) + l1 * (4)) + 1] = ca[1] - cb[1];
                    }
                }
                {
                    let mut ca: [f64; 2] = [0.0; 2];
                    let mut cb: [f64; 2] = [0.0; 2];
                    ca[0] = t0[0] + tw2r * t1[0] + tw1r * t2[0];
                    ca[1] = t0[1] + tw2r * t1[1] + tw1r * t2[1];
                    cb[1] = tw2i * t4[0] - tw1i * t3[0];
                    cb[0] = -(tw2i * t4[1] - tw1i * t3[1]);
                    {
                        ch[(0) + ido * ((k) + l1 * (2))] = ca[0] + cb[0];
                        ch[(0) + ido * ((k) + l1 * (2)) + 1] = ca[1] + cb[1];
                        ch[(0) + ido * ((k) + l1 * (3))] = ca[0] - cb[0];
                        ch[(0) + ido * ((k) + l1 * (3)) + 1] = ca[1] - cb[1];
                    }
                }
            }
            for i in 1..ido {
                let mut t0: [f64; 2] = [0.0; 2];
                t0[0] = cc[(i) + ido * ((0) + cdim * (k))];
                t0[1] = cc[(i) + ido * ((0) + cdim * (k)) + 1];
                let mut t1: [f64; 2] = [0.0; 2];
                let mut t2: [f64; 2] = [0.0; 2];
                let mut t3: [f64; 2] = [0.0; 2];
                let mut t4: [f64; 2] = [0.0; 2];
                {
                    t1[0] = cc[(i) + ido * ((1) + cdim * (k))] + cc[(i) + ido * ((4) + cdim * (k))];
                    t1[1] = cc[(i) + ido * ((1) + cdim * (k)) + 1]
                        + cc[(i) + ido * ((4) + cdim * (k)) + 1];
                    t4[0] = cc[(i) + ido * ((1) + cdim * (k))] - cc[(i) + ido * ((4) + cdim * (k))];
                    t4[1] = cc[(i) + ido * ((1) + cdim * (k)) + 1]
                        - cc[(i) + ido * ((4) + cdim * (k)) + 1];
                }
                {
                    t2[0] = cc[(i) + ido * ((2) + cdim * (k))] + cc[(i) + ido * ((3) + cdim * (k))];
                    t2[1] = cc[(i) + ido * ((2) + cdim * (k)) + 1]
                        + cc[(i) + ido * ((3) + cdim * (k)) + 1];
                    t3[0] = cc[(i) + ido * ((2) + cdim * (k))] - cc[(i) + ido * ((3) + cdim * (k))];
                    t3[1] = cc[(i) + ido * ((2) + cdim * (k)) + 1]
                        - cc[(i) + ido * ((3) + cdim * (k)) + 1];
                }
                ch[(i) + ido * ((k) + l1 * (0))] = t0[0] + t1[0] + t2[0];
                ch[(i) + ido * ((k) + l1 * (0)) + 1] = t0[1] + t1[1] + t2[1];
                {
                    let mut ca: [f64; 2] = [0.0; 2];
                    let mut cb: [f64; 2] = [0.0; 2];
                    let mut da: [f64; 2] = [0.0; 2];
                    let mut db: [f64; 2] = [0.0; 2];
                    ca[0] = t0[0] + tw1r * t1[0] + tw2r * t2[0];
                    ca[1] = t0[1] + tw1r * t1[1] + tw2r * t2[1];
                    cb[1] = tw1i * t4[0] + tw2i * t3[0];
                    cb[0] = -(tw1i * t4[1] + tw2i * t3[1]);
                    {
                        da[0] = ca[0] + cb[0];
                        da[1] = ca[1] + cb[1];
                        db[0] = ca[0] - cb[0];
                        db[1] = ca[1] - cb[1];
                    }
                    {
                        ch[(i) + ido * ((k) + l1 * (1))] = wa[(i) - 1 + (1 - 1) * (ido - 1)]
                            * da[0]
                            - wa[(i) - 1 + (1 - 1) * (ido - 1) + 1] * da[1];
                        ch[(i) + ido * ((k) + l1 * (1)) + 1] = wa[(i) - 1 + (1 - 1) * (ido - 1)]
                            * da[1]
                            + wa[(i) - 1 + (1 - 1) * (ido - 1) + 1] * da[0];
                    }
                    {
                        ch[(i) + ido * ((k) + l1 * (4))] = wa[(i) - 1 + (4 - 1) * (ido - 1)]
                            * db[0]
                            - wa[(i) - 1 + (4 - 1) * (ido - 1) + 1] * db[1];
                        ch[(i) + ido * ((k) + l1 * (4)) + 1] = wa[(i) - 1 + (4 - 1) * (ido - 1)]
                            * db[1]
                            + wa[(i) - 1 + (4 - 1) * (ido - 1) + 1] * db[0];
                    }
                }
                {
                    let mut ca: [f64; 2] = [0.0; 2];
                    let mut cb: [f64; 2] = [0.0; 2];
                    let mut da: [f64; 2] = [0.0; 2];
                    let mut db: [f64; 2] = [0.0; 2];
                    ca[0] = t0[0] + tw2r * t1[0] + tw1r * t2[0];
                    ca[1] = t0[1] + tw2r * t1[1] + tw1r * t2[1];
                    cb[1] = tw2i * t4[0] - tw1i * t3[0];
                    cb[0] = -(tw2i * t4[1] - tw1i * t3[1]);
                    {
                        da[0] = ca[0] + cb[0];
                        da[1] = ca[1] + cb[1];
                        db[0] = ca[0] - cb[0];
                        db[1] = ca[1] - cb[1];
                    }
                    {
                        ch[(i) + ido * ((k) + l1 * (2))] = wa[(i) - 1 + (2 - 1) * (ido - 1)]
                            * da[0]
                            - wa[(i) - 1 + (2 - 1) * (ido - 1) + 1] * da[1];
                        ch[(i) + ido * ((k) + l1 * (2)) + 1] = wa[(i) - 1 + (2 - 1) * (ido - 1)]
                            * da[1]
                            + wa[(i) - 1 + (2 - 1) * (ido - 1) + 1] * da[0];
                    }
                    {
                        ch[(i) + ido * ((k) + l1 * (3))] = wa[(i) - 1 + (3 - 1) * (ido - 1)]
                            * db[0]
                            - wa[(i) - 1 + (3 - 1) * (ido - 1) + 1] * db[1];
                        ch[(i) + ido * ((k) + l1 * (3)) + 1] = wa[(i) - 1 + (3 - 1) * (ido - 1)]
                            * db[1]
                            + wa[(i) - 1 + (3 - 1) * (ido - 1) + 1] * db[0];
                    }
                }
            }
        }
    }
}

fn pass5f(ido: usize, l1: usize, cc: &[f64], ch: &mut [f64], wa: &[f64]) {
    let cdim: usize = 5;
    let tw1r: f64 = 0.3090169943749474241;
    let tw1i: f64 = -0.95105651629515357212;
    let tw2r: f64 = -0.8090169943749474241;
    let tw2i: f64 = -0.58778525229247312917;

    if ido == 1 {
        for k in 0..l1 {
            let mut t0: [f64; 2] = [0.0; 2];
            t0[0] = cc[(0) + ido * ((0) + cdim * (k))];
            t0[1] = cc[(0) + ido * ((0) + cdim * (k)) + 1];
            let mut t1: [f64; 2] = [0.0; 2];
            let mut t2: [f64; 2] = [0.0; 2];
            let mut t3: [f64; 2] = [0.0; 2];
            let mut t4: [f64; 2] = [0.0; 2];
            {
                t1[0] = cc[(0) + ido * ((1) + cdim * (k))] + cc[(0) + ido * ((4) + cdim * (k))];
                t1[1] =
                    cc[(0) + ido * ((1) + cdim * (k)) + 1] + cc[(0) + ido * ((4) + cdim * (k)) + 1];
                t4[0] = cc[(0) + ido * ((1) + cdim * (k))] - cc[(0) + ido * ((4) + cdim * (k))];
                t4[1] =
                    cc[(0) + ido * ((1) + cdim * (k)) + 1] - cc[(0) + ido * ((4) + cdim * (k)) + 1];
            }
            {
                t2[0] = cc[(0) + ido * ((2) + cdim * (k))] + cc[(0) + ido * ((3) + cdim * (k))];
                t2[1] =
                    cc[(0) + ido * ((2) + cdim * (k)) + 1] + cc[(0) + ido * ((3) + cdim * (k)) + 1];
                t3[0] = cc[(0) + ido * ((2) + cdim * (k))] - cc[(0) + ido * ((3) + cdim * (k))];
                t3[1] =
                    cc[(0) + ido * ((2) + cdim * (k)) + 1] - cc[(0) + ido * ((3) + cdim * (k)) + 1];
            }
            ch[(0) + ido * ((k) + l1 * (0))] = t0[0] + t1[0] + t2[0];
            ch[(0) + ido * ((k) + l1 * (0)) + 1] = t0[1] + t1[1] + t2[1];
            {
                let mut ca: [f64; 2] = [0.0; 2];
                let mut cb: [f64; 2] = [0.0; 2];
                ca[0] = t0[0] + tw1r * t1[0] + tw2r * t2[0];
                ca[1] = t0[1] + tw1r * t1[1] + tw2r * t2[1];
                cb[1] = tw1i * t4[0] + tw2i * t3[0];
                cb[0] = -(tw1i * t4[1] + tw2i * t3[1]);
                {
                    ch[(0) + ido * ((k) + l1 * (1))] = ca[0] + cb[0];
                    ch[(0) + ido * ((k) + l1 * (1)) + 1] = ca[1] + cb[1];
                    ch[(0) + ido * ((k) + l1 * (4))] = ca[0] - cb[0];
                    ch[(0) + ido * ((k) + l1 * (4)) + 1] = ca[1] - cb[1];
                }
            }
            {
                let mut ca: [f64; 2] = [0.0; 2];
                let mut cb: [f64; 2] = [0.0; 2];
                ca[0] = t0[0] + tw2r * t1[0] + tw1r * t2[0];
                ca[1] = t0[1] + tw2r * t1[1] + tw1r * t2[1];
                cb[1] = tw2i * t4[0] - tw1i * t3[0];
                cb[0] = -(tw2i * t4[1] - tw1i * t3[1]);
                {
                    ch[(0) + ido * ((k) + l1 * (2))] = ca[0] + cb[0];
                    ch[(0) + ido * ((k) + l1 * (2)) + 1] = ca[1] + cb[1];
                    ch[(0) + ido * ((k) + l1 * (3))] = ca[0] - cb[0];
                    ch[(0) + ido * ((k) + l1 * (3)) + 1] = ca[1] - cb[1];
                }
            }
        }
    } else {
        for k in 0..l1 {
            {
                let mut t0: [f64; 2] = [0.0; 2];
                t0[0] = cc[(0) + ido * ((0) + cdim * (k))];
                t0[1] = cc[(0) + ido * ((0) + cdim * (k)) + 1];
                let mut t1: [f64; 2] = [0.0; 2];
                let mut t2: [f64; 2] = [0.0; 2];
                let mut t3: [f64; 2] = [0.0; 2];
                let mut t4: [f64; 2] = [0.0; 2];
                {
                    t1[0] = cc[(0) + ido * ((1) + cdim * (k))] + cc[(0) + ido * ((4) + cdim * (k))];
                    t1[1] = cc[(0) + ido * ((1) + cdim * (k)) + 1]
                        + cc[(0) + ido * ((4) + cdim * (k)) + 1];
                    t4[0] = cc[(0) + ido * ((1) + cdim * (k))] - cc[(0) + ido * ((4) + cdim * (k))];
                    t4[1] = cc[(0) + ido * ((1) + cdim * (k)) + 1]
                        - cc[(0) + ido * ((4) + cdim * (k)) + 1];
                }
                {
                    t2[0] = cc[(0) + ido * ((2) + cdim * (k))] + cc[(0) + ido * ((3) + cdim * (k))];
                    t2[1] = cc[(0) + ido * ((2) + cdim * (k)) + 1]
                        + cc[(0) + ido * ((3) + cdim * (k)) + 1];
                    t3[0] = cc[(0) + ido * ((2) + cdim * (k))] - cc[(0) + ido * ((3) + cdim * (k))];
                    t3[1] = cc[(0) + ido * ((2) + cdim * (k)) + 1]
                        - cc[(0) + ido * ((3) + cdim * (k)) + 1];
                }
                ch[(0) + ido * ((k) + l1 * (0))] = t0[0] + t1[0] + t2[0];
                ch[(0) + ido * ((k) + l1 * (0)) + 1] = t0[1] + t1[1] + t2[1];
                {
                    let mut ca: [f64; 2] = [0.0; 2];
                    let mut cb: [f64; 2] = [0.0; 2];
                    ca[0] = t0[0] + tw1r * t1[0] + tw2r * t2[0];
                    ca[1] = t0[1] + tw1r * t1[1] + tw2r * t2[1];
                    cb[1] = tw1i * t4[0] + tw2i * t3[0];
                    cb[0] = -(tw1i * t4[1] + tw2i * t3[1]);
                    {
                        ch[(0) + ido * ((k) + l1 * (1))] = ca[0] + cb[0];
                        ch[(0) + ido * ((k) + l1 * (1)) + 1] = ca[1] + cb[1];
                        ch[(0) + ido * ((k) + l1 * (4))] = ca[0] - cb[0];
                        ch[(0) + ido * ((k) + l1 * (4)) + 1] = ca[1] - cb[1];
                    }
                }
                {
                    let mut ca: [f64; 2] = [0.0; 2];
                    let mut cb: [f64; 2] = [0.0; 2];
                    ca[0] = t0[0] + tw2r * t1[0] + tw1r * t2[0];
                    ca[1] = t0[1] + tw2r * t1[1] + tw1r * t2[1];
                    cb[1] = tw2i * t4[0] - tw1i * t3[0];
                    cb[0] = -(tw2i * t4[1] - tw1i * t3[1]);
                    {
                        ch[(0) + ido * ((k) + l1 * (2))] = ca[0] + cb[0];
                        ch[(0) + ido * ((k) + l1 * (2)) + 1] = ca[1] + cb[1];
                        ch[(0) + ido * ((k) + l1 * (3))] = ca[0] - cb[0];
                        ch[(0) + ido * ((k) + l1 * (3)) + 1] = ca[1] - cb[1];
                    }
                }
            }
            for i in 1..ido {
                let mut t0: [f64; 2] = [0.0; 2];
                t0[0] = cc[(i) + ido * ((0) + cdim * (k))];
                t0[1] = cc[(i) + ido * ((0) + cdim * (k)) + 1];
                let mut t1: [f64; 2] = [0.0; 2];
                let mut t2: [f64; 2] = [0.0; 2];
                let mut t3: [f64; 2] = [0.0; 2];
                let mut t4: [f64; 2] = [0.0; 2];
                {
                    t1[0] = cc[(i) + ido * ((1) + cdim * (k))] + cc[(i) + ido * ((4) + cdim * (k))];
                    t1[1] = cc[(i) + ido * ((1) + cdim * (k)) + 1]
                        + cc[(i) + ido * ((4) + cdim * (k)) + 1];
                    t4[0] = cc[(i) + ido * ((1) + cdim * (k))] - cc[(i) + ido * ((4) + cdim * (k))];
                    t4[1] = cc[(i) + ido * ((1) + cdim * (k)) + 1]
                        - cc[(i) + ido * ((4) + cdim * (k)) + 1];
                }
                {
                    t2[0] = cc[(i) + ido * ((2) + cdim * (k))] + cc[(i) + ido * ((3) + cdim * (k))];
                    t2[1] = cc[(i) + ido * ((2) + cdim * (k)) + 1]
                        + cc[(i) + ido * ((3) + cdim * (k)) + 1];
                    t3[0] = cc[(i) + ido * ((2) + cdim * (k))] - cc[(i) + ido * ((3) + cdim * (k))];
                    t3[1] = cc[(i) + ido * ((2) + cdim * (k)) + 1]
                        - cc[(i) + ido * ((3) + cdim * (k)) + 1];
                }
                ch[(i) + ido * ((k) + l1 * (0))] = t0[0] + t1[0] + t2[0];
                ch[(i) + ido * ((k) + l1 * (0)) + 1] = t0[1] + t1[1] + t2[1];
                {
                    let mut ca: [f64; 2] = [0.0; 2];
                    let mut cb: [f64; 2] = [0.0; 2];
                    let mut da: [f64; 2] = [0.0; 2];
                    let mut db: [f64; 2] = [0.0; 2];
                    ca[0] = t0[0] + tw1r * t1[0] + tw2r * t2[0];
                    ca[1] = t0[1] + tw1r * t1[1] + tw2r * t2[1];
                    cb[1] = tw1i * t4[0] + tw2i * t3[0];
                    cb[0] = -(tw1i * t4[1] + tw2i * t3[1]);
                    {
                        da[0] = ca[0] + cb[0];
                        da[1] = ca[1] + cb[1];
                        db[0] = ca[0] - cb[0];
                        db[1] = ca[1] - cb[1];
                    }
                    {
                        ch[(i) + ido * ((k) + l1 * (1))] = wa[(i) - 1 + (1 - 1) * (ido - 1)]
                            * da[0]
                            + wa[(i) - 1 + (1 - 1) * (ido - 1) + 1] * da[1];
                        ch[(i) + ido * ((k) + l1 * (1)) + 1] = wa[(i) - 1 + (1 - 1) * (ido - 1)]
                            * da[1]
                            - wa[(i) - 1 + (1 - 1) * (ido - 1) + 1] * da[0];
                    }
                    {
                        ch[(i) + ido * ((k) + l1 * (4))] = wa[(i) - 1 + (4 - 1) * (ido - 1)]
                            * db[0]
                            + wa[(i) - 1 + (4 - 1) * (ido - 1) + 1] * db[1];
                        ch[(i) + ido * ((k) + l1 * (4)) + 1] = wa[(i) - 1 + (4 - 1) * (ido - 1)]
                            * db[1]
                            - wa[(i) - 1 + (4 - 1) * (ido - 1) + 1] * db[0];
                    }
                }
                {
                    let mut ca: [f64; 2] = [0.0; 2];
                    let mut cb: [f64; 2] = [0.0; 2];
                    let mut da: [f64; 2] = [0.0; 2];
                    let mut db: [f64; 2] = [0.0; 2];
                    ca[0] = t0[0] + tw2r * t1[0] + tw1r * t2[0];
                    ca[1] = t0[1] + tw2r * t1[1] + tw1r * t2[1];
                    cb[1] = tw2i * t4[0] - tw1i * t3[0];
                    cb[0] = -(tw2i * t4[1] - tw1i * t3[1]);
                    {
                        da[0] = ca[0] + cb[0];
                        da[1] = ca[1] + cb[i];
                        db[0] = ca[0] - cb[0];
                        db[1] = ca[1] - cb[1];
                    }
                    {
                        ch[(i) + ido * ((k) + l1 * (2))] = wa[(i) - 1 + (2 - 1) * (ido - 1)]
                            * da[0]
                            + wa[(i) - 1 + (2 - 1) * (ido - 1) + 1] * da[1];
                        ch[(i) + ido * ((k) + l1 * (2)) + 1] = wa[(i) - 1 + (2 - 1) * (ido - 1)]
                            * da[1]
                            - wa[(i) - 1 + (2 - 1) * (ido - 1) + 1] * da[0];
                    }
                    {
                        ch[(i) + ido * ((k) + l1 * (3))] = wa[(i) - 1 + (3 - 1) * (ido - 1)]
                            * db[0]
                            + wa[(i) - 1 + (3 - 1) * (ido - 1) + 1] * db[1];
                        ch[(i) + ido * ((k) + l1 * (3)) + 1] = wa[(i) - 1 + (3 - 1) * (ido - 1)]
                            * db[1]
                            - wa[(i) - 1 + (3 - 1) * (ido - 1) + 1] * db[0];
                    }
                }
            }
        }
    }
}

fn pass7(ido: usize, l1: usize, cc: &[f64], ch: &mut [f64], wa: &[f64], sign: i64) {
    let cdim: usize = 7;
    let tw1r: f64 = 0.623489801858733530525;
    let tw1i: f64 = (sign as f64) * 0.7818314824680298087084;
    let tw2r = -0.222520933956314404289;
    let tw2i: f64 = (sign as f64) * 0.9749279121818236070181;
    let tw3r: f64 = -0.9009688679024191262361;
    let tw3i: f64 = (sign as f64) * 0.4338837391175581204758;

    if ido == 1 {
        for k in 0..l1 {
            let mut t1: [f64; 2] = [0.0; 2];
            t1[0] = cc[(0) + ido * ((0) + cdim * (k))];
            t1[1] = cc[(0) + ido * ((0) + cdim * (k)) + 1];
            let mut t2: [f64; 2] = [0.0; 2];
            let mut t3: [f64; 2] = [0.0; 2];
            let mut t4: [f64; 2] = [0.0; 2];
            let mut t5: [f64; 2] = [0.0; 2];
            let mut t6: [f64; 2] = [0.0; 2];
            let mut t7: [f64; 2] = [0.0; 2];
            {
                t2[0] = cc[(0) + ido * ((1) + cdim * (k))] + cc[(0) + ido * ((6) + cdim * (k))];
                t2[1] =
                    cc[(0) + ido * ((1) + cdim * (k)) + 1] + cc[(0) + ido * ((6) + cdim * (k)) + 1];
                t7[0] = cc[(0) + ido * ((1) + cdim * (k))] - cc[(0) + ido * ((6) + cdim * (k))];
                t7[1] =
                    cc[(0) + ido * ((1) + cdim * (k)) + 1] - cc[(0) + ido * ((6) + cdim * (k)) + 1];
            }
            {
                t3[0] = cc[(0) + ido * ((2) + cdim * (k))] + cc[(0) + ido * ((5) + cdim * (k))];
                t3[1] =
                    cc[(0) + ido * ((2) + cdim * (k)) + 1] + cc[(0) + ido * ((5) + cdim * (k)) + 1];
                t6[0] = cc[(0) + ido * ((2) + cdim * (k))] - cc[(0) + ido * ((5) + cdim * (k))];
                t6[1] =
                    cc[(0) + ido * ((2) + cdim * (k)) + 1] - cc[(0) + ido * ((5) + cdim * (k)) + 1];
            }
            {
                t4[0] = cc[(0) + ido * ((3) + cdim * (k))] + cc[(0) + ido * ((4) + cdim * (k))];
                t4[1] =
                    cc[(0) + ido * ((3) + cdim * (k)) + 1] + cc[(0) + ido * ((4) + cdim * (k)) + 1];
                t5[0] = cc[(0) + ido * ((3) + cdim * (k))] - cc[(0) + ido * ((4) + cdim * (k))];
                t5[1] =
                    cc[(0) + ido * ((3) + cdim * (k)) + 1] - cc[(0) + ido * ((4) + cdim * (k)) + 1];
            }
            ch[(0) + ido * ((k) + l1 * (0))] = t1[0] + t2[0] + t3[0] + t4[0];
            ch[(0) + ido * ((k) + l1 * (0)) + 1] = t1[1] + t2[1] + t3[1] + t4[1];
            {
                let mut ca: [f64; 2] = [0.0; 2];
                let mut cb: [f64; 2] = [0.0; 2];
                ca[0] = t1[0] + tw1r * t2[0] + tw2r * t3[0] + tw3r * t4[0];
                ca[1] = t1[1] + tw1r * t2[1] + tw2r * t3[1] + tw3r * t4[1];
                cb[1] = tw1i * t7[0] + tw2i * t6[0] + tw3i * t5[0];
                cb[0] = -(tw1i * t7[1] + tw2i * t6[1] + tw3i * t5[1]);
                {
                    ch[(0) + ido * ((k) + l1 * (1))] = ca[0] + cb[0];
                    ch[(0) + ido * ((k) + l1 * (1)) + 1] = ca[1] + cb[1];
                    ch[(0) + ido * ((k) + l1 * (6))] = ca[0] - cb[0];
                    ch[(0) + ido * ((k) + l1 * (6)) + 1] = ca[1] - cb[1];
                }
            }
            {
                let mut ca: [f64; 2] = [0.0; 2];
                let mut cb: [f64; 2] = [0.0; 2];
                ca[0] = t1[0] + tw2r * t2[0] + tw3r * t3[0] + tw1r * t4[0];
                ca[1] = t1[1] + tw2r * t2[1] + tw3r * t3[1] + tw1r * t4[1];
                cb[1] = tw2i * t7[0] - tw3i * t6[0] - tw1i * t5[0];
                cb[0] = -(tw2i * t7[1] - tw3i * t6[1] - tw1i * t5[1]);
                {
                    ch[(0) + ido * ((k) + l1 * (2))] = ca[0] + cb[0];
                    ch[(0) + ido * ((k) + l1 * (2)) + 1] = ca[1] + cb[1];
                    ch[(0) + ido * ((k) + l1 * (5))] = ca[0] - cb[0];
                    ch[(0) + ido * ((k) + l1 * (5)) + 1] = ca[1] - cb[1];
                }
            }
            {
                let mut ca: [f64; 2] = [0.0; 2];
                let mut cb: [f64; 2] = [0.0; 2];
                ca[0] = t1[0] + tw3r * t2[0] + tw1r * t3[0] + tw2r * t4[0];
                ca[1] = t1[1] + tw3r * t2[1] + tw1r * t3[1] + tw2r * t4[1];
                cb[1] = tw3i * t7[0] - tw1i * t6[0] + tw2i * t5[0];
                cb[0] = -(tw3i * t7[1] - tw1i * t6[1] + tw2i * t5[1]);
                {
                    ch[(0) + ido * ((k) + l1 * (3))] = ca[0] + cb[0];
                    ch[(0) + ido * ((k) + l1 * (3)) + 1] = ca[1] + cb[1];
                    ch[(0) + ido * ((k) + l1 * (4))] = ca[0] - cb[0];
                    ch[(0) + ido * ((k) + l1 * (4)) + 1] = ca[1] - cb[1];
                }
            }
        }
    } else {
        for k in 0..l1 {
            {
                let mut t1: [f64; 2] = [0.0; 2];
                t1[0] = cc[(0) + ido * ((0) + cdim * (k))];
                t1[1] = cc[(0) + ido * ((0) + cdim * (k)) + 1];
                let mut t2: [f64; 2] = [0.0; 2];
                let mut t3: [f64; 2] = [0.0; 2];
                let mut t4: [f64; 2] = [0.0; 2];
                let mut t5: [f64; 2] = [0.0; 2];
                let mut t6: [f64; 2] = [0.0; 2];
                let mut t7: [f64; 2] = [0.0; 2];
                {
                    t2[0] = cc[(0) + ido * ((1) + cdim * (k))] + cc[(0) + ido * ((6) + cdim * (k))];
                    t2[1] = cc[(0) + ido * ((1) + cdim * (k)) + 1]
                        + cc[(0) + ido * ((6) + cdim * (k)) + 1];
                    t7[0] = cc[(0) + ido * ((1) + cdim * (k))] - cc[(0) + ido * ((6) + cdim * (k))];
                    t7[1] = cc[(0) + ido * ((1) + cdim * (k)) + 1]
                        - cc[(0) + ido * ((6) + cdim * (k)) + 1];
                }
                {
                    t3[0] = cc[(0) + ido * ((2) + cdim * (k))] + cc[(0) + ido * ((5) + cdim * (k))];
                    t3[1] = cc[(0) + ido * ((2) + cdim * (k)) + 1]
                        + cc[(0) + ido * ((5) + cdim * (k)) + 1];
                    t6[0] =
                        cc[(0) + ido * ((2) + cdim * (k)) + 1] - cc[(0) + ido * ((5) + cdim * (k))];
                    t6[1] = cc[(0) + ido * ((2) + cdim * (k)) + 1]
                        - cc[(0) + ido * ((5) + cdim * (k)) + 1];
                }
                {
                    t4[0] = cc[(0) + ido * ((3) + cdim * (k))] + cc[(0) + ido * ((4) + cdim * (k))];
                    t4[1] = cc[(0) + ido * ((3) + cdim * (k)) + 1]
                        + cc[(0) + ido * ((4) + cdim * (k)) + 1];
                    t5[0] = cc[(0) + ido * ((3) + cdim * (k))] - cc[(0) + ido * ((4) + cdim * (k))];
                    t5[1] = cc[(0) + ido * ((3) + cdim * (k)) + 1]
                        - cc[(0) + ido * ((4) + cdim * (k)) + 1];
                }
                ch[(0) + ido * ((k) + l1 * (0))] = t1[0] + t2[0] + t3[0] + t4[0];
                ch[(0) + ido * ((k) + l1 * (0)) + 1] = t1[1] + t2[1] + t3[1] + t4[1];
                {
                    let mut ca: [f64; 2] = [0.0; 2];
                    let mut cb: [f64; 2] = [0.0; 2];
                    ca[0] = t1[0] + tw1r * t2[0] + tw2r * t3[0] + tw3r * t4[0];
                    ca[1] = t1[1] + tw1r * t2[1] + tw2r * t3[1] + tw3r * t4[1];
                    cb[1] = tw1i * t7[0] + tw2i * t6[0] + tw3i * t5[0];
                    cb[0] = -(tw1i * t7[1] + tw2i * t6[1] + tw3i * t5[1]);
                    {
                        ch[(0) + ido * ((k) + l1 * (1))] = ca[0] + cb[0];
                        ch[(0) + ido * ((k) + l1 * (1)) + 1] = ca[1] + cb[1];
                        ch[(0) + ido * ((k) + l1 * (6))] = ca[0] - cb[0];
                        ch[(0) + ido * ((k) + l1 * (6)) + 1] = ca[1] - cb[1];
                    }
                }
                {
                    let mut ca: [f64; 2] = [0.0; 2];
                    let mut cb: [f64; 2] = [0.0; 2];
                    ca[0] = t1[0] + tw2r * t2[0] + tw3r * t3[0] + tw1r * t4[0];
                    ca[1] = t1[1] + tw2r * t2[1] + tw3r * t3[1] + tw1r * t4[1];
                    cb[1] = tw2i * t7[0] - tw3i * t6[0] - tw1i * t5[0];
                    cb[0] = -(tw2i * t7[1] - tw3i * t6[1] - tw1i * t5[1]);
                    {
                        ch[(0) + ido * ((k) + l1 * (2))] = ca[0] + cb[0];
                        ch[(0) + ido * ((k) + l1 * (2)) + 1] = ca[1] + cb[1];
                        ch[(0) + ido * ((k) + l1 * (5))] = ca[0] - cb[0];
                        ch[(0) + ido * ((k) + l1 * (5)) + 1] = ca[1] - cb[1];
                    }
                }
                {
                    let mut ca: [f64; 2] = [0.0; 2];
                    let mut cb: [f64; 2] = [0.0; 2];
                    ca[0] = t1[0] + tw3r * t2[0] + tw1r * t3[0] + tw2r * t4[0];
                    ca[1] = t1[1] + tw3r * t2[1] + tw1r * t3[1] + tw2r * t4[1];
                    cb[1] = tw3i * t7[0] - tw1i * t6[0] + tw2i * t5[0];
                    cb[0] = -(tw3i * t7[1] - tw1i * t6[1] + tw2i * t5[1]);
                    {
                        ch[(0) + ido * ((k) + l1 * (3))] = ca[0] + cb[0];
                        ch[(0) + ido * ((k) + l1 * (3)) + 1] = ca[1] + cb[1];
                        ch[(0) + ido * ((k) + l1 * (4))] = ca[0] - cb[0];
                        ch[(0) + ido * ((k) + l1 * (4)) + 1] = ca[1] - cb[1];
                    }
                }
            }
            for i in 1..ido {
                let mut t1: [f64; 2] = [0.0; 2];
                t1[0] = cc[(i) + ido * ((0) + cdim * (k))];
                t1[1] = cc[(i) + ido * ((0) + cdim * (k)) + 1];
                let mut t2: [f64; 2] = [0.0; 2];
                let mut t3: [f64; 2] = [0.0; 2];
                let mut t4: [f64; 2] = [0.0; 2];
                let mut t5: [f64; 2] = [0.0; 2];
                let mut t6: [f64; 2] = [0.0; 2];
                let mut t7: [f64; 2] = [0.0; 2];
                {
                    t2[0] = cc[(i) + ido * ((1) + cdim * (k))] + cc[(i) + ido * ((6) + cdim * (k))];
                    t2[1] = cc[(i) + ido * ((1) + cdim * (k)) + 1]
                        + cc[(i) + ido * ((6) + cdim * (k)) + 1];
                    t7[0] = cc[(i) + ido * ((1) + cdim * (k))] - cc[(i) + ido * ((6) + cdim * (k))];
                    t7[1] = cc[(i) + ido * ((1) + cdim * (k)) + 1]
                        - cc[(i) + ido * ((6) + cdim * (k)) + 1];
                }
                {
                    t3[0] = cc[(i) + ido * ((2) + cdim * (k))] + cc[(i) + ido * ((5) + cdim * (k))];
                    t3[1] = cc[(i) + ido * ((2) + cdim * (k)) + 1]
                        + cc[(i) + ido * ((5) + cdim * (k)) + 1];
                    t6[0] = cc[(i) + ido * ((2) + cdim * (k))] - cc[(i) + ido * ((5) + cdim * (k))];
                    t6[1] = cc[(i) + ido * ((2) + cdim * (k)) + 1]
                        - cc[(i) + ido * ((5) + cdim * (k)) + 1];
                }
                {
                    t4[0] = cc[(i) + ido * ((3) + cdim * (k))] + cc[(i) + ido * ((4) + cdim * (k))];
                    t4[1] = cc[(i) + ido * ((3) + cdim * (k)) + 1]
                        + cc[(i) + ido * ((4) + cdim * (k)) + 1];
                    t5[0] = cc[(i) + ido * ((3) + cdim * (k))] - cc[(i) + ido * ((4) + cdim * (k))];
                    t5[1] = cc[(i) + ido * ((3) + cdim * (k)) + 1]
                        - cc[(i) + ido * ((4) + cdim * (k)) + 1];
                }
                ch[(i) + ido * ((k) + l1 * (0))] = t1[0] + t2[0] + t3[0] + t4[0];
                ch[(i) + ido * ((k) + l1 * (0)) + 1] = t1[1] + t2[1] + t3[1] + t4[1];
                {
                    let mut da: [f64; 2] = [0.0; 2];
                    let mut db: [f64; 2] = [0.0; 2];
                    {
                        let mut ca: [f64; 2] = [0.0; 2];
                        let mut cb: [f64; 2] = [0.0; 2];
                        ca[0] = t1[0] + tw1r * t2[0] + tw2r * t3[0] + tw3r * t4[0];
                        ca[1] = t1[1] + tw1r * t2[1] + tw2r * t3[1] + tw3r * t4[1];
                        cb[1] = tw1i * t7[0] + tw2i * t6[0] + tw3i * t5[0];
                        cb[0] = -(tw1i * t7[1] + tw2i * t6[1] + tw3i * t5[1]);
                        {
                            da[0] = ca[0] + cb[0];
                            da[1] = ca[1] + cb[1];
                            db[0] = ca[0] - cb[0];
                            db[1] = ca[1] - cb[1];
                        }
                    }
                    {
                        ch[(i) + ido * ((k) + l1 * (1))] = wa[(i) - 1 + (1 - 1) * (ido - 1)]
                            * da[0]
                            - (sign as f64) * wa[(i) - 1 + (1 - 1) * (ido - 1) + 1] * da[1];
                        ch[(i) + ido * ((k) + l1 * (1)) + 1] = wa[(i) - 1 + (1 - 1) * (ido - 1)]
                            * da[1]
                            + (sign as f64) * wa[(i) - 1 + (1 - 1) * (ido - 1) + 1] * da[0];
                    }
                    {
                        ch[(i) + ido * ((k) + l1 * (6))] = wa[(i) - 1 + (6 - 1) * (ido - 1)]
                            * db[0]
                            - (sign as f64) * wa[(i) - 1 + (6 - 1) * (ido - 1) + 1] * db[1];
                        ch[(i) + ido * ((k) + l1 * (6)) + 1] = wa[(i) - 1 + (6 - 1) * (ido - 1)]
                            * db[1]
                            + (sign as f64) * wa[(i) - 1 + (6 - 1) * (ido - 1) + 1] * db[0];
                    }
                }
                {
                    let mut da: [f64; 2] = [0.0; 2];
                    let mut db: [f64; 2] = [0.0; 2];
                    {
                        let mut ca: [f64; 2] = [0.0; 2];
                        let mut cb: [f64; 2] = [0.0; 2];
                        ca[0] = t1[0] + tw2r * t2[0] + tw3r * t3[0] + tw1r * t4[0];
                        ca[1] = t1[1] + tw2r * t2[1] + tw3r * t3[1] + tw1r * t4[1];
                        cb[1] = tw2i * t7[0] - tw3i * t6[0] - tw1i * t5[0];
                        cb[0] = -(tw2i * t7[1] - tw3i * t6[1] - tw1i * t5[1]);
                        {
                            da[0] = ca[0] + cb[0];
                            da[1] = ca[1] + cb[1];
                            db[0] = ca[0] - cb[0];
                            db[1] = ca[1] - cb[1];
                        }
                    }
                    {
                        ch[(i) + ido * ((k) + l1 * (2))] = wa[(i) - 1 + (2 - 1) * (ido - 1)]
                            * da[0]
                            - (sign as f64) * wa[(i) - 1 + (2 - 1) * (ido - 1) + 1] * da[1];
                        ch[(i) + ido * ((k) + l1 * (2)) + 1] = wa[(i) - 1 + (2 - 1) * (ido - 1)]
                            * da[1]
                            + (sign as f64) * wa[(i) - 1 + (2 - 1) * (ido - 1) + 1] * da[0];
                    }
                    {
                        ch[(i) + ido * ((k) + l1 * (5))] = wa[(i) - 1 + (5 - 1) * (ido - 1)]
                            * db[0]
                            - (sign as f64) * wa[(i) - 1 + (5 - 1) * (ido - 1) + 1] * db[1];
                        ch[(i) + ido * ((k) + l1 * (5)) + 1] = wa[(i) - 1 + (5 - 1) * (ido - 1)]
                            * db[1]
                            + (sign as f64) * wa[(i) - 1 + (5 - 1) * (ido - 1) + 1] * db[0];
                    }
                }
                {
                    let mut da: [f64; 2] = [0.0; 2];
                    let mut db: [f64; 2] = [0.0; 2];
                    {
                        let mut ca: [f64; 2] = [0.0; 2];
                        let mut cb: [f64; 2] = [0.0; 2];
                        ca[0] = t1[0] + tw3r * t2[0] + tw1r * t3[0] + tw2r * t4[0];
                        ca[1] = t1[1] + tw3r * t2[1] + tw1r * t3[1] + tw2r * t4[1];
                        cb[1] = tw3i * t7[0] - tw1i * t6[0] + tw2i * t5[0];
                        cb[0] = -(tw3i * t7[1] - tw1i * t6[1] + tw2i * t5[1]);
                        {
                            da[0] = ca[0] + cb[0];
                            da[1] = ca[1] + cb[1];
                            db[0] = ca[0] - cb[0];
                            db[1] = ca[1] - cb[1];
                        }
                    }
                    {
                        ch[(i) + ido * ((k) + l1 * (3))] = wa[(i) - 1 + (3 - 1) * (ido - 1)]
                            * da[0]
                            - (sign as f64) * wa[(i) - 1 + (3 - 1) * (ido - 1) + 1] * da[1];
                        ch[(i) + ido * ((k) + l1 * (3)) + 1] = wa[(i) - 1 + (3 - 1) * (ido - 1)]
                            * da[1]
                            + (sign as f64) * wa[(i) - 1 + (3 - 1) * (ido - 1) + 1] * da[0];
                    }
                    {
                        ch[(i) + ido * ((k) + l1 * (4))] = wa[(i) - 1 + (4 - 1) * (ido - 1)]
                            * db[0]
                            - (sign as f64) * wa[(i) - 1 + (4 - 1) * (ido - 1) + 1] * db[1];
                        ch[(i) + ido * ((k) + l1 * (4)) + 1] = wa[(i) - 1 + (4 - 1) * (ido - 1)]
                            * db[1]
                            + (sign as f64) * wa[(i) - 1 + (4 - 1) * (ido - 1) + 1] * db[0];
                    }
                }
            }
        }
    }
}

fn pass11(ido: usize, l1: usize, cc: &[f64], ch: &mut [f64], wa: &[f64], sign: i64) {
    let cdim: usize = 11;
    let tw1r: f64 = 0.8412535328311811688618;
    let tw1i: f64 = (sign as f64) * 0.5406408174555975821076;
    let tw2r: f64 = 0.4154150130018864255293;
    let tw2i: f64 = (sign as f64) * 0.9096319953545183714117;
    let tw3r: f64 = -0.1423148382732851404438;
    let tw3i: f64 = (sign as f64) * 0.9898214418809327323761;
    let tw4r: f64 = -0.6548607339452850640569;
    let tw4i: f64 = (sign as f64) * 0.755749574354258283774;
    let tw5r: f64 = -0.9594929736144973898904;
    let tw5i: f64 = (sign as f64) * 0.2817325568414296977114;

    if ido == 1 {
        for k in 0..l1 {
            let mut t1: [f64; 2] = [0.0; 2];
            t1[0] = cc[(0) + ido * ((0) + cdim * (k))];
            t1[1] = cc[(0) + ido * ((0) + cdim * (k)) + 1];
            let mut t2: [f64; 2] = [0.0; 2];
            let mut t3: [f64; 2] = [0.0; 2];
            let mut t4: [f64; 2] = [0.0; 2];
            let mut t5: [f64; 2] = [0.0; 2];
            let mut t6: [f64; 2] = [0.0; 2];
            let mut t7: [f64; 2] = [0.0; 2];
            let mut t8: [f64; 2] = [0.0; 2];
            let mut t9: [f64; 2] = [0.0; 2];
            let mut t10: [f64; 2] = [0.0; 2];
            let mut t11: [f64; 2] = [0.0; 2];
            {
                t2[0] = &cc[(0) + ido * ((1) + cdim * (k))] + cc[(0) + ido * ((10) + cdim * (k))];
                t2[1] = &cc[(0) + ido * ((1) + cdim * (k)) + 1]
                    + cc[(0) + ido * ((10) + cdim * (k)) + 1];
                t11[0] = &cc[(0) + ido * ((1) + cdim * (k))] - cc[(0) + ido * ((10) + cdim * (k))];
                t11[1] = &cc[(0) + ido * ((1) + cdim * (k)) + 1]
                    - cc[(0) + ido * ((10) + cdim * (k)) + 1];
            }
            {
                t3[0] = &cc[(0) + ido * ((2) + cdim * (k))] + cc[(0) + ido * ((9) + cdim * (k))];
                t3[1] = &cc[(0) + ido * ((2) + cdim * (k)) + 1]
                    + cc[(0) + ido * ((9) + cdim * (k)) + 1];
                t10[0] = &cc[(0) + ido * ((2) + cdim * (k))] - cc[(0) + ido * ((9) + cdim * (k))];
                t10[1] = &cc[(0) + ido * ((2) + cdim * (k)) + 1]
                    - cc[(0) + ido * ((9) + cdim * (k)) + 1];
            }
            {
                t4[0] = &cc[(0) + ido * ((3) + cdim * (k))] + cc[(0) + ido * ((8) + cdim * (k))];
                t4[1] = &cc[(0) + ido * ((3) + cdim * (k)) + 1]
                    + cc[(0) + ido * ((8) + cdim * (k)) + 1];
                t9[0] = &cc[(0) + ido * ((3) + cdim * (k))] - cc[(0) + ido * ((8) + cdim * (k))];
                t9[1] = &cc[(0) + ido * ((3) + cdim * (k)) + 1]
                    - cc[(0) + ido * ((8) + cdim * (k)) + 1];
            }
            {
                t5[0] = &cc[(0) + ido * ((4) + cdim * (k))] + cc[(0) + ido * ((7) + cdim * (k))];
                t5[1] = &cc[(0) + ido * ((4) + cdim * (k)) + 1]
                    + cc[(0) + ido * ((7) + cdim * (k)) + 1];
                t8[0] = &cc[(0) + ido * ((4) + cdim * (k))] - cc[(0) + ido * ((7) + cdim * (k))];
                t8[1] = &cc[(0) + ido * ((4) + cdim * (k)) + 1]
                    - cc[(0) + ido * ((7) + cdim * (k)) + 1];
            }
            {
                t6[0] = cc[(0) + ido * ((5) + cdim * (k))] + cc[(0) + ido * ((6) + cdim * (k))];
                t6[1] =
                    cc[(0) + ido * ((5) + cdim * (k)) + 1] + cc[(0) + ido * ((6) + cdim * (k)) + 1];
                t7[0] = cc[(0) + ido * ((5) + cdim * (k))] - cc[(0) + ido * ((6) + cdim * (k))];
                t7[1] =
                    cc[(0) + ido * ((5) + cdim * (k)) + 1] - cc[(0) + ido * ((6) + cdim * (k)) + 1];
            }
            ch[(0) + ido * ((k) + l1 * (0))] = t1[0] + t2[0] + t3[0] + t4[0] + t5[0] + t6[0];
            ch[(0) + ido * ((k) + l1 * (0)) + 1] = t1[1] + t2[1] + t3[1] + t4[1] + t5[1] + t6[1];
            {
                let mut ca: [f64; 2] = [0.0; 2];
                let mut cb: [f64; 2] = [0.0; 2];
                ca[0] = t1[0]
                    + tw1r * t2[0]
                    + tw2r * t3[0]
                    + tw3r * t4[0]
                    + tw4r * t5[0]
                    + tw5r * t6[0];
                ca[1] = t1[1]
                    + tw1r * t2[1]
                    + tw2r * t3[1]
                    + tw3r * t4[1]
                    + tw4r * t5[1]
                    + tw5r * t6[1];
                cb[1] = tw1i * t11[0] + tw2i * t10[01] + tw3i * t9[0] + tw4i * t8[0] + tw5i * t7[0];
                cb[0] =
                    -(tw1i * t11[1] + tw2i * t10[1] + tw3i * t9[1] + tw4i * t8[1] + tw5i * t7[1]);
                {
                    ch[(0) + ido * ((k) + l1 * (1))] = ca[0] + cb[0];
                    ch[(0) + ido * ((k) + l1 * (1)) + 1] = ca[1] + cb[1];
                    ch[(0) + ido * ((k) + l1 * (10))] = ca[0] - cb[0];
                    ch[(0) + ido * ((k) + l1 * (10)) + 1] = ca[1] - cb[1];
                }
            }
            {
                let mut ca: [f64; 2] = [0.0; 2];
                let mut cb: [f64; 2] = [0.0; 2];
                ca[0] = t1[0]
                    + tw2r * t2[0]
                    + tw4r * t3[0]
                    + tw5r * t4[0]
                    + tw3r * t5[0]
                    + tw1r * t6[0];
                ca[1] = t1[1]
                    + tw2r * t2[1]
                    + tw4r * t3[1]
                    + tw5r * t4[1]
                    + tw3r * t5[1]
                    + tw1r * t6[1];
                cb[1] = tw2i * t11[0] + tw4i * t10[0] - tw5i * t9[0] - tw3i * t8[0] - tw1i * t7[0];
                cb[0] =
                    -(tw2i * t11[1] + tw4i * t10[1] - tw5i * t9[1] - tw3i * t8[1] - tw1i * t7[1]);
                {
                    ch[(0) + ido * ((k) + l1 * (2))] = ca[0] + cb[0];
                    ch[(0) + ido * ((k) + l1 * (2)) + 1] = ca[1] + cb[1];
                    ch[(0) + ido * ((k) + l1 * (9))] = ca[0] - cb[0];
                    ch[(0) + ido * ((k) + l1 * (9)) + 1] = ca[1] - cb[1];
                }
            }
            {
                let mut ca: [f64; 2] = [0.0; 2];
                let mut cb: [f64; 2] = [0.0; 2];
                ca[0] = t1[0]
                    + tw3r * t2[0]
                    + tw5r * t3[0]
                    + tw2r * t4[0]
                    + tw1r * t5[0]
                    + tw4r * t6[0];
                ca[1] = t1[1]
                    + tw3r * t2[1]
                    + tw5r * t3[1]
                    + tw2r * t4[1]
                    + tw1r * t5[1]
                    + tw4r * t6[1];
                cb[1] = tw3i * t11[0] - tw5i * t10[0] - tw2i * t9[0] + tw1i * t8[0] + tw4i * t7[0];
                cb[0] =
                    -(tw3i * t11[1] - tw5i * t10[1] - tw2i * t9[1] + tw1i * t8[1] + tw4i * t7[1]);
                {
                    ch[(0) + ido * ((k) + l1 * (3))] = ca[0] + cb[0];
                    ch[(0) + ido * ((k) + l1 * (3)) + 1] = ca[1] + cb[1];
                    ch[(0) + ido * ((k) + l1 * (8))] = ca[0] - cb[0];
                    ch[(0) + ido * ((k) + l1 * (8)) + 1] = ca[1] - cb[1];
                }
            }
            {
                let mut ca: [f64; 2] = [0.0; 2];
                let mut cb: [f64; 2] = [0.0; 2];
                ca[0] = t1[0]
                    + tw4r * t2[0]
                    + tw3r * t3[0]
                    + tw1r * t4[0]
                    + tw5r * t5[0]
                    + tw2r * t6[0];
                ca[1] = t1[1]
                    + tw4r * t2[1]
                    + tw3r * t3[1]
                    + tw1r * t4[1]
                    + tw5r * t5[1]
                    + tw2r * t6[1];
                cb[1] = tw4i * t11[0] - tw3i * t10[0] + tw1i * t9[0] + tw5i * t8[0] - tw2i * t7[0];
                cb[0] =
                    -(tw4i * t11[1] - tw3i * t10[1] + tw1i * t9[1] + tw5i * t8[1] - tw2i * t7[1]);
                {
                    ch[(0) + ido * ((k) + l1 * (4))] = ca[0] + cb[0];
                    ch[(0) + ido * ((k) + l1 * (4)) + 1] = ca[1] + cb[1];
                    ch[(0) + ido * ((k) + l1 * (7))] = ca[0] - cb[0];
                    ch[(0) + ido * ((k) + l1 * (7)) + 1] = ca[1] - cb[1];
                }
            }
            {
                let mut ca: [f64; 2] = [0.0; 2];
                let mut cb: [f64; 2] = [0.0; 2];
                ca[0] = t1[0]
                    + tw5r * t2[0]
                    + tw1r * t3[0]
                    + tw4r * t4[0]
                    + tw2r * t5[0]
                    + tw3r * t6[0];
                ca[1] = t1[1]
                    + tw5r * t2[1]
                    + tw1r * t3[1]
                    + tw4r * t4[1]
                    + tw2r * t5[1]
                    + tw3r * t6[1];
                cb[1] = tw5i * t11[0] - tw1i * t10[0] + tw4i * t9[0] - tw2i * t8[0] + tw3i * t7[0];
                cb[0] =
                    -(tw5i * t11[1] - tw1i * t10[1] + tw4i * t9[1] - tw2i * t8[1] + tw3i * t7[1]);
                {
                    ch[(0) + ido * ((k) + l1 * (5))] = ca[0] + cb[0];
                    ch[(0) + ido * ((k) + l1 * (5)) + 1] = ca[1] + cb[1];
                    ch[(0) + ido * ((k) + l1 * (6))] = ca[0] - cb[0];
                    ch[(0) + ido * ((k) + l1 * (6)) + 1] = ca[1] - cb[1];
                }
            }
        }
    } else {
        for k in 0..l1 {
            {
                let mut t1: [f64; 2] = [0.0; 2];
                t1[0] = cc[(0) + ido * ((0) + cdim * (k))];
                t1[1] = cc[(0) + ido * ((0) + cdim * (k)) + 1];
                let mut t2: [f64; 2] = [0.0; 2];
                let mut t3: [f64; 2] = [0.0; 2];
                let mut t4: [f64; 2] = [0.0; 2];
                let mut t5: [f64; 2] = [0.0; 2];
                let mut t6: [f64; 2] = [0.0; 2];
                let mut t7: [f64; 2] = [0.0; 2];
                let mut t8: [f64; 2] = [0.0; 2];
                let mut t9: [f64; 2] = [0.0; 2];
                let mut t10: [f64; 2] = [0.0; 2];
                let mut t11: [f64; 2] = [0.0; 2];
                {
                    t2[0] =
                        cc[(0) + ido * ((1) + cdim * (k))] + cc[(0) + ido * ((10) + cdim * (k))];
                    t2[1] = cc[(0) + ido * ((1) + cdim * (k)) + 1]
                        + cc[(0) + ido * ((10) + cdim * (k)) + 1];
                    t11[0] =
                        cc[(0) + ido * ((1) + cdim * (k))] - cc[(0) + ido * ((10) + cdim * (k))];
                    t11[1] = cc[(0) + ido * ((1) + cdim * (k)) + 1]
                        - cc[(0) + ido * ((10) + cdim * (k)) + 1];
                }
                {
                    t3[0] = cc[(0) + ido * ((2) + cdim * (k))] + cc[(0) + ido * ((9) + cdim * (k))];
                    t3[1] = cc[(0) + ido * ((2) + cdim * (k)) + 1]
                        + cc[(0) + ido * ((9) + cdim * (k)) + 1];
                    t10[0] =
                        cc[(0) + ido * ((2) + cdim * (k))] - cc[(0) + ido * ((9) + cdim * (k))];
                    t10[1] = cc[(0) + ido * ((2) + cdim * (k)) + 1]
                        - cc[(0) + ido * ((9) + cdim * (k)) + 1];
                }
                {
                    t4[0] = cc[(0) + ido * ((3) + cdim * (k))] + cc[(0) + ido * ((8) + cdim * (k))];
                    t4[1] = cc[(0) + ido * ((3) + cdim * (k)) + 1]
                        + cc[(0) + ido * ((8) + cdim * (k)) + 1];
                    t9[0] = cc[(0) + ido * ((3) + cdim * (k))] - cc[(0) + ido * ((8) + cdim * (k))];
                    t9[1] = cc[(0) + ido * ((3) + cdim * (k)) + 1]
                        - cc[(0) + ido * ((8) + cdim * (k)) + 1];
                }
                {
                    t5[0] = cc[(0) + ido * ((4) + cdim * (k))] + cc[(0) + ido * ((7) + cdim * (k))];
                    t5[1] = cc[(0) + ido * ((4) + cdim * (k)) + 1]
                        + cc[(0) + ido * ((7) + cdim * (k)) + 1];
                    t8[0] = cc[(0) + ido * ((4) + cdim * (k))] - cc[(0) + ido * ((7) + cdim * (k))];
                    t8[1] = cc[(0) + ido * ((4) + cdim * (k)) + 1]
                        - cc[(0) + ido * ((7) + cdim * (k)) + 1];
                }
                {
                    t6[0] = cc[(0) + ido * ((5) + cdim * (k))] + cc[(0) + ido * ((6) + cdim * (k))];
                    t6[1] = cc[(0) + ido * ((5) + cdim * (k)) + 1]
                        + cc[(0) + ido * ((6) + cdim * (k)) + 1];
                    t7[0] = cc[(0) + ido * ((5) + cdim * (k))] - cc[(0) + ido * ((6) + cdim * (k))];
                    t7[1] = cc[(0) + ido * ((5) + cdim * (k)) + 1]
                        - cc[(0) + ido * ((6) + cdim * (k)) + 1];
                }
                ch[(0) + ido * ((k) + l1 * (0))] = t1[0] + t2[0] + t3[0] + t4[0] + t5[0] + t6[0];
                ch[(0) + ido * ((k) + l1 * (0)) + 1] =
                    t1[1] + t2[1] + t3[1] + t4[1] + t5[1] + t6[1];
                {
                    let mut ca: [f64; 2] = [0.0; 2];
                    let mut cb: [f64; 2] = [0.0; 2];
                    ca[0] = t1[0]
                        + tw1r * t2[0]
                        + tw2r * t3[0]
                        + tw3r * t4[0]
                        + tw4r * t5[0]
                        + tw5r * t6[0];
                    ca[1] = t1[1]
                        + tw1r * t2[1]
                        + tw2r * t3[1]
                        + tw3r * t4[1]
                        + tw4r * t5[1]
                        + tw5r * t6[1];
                    cb[1] =
                        tw1i * t11[0] + tw2i * t10[0] + tw3i * t9[0] + tw4i * t8[0] + tw5i * t7[0];
                    cb[0] = -(tw1i * t11[1]
                        + tw2i * t10[1]
                        + tw3i * t9[1]
                        + tw4i * t8[1]
                        + tw5i * t7[1]);
                    {
                        ch[(0) + ido * ((k) + l1 * (1))] = ca[0] + cb[0];
                        ch[(0) + ido * ((k) + l1 * (1)) + 1] = ca[1] + cb[1];
                        ch[(0) + ido * ((k) + l1 * (10))] = ca[0] - cb[0];
                        ch[(0) + ido * ((k) + l1 * (10)) + 1] = ca[1] - cb[1];
                    }
                }
                {
                    let mut ca: [f64; 2] = [0.0; 2];
                    let mut cb: [f64; 2] = [0.0; 2];
                    ca[0] = t1[0]
                        + tw2r * t2[0]
                        + tw4r * t3[0]
                        + tw5r * t4[0]
                        + tw3r * t5[0]
                        + tw1r * t6[0];
                    ca[1] = t1[1]
                        + tw2r * t2[1]
                        + tw4r * t3[1]
                        + tw5r * t4[1]
                        + tw3r * t5[1]
                        + tw1r * t6[1];
                    cb[1] =
                        tw2i * t11[0] + tw4i * t10[0] - tw5i * t9[0] - tw3i * t8[0] - tw1i * t7[0];
                    cb[0] = -(tw2i * t11[1] + tw4i * t10[1]
                        - tw5i * t9[1]
                        - tw3i * t8[1]
                        - tw1i * t7[1]);
                    {
                        ch[(0) + ido * ((k) + l1 * (2))] = ca[0] + cb[0];
                        ch[(0) + ido * ((k) + l1 * (2)) + 1] = ca[1] + cb[1];
                        ch[(0) + ido * ((k) + l1 * (9))] = ca[0] - cb[0];
                        ch[(0) + ido * ((k) + l1 * (9)) + 1] = ca[1] - cb[1];
                    }
                }
                {
                    let mut ca: [f64; 2] = [0.0; 2];
                    let mut cb: [f64; 2] = [0.0; 2];
                    ca[0] = t1[0]
                        + tw3r * t2[0]
                        + tw5r * t3[0]
                        + tw2r * t4[0]
                        + tw1r * t5[0]
                        + tw4r * t6[0];
                    ca[1] = t1[1]
                        + tw3r * t2[1]
                        + tw5r * t3[1]
                        + tw2r * t4[1]
                        + tw1r * t5[1]
                        + tw4r * t6[1];
                    cb[1] =
                        tw3i * t11[0] - tw5i * t10[0] - tw2i * t9[0] + tw1i * t8[0] + tw4i * t7[0];
                    cb[0] = -(tw3i * t11[1] - tw5i * t10[1] - tw2i * t9[1]
                        + tw1i * t8[1]
                        + tw4i * t7[1]);
                    {
                        ch[(0) + ido * ((k) + l1 * (3))] = ca[0] + cb[0];
                        ch[(0) + ido * ((k) + l1 * (3)) + 1] = ca[1] + cb[1];
                        ch[(0) + ido * ((k) + l1 * (8))] = ca[0] - cb[0];
                        ch[(0) + ido * ((k) + l1 * (8)) + 1] = ca[1] - cb[1];
                    }
                }
                {
                    let mut ca: [f64; 2] = [0.0; 2];
                    let mut cb: [f64; 2] = [0.0; 2];
                    ca[0] = t1[0]
                        + tw4r * t2[0]
                        + tw3r * t3[0]
                        + tw1r * t4[0]
                        + tw5r * t5[0]
                        + tw2r * t6[0];
                    ca[1] = t1[1]
                        + tw4r * t2[1]
                        + tw3r * t3[1]
                        + tw1r * t4[1]
                        + tw5r * t5[1]
                        + tw2r * t6[1];
                    cb[1] =
                        tw4i * t11[0] - tw3i * t10[0] + tw1i * t9[0] + tw5i * t8[0] - tw2i * t7[0];
                    cb[0] = -(tw4i * t11[1] - tw3i * t10[1] + tw1i * t9[1] + tw5i * t8[1]
                        - tw2i * t7[1]);
                    {
                        ch[(0) + ido * ((k) + l1 * (4))] = ca[0] + cb[0];
                        ch[(0) + ido * ((k) + l1 * (4)) + 1] = ca[1] + cb[1];
                        ch[(0) + ido * ((k) + l1 * (7))] = ca[0] - cb[0];
                        ch[(0) + ido * ((k) + l1 * (7)) + 1] = ca[1] - cb[1];
                    }
                }
                {
                    let mut ca: [f64; 2] = [0.0; 2];
                    let mut cb: [f64; 2] = [0.0; 2];
                    ca[0] = t1[0]
                        + tw5r * t2[0]
                        + tw1r * t3[0]
                        + tw4r * t4[0]
                        + tw2r * t5[0]
                        + tw3r * t6[0];
                    ca[1] = t1[1]
                        + tw5r * t2[1]
                        + tw1r * t3[1]
                        + tw4r * t4[1]
                        + tw2r * t5[1]
                        + tw3r * t6[1];
                    cb[1] =
                        tw5i * t11[0] - tw1i * t10[0] + tw4i * t9[0] - tw2i * t8[0] + tw3i * t7[0];
                    cb[0] = -(tw5i * t11[1] - tw1i * t10[1] + tw4i * t9[1] - tw2i * t8[1]
                        + tw3i * t7[1]);
                    {
                        ch[(0) + ido * ((k) + l1 * (5))] = ca[0] + cb[0];
                        ch[(0) + ido * ((k) + l1 * (5)) + 1] = ca[1] + cb[1];
                        ch[(0) + ido * ((k) + l1 * (6))] = ca[0] - cb[0];
                        ch[(0) + ido * ((k) + l1 * (6)) + 1] = ca[1] - cb[1];
                    }
                }
            }
            for i in 1..ido {
                let mut t1: [f64; 2] = [0.0; 2];
                t1[0] = cc[(i) + ido * ((0) + cdim * (k))];
                t1[1] = cc[(i) + ido * ((0) + cdim * (k)) + 1];
                let mut t2: [f64; 2] = [0.0; 2];
                let mut t3: [f64; 2] = [0.0; 2];
                let mut t4: [f64; 2] = [0.0; 2];
                let mut t5: [f64; 2] = [0.0; 2];
                let mut t6: [f64; 2] = [0.0; 2];
                let mut t7: [f64; 2] = [0.0; 2];
                let mut t8: [f64; 2] = [0.0; 2];
                let mut t9: [f64; 2] = [0.0; 2];
                let mut t10: [f64; 2] = [0.0; 2];
                let mut t11: [f64; 2] = [0.0; 2];
                {
                    t2[0] =
                        cc[(i) + ido * ((1) + cdim * (k))] + cc[(i) + ido * ((10) + cdim * (k))];
                    t2[1] = cc[(i) + ido * ((1) + cdim * (k)) + 1]
                        + cc[(i) + ido * ((10) + cdim * (k)) + 1];
                    t11[0] =
                        cc[(i) + ido * ((1) + cdim * (k))] - cc[(i) + ido * ((10) + cdim * (k))];
                    t11[1] = cc[(i) + ido * ((1) + cdim * (k)) + 1]
                        - cc[(i) + ido * ((10) + cdim * (k)) + 1];
                }
                {
                    t3[0] = cc[(i) + ido * ((2) + cdim * (k))] + cc[(i) + ido * ((9) + cdim * (k))];
                    t3[1] = cc[(i) + ido * ((2) + cdim * (k)) + 1]
                        + cc[(i) + ido * ((9) + cdim * (k)) + 1];
                    t10[0] =
                        cc[(i) + ido * ((2) + cdim * (k))] - cc[(i) + ido * ((9) + cdim * (k))];
                    t10[1] = cc[(i) + ido * ((2) + cdim * (k)) + 1]
                        - cc[(i) + ido * ((9) + cdim * (k)) + 1];
                }
                {
                    t4[0] = cc[(i) + ido * ((3) + cdim * (k))] + cc[(i) + ido * ((8) + cdim * (k))];
                    t4[1] = cc[(i) + ido * ((3) + cdim * (k)) + 1]
                        + cc[(i) + ido * ((8) + cdim * (k)) + 1];
                    t9[0] = cc[(i) + ido * ((3) + cdim * (k))] - cc[(i) + ido * ((8) + cdim * (k))];
                    t9[1] = cc[(i) + ido * ((3) + cdim * (k)) + 1]
                        - cc[(i) + ido * ((8) + cdim * (k)) + 1];
                }
                {
                    t5[0] = cc[(i) + ido * ((4) + cdim * (k))] + cc[(i) + ido * ((7) + cdim * (k))];
                    t5[1] = cc[(i) + ido * ((4) + cdim * (k)) + 1]
                        + cc[(i) + ido * ((7) + cdim * (k)) + 1];
                    t8[0] = cc[(i) + ido * ((4) + cdim * (k))] - cc[(i) + ido * ((7) + cdim * (k))];
                    t8[1] = cc[(i) + ido * ((4) + cdim * (k)) + 1]
                        - cc[(i) + ido * ((7) + cdim * (k)) + 1];
                }
                {
                    t6[0] = cc[(i) + ido * ((5) + cdim * (k))] + cc[(i) + ido * ((6) + cdim * (k))];
                    t6[1] = cc[(i) + ido * ((5) + cdim * (k)) + 1]
                        + cc[(i) + ido * ((6) + cdim * (k)) + 1];
                    t7[0] = cc[(i) + ido * ((5) + cdim * (k))] - cc[(i) + ido * ((6) + cdim * (k))];
                    t7[1] = cc[(i) + ido * ((5) + cdim * (k)) + 1]
                        - cc[(i) + ido * ((6) + cdim * (k)) + 1];
                }
                ch[(i) + ido * ((k) + l1 * (0))] = t1[0] + t2[0] + t3[0] + t4[0] + t5[0] + t6[0];
                ch[(i) + ido * ((k) + l1 * (0)) + 1] =
                    t1[1] + t2[1] + t3[1] + t4[1] + t5[1] + t6[1];
                {
                    let mut da: [f64; 2] = [0.0; 2];
                    let mut db: [f64; 2] = [0.0; 2];
                    {
                        let mut ca: [f64; 2] = [0.0; 2];
                        let mut cb: [f64; 2] = [0.0; 2];
                        ca[0] = t1[0]
                            + tw1r * t2[0]
                            + tw2r * t3[0]
                            + tw3r * t4[0]
                            + tw4r * t5[0]
                            + tw5r * t6[0];
                        ca[1] = t1[1]
                            + tw1r * t2[1]
                            + tw2r * t3[1]
                            + tw3r * t4[1]
                            + tw4r * t5[1]
                            + tw5r * t6[1];
                        cb[1] = tw1i * t11[0]
                            + tw2i * t10[0]
                            + tw3i * t9[0]
                            + tw4i * t8[0]
                            + tw5i * t7[0];
                        cb[0] = -(tw1i * t11[1]
                            + tw2i * t10[1]
                            + tw3i * t9[1]
                            + tw4i * t8[1]
                            + tw5i * t7[1]);
                        {
                            da[0] = ca[0] + cb[0];
                            da[1] = ca[1] + cb[1];
                            db[0] = ca[0] - cb[0];
                            db[1] = ca[1] - cb[1];
                        }
                    }
                    {
                        ch[(i) + ido * ((k) + l1 * (1))] = wa[(i) - 1 + (1 - 1) * (ido - 1)]
                            * da[0]
                            - (sign as f64) * wa[(i) - 1 + (1 - 1) * (ido - 1) + 1] * da[1];
                        ch[(i) + ido * ((k) + l1 * (1)) + 1] = wa[(i) - 1 + (1 - 1) * (ido - 1)]
                            * da[1]
                            + (sign as f64) * wa[(i) - 1 + (1 - 1) * (ido - 1) + 1] * da[0];
                    }
                    {
                        ch[(i) + ido * ((k) + l1 * (10))] = wa[(i) - 1 + (10 - 1) * (ido - 1)]
                            * db[0]
                            - (sign as f64) * wa[(i) - 1 + (10 - 1) * (ido - 1) + 1] * db[1];
                        ch[(i) + ido * ((k) + l1 * (10)) + 1] = wa[(i) - 1 + (10 - 1) * (ido - 1)]
                            * db[1]
                            + (sign as f64) * wa[(i) - 1 + (10 - 1) * (ido - 1) + 1] * db[0];
                    }
                }
                {
                    let mut da: [f64; 2] = [0.0; 2];
                    let mut db: [f64; 2] = [0.0; 2];
                    {
                        let mut ca: [f64; 2] = [0.0; 2];
                        let mut cb: [f64; 2] = [0.0; 2];
                        ca[0] = t1[0]
                            + tw2r * t2[0]
                            + tw4r * t3[0]
                            + tw5r * t4[0]
                            + tw3r * t5[0]
                            + tw1r * t6[0];
                        ca[1] = t1[1]
                            + tw2r * t2[1]
                            + tw4r * t3[1]
                            + tw5r * t4[1]
                            + tw3r * t5[1]
                            + tw1r * t6[1];
                        cb[1] = tw2i * t11[0] + tw4i * t10[0]
                            - tw5i * t9[0]
                            - tw3i * t8[0]
                            - tw1i * t7[0];
                        cb[0] = -(tw2i * t11[1] + tw4i * t10[1]
                            - tw5i * t9[1]
                            - tw3i * t8[1]
                            - tw1i * t7[1]);
                        {
                            da[0] = ca[0] + cb[0];
                            da[1] = ca[1] + cb[1];
                            db[0] = ca[0] - cb[0];
                            db[1] = ca[1] - cb[1];
                        }
                    }
                    {
                        ch[(i) + ido * ((k) + l1 * (2))] = wa[(i) - 1 + (2 - 1) * (ido - 1)]
                            * da[0]
                            - (sign as f64) * wa[(i) - 1 + (2 - 1) * (ido - 1) + 1] * da[1];
                        ch[(i) + ido * ((k) + l1 * (2)) + 1] = wa[(i) - 1 + (2 - 1) * (ido - 1)]
                            * da[1]
                            + (sign as f64) * wa[(i) - 1 + (2 - 1) * (ido - 1) + 1] * da[0];
                    }
                    {
                        ch[(i) + ido * ((k) + l1 * (9))] = wa[(i) - 1 + (9 - 1) * (ido - 1)]
                            * db[0]
                            - (sign as f64) * wa[(i) - 1 + (9 - 1) * (ido - 1) + 1] * db[1];
                        ch[(i) + ido * ((k) + l1 * (9)) + 1] = wa[(i) - 1 + (9 - 1) * (ido - 1)]
                            * db[1]
                            + (sign as f64) * wa[(i) - 1 + (9 - 1) * (ido - 1) + 1] * db[0];
                    }
                }
                {
                    let mut da: [f64; 2] = [0.0; 2];
                    let mut db: [f64; 2] = [0.0; 2];
                    {
                        let mut ca: [f64; 2] = [0.0; 2];
                        let mut cb: [f64; 2] = [0.0; 2];
                        ca[0] = t1[0]
                            + tw3r * t2[0]
                            + tw5r * t3[0]
                            + tw2r * t4[0]
                            + tw1r * t5[0]
                            + tw4r * t6[0];
                        ca[1] = t1[1]
                            + tw3r * t2[1]
                            + tw5r * t3[1]
                            + tw2r * t4[1]
                            + tw1r * t5[1]
                            + tw4r * t6[1];
                        cb[1] = tw3i * t11[0] - tw5i * t10[0] - tw2i * t9[0]
                            + tw1i * t8[0]
                            + tw4i * t7[0];
                        cb[0] = -(tw3i * t11[1] - tw5i * t10[1] - tw2i * t9[1]
                            + tw1i * t8[1]
                            + tw4i * t7[1]);
                        {
                            da[0] = ca[0] + cb[0];
                            da[1] = ca[1] + cb[1];
                            db[0] = ca[0] - cb[0];
                            db[1] = ca[1] - cb[1];
                        }
                    }
                    {
                        ch[(i) + ido * ((k) + l1 * (3))] = wa[(i) - 1 + (3 - 1) * (ido - 1)]
                            * da[0]
                            - (sign as f64) * wa[(i) - 1 + (3 - 1) * (ido - 1) + 1] * da[1];
                        ch[(i) + ido * ((k) + l1 * (3)) + 1] = wa[(i) - 1 + (3 - 1) * (ido - 1)]
                            * da[1]
                            + (sign as f64) * wa[(i) - 1 + (3 - 1) * (ido - 1) + 1] * da[0];
                    }
                    {
                        ch[(i) + ido * ((k) + l1 * (8))] = wa[(i) - 1 + (8 - 1) * (ido - 1)]
                            * db[0]
                            - (sign as f64) * wa[(i) - 1 + (8 - 1) * (ido - 1) + 1] * db[1];
                        ch[(i) + ido * ((k) + l1 * (8)) + 1] = wa[(i) - 1 + (8 - 1) * (ido - 1)]
                            * db[1]
                            + (sign as f64) * wa[(i) - 1 + (8 - 1) * (ido - 1) + 1] * db[0];
                    }
                }
                {
                    let mut da: [f64; 2] = [0.0; 2];
                    let mut db: [f64; 2] = [0.0; 2];
                    {
                        let mut ca: [f64; 2] = [0.0; 2];
                        let mut cb: [f64; 2] = [0.0; 2];
                        ca[0] = t1[0]
                            + tw4r * t2[0]
                            + tw3r * t3[0]
                            + tw1r * t4[0]
                            + tw5r * t5[0]
                            + tw2r * t6[0];
                        ca[1] = t1[1]
                            + tw4r * t2[1]
                            + tw3r * t3[1]
                            + tw1r * t4[1]
                            + tw5r * t5[1]
                            + tw2r * t6[1];
                        cb[1] = tw4i * t11[0] - tw3i * t10[0] + tw1i * t9[0] + tw5i * t8[0]
                            - tw2i * t7[0];
                        cb[0] = -(tw4i * t11[1] - tw3i * t10[1] + tw1i * t9[1] + tw5i * t8[1]
                            - tw2i * t7[1]);
                        {
                            da[0] = ca[0] + cb[0];
                            da[1] = ca[1] + cb[1];
                            db[0] = ca[0] - cb[0];
                            db[1] = ca[1] - cb[1];
                        }
                    }
                    {
                        ch[(i) + ido * ((k) + l1 * (4))] = wa[(i) - 1 + (4 - 1) * (ido - 1)]
                            * da[0]
                            - (sign as f64) * wa[(i) - 1 + (4 - 1) * (ido - 1) + 1] * da[1];
                        ch[(i) + ido * ((k) + l1 * (4)) + 1] = wa[(i) - 1 + (4 - 1) * (ido - 1)]
                            * da[1]
                            + (sign as f64) * wa[(i) - 1 + (4 - 1) * (ido - 1) + 1] * da[0];
                    }
                    {
                        ch[(i) + ido * ((k) + l1 * (7))] = wa[(i) - 1 + (7 - 1) * (ido - 1)]
                            * db[0]
                            - (sign as f64) * wa[(i) - 1 + (7 - 1) * (ido - 1) + 1] * db[1];
                        ch[(i) + ido * ((k) + l1 * (7)) + 1] = wa[(i) - 1 + (7 - 1) * (ido - 1)]
                            * db[1]
                            + (sign as f64) * wa[(i) - 1 + (7 - 1) * (ido - 1) + 1] * db[0];
                    }
                }
                {
                    let mut da: [f64; 2] = [0.0; 2];
                    let mut db: [f64; 2] = [0.0; 2];
                    {
                        let mut ca: [f64; 2] = [0.0; 2];
                        let mut cb: [f64; 2] = [0.0; 2];
                        ca[0] = t1[0]
                            + tw5r * t2[0]
                            + tw1r * t3[0]
                            + tw4r * t4[0]
                            + tw2r * t5[0]
                            + tw3r * t6[0];
                        ca[1] = t1[1]
                            + tw5r * t2[1]
                            + tw1r * t3[1]
                            + tw4r * t4[1]
                            + tw2r * t5[1]
                            + tw3r * t6[1];
                        cb[1] = tw5i * t11[0] - tw1i * t10[0] + tw4i * t9[0] - tw2i * t8[0]
                            + tw3i * t7[0];
                        cb[0] = -(tw5i * t11[1] - tw1i * t10[1] + tw4i * t9[1] - tw2i * t8[1]
                            + tw3i * t7[1]);
                        {
                            da[0] = ca[0] + cb[0];
                            da[1] = ca[1] + cb[1];
                            db[0] = ca[0] - cb[0];
                            db[1] = ca[1] - cb[1];
                        }
                    }
                    {
                        ch[(i) + ido * ((k) + l1 * (5))] = wa[(i) - 1 + (5 - 1) * (ido - 1)]
                            * da[0]
                            - (sign as f64) * wa[(i) - 1 + (5 - 1) * (ido - 1) + 1] * da[1];
                        ch[(i) + ido * ((k) + l1 * (5)) + 1] = wa[(i) - 1 + (5 - 1) * (ido - 1)]
                            * da[1]
                            + (sign as f64) * wa[(i) - 1 + (5 - 1) * (ido - 1) + 1] * da[0];
                    }
                    {
                        ch[(i) + ido * ((k) + l1 * (6))] = wa[(i) - 1 + (6 - 1) * (ido - 1)]
                            * db[0]
                            - (sign as f64) * wa[(i) - 1 + (6 - 1) * (ido - 1) + 1] * db[1];
                        ch[(i) + ido * ((k) + l1 * (6)) + 1] = wa[(i) - 1 + (6 - 1) * (ido - 1)]
                            * db[1]
                            + (sign as f64) * wa[(i) - 1 + (6 - 1) * (ido - 1) + 1] * db[0];
                    }
                }
            }
        }
    }
}

fn passg(
    ido: usize,
    ip: usize,
    l1: usize,
    cc: &mut [f64],
    ch: &mut [f64],
    wa: &[f64],
    csarr: &[f64],
    sign: i64,
) {
    let cdim: usize = ip;
    let ipph: usize = ((ip + 1) / 2) as usize;
    let idl1: usize = ido * l1;

    //cmplx * restrict wal=((cmplx *)malloc((ip)*sizeof(cmplx)));
    let mut wal: Vec<f64> = Vec::new();
    //if (!wal) return -1;
    wal.push(0.0);
    wal.push(0.0);
    for i in 1..ip {
        wal.push(csarr[i]);
        wal.push((sign as f64) * csarr[i + 1]);
    }

    for k in 0..l1 {
        for i in 0..ido {
            let tmp_r = cc[(i) + ido * ((0) + cdim * (k))];
            let tmp_i = cc[(i) + ido * ((0) + cdim * (k)) + 1];
            ch[(i) + ido * ((k) + l1 * (0))] = tmp_r;
            ch[(i) + ido * ((k) + l1 * (0)) + 1] = tmp_i;
        }
    }
    let mut j = 1;
    let mut jc = ip - 1;
    while j < ipph {
        for k in 0..l1 {
            for i in 0..ido {
                ch[(i) + ido * ((k) + l1 * (j))] =
                    &cc[(i) + ido * ((j) + cdim * (k))] + cc[(i) + ido * ((jc) + cdim * (k))];
                ch[(i) + ido * ((k) + l1 * (j)) + 1] = &cc[(i) + ido * ((j) + cdim * (k)) + 1]
                    + cc[(i) + ido * ((jc) + cdim * (k)) + 1];
                ch[(i) + ido * ((k) + l1 * (jc))] =
                    &cc[(i) + ido * ((j) + cdim * (k))] - cc[(i) + ido * ((jc) + cdim * (k))];
                ch[(i) + ido * ((k) + l1 * (jc)) + 1] = &cc[(i) + ido * ((j) + cdim * (k)) + 1]
                    - cc[(i) + ido * ((jc) + cdim * (k)) + 1];
            }
        }
        j += 1;
        jc -= 1;
    }
    for k in 0..l1 {
        for i in 0..ido {
            let mut tmp_r = ch[(i) + ido * ((k) + l1 * (0))];
            let mut tmp_i = ch[(i) + ido * ((k) + l1 * (0)) + 1];
            //let mut tmp = &mut ch[(i)+ido*((k)+l1*(0))];
            for j in 1..ipph {
                //tmp.r=tmp.r+ch[(i)+ido*((k)+l1*(j))];
                //tmp.i=tmp.i+ch[(i)+ido*((k)+l1*(j))+1];
                tmp_r = tmp_r + ch[(i) + ido * ((k) + l1 * (j))];
                tmp_i = tmp_i + ch[(i) + ido * ((k) + l1 * (j)) + 1];
            }
            cc[(i) + ido * ((k) + l1 * (0))] = tmp_r;
            cc[(i) + ido * ((k) + l1 * (0)) + 1] = tmp_i;
        }
    }
    let mut l = 1;
    let mut lc = ip - 1;
    while l < ipph {
        for ik in 0..idl1 {
            cc[(ik) + idl1 * (l)] = ch[(ik) + idl1 * (0)]
                + wal[l] * ch[(ik) + idl1 * (1)]
                + wal[2 * l] * ch[(ik) + idl1 * (2)];
            cc[(ik) + idl1 * (l) + 1] = ch[(ik) + idl1 * (0) + 1]
                + wal[l] * ch[(ik) + idl1 * (1) + 1]
                + wal[2 * l] * ch[(ik) + idl1 * (2) + 1];
            cc[(ik) + idl1 * (lc)] = -wal[l + 1] * ch[(ik) + idl1 * (ip - 1) + 1]
                - wal[2 * l + 1] * ch[(ik) + idl1 * (ip - 2) + 1];
            cc[(ik) + idl1 * (lc) + 1] = wal[l + 1] * ch[(ik) + idl1 * (ip - 1)]
                + wal[2 * l + 1] * ch[(ik) + idl1 * (ip - 2)];
        }

        let mut iwal: usize = 2 * l;
        let mut j: usize = 3;
        let mut jc: usize = ip - 3;
        while j < ipph - 1 {
            iwal += l;
            if iwal > ip {
                iwal -= ip;
            }
            let mut xwal: [f64; 2] = [0.0; 2];
            xwal[0] = wal[iwal];
            xwal[1] = wal[iwal + 1];
            iwal += l;
            if iwal > ip {
                iwal -= ip;
            }
            let mut xwal2: [f64; 2] = [0.0; 2];
            xwal2[0] = wal[iwal];
            xwal2[1] = wal[iwal + 1];
            for ik in 0..idl1 {
                cc[(ik) + idl1 * (l)] +=
                    ch[(ik) + idl1 * (j)] * xwal[0] + ch[(ik) + idl1 * (j + 1)] * xwal2[0];
                cc[(ik) + idl1 * (l) + 1] +=
                    ch[(ik) + idl1 * (j) + 1] * xwal[0] + ch[(ik) + idl1 * (j + 1) + 1] * xwal2[0];
                cc[(ik) + idl1 * (lc)] -= ch[(ik) + idl1 * (jc) + 1] * xwal[1]
                    + ch[(ik) + idl1 * (jc - 1) + 1] * xwal2[1];
                cc[(ik) + idl1 * (lc) + 1] +=
                    ch[(ik) + idl1 * (jc)] * xwal[1] + ch[(ik) + idl1 * (jc - 1)] * xwal2[1];
            }
            j += 2;
            jc -= 2;
        }
        while j < ipph {
            iwal += l;
            if iwal > ip {
                iwal -= ip;
            }
            let mut xwal: [f64; 2] = [0.0; 2];
            xwal[0] = wal[iwal];
            xwal[1] = wal[iwal + 1];
            for ik in 0..idl1 {
                cc[(ik) + idl1 * (l)] += ch[(ik) + idl1 * (j)] * xwal[0];
                cc[(ik) + idl1 * (l) + 1] += ch[(ik) + idl1 * (j) + 1] * xwal[0];
                cc[(ik) + idl1 * (lc)] -= ch[(ik) + idl1 * (jc) + 1] * xwal[1];
                cc[(ik) + idl1 * (lc) + 1] += ch[(ik) + idl1 * (jc)] * xwal[1];
            }
            j += 1;
            jc -= 1;
        }
        l += 1;
        lc -= 1;
    }

    if ido == 1 {
        let mut j = 1;
        let mut jc = ip - 1;
        while j < ipph {
            for ik in 0..idl1 {
                let t1_r = cc[(ik) + idl1 * (j)];
                let t1_i = cc[(ik) + idl1 * (j) + 1];
                let t2_r = cc[(ik) + idl1 * (jc)];
                let t2_i = cc[(ik) + idl1 * (jc) + 1];
                {
                    cc[(ik) + idl1 * (j)] = t1_r + t2_r;
                    cc[(ik) + idl1 * (j) + 1] = t1_i + t2_i;
                    cc[(ik) + idl1 * (jc)] = t1_r - t2_r;
                    cc[(ik) + idl1 * (jc) + 1] = t1_i - t2_i;
                }
            }
            j += 1;
            jc -= 1;
        }
    } else {
        let mut j = 1;
        let mut jc = ip - 1;
        while j < ipph {
            for k in 0..l1 {
                //let t1=&cc[(0)+ido*((k)+l1*(j))];
                let t1_r = cc[(0) + ido * ((k) + l1 * (j))];
                let t1_i = cc[(0) + ido * ((k) + l1 * (j)) + 1];
                //let  t2=&cc[(0)+ido*((k)+l1*(jc))];
                let t2_r = cc[(0) + ido * ((k) + l1 * (jc))];
                let t2_i = cc[(0) + ido * ((k) + l1 * (jc)) + 1];
                {
                    /*cc[(0)+ido*((k)+l1*(j))]=t1.r+t2.r;
                    cc[(0)+ido*((k)+l1*(j))+1]=t1.i+t2.i;
                    cc[(0)+ido*((k)+l1*(jc))]=t1.r-t2.r;
                    cc[(0)+ido*((k)+l1*(jc))+1]=t1.i-t2.i;*/
                    cc[(0) + ido * ((k) + l1 * (j))] = t1_r + t2_r;
                    cc[(0) + ido * ((k) + l1 * (j)) + 1] = t1_i + t2_i;
                    cc[(0) + ido * ((k) + l1 * (jc))] = t1_r - t2_r;
                    cc[(0) + ido * ((k) + l1 * (jc)) + 1] = t1_i - t2_i;
                }
                for i in 1..ido {
                    let mut x1: [f64; 2] = [0.0; 2];
                    let mut x2: [f64; 2] = [0.0; 2];
                    {
                        x1[0] =
                            cc[(i) + ido * ((k) + l1 * (j))] + cc[(i) + ido * ((k) + l1 * (jc))];
                        x1[1] = cc[(i) + ido * ((k) + l1 * (j)) + 1]
                            + cc[(i) + ido * ((k) + l1 * (jc)) + 1];
                        x2[0] =
                            cc[(i) + ido * ((k) + l1 * (j))] - cc[(i) + ido * ((k) + l1 * (jc))];
                        x2[1] = cc[(i) + ido * ((k) + l1 * (j)) + 1]
                            - cc[(i) + ido * ((k) + l1 * (jc)) + 1];
                    }
                    let mut idij: usize = (j - 1) * (ido - 1) + i - 1;
                    {
                        cc[(i) + ido * ((k) + l1 * (j))] =
                            &wa[idij] * x1[0] - (sign as f64) * wa[idij + 1] * x1[1];
                        cc[(i) + ido * ((k) + l1 * (j)) + 1] =
                            &wa[idij] * x1[1] + (sign as f64) * wa[idij + 1] * x1[0];
                    }
                    idij = (jc - 1) * (ido - 1) + i - 1;
                    {
                        cc[(i) + ido * ((k) + l1 * (jc))] =
                            &wa[idij] * x2[0] - (sign as f64) * wa[idij + 1] * x2[1];
                        cc[(i) + ido * ((k) + l1 * (jc)) + 1] =
                            &wa[idij] * x2[1] + (sign as f64) * wa[idij + 1] * x2[0];
                    }
                }
            }
            j += 1;
            jc -= 1;
        }
    }
}

fn pass_all(plan: *mut cfftp_plan_i, c: &mut [f64], fct: f64, sign: i64) -> i32 {
    unsafe {
        //if (plan->length==1) return 0;
        let len: usize = (*plan).length;
        let mut l1: usize = 1;
        let nf: usize = (*plan).nfct;
        let mut ch: Vec<f64> = Vec::new();
        //if (!ch) return -1;
        //let mut p1=c;
        //let mut p2=ch;
        let mut ch_slice = ch.as_mut_slice();

        for k1 in 0..nf {
            let ip = (*plan).fct[k1].fct;
            let l2 = ip * l1;
            let ido = len / l2;
            if ip == 4 {
                if sign > 0 {
                    pass4b(ido, l1, c, ch_slice, (*plan).fct[k1].tw.as_slice())
                } else {
                    pass4f(ido, l1, &c, ch_slice, (*plan).fct[k1].tw.as_slice());
                }
            } else if ip == 2 {
                if sign > 0 {
                    pass2b(ido, l1, &c, ch_slice, (*plan).fct[k1].tw.as_slice())
                } else {
                    pass2f(ido, l1, &c, ch_slice, (*plan).fct[k1].tw.as_slice());
                }
            } else if ip == 3 {
                if sign > 0 {
                    pass3b(ido, l1, &c, ch_slice, (*plan).fct[k1].tw.as_slice())
                } else {
                    pass3f(ido, l1, &c, ch_slice, (*plan).fct[k1].tw.as_slice());
                }
            } else if ip == 5 {
                if sign > 0 {
                    pass5b(ido, l1, &c, ch_slice, (*plan).fct[k1].tw.as_slice())
                } else {
                    pass5f(ido, l1, &c, ch_slice, (*plan).fct[k1].tw.as_slice());
                }
            } else if ip == 7 {
                pass7(ido, l1, &c, ch_slice, (*plan).fct[k1].tw.as_slice(), sign);
            } else if ip == 11 {
                pass11(ido, l1, &c, ch_slice, (*plan).fct[k1].tw.as_slice(), sign);
            } else {
                passg(
                    ido,
                    ip,
                    l1,
                    c,
                    ch_slice,
                    (*plan).fct[k1].tw.as_slice(),
                    (*plan).fct[k1].tws.as_slice(),
                    sign,
                );
                //if ()        { return -1; }
                /*do {
                    cmplx * tmp_=(p1);
                    (p1)=(p2);
                    (p2)=tmp_;
                } while(0);
                }*/
                /*do {
                cmplx * tmp_=(p1);
                (p1)=(p2);
                (p2)=tmp_;
                } while(0);*/
                l1 = l2;
            }
            /*if (p1!=c) {
                if (fct!=1.){
                for i in 0..len {
                    c[i] = ch[i]*fct;
                    c[i+1] = ch[i+1]*fct;
                    }
                }
                else {
                memcpy (c,p1,len*sizeof(cmplx));
                }
                }
            else{
                if (fct!=1.) {
                for i in 0..len {
                    c[i] *= fct;
                    c[i+1] *= fct;
                    }
            }
                }*/
            //
        }
    }
    return 0;
}

fn cfftp_factorize(plan: &mut cfftp_plan_i) -> i32 {
    let mut length = plan.length;
    let mut nfct: usize = 0;
    while (length % 4) == 0 {
        if nfct >= NFCT {
            return -1;
        }
        plan.fct[nfct].fct = 4;
        nfct += 1;
        length >>= 2;
    }
    if (length % 2) == 0 {
        length >>= 1;

        if nfct >= NFCT {
            return -1;
        }
        plan.fct[nfct].fct = 2;
        nfct += 1;
        let tmp_ = plan.fct[0].fct;
        (plan.fct[0].fct) = plan.fct[nfct - 1].fct;
        (plan.fct[nfct - 1].fct) = tmp_;
    }
    let mut maxl: usize = ((length as f64).sqrt() as usize) + 1;
    let mut divisor: usize = 3;
    while (length > 1) && (divisor < maxl) {
        if (length % divisor) == 0 {
            while (length % divisor) == 0 {
                if nfct >= NFCT {
                    return -1;
                }
                plan.fct[nfct].fct = divisor;
                nfct += 1;
                length /= divisor;
            }
            maxl = ((length as f64).sqrt() as usize) + 1;
        }
        divisor += 2;
    }
    if length > 1 {
        plan.fct[nfct].fct = length;
    }
    nfct += 1;
    plan.nfct = nfct;
    return 0;
}

fn cfftp_twsize(plan: *mut cfftp_plan_i) -> usize {
    let mut twsize: usize = 0;
    let mut l1: usize = 1;
    unsafe {
        for k in 0..(*plan).nfct {
            let ip = (*plan).fct[k].fct;
            let ido = (*plan).length / (l1 * ip);
            twsize += (ip - 1) * (ido - 1);
            if ip > 11 {
                twsize += ip;
            }
            l1 *= ip;
        }
    }
    return twsize;
}

fn cfftp_comp_twiddle(plan: &mut cfftp_plan_i) -> i32 {
    let length: usize = unsafe { (*plan).length };
    let mut twid: Vec<f64> = Vec::with_capacity(2 * length);
    //if (!twid) {return -1;}
    sincos_2pibyn(length, twid.as_mut_slice());
    let mut l1: usize = 1;
    let mut memofs: usize = 0;

    for k in 0..plan.nfct {
        let ip = plan.fct[k].fct;
        let ido = length / (l1 * ip);
        //todo: plan.fct[k].tw = plan.mem + memofs;
        memofs += (ip - 1) * (ido - 1);
        for j in 1..ip {
            for i in 1..ido {
                plan.fct[k].tw[(j - 1) * (ido - 1) + i - 1] = twid[2 * j * l1 * i];
                plan.fct[k].tw[(j - 1) * (ido - 1) + i - 1 + 1] = twid[2 * j * l1 * i + 1];
            }
        }
        if ip > 11 {
            //todo: plan.fct[k].tws=plan.mem + memofs;
            memofs += ip;
            for j in 0..ip {
                plan.fct[k].tws[j] = twid[2 * j * l1 * ido];
                plan.fct[k].tws[j + 1] = twid[2 * j * l1 * ido + 1];
            }
        }
        l1 *= ip;
    }
    return 0;
}

fn make_cfftp_plan(len: usize) -> cfftp_plan {
    //if (len==0) return ((void *)0);
    //let plan: cfftp_plan = malloc(size_of::<cfftp_plan_i>()) as cfftp_plan;
    //if (!plan) return ((void *)0);

    let mut tmp_fct: Vec<cfftp_fctdata> = Vec::new();
    for i in 0..NFCT {
        tmp_fct.insert(i, cfftp_fctdata {
            fct: 0,
            tw: Vec::new(),
            tws: Vec::new(),
        });
    } //(cfftp_fctdata){0,0,0};
    //plan.mem=0;
    let tmp_cfftp_plan_i = cfftp_plan_i {
        length: len,
        nfct: 0,
        mem: Vec::new(),
        fct: tmp_fct,
    };
    let plan: cfftp_plan = Box::into_raw(Box::new(tmp_cfftp_plan_i));
    if len == 1 {
        return plan;
    }
    unsafe {
        if cfftp_factorize(&mut *plan) != 0 {
            Box::from_raw(plan);
            return null_mut();
        }
    }
    //let tws = cfftp_twsize(plan);
    //(*plan).mem = (malloc((tws) * size_of::<f64>())) as *mut f64;
    //if (!plan.mem) { do { free(plan); (plan)= ((void *)0); } while(0); return ((void *)0); }
    unsafe {
        if cfftp_comp_twiddle(&mut *plan) != 0 {
            //free((*plan).mem as *mut c_void);

            Box::from_raw(plan);

            return null_mut();
        }
    }
    return plan;
}

struct fftblue_plan_i {
    n: usize,
    n2: usize,
    plan: *mut cfftp_plan_i,
    mem: Vec<f64>,
    bk: Vec<f64>,
    bkf: Vec<f64>,
}

type fftblue_plan = *mut fftblue_plan_i;

fn make_fftblue_plan(length: usize) -> fftblue_plan {
    //let mut plan: fftblue_plan = malloc((1) * size_of::<fftblue_plan_i>()) as fftblue_plan;
    let tmp_n2 = good_size(length * 2 - 1);
    let tmp_mem: Vec<f64> = Vec::with_capacity(2*length+2*tmp_n2);
    let mut tmp_bk: Vec<f64> = Vec::with_capacity(2*length+2*tmp_n2);
    let mut tmp_bkf: Vec<f64> = Vec::with_capacity(2*length);
    //if (!plan) return ((void *)0) ;
    //(*plan).n = length;
    //(*plan).n2 = good_size((*plan).n * 2 - 1);
    //fixme: (*plan).mem = malloc((2*(*plan).n+2*(*plan).n2)*size_of::<f64>()) as *mut f64;
    //(*plan).mem = Vec::new();

    /*if (!plan->mem) {
      do {
          free(plan);
          (plan)= ((void *)0);
         } while(0);
          return ((void *)0);
    }*/
    //fixme: (*plan).bk = (*plan).mem.to_vec();
    //fixme: (*plan).bkf = (*plan).bk + 2 * (*plan).n;
    //fixme: let mut tmp: *mut f64 = (malloc((4 * (*plan).n) * size_of::<f64>())) as *mut f64;
    let mut tmp: Vec<f64> = Vec::with_capacity(4 * length);

    /*if (!tmp) { do { free(plan->mem); (plan->mem)=
    ((void *)0)
    ; } while(0); do { free(plan); (plan)=
                        ((void *)0)
                        ; } while(0); return
                                              ((void *)0)
                                                  ; }*/
    sincos_2pibyn(2 * length, tmp.as_mut_slice());
    
    tmp_bk.insert(0, 1.0);
    tmp_bk.insert(1, 0.0);


    let mut coeff: usize = 0;
    for m in 1..length {
        coeff += 2 * m - 1;
        if coeff >= 2 * length {
            coeff -= 2 * length;
        }
        tmp_bk.insert(2 * m, tmp[2 * coeff]);
        tmp_bk.insert(2 * m + 1, tmp[2 * coeff + 1]);
    }

    let xn2: f64 = 1.0 / tmp_n2 as f64;
    tmp_bkf.insert(0, tmp_bk[0] * xn2);
    tmp_bkf.insert(1, tmp_bk[1] * xn2);

    let mut m = 2;
    while m < 2 * length {
        tmp_bkf.insert(2 * tmp_n2 - m, tmp_bk[m] * xn2);
        tmp_bkf.insert(m, tmp_bkf[2 * tmp_n2 - m]);
        tmp_bkf.insert(2 * tmp_n2 - m + 1, tmp_bk[m + 1] * xn2);
        tmp_bkf.insert(m + 1, tmp_bkf[2 * tmp_n2 - m + 1]);
        m += 2;
    }

    m = 2 * length;
    while m <= (2 * tmp_n2 - 2 * length + 1) {
        tmp_bkf.push(0.0);
        m += 1;
    }
    /*if (cfftp_forward(plan->plan,plan->bkf,1.)!=0)
    { do { free(tmp); (tmp)=

     ((void *)0)

     ; } while(0); do { free(plan->mem); (plan->mem)=

                   ((void *)0)
                   ; } while(0); do { free(plan); (plan)=
                                       ((void *)0)
                                       ; } while(0); return
                                                             ((void *)0)
                                                                 ; }*/
    //do { free(tmp); (tmp)=((void *)0) ; } while(0);

    let tmp_fftblue_plan_i = fftblue_plan_i{
        n: length, 
        n2: tmp_n2, 
        plan: make_cfftp_plan(tmp_n2),
         mem: tmp_mem,
         bk: tmp_bk,
         bkf: tmp_bkf};
    let plan: fftblue_plan = Box::into_raw(Box::new(tmp_fftblue_plan_i));
    return plan;
}

fn destroy_cfftp_plan(plan: cfftp_plan) {
    unsafe {
        Box::from_raw(plan);
    }
}

fn destroy_fftblue_plan(fftblueplan: fftblue_plan) {
    unsafe {
        destroy_cfftp_plan((*fftblueplan).plan);
        Box::from_raw(fftblueplan);
    }
}

pub struct cfft_plan_i {
    packplan: cfftp_plan,
    blueplan: fftblue_plan,
}
type cfft_plan = *mut cfft_plan_i;

#[no_mangle]
pub unsafe extern "C" fn make_cfft_plan(length: usize) -> cfft_plan {
    if length == 0 {
        return null_mut() as cfft_plan;
    //return 0 as *mut c_void as *mut _ as cfft_plan;
    } else {
        //let mut plan: cfft_plan = malloc(size_of::<cfft_plan_i>()) as cfft_plan;
        let mut plan: cfft_plan = Box::into_raw(Box::new(cfft_plan_i{packplan: null_mut(), blueplan: null_mut()}));
        if ((length < 50) || (largest_prime_factor(length) <= ((length as f64).sqrt() as usize))) {
            (*plan).packplan = make_cfftp_plan(length);
            //forget(plan);
            return plan;
        }
        let comp1: f64 = cost_guess(length);
        let mut comp2: f64 = 2.0 * cost_guess(good_size(2 * length - 1));
        comp2 *= 1.5;
        if comp2 < comp1 {
            (*plan).blueplan = make_fftblue_plan(length);
        } else {
            (*plan).packplan = make_cfftp_plan(length);
        }
        //forget(plan);
        return plan;
    }
}

#[no_mangle]
pub unsafe extern "C" fn destroy_cfft_plan(plan: cfft_plan) {
    destroy_fftblue_plan((*plan).blueplan);
    destroy_cfftp_plan((*plan).packplan);
    Box::from_raw(plan as *mut c_void);
}

unsafe fn fftblue_fft(plan: fftblue_plan, c: &mut [f64], isign: i32, fct: f64) -> i32 {
    let n: usize = (*plan).n;
    let n2: usize = (*plan).n2;
    let bk = &(*plan).bk;
    let bkf = &(*plan).bkf;
    //let akf = ((double *)malloc((2*n2)*sizeof(double)));
    let mut akf: Vec<f64> = Vec::with_capacity(2 * n2);
    //if (!akf) return -1;

    if isign > 0 {
        let mut m = 0;
        while m < 2 * n {
            akf[m] = c[m] * bk[m] - c[m + 1] * bk[m + 1];
            akf[m + 1] = c[m] * bk[m + 1] + c[m + 1] * bk[m];
            m += 2;
        }
    } else {
        let mut m = 0;
        while m < 2 * n {
            akf[m] = c[m] * bk[m] + c[m + 1] * bk[m + 1];
            akf[m + 1] = -c[m] * bk[m + 1] + c[m + 1] * bk[m];
            m += 2;
        }
    }
    let mut m: usize = 2 * n;
    while m < 2 * n2 {
        akf[m] = 0.0;
        m += 1;
    }

    /*if (cfftp_forward (plan->plan,akf,fct)!=0){
       do { free(akf); (akf)=((void *)0); } while(0); return -1;
    }*/

    if isign > 0 {
        let mut m: usize = 0;
        while m < 2 * n2 {
            let im = -akf[m] * bkf[m + 1] + akf[m + 1] * bkf[m];
            akf[m] = akf[m] * bkf[m] + akf[m + 1] * bkf[m + 1];
            akf[m + 1] = im;
        }
        m += 2;
    } else {
        let mut m: usize = 0;
        while m < 2 * n2 {
            let im = akf[m] * bkf[m + 1] + akf[m + 1] * bkf[m];
            akf[m] = akf[m] * bkf[m] - akf[m + 1] * bkf[m + 1];
            akf[m + 1] = im;
            m += 2;
        }
    }

    /*if (cfftp_backward (plan->plan,akf,1.0)!=0){
       do { free(akf); (akf)=((void *)0); } while(0); return -1;
    }*/

    if isign > 0 {
        let mut m: usize = 0;
        while m < 2 * n {
            c[m] = bk[m] * akf[m] - bk[m + 1] * akf[m + 1];
            c[m + 1] = bk[m + 1] * akf[m] + bk[m] * akf[m + 1];
            m += 2;
        }
    } else {
        let mut m: usize = 0;
        while m < 2 * n {
            c[m] = bk[m] * akf[m] + bk[m + 1] * akf[m + 1];
            c[m + 1] = -bk[m + 1] * akf[m] + bk[m] * akf[m + 1];
            m += 2;
        }
    }
    //do { free(akf); (akf)= ((void *)0) ; } while(0);
    return 0;
}

unsafe fn cfftblue_backward(plan: fftblue_plan, c: &mut [f64], fct: f64) -> i32 {
    return fftblue_fft(plan, c, 1, fct);
}

fn cfftp_backward(plan: cfftp_plan, c: &mut [f64], fct: f64) -> i32 {
    unsafe {
        return pass_all(plan, c, fct, 1);
    };
}

#[no_mangle]
pub unsafe extern "C" fn cfft_backward(plan: cfft_plan, c: *mut f64, fct: f64) -> i32 {
    use std::slice::from_raw_parts_mut;
    if (*plan).packplan.is_null() == false {
        let tmp_packplan = (*plan).packplan;
        let len = (*tmp_packplan).length;
        let mut tmp_c = from_raw_parts_mut(c, len);
        return cfftp_backward((*plan).packplan, tmp_c, fct);
    }
    let tmp_blueplan = (*plan).blueplan;
    let mut len = (*tmp_blueplan).n * 2;
    let mut tmp_c = from_raw_parts_mut(c, len);
    return cfftblue_backward((*plan).blueplan, tmp_c, fct);
}

#[no_mangle]
pub unsafe extern "C" fn cfft_forward(plan: cfft_plan, c: *mut f64, fct: f64) -> i32 {
    if (*plan).packplan.is_null() == false {
        let tmp_packplan = (*plan).packplan;
        let len = (*tmp_packplan).length;
        let mut tmp_c = std::slice::from_raw_parts_mut(c, len);
        return cfftp_forward((*plan).packplan, tmp_c, fct);
    }
    // if (plan->blueplan)
    let tmp_blueplan = (*plan).blueplan;
    let mut len = (*tmp_blueplan).n * 2;
    let mut tmp_c = std::slice::from_raw_parts_mut(c, len);
    return cfftblue_forward((*plan).blueplan, tmp_c, fct);
}

fn cfftp_forward(plan: cfftp_plan, c: &mut [f64], fct: f64) -> i32 {
    return pass_all(plan, c, fct, -1);
}

fn cfftblue_forward(plan: fftblue_plan, c: &mut [f64], fct: f64) -> i32 {
    unsafe {
        return fftblue_fft(plan, c, -1, fct);
    }
}

#[no_mangle]
pub unsafe extern "C" fn cfft_length(plan: cfft_plan) -> usize {
    if (*plan).packplan.is_null() == false {
        let tmp_packplan = (*plan).packplan;
        return (*tmp_packplan).length;
    }
    let tmp_blueplan = (*plan).blueplan;
    return (*tmp_blueplan).n;
}

#[repr(C)]
struct rfftp_fctdata {
    fct: usize,
    tw: Vec<f64>,
    tws: Vec<f64>,
}

#[repr(C)]
struct rfftp_plan_i {
    length: usize,
    nfct: usize,
    mem: Vec<f64>,
    fct: [rfftp_fctdata; NFCT],
}
type rfftp_plan = *mut rfftp_plan_i;

pub struct rfft_plan_i {
    packplan: rfftp_plan,
    blueplan: fftblue_plan,
}
type rfft_plan = *mut rfft_plan_i;

#[no_mangle]
pub unsafe extern "C" fn make_rfft_plan(length: usize) -> rfft_plan {
    if length == 0 {
        return null_mut();
    }
    let mut plan: rfft_plan = (malloc((1) * size_of::<rfft_plan_i>())) as rfft_plan;
    if plan.is_null() {
        return null_mut();
    }
    (*plan).blueplan = null_mut();
    (*plan).packplan = null_mut() as rfftp_plan;
    let length_sqrt = (length as f64).sqrt() as usize;
    if (length < 50) || (largest_prime_factor(length) <= length_sqrt) {
        (*plan).packplan = make_rfftp_plan(length);
        if (*plan).packplan.is_null() {
            free(plan as *mut c_void);
            plan = null_mut();
            return null_mut();
        }
        return plan;
    }
    let comp1: f64 = 0.5 * cost_guess(length);
    let mut comp2: f64 = 2.0 * cost_guess(good_size(2 * length - 1));
    comp2 *= 1.5;
    if comp2 < comp1 {
        (*plan).blueplan = make_fftblue_plan(length);
        let tmp_blueplan = (*plan).blueplan;
        if tmp_blueplan.is_null() {
            free(plan as *mut c_void);
            plan = null_mut();
            return null_mut();
        }
    } else {
        (*plan).packplan = make_rfftp_plan(length);
        let tmp_packplan = (*plan).packplan;
        if tmp_packplan.is_null() {
            free(plan as *mut c_void);
            plan = null_mut();
            return null_mut();
        }
    }
    return plan;
}

unsafe fn make_rfftp_plan(length: usize) -> rfftp_plan {
    if length == 0 {
        return null_mut();
    };
    let mut plan: rfftp_plan = (malloc((1) * size_of::<rfftp_plan_i>())) as rfftp_plan;
    if plan.is_null() {
        return null_mut();
    }
    (*plan).length = length;
    (*plan).nfct = 0;
    for i in 0..NFCT {
        (*plan).fct[i] = rfftp_fctdata {
            fct: 0,
            tw: Vec::new(),
            tws: Vec::new(),
        };
    }
    if length == 1 {
        return plan;
    }
    if rfftp_factorize(plan) != 0 {
        free(plan as *mut c_void);
        plan = null_mut();
        return null_mut();
    }
    let tws: usize = rfftp_twsize(plan);
    //(*plan).mem=(malloc((tws)*size_of::<f64>())) as *mut f64;
    (*plan).mem = Vec::with_capacity(tws);
    if rfftp_comp_twiddle(plan) != 0 {
        free(plan as *mut c_void);
        plan = null_mut();
        return null_mut();
    }
    return plan;
}

unsafe fn rfftp_factorize(plan: rfftp_plan) -> i32 {
    let mut length = (*plan).length;
    let mut nfct: usize = 0;
    while (length % 4) == 0 {
        if nfct >= NFCT {
            return -1;
        }
        (*plan).fct[nfct].fct = 4;
        nfct += 1;
        length >>= 2;
    }
    if (length % 2) == 0 {
        length >>= 1;

        if nfct >= NFCT {
            return -1;
        }
        (*plan).fct[nfct].fct = 2;
        nfct += 1;

        let tmp_: usize = (*plan).fct[0].fct;
        (*plan).fct[0].fct = (*plan).fct[nfct - 1].fct;
        (*plan).fct[nfct - 1].fct = tmp_;
    }
    let mut maxl: usize = ((length as f64).sqrt() as usize) + 1;
    let mut divisor: usize = 3;
    while (length > 1) && (divisor < maxl) {
        if (length % divisor) == 0 {
            while (length % divisor) == 0 {
                if nfct >= NFCT {
                    return -1;
                }
                (*plan).fct[nfct].fct = divisor;
                nfct += 1;
                length /= divisor;
            }
            maxl = ((length as f64).sqrt() as usize) + 1;
        }
        divisor += 2;
    }

    if length > 1 {
        (*plan).fct[nfct].fct = length;
        nfct += 1;
    }
    (*plan).nfct = nfct;
    return 0;
}

unsafe fn rfftp_twsize(plan: rfftp_plan) -> usize {
    let mut twsize: usize = 0;
    let mut l1: usize = 1;
    for k in 0..(*plan).nfct {
        let ip = (*plan).fct[k].fct;
        let ido = (*plan).length / (l1 * ip);
        twsize += (ip - 1) * (ido - 1);
        if ip > 5 {
            twsize += 2 * ip;
        }
        l1 *= ip;
    }
    return twsize;
}

unsafe fn rfftp_comp_twiddle(plan: rfftp_plan) -> i32 {
    let length: usize = (*plan).length;
    //let twid: *mut f64 = (malloc((2*length)*size_of::<f64>())) as (*mut f64);
    let mut twid: Vec<f64> = Vec::with_capacity(2 * length);
    sincos_2pibyn_half(length, twid.as_mut_slice());
    let mut l1: usize = 1;
    //fixme:  let mut ptr = (*plan).mem;
    for k in 0..(*plan).nfct {
        let ip: usize = (*plan).fct[k].fct;
        let ido = length / (l1 * ip);
        if k < (*plan).nfct - 1 {
            //fixme:  (*plan).fct[k].tw = ptr;
            //fixme: ptr+=(ip-1)*(ido-1);
            for j in 1..ip {
                for i in 1..((ido - 1) / 2 + 1) {
                    (*plan).fct[k].tw[(j - 1) * (ido - 1) + 2 * i - 2] = twid[2 * j * l1 * i];
                    (*plan).fct[k].tw[(j - 1) * (ido - 1) + 2 * i - 1] = twid[2 * j * l1 * i + 1];
                }
            }
            if ip > 5 {
                //fixme:  (*plan).fct[k].tws = ptr;
                //fixme:  ptr+=2*ip;
                (*plan).fct[k].tws[0] = 1.0;
                (*plan).fct[k].tws[1] = 0.0;
                for i in 1..(ip >> 1) + 1 {
                    (*plan).fct[k].tws[2 * i] = twid[2 * i * (length / ip)];
                    (*plan).fct[k].tws[2 * i + 1] = twid[2 * i * (length / ip) + 1];
                    (*plan).fct[k].tws[2 * (ip - i)] = twid[2 * i * (length / ip)];
                    (*plan).fct[k].tws[2 * (ip - i) + 1] = -twid[2 * i * (length / ip) + 1];
                }
            }
            l1 *= ip;
        }
    }
    return 0;
}

#[no_mangle]
pub unsafe extern "C" fn destroy_rfft_plan(plan: rfft_plan) {
    if (*plan).blueplan.is_null() == false {
        destroy_fftblue_plan((*plan).blueplan);
    }
    if (*plan).packplan.is_null() == false {
        destroy_rfftp_plan((*plan).packplan);
    }
    free(plan as *mut c_void);
}

unsafe fn destroy_rfftp_plan(plan: rfftp_plan) {
    free(plan as *mut c_void);
}

#[no_mangle]
pub unsafe extern "C" fn rfft_backward(plan: rfft_plan, c: *mut f64, fct: f64) -> i32 {
    if (*plan).packplan.is_null() == false {
        let tmp_packplan = (*plan).packplan;
        let len = (*tmp_packplan).length;
        let tmp_c = from_raw_parts_mut(c, len);
        return rfftp_backward((*plan).packplan, tmp_c, fct);
    } else {
        let tmp_blueplan = (*plan).blueplan;
        let len = (*tmp_blueplan).n;// fixme: is that n or n2, or other field?
        let tmp_c = from_raw_parts_mut(c, len);
        return rfftblue_backward(tmp_blueplan, tmp_c, fct);
    }
}

unsafe fn rfftp_backward(plan: rfftp_plan, c: &mut [f64], fct: f64) -> i32 {
    if (*plan).length == 1 {
        return 0;
    }
    let n: usize = (*plan).length;
    let mut l1: usize = 1;
    let nf = (*plan).nfct;
    let ch: *mut f64 = (malloc((n) * size_of::<f64>())) as *mut f64;
    if ch.is_null() {
        return -1;
    }
    let mut p1: *mut f64 = c.as_mut_ptr();
    let mut p2: *mut f64 = ch;

    let mut k: usize = 0;
    while k < nf {
        let ip: usize = (*plan).fct[k].fct;
        let ido: usize = n / (ip * l1);
        if ip == 4 {
            radb4(ido, l1, p1, p2, (*plan).fct[k].tw.as_ptr());
        } else if ip == 2 {
            radb2(ido, l1, p1, p2, (*plan).fct[k].tw.as_ptr());
        } else if ip == 3 {
            radb3(ido, l1, p1, p2, (*plan).fct[k].tw.as_ptr());
        } else if ip == 5 {
            radb5(ido, l1, p1, p2, (*plan).fct[k].tw.as_ptr());
        } else {
            radbg(
                ido,
                ip,
                l1,
                p1,
                p2,
                (*plan).fct[k].tw.as_ptr(),
                (*plan).fct[k].tws.as_ptr(),
            );
        }
        let mut tmp_: *mut f64 = p1;
        p1 = p2;
        p2 = tmp_;
        l1 *= ip;
        k += 1;
    }
    copy_and_norm(c.as_mut_ptr(), p1, n, fct);
    free(ch as *mut c_void);

    return 0;
}

unsafe fn radb4(ido: usize, l1: usize, cc: *const f64, ch: *mut f64, wa: *const f64) {
    let cdim: usize = 4;
    let sqrt2 = 1.41421356237309504880;

    let tmp_cc = from_raw_parts(cc, l1); //fixme: check if l1 is actual length of cc
    let tmp_ch = from_raw_parts_mut(ch, l1); //fixme: check if l1 is actual length of ch
    let tmp_wa = from_raw_parts(wa, l1); //fixme: check if l1 is actual length of wa
    let mut k: usize = 0;
    while k < l1 {
        let mut tr1: f64 = 0.0;
        let mut tr2: f64 = 0.0;
        {
            tr2 = tmp_cc[(0) + ido * ((0) + cdim * (k))]
                + tmp_cc[(ido - 1) + ido * ((3) + cdim * (k))];
            tr1 = tmp_cc[(0) + ido * ((0) + cdim * (k))]
                - tmp_cc[(ido - 1) + ido * ((3) + cdim * (k))];
        }
        let mut tr3: f64 = 2.0 * tmp_cc[(ido - 1) + ido * ((1) + cdim * (k))];
        let mut tr4: f64 = 2.0 * tmp_cc[(0) + ido * ((2) + cdim * (k))];
        {
            tmp_ch[(0) + ido * ((k) + l1 * (0))] = tr2 + tr3;
            tmp_ch[(0) + ido * ((k) + l1 * (2))] = tr2 - tr3;
        }
        {
            tmp_ch[(0) + ido * ((k) + l1 * (3))] = tr1 + tr4;
            tmp_ch[(0) + ido * ((k) + l1 * (1))] = tr1 - tr4;
        }
        k += 1;
    }
    if (ido & 1) == 0 {
        let mut k: usize = 0;
        while k < l1 {
            let mut ti1: f64 = 0.0;
            let mut ti2: f64 = 0.0;
            let mut tr1: f64 = 0.0;
            let mut tr2: f64 = 0.0;
            {
                ti1 =
                    tmp_cc[(0) + ido * ((3) + cdim * (k))] + tmp_cc[(0) + ido * ((1) + cdim * (k))];
                ti2 =
                    tmp_cc[(0) + ido * ((3) + cdim * (k))] - tmp_cc[(0) + ido * ((1) + cdim * (k))];
            }
            {
                tr2 = tmp_cc[(ido - 1) + ido * ((0) + cdim * (k))]
                    + tmp_cc[(ido - 1) + ido * ((2) + cdim * (k))];
                tr1 = tmp_cc[(ido - 1) + ido * ((0) + cdim * (k))]
                    - tmp_cc[(ido - 1) + ido * ((2) + cdim * (k))];
            }
            tmp_ch[(ido - 1) + ido * ((k) + l1 * (0))] = tr2 + tr2;
            tmp_ch[(ido - 1) + ido * ((k) + l1 * (1))] = sqrt2 * (tr1 - ti1);
            tmp_ch[(ido - 1) + ido * ((k) + l1 * (2))] = ti2 + ti2;
            tmp_ch[(ido - 1) + ido * ((k) + l1 * (3))] = -sqrt2 * (tr1 + ti1);
            k += 1;
        }
    }
    if ido <= 2 {
        return;
    };
    let mut k: usize = 0;
    while k < l1 {
        let mut i: usize = 2;
        while i < ido {
            let mut ci2: f64 = 0.0;
            let mut ci3: f64 = 0.0;
            let mut ci4: f64 = 0.0;
            let mut cr2: f64 = 0.0;
            let mut cr3: f64 = 0.0;
            let mut cr4: f64 = 0.0;
            let mut ti1: f64 = 0.0;
            let mut ti2: f64 = 0.0;
            let mut ti3: f64 = 0.0;
            let mut ti4: f64 = 0.0;
            let mut tr1: f64 = 0.0;
            let mut tr2: f64 = 0.0;
            let mut tr3: f64 = 0.0;
            let mut tr4: f64 = 0.0;
            let ic: usize = ido - i;
            {
                tr2 = tmp_cc[(i - 1) + ido * ((0) + cdim * (k))]
                    + tmp_cc[(ic - 1) + ido * ((3) + cdim * (k))];
                tr1 = tmp_cc[(i - 1) + ido * ((0) + cdim * (k))]
                    - tmp_cc[(ic - 1) + ido * ((3) + cdim * (k))];
            }
            {
                ti1 = tmp_cc[(i) + ido * ((0) + cdim * (k))]
                    + tmp_cc[(ic) + ido * ((3) + cdim * (k))];
                ti2 = tmp_cc[(i) + ido * ((0) + cdim * (k))]
                    - tmp_cc[(ic) + ido * ((3) + cdim * (k))];
            }
            {
                tr4 = tmp_cc[(i) + ido * ((2) + cdim * (k))]
                    + tmp_cc[(ic) + ido * ((1) + cdim * (k))];
                ti3 = tmp_cc[(i) + ido * ((2) + cdim * (k))]
                    - tmp_cc[(ic) + ido * ((1) + cdim * (k))];
            }
            {
                tr3 = tmp_cc[(i - 1) + ido * ((2) + cdim * (k))]
                    + tmp_cc[(ic - 1) + ido * ((1) + cdim * (k))];
                ti4 = tmp_cc[(i - 1) + ido * ((2) + cdim * (k))]
                    - tmp_cc[(ic - 1) + ido * ((1) + cdim * (k))];
            }
            {
                tmp_ch[(i - 1) + ido * ((k) + l1 * (0))] = tr2 + tr3;
                cr3 = tr2 - tr3;
            }
            {
                tmp_ch[(i) + ido * ((k) + l1 * (0))] = ti2 + ti3;
                ci3 = ti2 - ti3;
            }
            {
                cr4 = tr1 + tr4;
                cr2 = tr1 - tr4;
            }
            {
                ci2 = ti1 + ti4;
                ci4 = ti1 - ti4;
            }
            {
                tmp_ch[(i) + ido * ((k) + l1 * (1))] = tmp_wa[(i - 2) + (0) * (ido - 1)] * ci2
                    + tmp_wa[(i - 1) + (0) * (ido - 1)] * cr2;
                tmp_ch[(i - 1) + ido * ((k) + l1 * (1))] = tmp_wa[(i - 2) + (0) * (ido - 1)] * cr2
                    - tmp_wa[(i - 1) + (0) * (ido - 1)] * ci2;
            }
            {
                tmp_ch[(i) + ido * ((k) + l1 * (2))] = tmp_wa[(i - 2) + (1) * (ido - 1)] * ci3
                    + tmp_wa[(i - 1) + (1) * (ido - 1)] * cr3;
                tmp_ch[(i - 1) + ido * ((k) + l1 * (2))] = tmp_wa[(i - 2) + (1) * (ido - 1)] * cr3
                    - tmp_wa[(i - 1) + (1) * (ido - 1)] * ci3;
            }
            {
                tmp_ch[(i) + ido * ((k) + l1 * (3))] = tmp_wa[(i - 2) + (2) * (ido - 1)] * ci4
                    + tmp_wa[(i - 1) + (2) * (ido - 1)] * cr4;
                tmp_ch[(i - 1) + ido * ((k) + l1 * (3))] = tmp_wa[(i - 2) + (2) * (ido - 1)] * cr4
                    - tmp_wa[(i - 1) + (2) * (ido - 1)] * ci4;
            }
            i += 2;
        }
        k += 1;
    }
}

unsafe fn radb2(ido: usize, l1: usize, cc: *const f64, ch: *mut f64, wa: *const f64) {
    let cdim: usize = 2;

    let tmp_cc = from_raw_parts(cc, l1); //fixme: check actual length of cc
    let tmp_ch = from_raw_parts_mut(ch, l1); //fixme: check actual length of ch
    let tmp_wa = from_raw_parts(wa, l1); //fixme: check actual length of wa
    let mut k: usize = 0;
    while k < l1 {
        tmp_ch[(0) + ido * ((k) + l1 * (0))] =
            tmp_cc[(0) + ido * ((0) + cdim * (k))] + tmp_cc[(ido - 1) + ido * ((1) + cdim * (k))];
        tmp_ch[(0) + ido * ((k) + l1 * (1))] =
            tmp_cc[(0) + ido * ((0) + cdim * (k))] - tmp_cc[(ido - 1) + ido * ((1) + cdim * (k))];
        k += 1;
    }
    if (ido & 1) == 0 {
        let mut k: usize = 0;
        while k < l1 {
            tmp_ch[(ido - 1) + ido * ((k) + l1 * (0))] =
                2.0 * tmp_cc[(ido - 1) + ido * ((0) + cdim * (k))];
            tmp_ch[(ido - 1) + ido * ((k) + l1 * (1))] =
                -2.0 * tmp_cc[(0) + ido * ((1) + cdim * (k))];
            k += 1;
        }
    }
    if ido <= 2 {
        return;
    }
    let mut k: usize = 0;
    while k < l1 {
        let mut i: usize = 2;
        while i < ido {
            let ic = ido - i;
            let mut ti2: f64 = 0.0;
            let mut tr2: f64 = 0.0;
            tr2;
            {
                tmp_ch[(i - 1) + ido * ((k) + l1 * (0))] = tmp_cc
                    [(i - 1) + ido * ((0) + cdim * (k))]
                    + tmp_cc[(ic - 1) + ido * ((1) + cdim * (k))];
                tr2 = tmp_cc[(i - 1) + ido * ((0) + cdim * (k))]
                    - tmp_cc[(ic - 1) + ido * ((1) + cdim * (k))];
            }
            {
                ti2 = tmp_cc[(i) + ido * ((0) + cdim * (k))]
                    + tmp_cc[(ic) + ido * ((1) + cdim * (k))];
                tmp_ch[(i) + ido * ((k) + l1 * (0))] = tmp_cc[(i) + ido * ((0) + cdim * (k))]
                    - tmp_cc[(ic) + ido * ((1) + cdim * (k))];
            }
            {
                tmp_ch[(i) + ido * ((k) + l1 * (1))] = tmp_wa[(i - 2) + (0) * (ido - 1)] * ti2
                    + tmp_wa[(i - 1) + (0) * (ido - 1)] * tr2;
                tmp_ch[(i - 1) + ido * ((k) + l1 * (1))] = tmp_wa[(i - 2) + (0) * (ido - 1)] * tr2
                    - tmp_wa[(i - 1) + (0) * (ido - 1)] * ti2;
            }
            i += 2;
        }
        k += 1;
    }
}

unsafe fn radb3(ido: usize, l1: usize, cc: *const f64, ch: *mut f64, wa: *const f64) {
    let cdim = 3;
    let taur: f64 = -0.5;
    let taui: f64 = 0.86602540378443864676;

    let tmp_cc = from_raw_parts(cc, l1); //fixme: check actual length of cc
    let tmp_ch = from_raw_parts_mut(ch, l1); //fixme: check actual length of ch
    let tmp_wa = from_raw_parts(wa, l1); //fixme: check actual length of wa
    let mut k: usize = 0;
    while k < l1 {
        let tr2: f64 = 2.0 * tmp_cc[(ido - 1) + ido * ((1) + cdim * (k))];
        let cr2: f64 = tmp_cc[(0) + ido * ((0) + cdim * (k))] + taur * tr2;
        tmp_ch[(0) + ido * ((k) + l1 * (0))] = tmp_cc[(0) + ido * ((0) + cdim * (k))] + tr2;
        let ci3: f64 = 2.0 * taui * tmp_cc[(0) + ido * ((2) + cdim * (k))];
        {
            tmp_ch[(0) + ido * ((k) + l1 * (2))] = cr2 + ci3;
            tmp_ch[(0) + ido * ((k) + l1 * (1))] = cr2 - ci3;
        };
        k += 1;
    }
    if ido == 1 {
        return;
    }
    let mut k: usize = 0;
    while k < l1 {
        let mut i: usize = 2;
        while i < ido {
            let ic = ido - i;
            let tr2: f64 = tmp_cc[(i - 1) + ido * ((2) + cdim * (k))]
                + tmp_cc[(ic - 1) + ido * ((1) + cdim * (k))];
            let ti2: f64 =
                tmp_cc[(i) + ido * ((2) + cdim * (k))] - tmp_cc[(ic) + ido * ((1) + cdim * (k))];
            let mut cr2: f64 = tmp_cc[(i - 1) + ido * ((0) + cdim * (k))] + taur * tr2;
            let mut ci2: f64 = tmp_cc[(i) + ido * ((0) + cdim * (k))] + taur * ti2;
            tmp_ch[(i - 1) + ido * ((k) + l1 * (0))] =
                tmp_cc[(i - 1) + ido * ((0) + cdim * (k))] + tr2;
            tmp_ch[(i) + ido * ((k) + l1 * (0))] = tmp_cc[(i) + ido * ((0) + cdim * (k))] + ti2;
            let cr3: f64 = taui
                * (tmp_cc[(i - 1) + ido * ((2) + cdim * (k))]
                    - tmp_cc[(ic - 1) + ido * ((1) + cdim * (k))]);
            let ci3: f64 = taui
                * (tmp_cc[(i) + ido * ((2) + cdim * (k))]
                    + tmp_cc[(ic) + ido * ((1) + cdim * (k))]);
            let mut di2: f64 = 0.0;
            let mut di3: f64 = 0.0;
            let mut dr2: f64 = 0.0;
            let mut dr3: f64 = 0.0;
            {
                dr3 = cr2 + ci3;
                dr2 = cr2 - ci3;
            }
            {
                di2 = ci2 + cr3;
                di3 = ci2 - cr3;
            }
            {
                tmp_ch[(i) + ido * ((k) + l1 * (1))] = tmp_wa[(i - 2) + (0) * (ido - 1)] * di2
                    + tmp_wa[(i - 1) + (0) * (ido - 1)] * dr2;
                tmp_ch[(i - 1) + ido * ((k) + l1 * (1))] = tmp_wa[(i - 2) + (0) * (ido - 1)] * dr2
                    - tmp_wa[(i - 1) + (0) * (ido - 1)] * di2;
            }
            {
                tmp_ch[(i) + ido * ((k) + l1 * (2))] = tmp_wa[(i - 2) + (1) * (ido - 1)] * di3
                    + tmp_wa[(i - 1) + (1) * (ido - 1)] * dr3;
                tmp_ch[(i - 1) + ido * ((k) + l1 * (2))] = tmp_wa[(i - 2) + (1) * (ido - 1)] * dr3
                    - tmp_wa[(i - 1) + (1) * (ido - 1)] * di3;
            }
            i += 2;
        }
        k += 1;
    }
}

unsafe fn radb5(ido: usize, l1: usize, cc: *const f64, ch: *mut f64, wa: *const f64) {
    let cdim: usize = 5;
    let tr11: f64 = 0.3090169943749474241;
    let ti11: f64 = 0.95105651629515357212;
    let tr12: f64 = -0.8090169943749474241;
    let ti12: f64 = 0.58778525229247312917;

    let tmp_cc = from_raw_parts(cc, l1); //fixme: check actual length of cc
    let tmp_ch = from_raw_parts_mut(ch, l1); //fixme: check actual length of ch
    let tmp_wa = from_raw_parts(wa, l1); //fixme: check actual length of wa

    let mut k: usize = 0;
    while k < l1 {
        let ti5 = tmp_cc[(0) + ido * ((2) + cdim * (k))] + tmp_cc[(0) + ido * ((2) + cdim * (k))];
        let ti4 = tmp_cc[(0) + ido * ((4) + cdim * (k))] + tmp_cc[(0) + ido * ((4) + cdim * (k))];
        let tr2 = tmp_cc[(ido - 1) + ido * ((1) + cdim * (k))]
            + tmp_cc[(ido - 1) + ido * ((1) + cdim * (k))];
        let tr3 = tmp_cc[(ido - 1) + ido * ((3) + cdim * (k))]
            + tmp_cc[(ido - 1) + ido * ((3) + cdim * (k))];
        tmp_ch[(0) + ido * ((k) + l1 * (0))] = tmp_cc[(0) + ido * ((0) + cdim * (k))] + tr2 + tr3;
        let cr2 = tmp_cc[(0) + ido * ((0) + cdim * (k))] + tr11 * tr2 + tr12 * tr3;
        let cr3 = tmp_cc[(0) + ido * ((0) + cdim * (k))] + tr12 * tr2 + tr11 * tr3;
        let mut ci4: f64 = 0.0;
        let mut ci5: f64 = 0.0;
        {
            ci5 = ti5 * ti11 + ti4 * ti12;
            ci4 = ti5 * ti12 - ti4 * ti11;
        }
        {
            tmp_ch[(0) + ido * ((k) + l1 * (4))] = cr2 + ci5;
            tmp_ch[(0) + ido * ((k) + l1 * (1))] = cr2 - ci5;
        }
        {
            tmp_ch[(0) + ido * ((k) + l1 * (3))] = cr3 + ci4;
            tmp_ch[(0) + ido * ((k) + l1 * (2))] = cr3 - ci4;
        }
        k += 1;
    }
    if ido == 1 {
        return;
    }
    let mut k: usize = 0;
    while k < l1 {
        let mut i: usize = 2;
        while i < ido {
            let mut ic: usize = ido - i;
            let mut tr2: f64 = 0.0;
            let mut tr3: f64 = 0.0;
            let mut tr4: f64 = 0.0;
            let mut tr5: f64 = 0.0;
            let mut ti2: f64 = 0.0;
            let mut ti3: f64 = 0.0;
            let mut ti4: f64 = 0.0;
            let mut ti5: f64 = 0.0;
            {
                tr2 = tmp_cc[(i - 1) + ido * ((2) + cdim * (k))]
                    + tmp_cc[(ic - 1) + ido * ((1) + cdim * (k))];
                tr5 = tmp_cc[(i - 1) + ido * ((2) + cdim * (k))]
                    - tmp_cc[(ic - 1) + ido * ((1) + cdim * (k))];
            }
            {
                ti5 = tmp_cc[(i) + ido * ((2) + cdim * (k))]
                    + tmp_cc[(ic) + ido * ((1) + cdim * (k))];
                ti2 = tmp_cc[(i) + ido * ((2) + cdim * (k))]
                    - tmp_cc[(ic) + ido * ((1) + cdim * (k))];
            }
            {
                tr3 = tmp_cc[(i - 1) + ido * ((4) + cdim * (k))]
                    + tmp_cc[(ic - 1) + ido * ((3) + cdim * (k))];
                tr4 = tmp_cc[(i - 1) + ido * ((4) + cdim * (k))]
                    - tmp_cc[(ic - 1) + ido * ((3) + cdim * (k))];
            }
            {
                ti4 = tmp_cc[(i) + ido * ((4) + cdim * (k))]
                    + tmp_cc[(ic) + ido * ((3) + cdim * (k))];
                ti3 = tmp_cc[(i) + ido * ((4) + cdim * (k))]
                    - tmp_cc[(ic) + ido * ((3) + cdim * (k))];
            }
            tmp_ch[(i - 1) + ido * ((k) + l1 * (0))] =
                tmp_cc[(i - 1) + ido * ((0) + cdim * (k))] + tr2 + tr3;
            tmp_ch[(i) + ido * ((k) + l1 * (0))] =
                tmp_cc[(i) + ido * ((0) + cdim * (k))] + ti2 + ti3;
            let cr2: f64 = tmp_cc[(i - 1) + ido * ((0) + cdim * (k))] + tr11 * tr2 + tr12 * tr3;
            let ci2: f64 = tmp_cc[(i) + ido * ((0) + cdim * (k))] + tr11 * ti2 + tr12 * ti3;
            let cr3: f64 = tmp_cc[(i - 1) + ido * ((0) + cdim * (k))] + tr12 * tr2 + tr11 * tr3;
            let ci3: f64 = tmp_cc[(i) + ido * ((0) + cdim * (k))] + tr12 * ti2 + tr11 * ti3;
            let mut ci4: f64 = 0.0;
            let mut ci5: f64 = 0.0;
            let mut cr5: f64 = 0.0;
            let mut cr4: f64 = 0.0;
            {
                cr5 = tr5 * ti11 + tr4 * ti12;
                cr4 = tr5 * ti12 - tr4 * ti11;
            }
            {
                ci5 = ti5 * ti11 + ti4 * ti12;
                ci4 = ti5 * ti12 - ti4 * ti11;
            }
            let mut dr2: f64 = 0.0;
            let mut dr3: f64 = 0.0;
            let mut dr4: f64 = 0.0;
            let mut dr5: f64 = 0.0;
            let mut di2: f64 = 0.0;
            let mut di3: f64 = 0.0;
            let mut di4: f64 = 0.0;
            let mut di5: f64 = 0.0;
            {
                dr4 = cr3 + ci4;
                dr3 = cr3 - ci4;
            }
            {
                di3 = ci3 + cr4;
                di4 = ci3 - cr4;
            }
            {
                dr5 = cr2 + ci5;
                dr2 = cr2 - ci5;
            }
            {
                di2 = ci2 + cr5;
                di5 = ci2 - cr5;
            }
            {
                tmp_ch[(i) + ido * ((k) + l1 * (1))] = tmp_wa[(i - 2) + (0) * (ido - 1)] * di2
                    + tmp_wa[(i - 1) + (0) * (ido - 1)] * dr2;
                tmp_ch[(i - 1) + ido * ((k) + l1 * (1))] = tmp_wa[(i - 2) + (0) * (ido - 1)] * dr2
                    - tmp_wa[(i - 1) + (0) * (ido - 1)] * di2;
            }
            {
                tmp_ch[(i) + ido * ((k) + l1 * (2))] = tmp_wa[(i - 2) + (1) * (ido - 1)] * di3
                    + tmp_wa[(i - 1) + (1) * (ido - 1)] * dr3;
                tmp_ch[(i - 1) + ido * ((k) + l1 * (2))] = tmp_wa[(i - 2) + (1) * (ido - 1)] * dr3
                    - tmp_wa[(i - 1) + (1) * (ido - 1)] * di3;
            }
            {
                tmp_ch[(i) + ido * ((k) + l1 * (3))] = tmp_wa[(i - 2) + (2) * (ido - 1)] * di4
                    + tmp_wa[(i - 1) + (2) * (ido - 1)] * dr4;
                tmp_ch[(i - 1) + ido * ((k) + l1 * (3))] = tmp_wa[(i - 2) + (2) * (ido - 1)] * dr4
                    - tmp_wa[(i - 1) + (2) * (ido - 1)] * di4;
            }
            {
                tmp_ch[(i) + ido * ((k) + l1 * (4))] = tmp_wa[(i - 2) + (3) * (ido - 1)] * di5
                    + tmp_wa[(i - 1) + (3) * (ido - 1)] * dr5;
                tmp_ch[(i - 1) + ido * ((k) + l1 * (4))] = tmp_wa[(i - 2) + (3) * (ido - 1)] * dr5
                    - tmp_wa[(i - 1) + (3) * (ido - 1)] * di5;
            }
            i += 2;
        }
        k += 1;
    }
}

unsafe fn radbg(
    ido: usize,
    ip: usize,
    l1: usize,
    cc: *mut f64,
    ch: *mut f64,
    wa: *const f64,
    csarr: *const f64,
) {
    let cdim: usize = ip;
    let ipph: usize = (ip + 1) / 2;
    let idl1: usize = ido * l1;

    let tmp_cc = from_raw_parts_mut(cc, l1); //fixme: check actual size of cc
    let tmp_ch = from_raw_parts_mut(ch, l1); //fixme: check actual size of ch
    let tmp_wa = from_raw_parts(wa, l1); //fixme: check actual size of wa
    let tmp_csarr = from_raw_parts(csarr, l1); //fixme: check actual size of csarr
    let mut k: usize = 0;
    while k < l1 {
        let mut i: usize = 0;
        while i < ido {
            tmp_ch[(i) + ido * ((k) + l1 * (0))] = tmp_cc[(i) + ido * ((0) + cdim * (k))];
            i += 1;
        }
        k += 1;
    }
    let mut j: usize = 1;
    let mut jc: usize = ip - 1;
    while j < ipph {
        let j2: usize = 2 * j - 1;
        let mut k: usize = 0;
        while k < l1 {
            tmp_ch[(0) + ido * ((k) + l1 * (j))] =
                2.0 * tmp_cc[(ido - 1) + ido * ((j2) + cdim * (k))];
            tmp_ch[(0) + ido * ((k) + l1 * (jc))] =
                2.0 * tmp_cc[(0) + ido * ((j2 + 1) + cdim * (k))];
            k += 1;
        }
        j += 1;
        jc -= 1;
    }

    if ido != 1 {
        let mut j: usize = 1;
        let mut jc: usize = ip - 1;
        while j < ipph {
            let j2: usize = 2 * j - 1;
            let mut k: usize = 0;
            while k < l1 {
                let mut i: usize = 1;
                let mut ic: usize = ido - i - 2;
                while i <= ido - 2 {
                    tmp_ch[(i) + ido * ((k) + l1 * (j))] = tmp_cc
                        [(i) + ido * ((j2 + 1) + cdim * (k))]
                        + tmp_cc[(ic) + ido * ((j2) + cdim * (k))];
                    tmp_ch[(i) + ido * ((k) + l1 * (jc))] = tmp_cc
                        [(i) + ido * ((j2 + 1) + cdim * (k))]
                        - tmp_cc[(ic) + ido * ((j2) + cdim * (k))];
                    tmp_ch[(i + 1) + ido * ((k) + l1 * (j))] = tmp_cc
                        [(i + 1) + ido * ((j2 + 1) + cdim * (k))]
                        - tmp_cc[(ic + 1) + ido * ((j2) + cdim * (k))];
                    tmp_ch[(i + 1) + ido * ((k) + l1 * (jc))] = tmp_cc
                        [(i + 1) + ido * ((j2 + 1) + cdim * (k))]
                        + tmp_cc[(ic + 1) + ido * ((j2) + cdim * (k))];
                    i += 2;
                    ic -= 2;
                }
                k += 1;
            }
            j += 1;
            jc -= 1;
        }
    }
    let mut l: usize = 1;
    let mut lc = ip - 1;
    while l < ipph {
        let mut ik: usize = 0;
        while ik < idl1 {
            tmp_cc[(ik) + idl1 * (l)] = tmp_ch[(ik) + idl1 * (0)]
                + tmp_csarr[2 * l] * tmp_ch[(ik) + idl1 * (1)]
                + tmp_csarr[4 * l] * tmp_ch[(ik) + idl1 * (2)];
            tmp_cc[(ik) + idl1 * (lc)] = tmp_csarr[2 * l + 1] * tmp_ch[(ik) + idl1 * (ip - 1)]
                + tmp_csarr[4 * l + 1] * tmp_ch[(ik) + idl1 * (ip - 2)];
            ik += 1;
        }
        let mut iang: usize = 2 * l;
        let mut j: usize = 3;
        let mut jc = ip - 3;
        while j < ipph - 3 {
            iang += l;
            if iang > ip {
                iang -= ip;
            }
            let ar1 = tmp_csarr[2 * iang];
            let ai1 = tmp_csarr[2 * iang + 1];
            iang += l;
            if iang > ip {
                iang -= ip;
            }
            let ar2 = tmp_csarr[2 * iang];
            let ai2 = tmp_csarr[2 * iang + 1];
            iang += l;
            if iang > ip {
                iang -= ip;
            }
            let ar3 = tmp_csarr[2 * iang];
            let ai3 = tmp_csarr[2 * iang + 1];
            iang += l;
            if iang > ip {
                iang -= ip;
            }
            let ar4 = tmp_csarr[2 * iang];
            let ai4 = tmp_csarr[2 * iang + 1];
            let mut ik: usize = 0;
            while ik < idl1 {
                tmp_cc[(ik) + idl1 * (l)] += ar1 * tmp_ch[(ik) + idl1 * (j)]
                    + ar2 * tmp_ch[(ik) + idl1 * (j + 1)]
                    + ar3 * tmp_ch[(ik) + idl1 * (j + 2)]
                    + ar4 * tmp_ch[(ik) + idl1 * (j + 3)];
                tmp_cc[(ik) + idl1 * (lc)] += ai1 * tmp_ch[(ik) + idl1 * (jc)]
                    + ai2 * tmp_ch[(ik) + idl1 * (jc - 1)]
                    + ai3 * tmp_ch[(ik) + idl1 * (jc - 2)]
                    + ai4 * tmp_ch[(ik) + idl1 * (jc - 3)];
                ik += 1;
            }
            j += 4;
            jc -= 4;
        }
        while j < ipph - 1 {
            iang += l;
            if iang > ip {
                iang -= ip;
            }
            let ar1 = tmp_csarr[2 * iang];
            let ai1 = tmp_csarr[2 * iang + 1];
            iang += l;
            if iang > ip {
                iang -= ip;
            }
            let ar2 = tmp_csarr[2 * iang];
            let ai2 = tmp_csarr[2 * iang + 1];
            let mut ik: usize = 0;
            while ik < idl1 {
                tmp_cc[(ik) + idl1 * (l)] +=
                    ar1 * tmp_ch[(ik) + idl1 * (j)] + ar2 * tmp_ch[(ik) + idl1 * (j + 1)];
                tmp_cc[(ik) + idl1 * (lc)] +=
                    ai1 * tmp_ch[(ik) + idl1 * (jc)] + ai2 * tmp_ch[(ik) + idl1 * (jc - 1)];
                ik += 1;
            }
            j += 2;
            jc -= 2;
        }
        while j < ipph {
            iang += l;
            if iang > ip {
                iang -= ip;
            }
            let war = tmp_csarr[2 * iang];
            let wai = tmp_csarr[2 * iang + 1];
            let mut ik: usize = 0;
            while ik < idl1 {
                tmp_cc[(ik) + idl1 * (l)] += war * tmp_ch[(ik) + idl1 * (j)];
                tmp_cc[(ik) + idl1 * (lc)] += wai * tmp_ch[(ik) + idl1 * (jc)];
                ik += 1;
            }
            j += 1;
            jc -= 1;
        }
        l += 1;
        lc -= 1;
    }
    let mut j: usize = 1;
    while j < ipph {
        let mut ik: usize = 0;
        while ik < idl1 {
            tmp_ch[(ik) + idl1 * (0)] += tmp_ch[(ik) + idl1 * (j)];
            ik += 1;
        }
        j += 1;
    }
    let mut j: usize = 1;
    let mut jc: usize = ip - 1;
    while j < ipph {
        let mut k: usize = 0;
        while k < l1 {
            tmp_ch[(0) + ido * ((k) + l1 * (j))] =
                tmp_cc[(0) + ido * ((k) + l1 * (j))] - tmp_cc[(0) + ido * ((k) + l1 * (jc))];
            tmp_ch[(0) + ido * ((k) + l1 * (jc))] =
                tmp_cc[(0) + ido * ((k) + l1 * (j))] + tmp_cc[(0) + ido * ((k) + l1 * (jc))];
            k += 1;
        }
        j += 1;
        jc -= 1;
    }
    if ido == 1 {
        return;
    }

    let mut j: usize = 1;
    let mut jc: usize = ip - 1;
    while j < ipph {
        let mut k: usize = 0;
        while k < l1 {
            let mut i: usize = 1;
            while i <= ido - 2 {
                tmp_ch[(i) + ido * ((k) + l1 * (j))] = tmp_cc[(i) + ido * ((k) + l1 * (j))]
                    - tmp_cc[(i + 1) + ido * ((k) + l1 * (jc))];
                tmp_ch[(i) + ido * ((k) + l1 * (jc))] = tmp_cc[(i) + ido * ((k) + l1 * (j))]
                    + tmp_cc[(i + 1) + ido * ((k) + l1 * (jc))];
                tmp_ch[(i + 1) + ido * ((k) + l1 * (j))] = tmp_cc[(i + 1) + ido * ((k) + l1 * (j))]
                    + tmp_cc[(i) + ido * ((k) + l1 * (jc))];
                tmp_ch[(i + 1) + ido * ((k) + l1 * (jc))] = tmp_cc
                    [(i + 1) + ido * ((k) + l1 * (j))]
                    - tmp_cc[(i) + ido * ((k) + l1 * (jc))];
                i += 2;
            }
            k += 1;
        }
        j += 1;
        jc -= 1;
    }

    let mut j: usize = 1;
    while j < ip {
        let is: usize = (j - 1) * (ido - 1);
        let mut k: usize = 0;
        while k < l1 {
            let mut idij: usize = is;
            let mut i: usize = 1;
            while i <= ido - 2 {
                let t1 = tmp_ch[(i) + ido * ((k) + l1 * (j))];
                let t2 = tmp_ch[(i + 1) + ido * ((k) + l1 * (j))];
                tmp_ch[(i) + ido * ((k) + l1 * (j))] = tmp_wa[idij] * t1 - tmp_wa[idij + 1] * t2;
                tmp_ch[(i + 1) + ido * ((k) + l1 * (j))] =
                    tmp_wa[idij] * t2 + tmp_wa[idij + 1] * t1;
                idij += 2;
                i += 2;
            }
            k += 1;
        }
        j += 1;
    }
}

unsafe fn rfftblue_backward(plan: fftblue_plan, c: &mut [f64], fct: f64) -> i32 {
    let n: usize = (*plan).n;
    let tmp: *mut f64 = (malloc((2 * n) * size_of::<f64>())) as *mut f64;
    if tmp.is_null() {
        return -1;
    }
    *tmp.offset(0) = c[0];
    *tmp.offset(1) = 0.0;
    memcpy(
        tmp.offset(2) as *mut c_void,
        c.as_ptr().offset(1) as *const c_void,
        (n - 1) * size_of::<f64>(),
    );
    if (n & 1) == 0 {
        *tmp.offset((n as isize) + 1) = 0.0;
    }
    let mut m: usize = 2;
    while m < n {
        *tmp.offset(2 * (n as isize) - (m as isize)) = *tmp.offset(m as isize);
        *tmp.offset(2 * (n as isize) - (m as isize) + 1) = *tmp.offset((m as isize) + 1) * (-1.0);
        m += 2;
    }
    if fftblue_fft(plan, from_raw_parts_mut(tmp, n * 2), 1, fct) != 0 {
        free(tmp as *mut c_void);
        return -1;
    }
    m = 0;
    while m < n {
        c[m] = *tmp.offset(2 * (m as isize));
        m += 1;
    }
    free(tmp as *mut c_void);
    return 0;
}

#[no_mangle]
pub unsafe extern "C" fn rfft_forward(plan: rfft_plan, c: *mut f64, fct: f64) -> i32 {
    if (*plan).packplan.is_null() == false {
        return rfftp_forward((*plan).packplan, c, fct);
    } else {
        return rfftblue_forward((*plan).blueplan, c, fct);
    }
}

unsafe fn rfftp_forward(plan: rfftp_plan, c: *mut f64, fct: f64) -> i32 {
    if (*plan).length == 1 {
        return 0;
    }
    let n: usize = (*plan).length;
    let mut l1: usize = n;
    let nf = (*plan).nfct;
    let ch: *mut f64 = (malloc((n) * size_of::<f64>())) as *mut f64;
    if ch.is_null() {
        return -1;
    }
    let mut p1: *mut f64 = c;
    let mut p2: *mut f64 = ch;

    let mut k1: usize = 0;
    while k1 < nf {
        let mut k: usize = nf - k1 - 1;
        let mut ip: usize = (*plan).fct[k].fct;
        let mut ido = n / l1;
        l1 /= ip;
        if ip == 4 {
            radf4(ido, l1, p1, p2, (*plan).fct[k].tw.as_ptr());
        } else if ip == 2 {
            radf2(ido, l1, p1, p2, (*plan).fct[k].tw.as_ptr());
        } else if ip == 3 {
            radf3(ido, l1, p1, p2, (*plan).fct[k].tw.as_ptr());
        } else if ip == 5 {
            radf5(ido, l1, p1, p2, (*plan).fct[k].tw.as_ptr());
        } else {
            radfg(
                ido,
                ip,
                l1,
                p1,
                p2,
                (*plan).fct[k].tw.as_ptr(),
                (*plan).fct[k].tws.as_ptr(),
            );
            let mut tmp_: *mut f64 = p1;
            p1 = p2;
            p2 = tmp_;
        }
        let mut tmp_: *mut f64 = p1;
        p1 = p2;
        p2 = tmp_;
        k1 += 1;
    }
    copy_and_norm(c, p1, n, fct);
    free(ch as *mut c_void);
    return 0;
}
unsafe fn radf4(ido: usize, l1: usize, cc: *const f64, ch: *mut f64, wa: *const f64) {
    let cdim: usize = 4;
    let hsqt2: f64 = 0.70710678118654752440;
    let mut k: usize = 0;
    let tmp_cc = from_raw_parts(cc, l1); //fixme: check actual size of cc
    let tmp_ch = from_raw_parts_mut(ch, l1); //fixme: check actual size of cc
    let tmp_wa = from_raw_parts(wa, l1); //fixme:check actual size of wa
    while k < l1 {
        let mut tr1: f64 = 0.0;
        let mut tr2: f64 = 0.0;
        {
            tr1 = tmp_cc[(0) + ido * ((k) + l1 * (3))] + tmp_cc[(0) + ido * ((k) + l1 * (1))];
            tmp_ch[(0) + ido * ((2) + cdim * (k))] =
                tmp_cc[(0) + ido * ((k) + l1 * (3))] - tmp_cc[(0) + ido * ((k) + l1 * (1))];
        }
        {
            tr2 = tmp_cc[(0) + ido * ((k) + l1 * (0))] + tmp_cc[(0) + ido * ((k) + l1 * (2))];
            tmp_ch[(ido - 1) + ido * ((1) + cdim * (k))] =
                tmp_cc[(0) + ido * ((k) + l1 * (0))] - tmp_cc[(0) + ido * ((k) + l1 * (2))];
        }
        {
            tmp_ch[(0) + ido * ((0) + cdim * (k))] = tr2 + tr1;
            tmp_ch[(ido - 1) + ido * ((3) + cdim * (k))] = tr2 - tr1;
        }
        k += 1;
    }
    if (ido & 1) == 0 {
        let mut k: usize = 0;
        while k < l1 {
            let ti1: f64 = -hsqt2
                * (tmp_cc[(ido - 1) + ido * ((k) + l1 * (1))]
                    + tmp_cc[(ido - 1) + ido * ((k) + l1 * (3))]);
            let tr1: f64 = hsqt2
                * (tmp_cc[(ido - 1) + ido * ((k) + l1 * (1))]
                    - tmp_cc[(ido - 1) + ido * ((k) + l1 * (3))]);
            {
                tmp_ch[(ido - 1) + ido * ((0) + cdim * (k))] =
                    tmp_cc[(ido - 1) + ido * ((k) + l1 * (0))] + tr1;
                tmp_ch[(ido - 1) + ido * ((2) + cdim * (k))] =
                    tmp_cc[(ido - 1) + ido * ((k) + l1 * (0))] - tr1;
            }
            {
                tmp_ch[(0) + ido * ((3) + cdim * (k))] =
                    ti1 + tmp_cc[(ido - 1) + ido * ((k) + l1 * (2))];
                tmp_ch[(0) + ido * ((1) + cdim * (k))] =
                    ti1 - tmp_cc[(ido - 1) + ido * ((k) + l1 * (2))];
            }
            k += 1;
        }
    }
    if ido <= 2 {
        return;
    }
    let mut k: usize = 0;
    while k < l1 {
        let mut i: usize = 2;
        while i < ido {
            let ic: usize = ido - i;
            let mut ci2: f64 = 0.0;
            let mut ci3: f64 = 0.0;
            let mut ci4: f64 = 0.0;
            let mut cr2: f64 = 0.0;
            let mut cr3: f64 = 0.0;
            let mut cr4: f64 = 0.0;
            let mut ti1: f64 = 0.0;
            let mut ti2: f64 = 0.0;
            let mut ti3: f64 = 0.0;
            let mut ti4: f64 = 0.0;
            let mut tr1: f64 = 0.0;
            let mut tr2: f64 = 0.0;
            let mut tr3: f64 = 0.0;
            let mut tr4: f64 = 0.0;
            {
                cr2 = tmp_wa[(i - 2) + (0) * (ido - 1)] * tmp_cc[(i - 1) + ido * ((k) + l1 * (1))]
                    + tmp_wa[(i - 1) + (0) * (ido - 1)] * tmp_cc[(i) + ido * ((k) + l1 * (1))];
                ci2 = tmp_wa[(i - 2) + (0) * (ido - 1)] * tmp_cc[(i) + ido * ((k) + l1 * (1))]
                    - tmp_wa[(i - 1) + (0) * (ido - 1)] * tmp_cc[(i - 1) + ido * ((k) + l1 * (1))];
            }
            {
                cr3 = tmp_wa[(i - 2) + (1) * (ido - 1)] * tmp_cc[(i - 1) + ido * ((k) + l1 * (2))]
                    + tmp_wa[(i - 1) + (1) * (ido - 1)] * tmp_cc[(i) + ido * ((k) + l1 * (2))];
                ci3 = tmp_wa[(i - 2) + (1) * (ido - 1)] * tmp_cc[(i) + ido * ((k) + l1 * (2))]
                    - tmp_wa[(i - 1) + (1) * (ido - 1)] * tmp_cc[(i - 1) + ido * ((k) + l1 * (2))];
            }
            {
                cr4 = tmp_wa[(i - 2) + (2) * (ido - 1)] * tmp_cc[(i - 1) + ido * ((k) + l1 * (3))]
                    + tmp_wa[(i - 1) + (2) * (ido - 1)] * tmp_cc[(i) + ido * ((k) + l1 * (3))];
                ci4 = tmp_wa[(i - 2) + (2) * (ido - 1)] * tmp_cc[(i) + ido * ((k) + l1 * (3))]
                    - tmp_wa[(i - 1) + (2) * (ido - 1)] * tmp_cc[(i - 1) + ido * ((k) + l1 * (3))];
            }
            {
                tr1 = cr4 + cr2;
                tr4 = cr4 - cr2;
            }
            {
                ti1 = ci2 + ci4;
                ti4 = ci2 - ci4;
            }
            {
                tr2 = tmp_cc[(i - 1) + ido * ((k) + l1 * (0))] + cr3;
                tr3 = tmp_cc[(i - 1) + ido * ((k) + l1 * (0))] - cr3;
            }
            {
                ti2 = tmp_cc[(i) + ido * ((k) + l1 * (0))] + ci3;
                ti3 = tmp_cc[(i) + ido * ((k) + l1 * (0))] - ci3;
            }
            {
                tmp_ch[(i - 1) + ido * ((0) + cdim * (k))] = tr2 + tr1;
                tmp_ch[(ic - 1) + ido * ((3) + cdim * (k))] = tr2 - tr1;
            }
            {
                tmp_ch[(i) + ido * ((0) + cdim * (k))] = ti1 + ti2;
                tmp_ch[(ic) + ido * ((3) + cdim * (k))] = ti1 - ti2;
            }
            {
                tmp_ch[(i - 1) + ido * ((2) + cdim * (k))] = tr3 + ti4;
                tmp_ch[(ic - 1) + ido * ((1) + cdim * (k))] = tr3 - ti4;
            }
            {
                tmp_ch[(i) + ido * ((2) + cdim * (k))] = tr4 + ti3;
                tmp_ch[(ic) + ido * ((1) + cdim * (k))] = tr4 - ti3;
            }
            i += 2;
        }
        k += 1;
    }
}

unsafe fn radf2(ido: usize, l1: usize, cc: *const f64, ch: *mut f64, wa: *const f64) {
    let cdim: usize = 2;

    let tmp_cc = from_raw_parts(cc, l1); //fixme: check actual size of cc
    let tmp_ch = from_raw_parts_mut(ch, l1); //fixme: check actual size of ch
    let tmp_wa = from_raw_parts(wa, l1); //fixme: check actual size of cc
    let mut k: usize = 0;
    while k < l1 {
        tmp_ch[(0) + ido * ((0) + cdim * (k))] =
            tmp_cc[(0) + ido * ((k) + l1 * (0))] + tmp_cc[(0) + ido * ((k) + l1 * (1))];
        tmp_ch[(ido - 1) + ido * ((1) + cdim * (k))] =
            tmp_cc[(0) + ido * ((k) + l1 * (0))] - tmp_cc[(0) + ido * ((k) + l1 * (1))];
        k += 1;
    }
    if (ido & 1) == 0 {
        let mut k: usize = 0;
        while k < l1 {
            tmp_ch[(0) + ido * ((1) + cdim * (k))] = -tmp_cc[(ido - 1) + ido * ((k) + l1 * (1))];
            tmp_ch[(ido - 1) + ido * ((0) + cdim * (k))] =
                tmp_cc[(ido - 1) + ido * ((k) + l1 * (0))];
            k += 1;
        }
    }
    if ido <= 2 {
        return;
    }
    let mut k: usize = 0;
    while k < l1 {
        let mut i: usize = 2;
        while i < ido {
            let ic: usize = ido - i;
            let mut tr2: f64 = 0.0;
            let mut ti2: f64 = 0.0;
            {
                tr2 = tmp_wa[(i - 2) + (0) * (ido - 1)] * tmp_cc[(i - 1) + ido * ((k) + l1 * (1))]
                    + tmp_wa[(i - 1) + (0) * (ido - 1)] * tmp_cc[(i) + ido * ((k) + l1 * (1))];
                ti2 = tmp_wa[(i - 2) + (0) * (ido - 1)] * tmp_cc[(i) + ido * ((k) + l1 * (1))]
                    - tmp_wa[(i - 1) + (0) * (ido - 1)] * tmp_cc[(i - 1) + ido * ((k) + l1 * (1))];
            }
            {
                tmp_ch[(i - 1) + ido * ((0) + cdim * (k))] =
                    tmp_cc[(i - 1) + ido * ((k) + l1 * (0))] + tr2;
                tmp_ch[(ic - 1) + ido * ((1) + cdim * (k))] =
                    tmp_cc[(i - 1) + ido * ((k) + l1 * (0))] - tr2;
            }
            {
                tmp_ch[(i) + ido * ((0) + cdim * (k))] = ti2 + tmp_cc[(i) + ido * ((k) + l1 * (0))];
                tmp_ch[(ic) + ido * ((1) + cdim * (k))] =
                    ti2 - tmp_cc[(i) + ido * ((k) + l1 * (0))];
            }
            i += 2;
        }
        k += 1;
    }
}

unsafe fn radf3(ido: usize, l1: usize, cc: *const f64, ch: *mut f64, wa: *const f64) {
    let cdim: usize = 3;
    let taur: f64 = -0.5;
    let taui: f64 = 0.86602540378443864676;

    let tmp_cc = from_raw_parts(cc, l1); //fixme: check actual size of cc
    let tmp_ch = from_raw_parts_mut(ch, l1); //fixme: check actual size of ch
    let tmp_wa = from_raw_parts(wa, l1); //fixme: check actual size of cc

    let mut k: usize = 0;
    while k < l1 {
        let cr2: f64 = tmp_cc[(0) + ido * ((k) + l1 * (1))] + tmp_cc[(0) + ido * ((k) + l1 * (2))];
        tmp_ch[(0) + ido * ((0) + cdim * (k))] = tmp_cc[(0) + ido * ((k) + l1 * (0))] + cr2;
        tmp_ch[(0) + ido * ((2) + cdim * (k))] =
            taui * (tmp_cc[(0) + ido * ((k) + l1 * (2))] - tmp_cc[(0) + ido * ((k) + l1 * (1))]);
        tmp_ch[(ido - 1) + ido * ((1) + cdim * (k))] =
            tmp_cc[(0) + ido * ((k) + l1 * (0))] + taur * cr2;
        k += 1;
    }
    if ido == 1 {
        return;
    }
    let mut k: usize = 0;
    while k < l1 {
        let mut i: usize = 2;
        while i < ido {
            let ic: usize = ido - i;
            let mut di2: f64 = 0.0;
            let mut di3: f64 = 0.0;
            let mut dr2: f64 = 0.0;
            let mut dr3: f64 = 0.0;
            {
                dr2 = tmp_wa[(i - 2) + (0) * (ido - 1)] * tmp_cc[(i - 1) + ido * ((k) + l1 * (1))]
                    + tmp_wa[(i - 1) + (0) * (ido - 1)] * tmp_cc[(i) + ido * ((k) + l1 * (1))];
                di2 = tmp_wa[(i - 2) + (0) * (ido - 1)] * tmp_cc[(i) + ido * ((k) + l1 * (1))]
                    - tmp_wa[(i - 1) + (0) * (ido - 1)] * tmp_cc[(i - 1) + ido * ((k) + l1 * (1))];
            }
            {
                dr3 = tmp_wa[(i - 2) + (1) * (ido - 1)] * tmp_cc[(i - 1) + ido * ((k) + l1 * (2))]
                    + tmp_wa[(i - 1) + (1) * (ido - 1)] * tmp_cc[(i) + ido * ((k) + l1 * (2))];
                di3 = tmp_wa[(i - 2) + (1) * (ido - 1)] * tmp_cc[(i) + ido * ((k) + l1 * (2))]
                    - tmp_wa[(i - 1) + (1) * (ido - 1)] * tmp_cc[(i - 1) + ido * ((k) + l1 * (2))];
            }
            let cr2: f64 = dr2 + dr3;
            let ci2: f64 = di2 + di3;
            tmp_ch[(i - 1) + ido * ((0) + cdim * (k))] =
                tmp_cc[(i - 1) + ido * ((k) + l1 * (0))] + cr2;
            tmp_ch[(i) + ido * ((0) + cdim * (k))] = tmp_cc[(i) + ido * ((k) + l1 * (0))] + ci2;
            let tr2: f64 = tmp_cc[(i - 1) + ido * ((k) + l1 * (0))] + taur * cr2;
            let ti2: f64 = tmp_cc[(i) + ido * ((k) + l1 * (0))] + taur * ci2;
            let tr3: f64 = taui * (di2 - di3);
            let ti3: f64 = taui * (dr3 - dr2);
            {
                tmp_ch[(i - 1) + ido * ((2) + cdim * (k))] = tr2 + tr3;
                tmp_ch[(ic - 1) + ido * ((1) + cdim * (k))] = tr2 - tr3;
            }
            {
                tmp_ch[(i) + ido * ((2) + cdim * (k))] = ti3 + ti2;
                tmp_ch[(ic) + ido * ((1) + cdim * (k))] = ti3 - ti2;
            }
            i += 2;
        }
        k += 1;
    }
}

unsafe fn radf5(ido: usize, l1: usize, cc: *const f64, ch: *mut f64, wa: *const f64) {
    let cdim: usize = 5;
    let tr11: f64 = 0.3090169943749474241;
    let ti11: f64 = 0.95105651629515357212;
    let tr12: f64 = -0.8090169943749474241;
    let ti12: f64 = 0.58778525229247312917;

    let tmp_cc = from_raw_parts(cc, l1); //fixme: check actual size of cc
    let tmp_ch = from_raw_parts_mut(ch, l1); //fixme: check actual size of ch
    let tmp_wa = from_raw_parts(wa, l1); //fixme: check actual size of cc
    let mut k: usize = 0;
    while k < l1 {
        let mut cr2: f64 = 0.0;
        let mut cr3: f64 = 0.0;
        let mut ci4: f64 = 0.0;
        let mut ci5: f64 = 0.0;
        {
            cr2 = tmp_cc[(0) + ido * ((k) + l1 * (4))] + tmp_cc[(0) + ido * ((k) + l1 * (1))];
            ci5 = tmp_cc[(0) + ido * ((k) + l1 * (4))] - tmp_cc[(0) + ido * ((k) + l1 * (1))];
        }
        {
            cr3 = tmp_cc[(0) + ido * ((k) + l1 * (3))] + tmp_cc[(0) + ido * ((k) + l1 * (2))];
            ci4 = tmp_cc[(0) + ido * ((k) + l1 * (3))] - tmp_cc[(0) + ido * ((k) + l1 * (2))];
        }
        tmp_ch[(0) + ido * ((0) + cdim * (k))] = tmp_cc[(0) + ido * ((k) + l1 * (0))] + cr2 + cr3;
        tmp_ch[(ido - 1) + ido * ((1) + cdim * (k))] =
            tmp_cc[(0) + ido * ((k) + l1 * (0))] + tr11 * cr2 + tr12 * cr3;
        tmp_ch[(0) + ido * ((2) + cdim * (k))] = ti11 * ci5 + ti12 * ci4;
        tmp_ch[(ido - 1) + ido * ((3) + cdim * (k))] =
            tmp_cc[(0) + ido * ((k) + l1 * (0))] + tr12 * cr2 + tr11 * cr3;
        tmp_ch[(0) + ido * ((4) + cdim * (k))] = ti12 * ci5 - ti11 * ci4;
        k += 1;
    }
    if ido == 1 {
        return;
    }
    let mut k: usize = 0;
    while k < l1 {
        let mut i: usize = 2;
        while i < ido {
            let mut ci2: f64 = 0.0;
            let mut di2: f64 = 0.0;
            let mut ci4: f64 = 0.0;
            let mut ci5: f64 = 0.0;
            let mut di3: f64 = 0.0;
            let mut di4: f64 = 0.0;
            let mut di5: f64 = 0.0;
            let mut ci3: f64 = 0.0;
            let mut cr2: f64 = 0.0;
            let mut cr3: f64 = 0.0;
            let mut dr2: f64 = 0.0;
            let mut dr3: f64 = 0.0;
            let mut dr4: f64 = 0.0;
            let mut dr5: f64 = 0.0;
            let mut cr5: f64 = 0.0;
            let mut cr4: f64 = 0.0;
            let mut ti2: f64 = 0.0;
            let mut ti3: f64 = 0.0;
            let mut ti5: f64 = 0.0;
            let mut ti4: f64 = 0.0;
            let mut tr2: f64 = 0.0;
            let mut tr3: f64 = 0.0;
            let mut tr4: f64 = 0.0;
            let mut tr5: f64 = 0.0;
            let mut ic: usize = ido - i;
            {
                dr2 = tmp_wa[(i - 2) + (0) * (ido - 1)] * tmp_cc[(i - 1) + ido * ((k) + l1 * (1))]
                    + tmp_wa[(i - 1) + (0) * (ido - 1)] * tmp_cc[(i) + ido * ((k) + l1 * (1))];
                di2 = tmp_wa[(i - 2) + (0) * (ido - 1)] * tmp_cc[(i) + ido * ((k) + l1 * (1))]
                    - tmp_wa[(i - 1) + (0) * (ido - 1)] * tmp_cc[(i - 1) + ido * ((k) + l1 * (1))];
            }
            {
                dr3 = tmp_wa[(i - 2) + (1) * (ido - 1)] * tmp_cc[(i - 1) + ido * ((k) + l1 * (2))]
                    + tmp_wa[(i - 1) + (1) * (ido - 1)] * tmp_cc[(i) + ido * ((k) + l1 * (2))];
                di3 = tmp_wa[(i - 2) + (1) * (ido - 1)] * tmp_cc[(i) + ido * ((k) + l1 * (2))]
                    - tmp_wa[(i - 1) + (1) * (ido - 1)] * tmp_cc[(i - 1) + ido * ((k) + l1 * (2))];
            }
            {
                dr4 = tmp_wa[(i - 2) + (2) * (ido - 1)] * tmp_cc[(i - 1) + ido * ((k) + l1 * (3))]
                    + tmp_wa[(i - 1) + (2) * (ido - 1)] * tmp_cc[(i) + ido * ((k) + l1 * (3))];
                di4 = tmp_wa[(i - 2) + (2) * (ido - 1)] * tmp_cc[(i) + ido * ((k) + l1 * (3))]
                    - tmp_wa[(i - 1) + (2) * (ido - 1)] * tmp_cc[(i - 1) + ido * ((k) + l1 * (3))];
            }
            {
                dr5 = tmp_wa[(i - 2) + (3) * (ido - 1)] * tmp_cc[(i - 1) + ido * ((k) + l1 * (4))]
                    + tmp_wa[(i - 1) + (3) * (ido - 1)] * tmp_cc[(i) + ido * ((k) + l1 * (4))];
                di5 = tmp_wa[(i - 2) + (3) * (ido - 1)] * tmp_cc[(i) + ido * ((k) + l1 * (4))]
                    - tmp_wa[(i - 1) + (3) * (ido - 1)] * tmp_cc[(i - 1) + ido * ((k) + l1 * (4))];
            }
            {
                cr2 = dr5 + dr2;
                ci5 = dr5 - dr2;
            }
            {
                ci2 = di2 + di5;
                cr5 = di2 - di5;
            }
            {
                cr3 = dr4 + dr3;
                ci4 = dr4 - dr3;
            }
            {
                ci3 = di3 + di4;
                cr4 = di3 - di4;
            }
            tmp_ch[(i - 1) + ido * ((0) + cdim * (k))] =
                tmp_cc[(i - 1) + ido * ((k) + l1 * (0))] + cr2 + cr3;
            tmp_ch[(i) + ido * ((0) + cdim * (k))] =
                tmp_cc[(i) + ido * ((k) + l1 * (0))] + ci2 + ci3;
            tr2 = tmp_cc[(i - 1) + ido * ((k) + l1 * (0))] + tr11 * cr2 + tr12 * cr3;
            ti2 = tmp_cc[(i) + ido * ((k) + l1 * (0))] + tr11 * ci2 + tr12 * ci3;
            tr3 = tmp_cc[(i - 1) + ido * ((k) + l1 * (0))] + tr12 * cr2 + tr11 * cr3;
            ti3 = tmp_cc[(i) + ido * ((k) + l1 * (0))] + tr12 * ci2 + tr11 * ci3;
            {
                tr5 = cr5 * ti11 + cr4 * ti12;
                tr4 = cr5 * ti12 - cr4 * ti11;
            }
            {
                ti5 = ci5 * ti11 + ci4 * ti12;
                ti4 = ci5 * ti12 - ci4 * ti11;
            }
            {
                tmp_ch[(i - 1) + ido * ((2) + cdim * (k))] = tr2 + tr5;
                tmp_ch[(ic - 1) + ido * ((1) + cdim * (k))] = tr2 - tr5;
            }
            {
                tmp_ch[(i) + ido * ((2) + cdim * (k))] = ti5 + ti2;
                tmp_ch[(ic) + ido * ((1) + cdim * (k))] = ti5 - ti2;
            }
            {
                tmp_ch[(i - 1) + ido * ((4) + cdim * (k))] = tr3 + tr4;
                tmp_ch[(ic - 1) + ido * ((3) + cdim * (k))] = tr3 - tr4;
            }
            {
                tmp_ch[(i) + ido * ((4) + cdim * (k))] = ti4 + ti3;
                tmp_ch[(ic) + ido * ((3) + cdim * (k))] = ti4 - ti3;
            }
            i += 2;
        }
        k += 1;
    }
}

unsafe fn radfg(
    ido: usize,
    ip: usize,
    l1: usize,
    cc: *mut f64,
    ch: *mut f64,
    wa: *const f64,
    csarr: *const f64,
) {
    let cdim: usize = ip;
    let ipph: usize = (ip + 1) / 2;
    let idl1: usize = ido * l1;

    let tmp_cc = from_raw_parts_mut(cc, l1); //fixme: check actual size of cc
    let tmp_ch = from_raw_parts_mut(ch, l1); //fixme: check actual size of ch
    let tmp_wa = from_raw_parts(wa, l1); //fixme: check actual size of wa
    let tmp_csarr = from_raw_parts(csarr, l1); //fixme: check actual size of csarr

    if ido > 1 {
        let mut j: usize = 1;
        let mut jc: usize = ip - 1;
        while j < ipph {
            let is: usize = (j - 1) * (ido - 1);
            let is2: usize = (jc - 1) * (ido - 1);
            let mut k: usize = 0;
            while k < l1 {
                let mut idij: usize = is;
                let mut idij2: usize = is2;
                let mut i: usize = 1;
                while i <= ido - 2 {
                    let t1: f64 = tmp_cc[(i) + ido * ((k) + l1 * (j))];
                    let t2: f64 = tmp_cc[(i + 1) + ido * ((k) + l1 * (j))];
                    let t3: f64 = tmp_cc[(i) + ido * ((k) + l1 * (jc))];
                    let t4: f64 = tmp_cc[(i + 1) + ido * ((k) + l1 * (jc))];
                    let x1: f64 = tmp_wa[idij] * t1 + tmp_wa[idij + 1] * t2;
                    let x2: f64 = tmp_wa[idij] * t2 - tmp_wa[idij + 1] * t1;
                    let x3: f64 = tmp_wa[idij2] * t3 + tmp_wa[idij2 + 1] * t4;
                    let x4: f64 = tmp_wa[idij2] * t4 - tmp_wa[idij2 + 1] * t3;
                    tmp_cc[(i) + ido * ((k) + l1 * (j))] = x1 + x3;
                    tmp_cc[(i) + ido * ((k) + l1 * (jc))] = x2 - x4;
                    tmp_cc[(i + 1) + ido * ((k) + l1 * (j))] = x2 + x4;
                    tmp_cc[(i + 1) + ido * ((k) + l1 * (jc))] = x3 - x1;
                    idij += 2;
                    idij2 += 2;
                    i += 2;
                }
                k += 1;
            }
            j += 1;
            jc -= 1;
        }
    }

    let mut j: usize = 1;
    let mut jc: usize = ip - 1;
    while j < ipph {
        let mut k: usize = 0;
        while k < l1 {
            let t1: f64 = tmp_cc[(0) + ido * ((k) + l1 * (j))];
            let t2 = tmp_cc[(0) + ido * ((k) + l1 * (jc))];
            tmp_cc[(0) + ido * ((k) + l1 * (j))] = t1 + t2;
            tmp_cc[(0) + ido * ((k) + l1 * (jc))] = t2 - t1;
            k += 1;
        }
        j += 1;
        jc -= 1;
    }

    let mut l: usize = 1;
    let mut lc: usize = ip - 1;
    while l < ipph {
        let mut ik: usize = 0;
        while ik < idl1 {
            tmp_ch[(ik) + idl1 * (l)] = tmp_cc[(ik) + idl1 * (0)]
                + tmp_csarr[2 * l] * tmp_cc[(ik) + idl1 * (1)]
                + tmp_csarr[4 * l] * tmp_cc[(ik) + idl1 * (2)];
            tmp_ch[(ik) + idl1 * (lc)] = tmp_csarr[2 * l + 1] * tmp_cc[(ik) + idl1 * (ip - 1)]
                + tmp_csarr[4 * l + 1] * tmp_cc[(ik) + idl1 * (ip - 2)];
            ik += 1;
        }
        let mut iang: usize = 2 * l;
        let mut j: usize = 3;
        let mut jc = ip - 3;
        while j < ipph - 3 {
            iang += l;
            if iang >= ip {
                iang -= ip;
            }
            let ar1: f64 = tmp_csarr[2 * iang];
            let ai1: f64 = tmp_csarr[2 * iang + 1];
            iang += l;
            if iang >= ip {
                iang -= ip;
            }
            let ar2: f64 = tmp_csarr[2 * iang];
            let ai2: f64 = tmp_csarr[2 * iang + 1];
            iang += l;
            if iang >= ip {
                iang -= ip;
            }
            let ar3: f64 = tmp_csarr[2 * iang];
            let ai3: f64 = tmp_csarr[2 * iang + 1];
            iang += l;
            if iang >= ip {
                iang -= ip;
            }
            let ar4: f64 = tmp_csarr[2 * iang];
            let ai4: f64 = tmp_csarr[2 * iang + 1];
            let mut ik: usize = 0;
            while ik < idl1 {
                tmp_ch[(ik) + idl1 * (l)] += ar1 * tmp_cc[(ik) + idl1 * (j)]
                    + ar2 * tmp_cc[(ik) + idl1 * (j + 1)]
                    + ar3 * tmp_cc[(ik) + idl1 * (j + 2)]
                    + ar4 * tmp_cc[(ik) + idl1 * (j + 3)];
                tmp_ch[(ik) + idl1 * (lc)] += ai1 * tmp_cc[(ik) + idl1 * (jc)]
                    + ai2 * tmp_cc[(ik) + idl1 * (jc - 1)]
                    + ai3 * tmp_cc[(ik) + idl1 * (jc - 2)]
                    + ai4 * tmp_cc[(ik) + idl1 * (jc - 3)];
                ik += 1;
            }
            j += 4;
            jc -= 4;
        }
        while j < ipph - 1 {
            iang += l;
            if iang >= ip {
                iang -= ip;
            }
            let ar1: f64 = tmp_csarr[2 * iang];
            let ai1: f64 = tmp_csarr[2 * iang + 1];
            iang += l;
            if iang >= ip {
                iang -= ip;
            }
            let ar2: f64 = tmp_csarr[2 * iang];
            let ai2: f64 = tmp_csarr[2 * iang + 1];
            let mut ik: usize = 0;
            while ik < idl1 {
                tmp_ch[(ik) + idl1 * (l)] +=
                    ar1 * tmp_cc[(ik) + idl1 * (j)] + ar2 * tmp_cc[(ik) + idl1 * (j + 1)];
                tmp_ch[(ik) + idl1 * (lc)] +=
                    ai1 * tmp_cc[(ik) + idl1 * (jc)] + ai2 * tmp_cc[(ik) + idl1 * (jc - 1)];
                ik += 1;
            }
            j += 2;
            jc -= 2;
        }
        while j < ipph {
            iang += l;
            if iang >= ip {
                iang -= ip;
            }
            let ar: f64 = tmp_csarr[2 * iang];
            let ai: f64 = tmp_csarr[2 * iang + 1];
            let mut ik: usize = 0;
            while ik < idl1 {
                tmp_ch[(ik) + idl1 * (l)] += ar * tmp_cc[(ik) + idl1 * (j)];
                tmp_ch[(ik) + idl1 * (lc)] += ai * tmp_cc[(ik) + idl1 * (jc)];
                ik += 1;
            }
            j += 1;
            jc -= 1;
        }
        l += 1;
        lc -= 1;
    }
    let mut ik: usize = 0;
    while ik < idl1 {
        tmp_ch[(ik) + idl1 * (0)] = tmp_cc[(ik) + idl1 * (0)];
        ik += 1;
    }
    let mut j: usize = 1;
    while j < ipph {
        let mut ik: usize = 0;
        while ik < idl1 {
            tmp_ch[(ik) + idl1 * (0)] += tmp_cc[(ik) + idl1 * (j)];
            ik += 1;
        }
        j += 1;
    }

    let mut k: usize = 0;
    while k < l1 {
        let mut i: usize = 0;
        while i < ido {
            tmp_cc[(i) + ido * ((0) + cdim * (k))] = tmp_ch[(i) + ido * ((k) + l1 * (0))];
            i += 1;
        }
        k += 1;
    }

    let mut j: usize = 1;
    let mut jc: usize = ip - 1;
    while j < ipph {
        let j2: usize = 2 * j - 1;
        let mut k: usize = 0;
        while k < l1 {
            tmp_cc[(ido - 1) + ido * ((j2) + cdim * (k))] = tmp_ch[(0) + ido * ((k) + l1 * (j))];
            tmp_cc[(0) + ido * ((j2 + 1) + cdim * (k))] = tmp_ch[(0) + ido * ((k) + l1 * (jc))];
            k += 1;
        }
        j += 1;
        jc -= 1;
    }

    if ido == 1 {
        return;
    }

    let mut j: usize = 1;
    let mut jc: usize = ip - 1;
    while j < ipph {
        let j2: usize = 2 * j - 1;
        let mut k: usize = 0;
        while k < l1 {
            let mut i: usize = 1;
            let mut ic: usize = ido - i - 2;
            while i <= ido - 2 {
                tmp_cc[(i) + ido * ((j2 + 1) + cdim * (k))] =
                    tmp_ch[(i) + ido * ((k) + l1 * (j))] + tmp_ch[(i) + ido * ((k) + l1 * (jc))];
                tmp_cc[(ic) + ido * ((j2) + cdim * (k))] =
                    tmp_ch[(i) + ido * ((k) + l1 * (j))] - tmp_ch[(i) + ido * ((k) + l1 * (jc))];
                tmp_cc[(i + 1) + ido * ((j2 + 1) + cdim * (k))] = tmp_ch
                    [(i + 1) + ido * ((k) + l1 * (j))]
                    + tmp_ch[(i + 1) + ido * ((k) + l1 * (jc))];
                tmp_cc[(ic + 1) + ido * ((j2) + cdim * (k))] = tmp_ch
                    [(i + 1) + ido * ((k) + l1 * (jc))]
                    - tmp_ch[(i + 1) + ido * ((k) + l1 * (j))];
                i += 2;
                ic -= 2;
            }
            k += 1;
        }
        j += 1;
        jc -= 1;
    }
}

unsafe fn rfftblue_forward(plan: fftblue_plan, c: *mut f64, fct: f64) -> i32 {
    let n: usize = (*plan).n;
    let tmp_len = (2 * n) * size_of::<f64>();
    let tmp: *mut f64 = (malloc(tmp_len)) as *mut f64;
    if tmp.is_null() {
        return -1;
    }
    let mut m: usize = 0;
    while m < n {
        *tmp.offset(2 * (m as isize)) = *c.offset(m as isize);
        *tmp.offset(2 * (m as isize) + 1) = 0.0;
        m += 1;
    }

    let tmp_ref_mut = from_raw_parts_mut(tmp, tmp_len);
    let res = fftblue_fft(plan, tmp_ref_mut, -1, fct);
    if res != 0 {
        free(tmp as *mut c_void);
        return -1;
    }
    *c.offset(0) = *tmp.offset(0);
    memcpy(
        c.offset(1) as *mut c_void,
        tmp.offset(2) as *const c_void,
        (n - 1) * size_of::<f64>(),
    );
    free(tmp as *mut c_void);
    return 0;
}

unsafe fn copy_and_norm(c: *mut f64, p1: *mut f64, n: usize, fct: f64) {
    if p1 != c {
        if fct != 1.0 {
            let mut i: usize = 0;
            while i < n {
                *c.offset(i as isize) = *p1.offset(i as isize) * fct;
                i += 1;
            }
        } else {
            memcpy(c as *mut c_void, p1 as *const c_void, n * size_of::<f64>());
        }
    } else if fct != 1.0 {
        let mut i: usize = 0;
        while i < n {
            *c.offset(i as isize) *= fct;
            i += 1;
        }
    }
}

#[no_mangle]
pub unsafe extern "C" fn rfft_length(plan: rfft_plan)->usize{
  if (*plan).packplan.is_null() {
      let tmp_packplan = (*plan).packplan;
      return (*tmp_packplan).length;
    }
    let tmp_blueplan = (*plan).blueplan;
  return (*tmp_blueplan).n;
  }
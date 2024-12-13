//! Optimized multiexponentiation with precomputations
//!
//! Many ZK proofs often require computing `s^x t^y mod N` with s, t, and N being known in advance.
//! This module provides [`MultiexpTable`] that can compute multiexponent faster.

#![allow(non_snake_case)]

use rug::{Complete, Integer};

/// Precomputed table for performing faster multiexponentiation
#[derive(Debug, Clone)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub struct MultiexpTable {
    s: Vec<Integer>,
    ell_x: Integer,
    s_to_ell_x: Integer,
    t: Vec<Integer>,
    ell_y: Integer,
    t_to_ell_y: Integer,
    N: Integer,
}

impl MultiexpTable {
    /// Builds a multiexponentiation table to perform `s^x t^y mod N` faster
    /// where `x` and `y` are up to `x_bits` and `y_bits`
    ///
    /// Returns `None` is `s` or `t` are non-positive or if any of them are not co-prime to `N` or
    /// if `N` is less than 2.
    pub fn build(s: &Integer, t: &Integer, x_bits: u32, y_bits: u32, N: Integer) -> Option<Self> {
        if s.cmp0().is_le()
            || t.cmp0().is_le()
            || N <= *Integer::ONE
            || s.gcd_ref(&N).complete() != *Integer::ONE
            || t.gcd_ref(&N).complete() != *Integer::ONE
        {
            return None;
        }
        let k_x = x_bits / 8 + 1;
        let k_y = y_bits / 8 + 1;
        let mut s_table = Vec::with_capacity(k_x.try_into().ok()?);
        let mut t_table = Vec::with_capacity(k_y.try_into().ok()?);

        let B: u32 = 256;
        for i in 0..k_x {
            let B_to_i = Integer::u_pow_u(B, i).complete();
            s_table.push(s.clone().pow_mod(&B_to_i, &N).ok()?);
        }
        for i in 0..k_y {
            let B_to_i = Integer::u_pow_u(B, i).complete();
            t_table.push(t.clone().pow_mod(&B_to_i, &N).ok()?);
        }

        // smallest negative value possible for `x`
        let ell_x = -(Integer::ONE.clone() << (k_x * 8)) + 1;
        let s_to_ell_x = s.pow_mod_ref(&ell_x, &N)?.into();
        // smallest negative value possible for `y`
        let ell_y = -(Integer::ONE.clone() << (k_y * 8)) + 1;
        let t_to_ell_y = t.pow_mod_ref(&ell_y, &N)?.into();

        Some(Self {
            s: s_table,
            ell_x,
            s_to_ell_x,
            t: t_table,
            ell_y,
            t_to_ell_y,
            N,
        })
    }

    /// Calculates `s^x t^y mod N`
    ///
    /// Returns `None` if either `x` or `y` do not fit into `x_bits` or `y_bits` provided in [`MultiexpTable::build`].
    pub fn prod_exp(&self, x: &Integer, y: &Integer) -> Option<Integer> {
        let order = rug::integer::Order::Lsf;

        let x_is_neg = x.cmp0().is_lt();
        // `x_digits` correspond to digits of `x` is it's non-negative, and `x - ell_x` otherwise
        let x_digits = if !x_is_neg {
            x.to_digits::<u8>(order)
        } else {
            let x = (x - &self.ell_x).complete();
            if x.cmp0().is_lt() {
                // `x` is less than lower bound
                return None;
            }
            x.to_digits::<u8>(order)
        };

        let y_is_neg = y.cmp0().is_lt();
        // `y_digits` correspond to digits of `y` is it's non-negative, and `y - ell_y` otherwise
        let y_digits = if !y_is_neg {
            y.to_digits::<u8>(order)
        } else {
            let y = (y - &self.ell_y).complete();
            if y.cmp0().is_lt() {
                // `y` is less than lower bound
                return None;
            }
            y.to_digits::<u8>(order)
        };

        if x_digits.len() > self.s.len() || y_digits.len() > self.t.len() {
            // `x` or `y` are higher than upper bound
            return None;
        }

        let mut digits_table = [(); 255].map(|_| None);
        build_digits_table(&mut digits_table, &self.s, &x_digits, &self.N);
        build_digits_table(&mut digits_table, &self.t, &y_digits, &self.N);

        let mut res = Integer::ONE.clone();
        let mut acc = Integer::ONE.clone();
        for d in digits_table.iter().rev() {
            if let Some(d) = d {
                acc = (acc * d) % &self.N;
            }
            res = (res * &acc) % &self.N;
        }

        if x_is_neg {
            res = (res * &self.s_to_ell_x) % &self.N;
        }
        if y_is_neg {
            res = (res * &self.t_to_ell_y) % &self.N;
        }

        Some(res)
    }

    /// Returns max size of exponents (in bits) that can be computed
    ///
    /// Max exponent size is guaranteed to be equal or greater than `x_bits` and `y_bits`
    /// provided in [MultiexpTable::build]
    pub fn max_exponents_size(&self) -> (usize, usize) {
        (self.s.len() * 8, self.t.len() * 8)
    }

    /// Estimates size of the table in RAM in bytes
    pub fn size_in_bytes(&self) -> usize {
        let Self {
            s,
            ell_x,
            s_to_ell_x,
            t,
            ell_y,
            t_to_ell_y,
            N,
        } = self;

        // A few bytes to encode length of Vec `s` and `t`
        let vec_len = 2 * (usize::BITS as usize / 8);
        // And a few bytes more to encode length of each integer
        let int_len = (5 + s.len() + t.len()) * (usize::BITS as usize / 8);

        type Limb = u32;
        let s: usize = s.iter().map(|s_i| s_i.significant_digits::<Limb>()).sum();
        let ell_x = ell_x.significant_digits::<Limb>();
        let s_to_ell_x = s_to_ell_x.significant_digits::<Limb>();
        let t: usize = t.iter().map(|t_i| t_i.significant_digits::<Limb>()).sum();
        let ell_y = ell_y.significant_digits::<Limb>();
        let t_to_ell_y = t_to_ell_y.significant_digits::<Limb>();
        let N = N.significant_digits::<Limb>();

        let limbs_bytes =
            (Limb::BITS as usize / 8) * (s + ell_x + s_to_ell_x + t + ell_y + t_to_ell_y + N);

        vec_len + int_len + limbs_bytes
    }
}

fn build_digits_table(
    table: &mut [Option<Integer>; 255],
    base: &[Integer],
    digits: &[u8],
    N: &Integer,
) {
    for (i, digit) in digits.iter().copied().enumerate() {
        if digit != 0 {
            match &mut table[usize::from(digit - 1)] {
                Some(out) => {
                    *out *= &base[i];
                    *out %= N;
                }
                out @ None => *out = Some(base[i].clone()),
            }
        }
    }
}

#[cfg(test)]
mod test {
    use rug::Integer;

    use super::MultiexpTable;

    #[test]
    fn multiexp_works() {
        let N = Integer::from(100000);
        let s = Integer::from(3);
        let t = Integer::from(7);

        let x_bits = 48;
        let y_bits = 32;

        let table = MultiexpTable::build(&s, &t, x_bits, y_bits, N.clone()).unwrap();

        let mut rng = rug::rand::RandState::new_mersenne_twister();

        for _ in 0..100 {
            let mut x = Integer::from(Integer::random_bits(x_bits, &mut rng));
            if rng.bits(1) == 1 {
                x = -x
            }

            let mut y = Integer::from(Integer::random_bits(y_bits, &mut rng));
            if rng.bits(1) == 1 {
                y = -y
            }
            println!("x={x} y={y}");

            let actual = table.prod_exp(&x, &y).unwrap();
            let expected =
                (s.clone().pow_mod(&x, &N).unwrap() * t.clone().pow_mod(&y, &N).unwrap()) % &N;
            assert_eq!(actual, expected);
        }
    }
}

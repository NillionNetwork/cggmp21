//! Single operation benchmarks

use generic_ec::{Point, SecretScalar};
use generic_ec::curves::Secp256k1;
use crate::security_level::SecurityLevel128;
use crate::utils;
use std::time::Instant;
use rand::rngs::OsRng;
use paillier_zk::{
    fast_paillier,
    rug::{Complete, Integer},
};

fn computation_cost_curve() {
    let mut rng = OsRng;
    let mut total_duration = std::time::Duration::ZERO;
    
    for _ in 0..10 {
        let exponent = SecretScalar::<Secp256k1>::random(&mut rng);

        let start = Instant::now();
        let _group_exponentiation = Point::<Secp256k1>::generator() * &exponent;
        total_duration += start.elapsed();
    }

    let average_duration = total_duration / 10;
    println!("Average time elapsed for group exponentiation: {:?}", average_duration);
}

fn computation_cost_ring(n: Integer) {
    let mut total_duration = std::time::Duration::ZERO;

    for _ in 0..10 {
        let s: Integer = Integer::from(100);
        let mut rng = OsRng;
        let a = SecretScalar::<Secp256k1>::random(&mut rng);
        let start = Instant::now();
        let _ring_exponentiation = s.pow_mod_ref(&crate::utils::scalar_to_bignumber(&a), &n);
        total_duration += start.elapsed();
    }

    let average_duration = total_duration / 10;
    println!("Average time elapsed for ring exponentiation: {:?}", average_duration);
}

fn computation_cost_ring_squared(n: Integer) {
    let mut total_duration = std::time::Duration::ZERO;

    for _ in 0..10 {
        let s = Integer::from(100);
        let beta = Integer::from(100);
        let mut rng = OsRng;
        let exponent = SecretScalar::<Secp256k1>::random(&mut rng);
        let enc_j = fast_paillier::EncryptionKey::from_n(n.clone());
        let cipher = enc_j.encrypt_with(&beta, &s).unwrap();

        let start = Instant::now();
        let _ring_exponentiation_square = enc_j.omul(&utils::scalar_to_bignumber(&exponent), &cipher);
        total_duration += start.elapsed();
    }

    let average_duration = total_duration / 10;
    println!("Average time elapsed for ring squared exponentiation: {:?}", average_duration);
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_computation_cost() {
        let mut rng = OsRng;
        let pregenerated_primes = crate::PregeneratedPrimes::<SecurityLevel128>::generate(&mut rng);
        let (p, q) = pregenerated_primes.split();
        let n: Integer = (&p * &q).complete();
        computation_cost_curve();
        computation_cost_ring(n.clone());
        computation_cost_ring_squared(n.clone());
    }
}
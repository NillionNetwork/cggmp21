use cggmp21::{
    fast_paillier::utils::external_rand,
    rug::{self, Complete},
    security_level::{SecurityLevel, SecurityLevel128},
};
use generic_ec::{curves::Secp256k1 as E, Point, Scalar};

fn criterion_benchmark(c: &mut criterion::Criterion) {
    let mut rng = rand_dev::DevRng::new();
    c.bench_function("scalar at point mult (secp256k1)", |b| {
        b.iter_batched(
            || {
                let base = Point::generator() * Scalar::<E>::random(&mut rng);
                let exp = Scalar::<E>::random(&mut rng);
                (base, exp)
            },
            |(base, exp)| base * exp,
            criterion::BatchSize::SmallInput,
        )
    });

    let primes = cggmp21_tests::CACHED_PRIMES
        .iter::<SecurityLevel128>()
        .next()
        .unwrap();
    let (p, q) = primes.split();
    let n = p * q;

    let bits = [
        256, // something close to the curve order
        SecurityLevel128::ELL + SecurityLevel128::EPSILON,
        SecurityLevel128::ELL_PRIME + SecurityLevel128::EPSILON,
    ];
    let mut rand = external_rand(&mut rng);
    for bits in bits {
        let bits: u32 = bits.try_into().unwrap();
        c.bench_function(&format!("x^e mod N, |e| = {bits}"), |b| {
            b.iter_batched(
                || {
                    let base = rug::Integer::random_below_ref(&n, &mut rand).into();
                    let exp = rug::Integer::random_bits(bits, &mut rand).into();
                    (base, exp)
                },
                |(base, exp): (rug::Integer, rug::Integer)| base.pow_mod(&exp, &n).unwrap(),
                criterion::BatchSize::SmallInput,
            )
        });
    }

    let bits = [
        256, // something close to the curve order
        // SecurityLevel128::ELL, -- don't need it, ell = 256 is a repetition
        SecurityLevel128::ELL + SecurityLevel128::EPSILON,
        n.significant_bits().try_into().unwrap(),
    ];
    let nn = (&n * &n).complete();
    for bits in bits {
        let bits: u32 = bits.try_into().unwrap();
        c.bench_function(&format!("x^e mod N^2, |e| = {bits}"), |b| {
            b.iter_batched(
                || {
                    let x = rug::Integer::random_below_ref(&nn, &mut rand).into();
                    let e = rug::Integer::random_bits(bits, &mut rand).into();
                    (x, e)
                },
                |(x, e): (rug::Integer, rug::Integer)| x.pow_mod(&e, &nn).unwrap(),
                criterion::BatchSize::SmallInput,
            )
        });
    }
}

criterion::criterion_group!(benches, criterion_benchmark);
criterion::criterion_main!(benches);

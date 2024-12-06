use std::iter;

use generic_ec::{Curve, Point};
use rand::{seq::SliceRandom, Rng};
use rand_dev::DevRng;

use cggmp21::{key_share::reconstruct_secret_key, ExecutionId};

cggmp21_tests::test_suite! {
    test = keygen_works,
    generics = all_curves,
    suites = [
        n3(3, false, false),
        n5(5, false, false),
        n7(7, false, false),
        n10(10, false, false),
        n10_reliable(10, true, false),
        #[cfg(feature = "hd-wallet")]
        n3_hd(3, false, true),
        #[cfg(feature = "hd-wallet")]
        n5_hd(5, false, true),
        #[cfg(feature = "hd-wallet")]
        n7_hd(7, false, true),
        #[cfg(feature = "hd-wallet")]
        n10_hd(10, false, true),
    ]
}
fn keygen_works<E: Curve>(n: u16, reliable_broadcast: bool, hd_wallet: bool) {
    #[cfg(not(feature = "hd-wallet"))]
    assert!(!hd_wallet);

    let mut rng = DevRng::new();

    let eid: [u8; 32] = rng.gen();
    let eid = ExecutionId::new(&eid);

    let key_shares = round_based::sim::run(n, |i, party| {
        let mut party_rng = rng.fork();

        async move {
            let keygen =
                cggmp21::keygen::<E>(eid, i, n).enforce_reliable_broadcast(reliable_broadcast);

            #[cfg(feature = "hd-wallet")]
            let keygen = keygen.hd_wallet(hd_wallet);

            keygen.start(&mut party_rng, party).await
        }
    })
    .unwrap()
    .expect_ok()
    .into_vec();

    validate_keygen_output(&mut rng, &key_shares, hd_wallet);
}

cggmp21_tests::test_suite! {
    test = threshold_keygen_works,
    generics = all_curves,
    suites = [
        t2n3(2, 3, false, false),
        t3n5(3, 5, false, false),
        t3n5_reliable(3, 5, true, false),
        #[cfg(feature = "hd-wallet")]
        t2n3_hd(2, 3, false, true),
        #[cfg(feature = "hd-wallet")]
        t3n5_hd(3, 5, false, true),
    ]
}
fn threshold_keygen_works<E: Curve>(t: u16, n: u16, reliable_broadcast: bool, hd_wallet: bool) {
    #[cfg(not(feature = "hd-wallet"))]
    assert!(!hd_wallet);

    let mut rng = DevRng::new();

    let eid: [u8; 32] = rng.gen();
    let eid = ExecutionId::new(&eid);

    let key_shares = round_based::sim::run(n, |i, party| {
        let mut party_rng = rng.fork();

        async move {
            let keygen = cggmp21::keygen::<E>(eid, i, n)
                .enforce_reliable_broadcast(reliable_broadcast)
                .set_threshold(t);

            #[cfg(feature = "hd-wallet")]
            let keygen = keygen.hd_wallet(hd_wallet);

            keygen.start(&mut party_rng, party).await
        }
    })
    .unwrap()
    .expect_ok()
    .into_vec();

    validate_keygen_output(&mut rng, &key_shares, hd_wallet);
}

cggmp21_tests::test_suite! {
    test = threshold_keygen_sync_works,
    generics = all_curves,
    suites = [
        t3n5(3, 5, false),
        #[cfg(feature = "hd-wallet")]
        t3n5_hd(3, 5, true),
    ]
}
fn threshold_keygen_sync_works<E: Curve>(t: u16, n: u16, hd_wallet: bool) {
    #[cfg(not(feature = "hd-wallet"))]
    assert!(!hd_wallet);

    let mut rng = DevRng::new();

    let eid: [u8; 32] = rng.gen();
    let eid = ExecutionId::new(&eid);

    let mut party_rng = iter::repeat_with(|| rng.fork())
        .take(n.into())
        .collect::<Vec<_>>();

    let mut simulation = round_based::sim::Simulation::with_capacity(n);
    for (i, party_rng) in (0..).zip(&mut party_rng) {
        simulation.add_party({
            let keygen = cggmp21::keygen::<E>(eid, i, n).set_threshold(t);

            #[cfg(feature = "hd-wallet")]
            let keygen = keygen.hd_wallet(hd_wallet);

            keygen.into_state_machine(party_rng)
        })
    }
    let key_shares = simulation
        .run()
        .unwrap()
        .into_vec()
        .into_iter()
        .collect::<Result<Vec<_>, _>>()
        .unwrap();
    validate_keygen_output(&mut rng, &key_shares, hd_wallet);
}

fn validate_keygen_output<E: generic_ec::Curve>(
    rng: &mut impl rand::RngCore,
    key_shares: &[cggmp21::IncompleteKeyShare<E>],
    hd_wallet: bool,
) {
    #[cfg(not(feature = "hd-wallet"))]
    assert!(!hd_wallet);

    for (i, key_share) in (0u16..).zip(key_shares) {
        assert_eq!(key_share.i, i);
        assert_eq!(key_share.shared_public_key, key_shares[0].shared_public_key);
        assert_eq!(key_share.public_shares, key_shares[0].public_shares);
        assert_eq!(
            Point::<E>::generator() * &key_share.x,
            key_share.public_shares[usize::from(i)]
        );
    }

    #[cfg(feature = "hd-wallet")]
    if hd_wallet {
        assert!(key_shares[0].chain_code.is_some());
        for key_share in &key_shares[1..] {
            assert_eq!(key_share.chain_code, key_shares[0].chain_code);
        }
    } else {
        for key_share in key_shares {
            assert_eq!(key_share.chain_code, None);
        }
    }

    // Choose `t` random key shares and reconstruct a secret key
    let t = key_shares[0].min_signers();
    let t_shares = key_shares
        .choose_multiple(rng, t.into())
        .cloned()
        .collect::<Vec<_>>();

    let sk = reconstruct_secret_key(&t_shares).unwrap();
    assert_eq!(Point::generator() * sk, key_shares[0].shared_public_key);
}

use cggmp21::{key_share::AnyKeyShare, security_level::SecurityLevel128};
use cggmp21_tests::{convert_from_stark_scalar, convert_stark_scalar};
use generic_ec::{coords::HasAffineX, curves::Stark};
use rand::{seq::SliceRandom, Rng};
use rand_dev::DevRng;

#[test]
fn sign_transaction() {
    let mut rng = DevRng::new();
    let t = Some(2);
    let n = 3;

    let shares = cggmp21_tests::CACHED_SHARES
        .get_shares::<Stark, SecurityLevel128>(t, n, false)
        .expect("retrieve cached shares");

    let eid: [u8; 32] = rng.gen();
    let eid = cggmp21::ExecutionId::new(&eid);

    let fe = |hex| starknet_crypto::FieldElement::from_hex_be(hex).unwrap();
    let sep = |hex, idx| starknet_core::types::SierraEntryPoint {
        selector: fe(hex),
        function_idx: idx,
    };

    let account = starknet_accounts::single_owner::SingleOwnerAccount::new(
        starknet_providers::sequencer::SequencerGatewayProvider::new(
            url::Url::parse("http://example.com/gateway").unwrap(),
            url::Url::parse("http://example.com/feeder").unwrap(),
            fe("c6a14"),
        ),
        starknet_signers::local_wallet::LocalWallet::from_signing_key(
            starknet_signers::SigningKey::from_random(),
        ),
        fe("add4355"),
        fe("c4a1"),
        starknet_accounts::single_owner::ExecutionEncoding::New,
    );
    let declaration = starknet_accounts::Declaration::new(
        std::sync::Arc::new(starknet_core::types::FlattenedSierraClass {
            sierra_program: vec![fe("c0a1")],
            contract_class_version: "v1.3.3.7".to_owned(),
            entry_points_by_type: starknet_core::types::EntryPointsByType {
                constructor: vec![sep("01", 1)],
                external: vec![sep("02", 2)],
                l1_handler: vec![sep("03", 3)],
            },
            abi: "cdecl_for_web3".to_owned(),
        }),
        fe("deadf00d"),
        &account,
    )
    .nonce(fe("404ce"))
    .max_fee(fe("f333"));
    let declaration = declaration.prepared().unwrap();
    let transaction_hash = declaration.transaction_hash();

    // convert to cggmp scalar multiple ways to sanity check
    let bytes = transaction_hash.to_bytes_be();
    let s1 = cggmp21::generic_ec::Scalar::from_be_bytes_mod_order(bytes);
    let s2 = convert_from_stark_scalar(&transaction_hash).unwrap();
    assert_eq!(s1, s2);
    let cggmp_transaction_hash = cggmp21::DataToSign::from_scalar(s2);

    // Choose `t` signers to perform signing
    let t = shares[0].min_signers();
    let mut participants = (0..n).collect::<Vec<_>>();
    participants.shuffle(&mut rng);
    let participants = &participants[..usize::from(t)];
    println!("Signers: {participants:?}");
    let participants_shares = participants.iter().map(|i| &shares[usize::from(*i)]);

    let sig = round_based::simulation::run_with_setup(participants_shares, |i, party, share| {
        let mut party_rng = rng.fork();

        async move {
            cggmp21::signing(eid, i, participants, share)
                .sign(&mut party_rng, party, cggmp_transaction_hash)
                .await
        }
    })
    .unwrap()
    .expect_ok()
    .expect_eq();

    // verify with our lib
    sig.verify(&shares[0].core.shared_public_key, &cggmp_transaction_hash)
        .expect("signature is not valid");

    // verify with starknet lib
    let public_key_x = shares[0].core.shared_public_key.x().unwrap().to_scalar();
    let public_key = convert_stark_scalar(&public_key_x).unwrap();
    let r = convert_stark_scalar(&sig.r).unwrap();
    let s = convert_stark_scalar(&sig.s).unwrap();

    let r = starknet_crypto::verify(&public_key, &transaction_hash, &r, &s).unwrap();
    assert!(r, "failed to verify signature");
}

//! ZK-proof of paillier encryption in range with El-Gamal commitment.
//! Called Пenc-elg or Renc-elg in the CGGMP24
//! paper.
//!
//! ## Description
//!
//! A party P has `key` - public key in paillier
//! cryptosystem. P also has `plaintext`, `nonce`,`a`,`b` and
//! `ciphertext = key.encrypt_with(plaintext, nonce)`,
//! `g_to_a = Point::<E>::generator() * a.to_scalar()`,
//! `g_to_b = Point::<E>::generator() * b.to_scalar()`
//! `g_to_ab_plus_x = Point::<E>::generator() * ((&a * &b + &plaintext).complete()).to_scalar()`
//!
//! P wants to prove that `plaintext` is at most `l` bits, without disclosing
//! it, the `nonce`,`a`, and `b`

//! ## Example
//!
//! ```
//! use paillier_zk::{paillier_encryption_in_range_with_el_gamal as p, IntegerExt};
//! use rug::{Integer, Complete};
//! use generic_ec::{Point, curves::Secp256k1 as E};
//! # mod pregenerated {
//! #     use super::*;
//! #     paillier_zk::load_pregenerated_data!(
//! #         verifier_aux: p::Aux,
//! #         prover_decryption_key: fast_paillier::DecryptionKey,
//! #     );
//! # }
//! # fn main() -> Result<(), Box<dyn std::error::Error>> {
//!
//! let shared_state = "some shared state";
//!
//! let mut rng = rand_core::OsRng;
//! # let mut rng = rand_dev::DevRng::new();
//!
//! // 0. Setup: prover and verifier share common Ring-Pedersen parameters:
//!
//! let aux: p::Aux = pregenerated::verifier_aux();
//! let security = p::SecurityParams {
//!     l: 1024,
//!     epsilon: 128,
//!     q: (Integer::ONE << 128_u32).into(),
//! };
//!
//! // 1. Setup: prover prepares the paillier keys
//!
//! let private_key: fast_paillier::DecryptionKey =
//!     pregenerated::prover_decryption_key();
//! let key = private_key.encryption_key();
//!
//! // 2. Setup: prover has some plaintext and encrypts it
//!
//! let plaintext = Integer::from_rng_pm(&Integer::curve_order::<E>(), &mut rng);
//! let (ciphertext, nonce) = key.encrypt_with_random(&mut rng, &plaintext)?;
//!
//! let a = Integer::from_rng_pm(&Integer::curve_order::<E>(), &mut rng);
//! let b = Integer::from_rng_pm(&Integer::curve_order::<E>(), &mut rng);
//!
//! let g_to_a = Point::<E>::generator() * a.to_scalar();
//! let g_to_b = Point::<E>::generator() * b.to_scalar();
//! let exponent = (&a * &b + &plaintext).complete();
//! let g_to_ab_plus_x = Point::<E>::generator() * exponent.to_scalar();
//!
//! // 3. Prover computes a non-interactive proof that plaintext is at most 1024 bits:
//!
//! let data = p::Data { key, ciphertext: &ciphertext , g_to_a: &g_to_a,
//!     g_to_b: &g_to_b ,g_to_ab_plus_x: &g_to_ab_plus_x};
//! let (commitment, proof) = p::non_interactive::prove::<E,sha2::Sha256>(
//!     &shared_state,
//!     &aux,
//!     data,
//!     p::PrivateData {
//!         plaintext: &plaintext,
//!         nonce: &nonce,
//!         a: &a,
//!         b:&b
//!     },
//!     &security,
//!     &mut rng,
//! )?;
//!
//! // 4. Prover sends this data to verifier
//!
//! # use generic_ec::Curve;
//! # fn send<E: Curve>(_: &p::Data<E>, _: &p::Commitment<E>, _: &p::Proof) {  }
//! send(&data, &commitment, &proof);
//!
//! // 5. Verifier receives the data and the proof and verifies it
//!
//! # let recv = || (data, commitment, proof);
//! let (data, commitment, proof) = recv();
//! p::non_interactive::verify::<E, sha2::Sha256>(
//!     &shared_state,
//!     &aux,
//!     data,
//!     &commitment,
//!     &security,
//!     &proof,
//! );
//! # Ok(()) }
//! ```
//!
//! If the verification succeeded, verifier can continue communication with prover

use fast_paillier::{AnyEncryptionKey, Ciphertext, Nonce};
use generic_ec::{Curve, Point};
use rug::Integer;

#[cfg(feature = "serde")]
use serde::{Deserialize, Serialize};

pub use crate::common::Aux;
pub use crate::common::InvalidProof;

/// Security parameters for proof. Choosing the values is a tradeoff between
/// speed and chance of rejecting a valid proof or accepting an invalid proof
#[derive(Debug, Clone, udigest::Digestable)]
#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
pub struct SecurityParams {
    /// l in paper, security parameter for bit size of plaintext: it needs to
    /// be in range [-2^l; 2^l] or equivalently 2^l
    pub l: usize,
    /// Epsilon in paper, slackness parameter
    pub epsilon: usize,
    /// q in paper. Security parameter for challenge
    #[udigest(as = crate::common::encoding::Integer)]
    pub q: Integer,
}

/// Public data that both parties know
#[derive(Debug, Clone, Copy, udigest::Digestable)]
#[udigest(bound = "")]
pub struct Data<'a, C: Curve> {
    /// N0 in paper, public key that x -> C was encrypted on
    #[udigest(as = crate::common::encoding::AnyEncryptionKey)]
    pub key: &'a dyn AnyEncryptionKey,
    /// C in paper
    #[udigest(as = &crate::common::encoding::Integer)]
    pub ciphertext: &'a Ciphertext,
    /// A=g^a in paper
    pub g_to_a: &'a Point<C>,
    /// B=g^b in paper
    pub g_to_b: &'a Point<C>,
    /// X=g^{ab+x} in paper
    pub g_to_ab_plus_x: &'a Point<C>,
}

/// Private data of prover
#[derive(Clone, Copy)]
pub struct PrivateData<'a> {
    /// x in paper, plaintext of C
    pub plaintext: &'a Integer,
    /// rho in paper, nonce of encryption x -> C
    pub nonce: &'a Nonce,
    /// a in paper, preimage of A
    pub a: &'a Integer,
    /// b in paper, preimage of B
    pub b: &'a Integer,
}

// As described in cggmp24 at page 58
/// Prover's first message, obtained by [`interactive::commit`]
#[derive(Debug, Clone, udigest::Digestable)]
#[udigest(bound = "")]
#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
pub struct Commitment<C: Curve> {
    #[udigest(as = crate::common::encoding::Integer)]
    pub s: Integer,
    #[udigest(as = crate::common::encoding::Integer)]
    pub t: Integer,
    #[udigest(as = crate::common::encoding::Integer)]
    pub d: Integer,
    pub y: Point<C>,
    pub z: Point<C>,
}

/// Prover's data accompanying the commitment. Kept as state between rounds in
/// the interactive protocol.
#[derive(Clone)]
pub struct PrivateCommitment {
    pub alpha: Integer,
    pub mu: Integer,
    pub r: Integer,
    pub beta: Integer,
    pub gamma: Integer,
}

/// Verifier's challenge to prover. Can be obtained deterministically by
/// [`non_interactive::challenge`] or randomly by [`interactive::challenge`]
pub type Challenge = Integer;

// As described in cggmp24 at page 58
/// The ZK proof. Computed by [`interactive::prove`] or
/// [`non_interactive::prove`]
#[derive(Debug, Clone)]
#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
pub struct Proof {
    pub z1: Integer,
    pub z2: Integer,
    pub z3: Integer,
    pub w: Integer,
}

/// The interactive version of the ZK proof. Should be completed in 3 rounds:
/// prover commits to data, verifier responds with a random challenge, and
/// prover gives proof with commitment and challenge.
pub mod interactive {
    use generic_ec::{Curve, Point};
    use rand_core::RngCore;
    use rug::{Complete, Integer};

    use crate::{
        common::{fail_if, fail_if_ne, InvalidProofReason},
        BadExponent, Error,
    };

    use crate::common::{IntegerExt, InvalidProof};

    use super::{
        Aux, Challenge, Commitment, Data, PrivateCommitment, PrivateData, Proof, SecurityParams,
    };

    /// Create random commitment
    pub fn commit<C: Curve, R: RngCore>(
        aux: &Aux,
        data: Data<C>,
        pdata: PrivateData,
        security: &SecurityParams,
        rng: &mut R,
    ) -> Result<(Commitment<C>, PrivateCommitment), Error> {
        let two_to_l_plus_e = (Integer::ONE << (security.l + security.epsilon)).complete();
        let hat_n_at_two_to_l = (Integer::ONE << security.l).complete() * &aux.rsa_modulo;
        let hat_n_at_two_to_l_plus_e =
            (Integer::ONE << (security.l + security.epsilon)).complete() * &aux.rsa_modulo;

        let alpha = Integer::from_rng_pm(&two_to_l_plus_e, rng);
        let mu = Integer::from_rng_pm(&hat_n_at_two_to_l, rng);
        let r = Integer::gen_invertible(data.key.n(), rng);
        let beta = Integer::gen_invertible(&Integer::curve_order::<C>(), rng);
        let gamma = Integer::from_rng_pm(&hat_n_at_two_to_l_plus_e, rng);

        let s = aux.combine(pdata.plaintext, &mu)?;
        let t = aux.combine(&alpha, &gamma)?;
        let d = data.key.encrypt_with(&alpha, &r)?;
        let y = data.g_to_a * beta.to_scalar() + Point::<C>::generator() * alpha.to_scalar();
        let z = Point::<C>::generator() * beta.to_scalar();

        Ok((
            Commitment { s, t, d, y, z },
            PrivateCommitment {
                alpha,
                mu,
                r,
                beta,
                gamma,
            },
        ))
    }

    /// Compute proof for given data and prior protocol values
    pub fn prove<C: Curve>(
        data: Data<C>,
        pdata: PrivateData,
        private_commitment: &PrivateCommitment,
        challenge: &Challenge,
    ) -> Result<Proof, Error> {
        let z1 = (&private_commitment.alpha + (challenge * pdata.plaintext)).complete();
        let w = ((&private_commitment.beta + (challenge * pdata.b)).complete())
            .modulo(&Integer::curve_order::<C>());
        let nonce_to_challenge_mod_n: Integer = pdata
            .nonce
            .pow_mod_ref(challenge, data.key.n())
            .ok_or(BadExponent::undefined())?
            .into();
        let z2 = (&private_commitment.r * nonce_to_challenge_mod_n).modulo(data.key.n());
        let z3 = (&private_commitment.gamma + (challenge * &private_commitment.mu)).complete();
        Ok(Proof { z1, z2, z3, w })
    }

    /// Verify the proof
    pub fn verify<C: Curve>(
        aux: &Aux,
        data: Data<C>,
        commitment: &Commitment<C>,
        security: &SecurityParams,
        challenge: &Challenge,
        proof: &Proof,
    ) -> Result<(), InvalidProof> {
        {
            let lhs = data
                .key
                .encrypt_with(&proof.z1, &proof.z2)
                .map_err(|_| InvalidProofReason::PaillierEnc)?;
            let rhs = {
                let e_at_c = data
                    .key
                    .omul(challenge, data.ciphertext)
                    .map_err(|_| InvalidProofReason::PaillierOp)?;
                data.key
                    .oadd(&commitment.d, &e_at_c)
                    .map_err(|_| InvalidProofReason::PaillierOp)?
            };
            fail_if_ne(InvalidProofReason::EqualityCheck(1), lhs, rhs)?;
        }
        {
            let lhs =
                data.g_to_a * proof.w.to_scalar() + Point::<C>::generator() * proof.z1.to_scalar();
            let rhs = commitment.y + data.g_to_ab_plus_x * challenge.to_scalar();
            fail_if_ne(InvalidProofReason::EqualityCheck(2), lhs, rhs)?;
        }
        {
            let lhs = Point::<C>::generator() * proof.w.to_scalar();
            let rhs = commitment.z + data.g_to_b * challenge.to_scalar();
            fail_if_ne(InvalidProofReason::EqualityCheck(3), lhs, rhs)?;
        }
        {
            let lhs = aux.combine(&proof.z1, &proof.z3)?;
            let s_to_e = aux.pow_mod(&commitment.s, challenge)?;
            let rhs = (&commitment.t * s_to_e).modulo(&aux.rsa_modulo);
            fail_if_ne(InvalidProofReason::EqualityCheck(4), lhs, rhs)?;
        }

        fail_if(
            InvalidProofReason::RangeCheck(5),
            proof
                .z1
                .is_in_pm(&(Integer::ONE << (security.l + security.epsilon)).complete()),
        )?;

        Ok(())
    }

    /// Generate random challenge
    ///
    /// `security` parameter is used to generate challenge in correct range
    pub fn challenge<R: RngCore>(security: &SecurityParams, rng: &mut R) -> Challenge {
        Integer::from_rng_pm(&security.q, rng)
    }
}

/// The non-interactive version of proof. Completed in one round, for example
/// see the documentation of parent module.
pub mod non_interactive {
    use digest::Digest;
    use generic_ec::Curve;

    use crate::{Error, InvalidProof};

    use super::{Aux, Challenge, Commitment, Data, PrivateData, Proof, SecurityParams};

    /// Compute proof for the given data, producing random commitment and
    /// deriving determenistic challenge.
    ///
    /// Obtained from the above interactive proof via Fiat-Shamir heuristic.
    pub fn prove<C: Curve, D: Digest>(
        shared_state: &impl udigest::Digestable,
        aux: &Aux,
        data: Data<C>,
        pdata: PrivateData,
        security: &SecurityParams,
        rng: &mut impl rand_core::RngCore,
    ) -> Result<(Commitment<C>, Proof), Error> {
        let (comm, pcomm) = super::interactive::commit(aux, data, pdata, security, rng)?;
        let challenge = challenge::<C, D>(shared_state, aux, data, &comm, security);
        let proof = super::interactive::prove(data, pdata, &pcomm, &challenge)?;
        Ok((comm, proof))
    }

    /// Verify the proof, deriving challenge independently from same data
    pub fn verify<C: Curve, D: Digest>(
        shared_state: &impl udigest::Digestable,
        aux: &Aux,
        data: Data<C>,
        commitment: &Commitment<C>,
        security: &SecurityParams,
        proof: &Proof,
    ) -> Result<(), InvalidProof> {
        let challenge = challenge::<C, D>(shared_state, aux, data, commitment, security);
        super::interactive::verify(aux, data, commitment, security, &challenge, proof)
    }

    /// Deterministically compute challenge based on prior known values in protocol
    pub fn challenge<C: Curve, D: Digest>(
        shared_state: &impl udigest::Digestable,
        aux: &Aux,
        data: Data<C>,
        commitment: &Commitment<C>,
        security: &SecurityParams,
    ) -> Challenge {
        let tag = "paillier_zk.encryption_in_range_with_el_gamal.ni_challenge";
        let aux = aux.digest_public_data();
        let seed = udigest::inline_struct!(tag {
            shared_state,
            aux,
            security,
            data,
            commitment,
        });
        let mut rng = rand_hash::HashRng::<D, _>::from_seed(seed);
        super::interactive::challenge(security, &mut rng)
    }
}

#[cfg(test)]
mod test {
    use generic_ec::{Curve, Point};
    use rug::{Complete, Integer};
    use sha2::Digest;

    use crate::common::{IntegerExt, InvalidProofReason};

    fn run_with<C: Curve, D: Digest>(
        mut rng: &mut impl rand_core::CryptoRngCore,
        security: super::SecurityParams,
        plaintext: Integer,
        a: Integer,
        b: Integer,
    ) -> Result<(), crate::common::InvalidProof> {
        let aux = crate::common::test::aux(&mut rng);
        let private_key = crate::common::test::random_key(&mut rng).unwrap();
        let key = private_key.encryption_key();
        let (ciphertext, nonce) = key.encrypt_with_random(&mut rng, &plaintext).unwrap();
        let g_to_a = Point::<C>::generator() * a.to_scalar();
        let g_to_b = Point::<C>::generator() * b.to_scalar();
        let exponent = (&a * &b + &plaintext).complete();
        let g_to_exponent = Point::<C>::generator() * exponent.to_scalar();
        let data = super::Data {
            key,
            ciphertext: &ciphertext,
            g_to_a: &g_to_a,
            g_to_b: &g_to_b,
            g_to_ab_plus_x: &g_to_exponent,
        };
        let pdata = super::PrivateData {
            plaintext: &plaintext,
            nonce: &nonce,
            a: &a,
            b: &b,
        };

        let shared_state = "shared state";
        let (commitment, proof) =
            super::non_interactive::prove::<C, D>(&shared_state, &aux, data, pdata, &security, rng)
                .unwrap();
        super::non_interactive::verify::<C, D>(
            &shared_state,
            &aux,
            data,
            &commitment,
            &security,
            &proof,
        )
    }

    fn passing_test<C: Curve, D: Digest>() {
        let mut rng = rand_dev::DevRng::new();
        let security = super::SecurityParams {
            l: 1024,
            epsilon: 300,
            q: (Integer::ONE << 128_u32).into(),
        };
        let plaintext = Integer::from_rng_pm(&Integer::curve_order::<C>(), &mut rng);
        let a = Integer::from_rng_pm(&Integer::curve_order::<C>(), &mut rng);
        let b = Integer::from_rng_pm(&Integer::curve_order::<C>(), &mut rng);
        run_with::<C, D>(&mut rng, security, plaintext, a, b).expect("proof failed");
    }

    fn failing_test<C: Curve, D: Digest>() {
        let mut rng = rand_dev::DevRng::new();
        let security = super::SecurityParams {
            l: 1024,
            epsilon: 300,
            q: (Integer::ONE << 128_u32).complete(),
        };
        let plaintext = (Integer::ONE << (security.l + security.epsilon)).complete() + 1;
        let a = Integer::from_rng_pm(&Integer::curve_order::<C>(), &mut rng);
        let b = Integer::from_rng_pm(&Integer::curve_order::<C>(), &mut rng);
        let r = run_with::<C, D>(&mut rng, security, plaintext, a, b)
            .expect_err("proof should not pass");
        match r.reason() {
            InvalidProofReason::RangeCheck(5) => (),
            e => panic!("proof should not fail with: {e:?}"),
        }
    }

    #[test]
    fn passing_p256() {
        passing_test::<generic_ec::curves::Secp256r1, sha2::Sha256>()
    }
    #[test]
    fn failing_p256_add() {
        failing_test::<generic_ec::curves::Secp256r1, sha2::Sha256>()
    }

    #[test]
    fn passing_million() {
        passing_test::<crate::curve::C, sha2::Sha256>()
    }
    #[test]
    fn failing_million_add() {
        failing_test::<crate::curve::C, sha2::Sha256>()
    }
}

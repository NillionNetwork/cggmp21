//! ZK-proof, called Пlog* or Rlog* in the CGGMP21 paper.
//!
//! ## Description
//!
//! A party P has a number `X = x G`, with G being a generator of
//! curve `E`. P has encrypted x as C. P shares X and C with V and
//! wants to prove that the logarithm of X is the plaintext of C, and that the
//! plaintext (i.e. x) is at most l bits.
//!
//! Given:
//! - `key0`, `pkey0` - pair of public and private keys in paillier cryptosystem
//! - Curve `E`
//! - `X = x G` and `C = key0.encrypt(x)` - data to obtain proof about
//!
//! Prove:
//! - `decrypt(C) = log X`
//! - `bitsize(x) <= l`
//!
//! Disclosing only: `key0`, `C`, `X`
//!
//! ## Example
//!
//! ```rust
//! use rug::{Integer, Complete};
//! use generic_ec::{Point, curves::Secp256k1 as E};
//! use paillier_zk::{group_element_vs_paillier_encryption_in_range as p, IntegerExt};
//! # mod pregenerated {
//! #     use super::*;
//! #     paillier_zk::load_pregenerated_data!(
//! #         verifier_aux: p::Aux,
//! #         prover_decryption_key: fast_paillier::DecryptionKey,
//! #     );
//! # }
//!
//! # fn main() -> Result<(), Box<dyn std::error::Error>> {
//! // Prover and verifier have a shared protocol state
//! let shared_state = "some shared state";
//! let mut rng = rand_core::OsRng;
//! # let mut rng = rand_dev::DevRng::new();
//!
//! // 0. Setup: prover and verifier share common Ring-Pedersen parameters:
//!
//! let aux: p::Aux = pregenerated::verifier_aux();
//! let security = p::SecurityParams {
//!     l: 1024,
//!     epsilon: 300,
//!     q: (Integer::ONE << 128_u32).complete(),
//! };
//!
//! // 1. Setup: prover prepares the paillier keys
//!
//! let private_key: fast_paillier::DecryptionKey =
//!     pregenerated::prover_decryption_key();
//! let key0 = private_key.encryption_key();
//!
//! // 2. Setup: prover has some plaintext `x`, encrypts it and obtains `C`, and computes `X`
//!
//! let x = Integer::from_rng_pm(&(Integer::ONE << security.l).complete(), &mut rng);
//! let (C, nonce) = key0.encrypt_with_random(&mut rng, &x)?;
//! let X = Point::<E>::generator() * x.to_scalar();
//!
//! // 3. Prover computes a non-interactive proof that plaintext is at most `l` bits:
//!
//! let data = p::Data {
//!     key0,
//!     c: &C,
//!     x: &X,
//!     b: &Point::<E>::generator().into(),
//! };
//! let (commitment, proof) =
//!     p::non_interactive::prove::<E, sha2::Sha256>(
//!         &shared_state,
//!         &aux,
//!         data,
//!         p::PrivateData { x: &x, nonce: &nonce },
//!         &security,
//!         &mut rng,
//!     )?;
//!
//! // 4. Prover sends this data to verifier
//!
//! # fn send(_: &p::Data<E>, _: &p::Commitment<E>, _: &p::Proof) {  }
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
//! )?;
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

/// Security parameters for proof. Choosing the values is a tradeoff between
/// speed and chance of rejecting a valid proof or accepting an invalid proof
#[derive(Debug, Clone, udigest::Digestable)]
#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
pub struct SecurityParams {
    /// l in paper, bit size of +-plaintext
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
    /// N0 in paper, public key that C was encrypted on
    #[udigest(as = crate::common::encoding::AnyEncryptionKey)]
    pub key0: &'a dyn AnyEncryptionKey,
    /// C in paper, logarithm of X encrypted on N0
    #[udigest(as = &crate::common::encoding::Integer)]
    pub c: &'a Ciphertext,
    /// A basepoint, generator in group
    pub b: &'a Point<C>,
    /// X in paper, exponent of plaintext of C
    pub x: &'a Point<C>,
}

/// Private data of prover
#[derive(Clone, Copy)]
pub struct PrivateData<'a> {
    /// x in paper, logarithm of X and plaintext of C
    pub x: &'a Integer,
    /// rho in paper, nonce in encryption x -> C
    pub nonce: &'a Nonce,
}

/// Prover's first message, obtained by [`interactive::commit`]
#[derive(Debug, Clone, udigest::Digestable)]
#[udigest(bound = "")]
#[cfg_attr(feature = "serde", derive(Serialize, Deserialize), serde(bound = ""))]
pub struct Commitment<C: Curve> {
    #[udigest(as = crate::common::encoding::Integer)]
    pub s: Integer,
    #[udigest(as = crate::common::encoding::Integer)]
    pub a: Ciphertext,
    pub y: Point<C>,
    #[udigest(as = crate::common::encoding::Integer)]
    pub d: Integer,
}

/// Prover's data accompanying the commitment. Kept as state between rounds in
/// the interactive protocol.
#[derive(Clone)]
pub struct PrivateCommitment {
    pub alpha: Integer,
    pub mu: Integer,
    pub r: Nonce,
    pub gamma: Integer,
}

/// Verifier's challenge to prover. Can be obtained deterministically by
/// [`non_interactive::challenge`] or randomly by [`interactive::challenge`]
pub type Challenge = Integer;

/// The ZK proof. Computed by [`interactive::prove`] or
/// [`non_interactive::prove`]
#[derive(Debug, Clone)]
#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
pub struct Proof {
    pub z1: Integer,
    pub z2: Integer,
    pub z3: Integer,
}

/// The interactive version of the ZK proof. Should be completed in 3 rounds:
/// prover commits to data, verifier responds with a random challenge, and
/// prover gives proof with commitment and challenge.
pub mod interactive {
    use generic_ec::Curve;
    use rand_core::RngCore;
    use rug::{Complete, Integer};

    use crate::common::{fail_if, fail_if_ne, IntegerExt, InvalidProofReason};
    use crate::{Error, InvalidProof};

    use super::{
        Aux, Challenge, Commitment, Data, PrivateCommitment, PrivateData, Proof, SecurityParams,
    };

    /// Create random commitment
    pub fn commit<C: Curve, R: RngCore>(
        aux: &Aux,
        data: Data<C>,
        pdata: PrivateData,
        security: &SecurityParams,
        mut rng: R,
    ) -> Result<(Commitment<C>, PrivateCommitment), Error> {
        let two_to_l_e = (Integer::ONE << (security.l + security.epsilon)).complete();
        let hat_n_at_two_to_l = &aux.rsa_modulo * (Integer::ONE << security.l).complete();
        let hat_n_at_two_to_l_e = (&aux.rsa_modulo * &two_to_l_e).complete();

        let alpha = Integer::from_rng_pm(&two_to_l_e, &mut rng);
        let mu = Integer::from_rng_pm(&hat_n_at_two_to_l, &mut rng);
        let r = Integer::gen_invertible(data.key0.n(), &mut rng);
        let gamma = Integer::from_rng_pm(&hat_n_at_two_to_l_e, &mut rng);

        let commitment = Commitment {
            s: aux.combine(pdata.x, &mu)?,
            a: data.key0.encrypt_with(&alpha, &r)?,
            y: data.b * alpha.to_scalar(),
            d: aux.combine(&alpha, &gamma)?,
        };
        let private_commitment = PrivateCommitment {
            alpha,
            mu,
            r,
            gamma,
        };
        Ok((commitment, private_commitment))
    }

    /// Compute proof for given data and prior protocol values
    pub fn prove<C: Curve>(
        data: Data<C>,
        pdata: PrivateData,
        pcomm: &PrivateCommitment,
        challenge: &Challenge,
    ) -> Result<Proof, Error> {
        Ok(Proof {
            z1: (&pcomm.alpha + challenge * pdata.x).complete(),
            z2: data
                .key0
                .n()
                .combine(&pcomm.r, Integer::ONE, pdata.nonce, challenge)?,
            z3: (&pcomm.gamma + challenge * &pcomm.mu).complete(),
        })
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
                .key0
                .encrypt_with(&proof.z1, &proof.z2)
                .map_err(|_| InvalidProofReason::PaillierEnc)?;
            let rhs = {
                let e_at_c = data
                    .key0
                    .omul(challenge, data.c)
                    .map_err(|_| InvalidProofReason::PaillierOp)?;
                data.key0
                    .oadd(&commitment.a, &e_at_c)
                    .map_err(|_| InvalidProofReason::PaillierOp)?
            };
            fail_if_ne(InvalidProofReason::EqualityCheck(1), lhs, rhs)?;
        }
        {
            let lhs = data.b * proof.z1.to_scalar();
            let rhs = commitment.y + data.x * challenge.to_scalar();
            fail_if_ne(InvalidProofReason::EqualityCheck(2), lhs, rhs)?;
        }
        {
            let lhs = aux.combine(&proof.z1, &proof.z3)?;
            let s_to_e = aux.pow_mod(&commitment.s, challenge)?;
            let rhs = (&commitment.d * s_to_e).modulo(&aux.rsa_modulo);
            fail_if_ne(InvalidProofReason::EqualityCheck(3), lhs, rhs)?;
        }
        fail_if(
            InvalidProofReason::RangeCheck(4),
            proof
                .z1
                .is_in_pm(&(Integer::ONE << (security.l + security.epsilon)).complete()),
        )?;

        Ok(())
    }

    /// Generate random challenge
    ///
    /// `data` parameter is used to generate challenge in correct range
    pub fn challenge<R>(security: &SecurityParams, rng: &mut R) -> Integer
    where
        R: RngCore,
    {
        Integer::from_rng_pm(&security.q, rng)
    }
}

/// The non-interactive version of proof. Completed in one round, for example
/// see the documentation of parent module.
pub mod non_interactive {
    use digest::Digest;
    use generic_ec::Curve;
    use rand_core::RngCore;

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
        rng: &mut impl RngCore,
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

    /// Internal function for deriving challenge from protocol values
    /// deterministically
    pub fn challenge<C: Curve, D: Digest>(
        shared_state: &impl udigest::Digestable,
        aux: &Aux,
        data: Data<C>,
        commitment: &Commitment<C>,
        security: &SecurityParams,
    ) -> Challenge {
        let tag = "paillier-zk.group_element_vs_paillier_encryption_in_range.ni-challenge";
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
    use generic_ec::{Curve, Point, Scalar};
    use rug::{Complete, Integer};
    use sha2::Digest;

    use crate::common::test::random_key;
    use crate::common::{IntegerExt, InvalidProofReason};

    fn run<R: rand_core::RngCore + rand_core::CryptoRng, C: Curve, D: Digest>(
        mut rng: R,
        security: super::SecurityParams,
        plaintext: Integer,
    ) -> Result<(), crate::common::InvalidProof> {
        let private_key0 = random_key(&mut rng).unwrap();
        let key0 = private_key0.encryption_key().clone();

        let (ciphertext, nonce) = key0.encrypt_with_random(&mut rng, &plaintext).unwrap();
        let b = Point::<C>::generator() * Scalar::random(&mut rng);
        let x = b * plaintext.to_scalar();

        let data = super::Data {
            key0: &key0,
            c: &ciphertext,
            x: &x,
            b: &b,
        };
        let pdata = super::PrivateData {
            x: &plaintext,
            nonce: &nonce,
        };

        let aux = crate::common::test::aux(&mut rng);

        let shared_state = "shared state";

        let (commitment, proof) = super::non_interactive::prove::<C, D>(
            &shared_state,
            &aux,
            data,
            pdata,
            &security,
            &mut rng,
        )
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
            q: (Integer::ONE << 128_u32).complete(),
        };
        let plaintext = Integer::from_rng_pm(&(Integer::ONE << security.l).complete(), &mut rng);
        run::<_, C, D>(rng, security, plaintext).expect("proof failed");
    }

    fn failing_test<C: Curve, D: Digest>() {
        let rng = rand_dev::DevRng::new();
        let security = super::SecurityParams {
            l: 1024,
            epsilon: 300,
            q: (Integer::ONE << 128_u32).complete(),
        };
        let plaintext = (Integer::ONE << (security.l + security.epsilon + 1)).complete();
        let r = run::<_, C, D>(rng, security, plaintext).expect_err("proof should not pass");
        match r.reason() {
            InvalidProofReason::RangeCheck(_) => (),
            e => panic!("proof should not fail with: {e:?}"),
        }
    }

    #[test]
    fn passing_p256() {
        passing_test::<generic_ec::curves::Secp256r1, sha2::Sha256>()
    }
    #[test]
    fn failing_p256() {
        failing_test::<generic_ec::curves::Secp256r1, sha2::Sha256>()
    }

    #[test]
    fn passing_million() {
        passing_test::<crate::curve::C, sha2::Sha256>()
    }
    #[test]
    fn failing_million() {
        failing_test::<crate::curve::C, sha2::Sha256>()
    }
}

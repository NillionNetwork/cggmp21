//! ZK-proof of paillier encryption in range. Called Пenc or Renc in the CGGMP21
//! paper.
//!
//! ## Description
//!
//! A party P has `key`, `pkey` - public and private keys in paillier
//! cryptosystem. P also has `plaintext`, `nonce`, and
//! `ciphertext = key.encrypt_with(plaintext, nonce)`.
//!
//! P wants to prove that `plaintext` is at most `l` bits, without disclosing
//! it, the `pkey`, and `nonce`

//! ## Example
//!
//! ```
//! use paillier_zk::{paillier_encryption_in_range as p, IntegerExt};
//! use rug::{Integer, Complete};
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
//! let plaintext = Integer::from_rng_pm(&(Integer::ONE << security.l).complete(), &mut rng);
//! let (ciphertext, nonce) = key.encrypt_with_random(&mut rng, &plaintext)?;
//!
//! // 3. Prover computes a non-interactive proof that plaintext is at most 1024 bits:
//!
//! let data = p::Data { key, ciphertext: &ciphertext };
//! let (commitment, proof) = p::non_interactive::prove::<sha2::Sha256>(
//!     &shared_state,
//!     &aux,
//!     data,
//!     p::PrivateData {
//!         plaintext: &plaintext,
//!         nonce: &nonce,
//!     },
//!     &security,
//!     &mut rng,
//! )?;
//!
//! // 4. Prover sends this data to verifier
//!
//! # fn send(_: &p::Data, _: &p::Commitment, _: &p::Proof) {  }
//! send(&data, &commitment, &proof);
//!
//! // 5. Verifier receives the data and the proof and verifies it
//!
//! # let recv = || (data, commitment, proof);
//! let (data, commitment, proof) = recv();
//! p::non_interactive::verify::<sha2::Sha256>(
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
use rug::Integer;

#[cfg(feature = "serde")]
use serde::{Deserialize, Serialize};

pub use crate::common::Aux;
pub use crate::common::InvalidProof;

/// Security parameters for proof. Choosing the values is a tradeoff between
/// speed and chance of rejecting a valid proof or accepting an invalid proof
#[derive(Debug, Clone)]
#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
pub struct SecurityParams {
    /// l in paper, security parameter for bit size of plaintext: it needs to
    /// be in range [-2^l; 2^l] or equivalently 2^l
    pub l: usize,
    /// Epsilon in paper, slackness parameter
    pub epsilon: usize,
    /// q in paper. Security parameter for challenge
    pub q: Integer,
}

/// Public data that both parties know
#[derive(Debug, Clone, Copy, udigest::Digestable)]
pub struct Data<'a> {
    /// N0 in paper, public key that k -> K was encrypted on
    #[udigest(as = crate::common::encoding::AnyEncryptionKey)]
    pub key: &'a dyn AnyEncryptionKey,
    /// K in paper
    #[udigest(as = &crate::common::encoding::Integer)]
    pub ciphertext: &'a Ciphertext,
}

/// Private data of prover
#[derive(Clone, Copy)]
pub struct PrivateData<'a> {
    /// k in paper, plaintext of K
    pub plaintext: &'a Integer,
    /// rho in paper, nonce of encryption k -> K
    pub nonce: &'a Nonce,
}

// As described in cggmp21 at page 33
/// Prover's first message, obtained by [`interactive::commit`]
#[derive(Debug, Clone, udigest::Digestable)]
#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
pub struct Commitment {
    #[udigest(as = crate::common::encoding::Integer)]
    pub s: Integer,
    #[udigest(as = crate::common::encoding::Integer)]
    pub a: Integer,
    #[udigest(as = crate::common::encoding::Integer)]
    pub c: Integer,
}

/// Prover's data accompanying the commitment. Kept as state between rounds in
/// the interactive protocol.
#[derive(Clone)]
pub struct PrivateCommitment {
    pub alpha: Integer,
    pub mu: Integer,
    pub r: Integer,
    pub gamma: Integer,
}

/// Verifier's challenge to prover. Can be obtained deterministically by
/// [`non_interactive::challenge`] or randomly by [`interactive::challenge`]
pub type Challenge = Integer;

// As described in cggmp21 at page 33
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
    pub fn commit<R: RngCore>(
        aux: &Aux,
        data: Data,
        pdata: PrivateData,
        security: &SecurityParams,
        rng: &mut R,
    ) -> Result<(Commitment, PrivateCommitment), Error> {
        let two_to_l_plus_e = (Integer::ONE << (security.l + security.epsilon)).complete();
        let hat_n_at_two_to_l = (Integer::ONE << security.l).complete() * &aux.rsa_modulo;
        let hat_n_at_two_to_l_plus_e =
            (Integer::ONE << (security.l + security.epsilon)).complete() * &aux.rsa_modulo;

        let alpha = Integer::from_rng_pm(&two_to_l_plus_e, rng);
        let mu = Integer::from_rng_pm(&hat_n_at_two_to_l, rng);
        let r = Integer::gen_invertible(data.key.n(), rng);
        let gamma = Integer::from_rng_pm(&hat_n_at_two_to_l_plus_e, rng);

        let s = aux.combine(pdata.plaintext, &mu)?;
        let a = data.key.encrypt_with(&alpha, &r)?;
        let c = aux.combine(&alpha, &gamma)?;

        Ok((
            Commitment { s, a, c },
            PrivateCommitment {
                alpha,
                mu,
                r,
                gamma,
            },
        ))
    }

    /// Compute proof for given data and prior protocol values
    pub fn prove(
        data: Data,
        pdata: PrivateData,
        private_commitment: &PrivateCommitment,
        challenge: &Challenge,
    ) -> Result<Proof, Error> {
        let z1 = (&private_commitment.alpha + (challenge * pdata.plaintext)).complete();
        let nonce_to_challenge_mod_n: Integer = pdata
            .nonce
            .pow_mod_ref(challenge, data.key.n())
            .ok_or(BadExponent::undefined())?
            .into();
        let z2 = (&private_commitment.r * nonce_to_challenge_mod_n).modulo(data.key.n());
        let z3 = (&private_commitment.gamma + (challenge * &private_commitment.mu)).complete();
        Ok(Proof { z1, z2, z3 })
    }

    /// Verify the proof
    pub fn verify(
        aux: &Aux,
        data: Data,
        commitment: &Commitment,
        security: &SecurityParams,
        challenge: &Challenge,
        proof: &Proof,
    ) -> Result<(), InvalidProof> {
        {
            fail_if_ne(
                InvalidProofReason::EqualityCheck(1),
                &data.ciphertext.gcd_ref(data.key.n()).complete(),
                Integer::ONE,
            )?;
        }
        {
            let lhs = data
                .key
                .encrypt_with(&proof.z1, &proof.z2)
                .map_err(|_| InvalidProofReason::PaillierEnc)?;
            let rhs = {
                let e_at_k = data
                    .key
                    .omul(challenge, data.ciphertext)
                    .map_err(|_| InvalidProofReason::PaillierOp)?;
                data.key
                    .oadd(&commitment.a, &e_at_k)
                    .map_err(|_| InvalidProofReason::PaillierOp)?
            };
            fail_if_ne(InvalidProofReason::EqualityCheck(2), lhs, rhs)?;
        }

        {
            let lhs = aux.combine(&proof.z1, &proof.z3)?;
            let s_to_e = aux.pow_mod(&commitment.s, challenge)?;
            let rhs = (&commitment.c * s_to_e).modulo(&aux.rsa_modulo);
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
    /// `security` parameter is used to generate challenge in correct range
    pub fn challenge<R: RngCore>(security: &SecurityParams, rng: &mut R) -> Challenge {
        Integer::from_rng_pm(&security.q, rng)
    }
}

/// The non-interactive version of proof. Completed in one round, for example
/// see the documentation of parent module.
pub mod non_interactive {
    use digest::Digest;

    use crate::{Error, InvalidProof};

    use super::{Aux, Challenge, Commitment, Data, PrivateData, Proof, SecurityParams};

    /// Compute proof for the given data, producing random commitment and
    /// deriving determenistic challenge.
    ///
    /// Obtained from the above interactive proof via Fiat-Shamir heuristic.
    pub fn prove<D: Digest>(
        shared_state: &impl udigest::Digestable,
        aux: &Aux,
        data: Data,
        pdata: PrivateData,
        security: &SecurityParams,
        rng: &mut impl rand_core::RngCore,
    ) -> Result<(Commitment, Proof), Error> {
        let (comm, pcomm) = super::interactive::commit(aux, data, pdata, security, rng)?;
        let challenge = challenge::<D>(shared_state, aux, data, &comm, security);
        let proof = super::interactive::prove(data, pdata, &pcomm, &challenge)?;
        Ok((comm, proof))
    }

    /// Deterministically compute challenge based on prior known values in protocol
    pub fn challenge<D: Digest>(
        shared_state: &impl udigest::Digestable,
        aux: &Aux,
        data: Data,
        commitment: &Commitment,
        security: &SecurityParams,
    ) -> Challenge {
        let tag = "paillier_zk.encryption_in_range.ni_challenge";
        let seed = udigest::inline_struct!(tag {
            shared_state,
            aux: aux.digest_public_data(),
            data,
            commitment,
        });
        let mut rng = rand_hash::HashRng::<D, _>::from_seed(seed);
        super::interactive::challenge(security, &mut rng)
    }

    /// Verify the proof, deriving challenge independently from same data
    pub fn verify<D: Digest>(
        shared_state: &impl udigest::Digestable,
        aux: &Aux,
        data: Data,
        commitment: &Commitment,
        security: &SecurityParams,
        proof: &Proof,
    ) -> Result<(), InvalidProof> {
        let challenge = challenge::<D>(shared_state, aux, data, commitment, security);
        super::interactive::verify(aux, data, commitment, security, &challenge, proof)
    }
}

#[cfg(test)]
mod test {
    use rug::{Complete, Integer};
    use sha2::Digest;

    use crate::common::{IntegerExt, InvalidProofReason};

    fn run_with<D: Digest>(
        mut rng: &mut impl rand_core::CryptoRngCore,
        security: super::SecurityParams,
        plaintext: Integer,
    ) -> Result<(), crate::common::InvalidProof> {
        let aux = crate::common::test::aux(&mut rng);
        let private_key = crate::common::test::random_key(&mut rng).unwrap();
        let key = private_key.encryption_key();
        let (ciphertext, nonce) = key.encrypt_with_random(&mut rng, &plaintext).unwrap();
        let data = super::Data {
            key,
            ciphertext: &ciphertext,
        };
        let pdata = super::PrivateData {
            plaintext: &plaintext,
            nonce: &nonce,
        };

        let shared_state = "shared state";
        let (commitment, proof) =
            super::non_interactive::prove::<D>(&shared_state, &aux, data, pdata, &security, rng)
                .unwrap();
        super::non_interactive::verify::<D>(
            &shared_state,
            &aux,
            data,
            &commitment,
            &security,
            &proof,
        )
    }

    #[test]
    fn passing() {
        let mut rng = rand_dev::DevRng::new();
        let security = super::SecurityParams {
            l: 1024,
            epsilon: 256,
            q: (Integer::ONE << 128_u32).complete() - 1,
        };
        let plaintext = Integer::from_rng_pm(&(Integer::ONE << security.l).complete(), &mut rng);
        let r = run_with::<sha2::Sha256>(&mut rng, security, plaintext);
        match r {
            Ok(()) => (),
            Err(e) => panic!("{e:?}"),
        }
    }
    #[test]
    fn failing() {
        let mut rng = rand_dev::DevRng::new();
        let security = super::SecurityParams {
            l: 1024,
            epsilon: 256,
            q: (Integer::ONE << 128_u32).complete() - 1,
        };
        let plaintext = (Integer::ONE << (security.l + security.epsilon)).complete() + 1;
        let r = run_with::<sha2::Sha256>(&mut rng, security, plaintext);
        match r.map_err(|e| e.reason()) {
            Ok(()) => panic!("proof should not pass"),
            Err(InvalidProofReason::RangeCheck(_)) => (),
            Err(e) => panic!("proof should not fail with {e:?}"),
        }
    }
}

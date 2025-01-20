//! ZK-proof of discrete log with El-Gamal commitment.
//! Called Пelog or Relog in the CGGMP21/CGGMP24 papers.
//!
//! //! ## Description
//!
//! A party P has `L = g ^ lambda`, `M = (g ^ y) * (X ^ lambda)`, and `Y = h ^ y`,
//! with g being a generator of curve `E`, h is a point on the curve
//! and X is a public key (and a point on the curve).
//! P shares L, M, Y, X, and h with V and wants to prove that the
//! logarithm base h of Y is the discrete logarithm base g of the El-Gamal
//! plaintext associated with the ciphertext (L,M) and public key X.
//!
//! Given:
//! - Curve `E`
//! - `X` - public key, point on the curve
//! - `L = g * lambda`, `M = (g ^ y) * (X ^ lambda)`, and `Y = h ^ y` - data to obtain proof about
//!
//! Prove:
//! - `logarithm base h of Y= y`
//!
//! Disclosing only: `g`, `L`, `M`, `X`, `Y`, `h`
//!
//! ## Example
//!
//! ```rust
//! use paillier_zk::{dlog_with_el_gamal_commitment as p, IntegerExt};
//! use rug::{Integer, Complete};
//! use generic_ec::{Point, curves::Secp256k1 as E, Scalar};
//! # mod pregenerated {
//! #     use super::*;
//! #     paillier_zk::load_pregenerated_data!(
//! #         verifier_aux: p::Aux,
//! #     );
//! # }
//!
//! # fn main() -> Result<(), Box<dyn std::error::Error>> {
//! // Prover and verifier have a shared protocol state
//! let shared_state = "some shared state";
//!
//! let mut rng = rand_core::OsRng;
//! # let mut rng = rand_dev::DevRng::new();
//!
//! // 0. Setup: prover and verifier share common Ring-Pedersen parameters:
//!
//! let aux: p::Aux = pregenerated::verifier_aux();
//! let security = p::SecurityParams {
//!     q: (Integer::ONE << 128_u32).complete(),
//! };
//!
//! // 1. Setup: prover prepares the public key X
//!
//! // X in paper is a point on the Curve E
//! let x = Point::<E>::generator() * Scalar::random(&mut rng);
//!
//! // h in paper is a point on the Curve E
//! let h = Point::<E>::generator() * Scalar::random(&mut rng);
//!
//! // 2. Setup: prover prepares all plaintexts
//!
//! // y in paper
//! let plaintext_y = Integer::from_rng_pm(&Integer::curve_order::<E>(), &mut rng);
//! // lambda in paper
//! let plaintext_lambda = Integer::from_rng_pm(&Integer::curve_order::<E>(), &mut rng);
//!
//! // 3. Setup: prover encrypts everything on correct keys
//!
//! // L in paper
//! let ciphertext_l = Point::<E>::generator() * plaintext_lambda.to_scalar();
//! // M in paper
//! let ciphertext_m = Point::<E>::generator() * plaintext_y.to_scalar() + x * plaintext_lambda.to_scalar();
//! // Y in paper
//! let ciphertext_h_to_y = h * plaintext_y.to_scalar();
//!
//! // 4. Prover computes a non-interactive proof that logarithm base h of Y
//! //    and lambda are the same
//!
//! let data = p::Data {
//!     l: &ciphertext_l,
//!     m: &ciphertext_m,
//!     x: &x,
//!     h_to_y: &ciphertext_h_to_y,
//!     h: &h,
//! };
//! let pdata = p::PrivateData {
//!     y: &plaintext_y,
//!     lambda: &plaintext_lambda,
//! };
//! let (commitment, proof) =
//!     p::non_interactive::prove::<E, sha2::Sha256>(
//!         &shared_state,
//!         &aux,
//!         data,
//!         pdata,
//!         &security,
//!         &mut rng,
//!     )?;
//!
//! // 5. Prover sends this data to verifier
//!
//! # use generic_ec::Curve;
//! # fn send<E: Curve>(_: &p::Data<E>, _: &p::Commitment<E>, _: &p::Proof) {  }
//! send(&data, &commitment, &proof);
//!
//! // 6. Verifier receives the data and the proof and verifies it
//!
//! # let recv = || (data, commitment, proof);
//! let (data, commitment, proof) = recv();
//! let r = p::non_interactive::verify::<E, sha2::Sha256>(
//!     &shared_state,
//!     &aux,
//!     data,
//!     &commitment,
//!     &security,
//!     &proof,
//! )?;
//! #
//! # Ok(()) }
//! ```
//!
//! If the verification succeeded, verifier can continue communication with prover

use generic_ec::{Curve, Point};
use rug::Integer;

#[cfg(feature = "serde")]
use serde::{Deserialize, Serialize};

pub use crate::common::{Aux, InvalidProof};

/// Security parameters for proof. Choosing the values is a tradeoff between
/// speed and chance of rejecting a valid proof or accepting an invalid proof
#[derive(Debug, Clone, udigest::Digestable)]
#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
pub struct SecurityParams {
    /// q in paper. Security parameter for challenge
    #[udigest(as = crate::common::encoding::Integer)]
    pub q: Integer,
}

/// Public data that both parties know
#[derive(Debug, Clone, Copy, udigest::Digestable)]
#[udigest(bound = "")]
pub struct Data<'a, C: Curve> {
    /// L in paper, obtained as g^\lambda
    pub l: &'a Point<C>,
    /// M in paper, obtained as g^y X^\lambda
    pub m: &'a Point<C>,
    /// X in paper
    pub x: &'a Point<C>,
    /// Y in paper, obtained as h^y
    pub h_to_y: &'a Point<C>,
    /// h in paper
    pub h: &'a Point<C>,
}

/// Private data of prover
#[derive(Clone, Copy)]
pub struct PrivateData<'a> {
    /// y or epsilon in paper, log of Y base h
    pub y: &'a Integer,
    /// lambda in paper, preimage of L
    pub lambda: &'a Integer,
}

// As described in cggmp24 at page 57
/// Prover's first message, obtained by [`interactive::commit`]
#[derive(Debug, Clone, udigest::Digestable)]
#[udigest(bound = "")]
#[cfg_attr(feature = "serde", derive(Serialize, Deserialize), serde(bound = ""))]
pub struct Commitment<C: Curve> {
    pub a: Point<C>,
    pub cap_n: Point<C>,
    pub b: Point<C>,
}

/// Prover's data accompanying the commitment. Kept as state between rounds in
/// the interactive protocol.
#[derive(Clone)]
pub struct PrivateCommitment {
    pub alpha: Integer,
    pub m: Integer,
}

/// Verifier's challenge to prover. Can be obtained deterministically by
/// [`non_interactive::challenge`] or randomly by [`interactive::challenge`]
pub type Challenge = Integer;

/// The ZK proof. Computed by [`interactive::prove`] or
/// [`non_interactive::prove`]
#[derive(Debug, Clone)]
#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
pub struct Proof {
    pub z: Integer,
    pub u: Integer,
}

/// The interactive version of the ZK proof. Should be completed in 3 rounds:
/// prover commits to data, verifier responds with a random challenge, and
/// prover gives proof with commitment and challenge.
pub mod interactive {
    use generic_ec::{Curve, Point};
    use rand_core::RngCore;
    use rug::{Complete, Integer};

    use crate::common::{fail_if_ne, IntegerExt, InvalidProof, InvalidProofReason};
    use crate::Error;

    use super::*;

    /// Create random commitment
    pub fn commit<C: Curve, R: RngCore>(
        data: Data<C>,
        mut rng: R,
    ) -> Result<(Commitment<C>, PrivateCommitment), Error> {
        let alpha = Integer::gen_invertible(&Integer::curve_order::<C>(), &mut rng);
        let m = Integer::gen_invertible(&Integer::curve_order::<C>(), &mut rng);

        let a = Point::<C>::generator() * alpha.to_scalar();
        let enne = Point::<C>::generator() * m.to_scalar() + data.x * alpha.to_scalar();
        let b = data.h * m.to_scalar();

        let commitment = Commitment { a, cap_n: enne, b };
        let private_commitment = PrivateCommitment { alpha, m };
        Ok((commitment, private_commitment))
    }

    /// Compute proof for given data and prior protocol values
    pub fn prove<C: Curve>(
        pdata: PrivateData,
        pcomm: &PrivateCommitment,
        challenge: &Challenge,
    ) -> Result<Proof, Error> {
        let z = ((&pcomm.alpha + challenge * pdata.lambda).complete())
            .modulo(&Integer::curve_order::<C>());
        let u = ((&pcomm.m + challenge * pdata.y).complete()).modulo(&Integer::curve_order::<C>());
        Ok(Proof { z, u })
    }

    /// Verify the proof
    pub fn verify<C: Curve>(
        data: Data<C>,
        commitment: &Commitment<C>,
        challenge: &Challenge,
        proof: &Proof,
    ) -> Result<(), InvalidProof> {
        // Three equality checks
        {
            let lhs = Point::<C>::generator() * proof.z.to_scalar();
            let rhs = commitment.a + data.l * challenge.to_scalar();
            fail_if_ne(InvalidProofReason::EqualityCheck(1), lhs, rhs)?;
        }
        {
            let lhs = Point::<C>::generator() * proof.u.to_scalar() + data.x * proof.z.to_scalar();
            let rhs = commitment.cap_n + data.m * challenge.to_scalar();
            fail_if_ne(InvalidProofReason::EqualityCheck(2), lhs, rhs)?;
        }
        {
            let lhs = data.h * proof.u.to_scalar();
            let rhs = commitment.b + data.h_to_y * challenge.to_scalar();
            fail_if_ne(InvalidProofReason::EqualityCheck(3), lhs, rhs)?;
        }

        Ok(())
    }

    /// Generate random challenge
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
        let (comm, pcomm) = super::interactive::commit(data, rng)?;
        let challenge = challenge::<C, D>(shared_state, aux, data, &comm, security);
        let proof = super::interactive::prove::<C>(pdata, &pcomm, &challenge)?;
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
        super::interactive::verify::<C>(data, commitment, &challenge, proof)
    }

    /// Deterministically compute challenge based on prior known values in protocol
    pub fn challenge<C: Curve, D: Digest>(
        shared_state: &impl udigest::Digestable,
        aux: &Aux,
        data: Data<C>,
        commitment: &Commitment<C>,
        security: &SecurityParams,
    ) -> Challenge {
        let tag = "paillier_zk.dlog_with_el_gamal.ni_challenge";
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
    use rug::Integer;
    use sha2::Digest;

    use crate::common::{IntegerExt, InvalidProofReason};

    fn run<R: rand_core::RngCore + rand_core::CryptoRng, C: Curve, D: Digest>(
        rng: &mut R,
        security: super::SecurityParams,
        data: super::Data<C>,
        pdata: super::PrivateData,
    ) -> Result<(), crate::common::InvalidProof> {
        let aux = crate::common::test::aux(rng);

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
            q: (Integer::ONE << 128_u32).into(),
        };
        let y = Integer::from_rng_pm(&Integer::curve_order::<C>(), &mut rng);
        let lambda = Integer::from_rng_pm(&Integer::curve_order::<C>(), &mut rng);
        let x = Point::<C>::generator() * Scalar::random(&mut rng);
        let h = Point::<C>::generator() * Scalar::random(&mut rng);

        let l = Point::<C>::generator() * lambda.to_scalar();
        let m = Point::<C>::generator() * y.to_scalar() + x * lambda.to_scalar();
        let h_to_y = h * y.to_scalar();

        let data = super::Data {
            l: &l,
            m: &m,
            x: &x,
            h_to_y: &h_to_y,
            h: &h,
        };
        let pdata = super::PrivateData {
            y: &y,
            lambda: &lambda,
        };
        run::<_, C, D>(&mut rng, security, data, pdata).expect("proof failed");
    }

    fn failing_check_lambda_<C: Curve, D: Digest>() {
        // Scenario where the prover P does not know lambda
        let mut rng = rand_dev::DevRng::new();
        let security = super::SecurityParams {
            q: (Integer::ONE << 128_u32).into(),
        };
        let y = Integer::from_rng_pm(&Integer::curve_order::<C>(), &mut rng);
        let lambda = Integer::from_rng_pm(&Integer::curve_order::<C>(), &mut rng);
        let false_lambda = Integer::from_rng_pm(&Integer::curve_order::<C>(), &mut rng);

        let x = Point::<C>::generator() * Scalar::random(&mut rng);
        let h = Point::<C>::generator() * Scalar::random(&mut rng);

        let l = Point::<C>::generator() * lambda.to_scalar();
        let m = Point::<C>::generator() * y.to_scalar() + x * lambda.to_scalar();
        let h_to_y = h * y.to_scalar();

        let data = super::Data {
            l: &l,
            m: &m,
            x: &x,
            h_to_y: &h_to_y,
            h: &h,
        };
        let pdata = super::PrivateData {
            y: &y,
            lambda: &false_lambda,
        };
        let r = run::<_, C, D>(&mut rng, security, data, pdata).expect_err("proof should not pass");
        match r.reason() {
            InvalidProofReason::EqualityCheck(1) => (),
            e => panic!("proof should not fail with {e:?}"),
        }
    }

    fn failing_check_y_<C: Curve, D: Digest>() {
        // Scenario where the prover P does not know y
        let mut rng = rand_dev::DevRng::new();
        let security = super::SecurityParams {
            q: (Integer::ONE << 128_u32).into(),
        };
        let y = Integer::from_rng_pm(&Integer::curve_order::<C>(), &mut rng);
        let lambda = Integer::from_rng_pm(&Integer::curve_order::<C>(), &mut rng);
        let false_y = Integer::from_rng_pm(&Integer::curve_order::<C>(), &mut rng);

        let x = Point::<C>::generator() * Scalar::random(&mut rng);
        let h = Point::<C>::generator() * Scalar::random(&mut rng);

        let l = Point::<C>::generator() * lambda.to_scalar();
        let m = Point::<C>::generator() * y.to_scalar() + x * lambda.to_scalar();
        let h_to_y = h * y.to_scalar();

        let data = super::Data {
            l: &l,
            m: &m,
            x: &x,
            h_to_y: &h_to_y,
            h: &h,
        };
        let pdata = super::PrivateData {
            y: &false_y,
            lambda: &lambda,
        };
        let r = run::<_, C, D>(&mut rng, security, data, pdata).expect_err("proof should not pass");
        match r.reason() {
            InvalidProofReason::EqualityCheck(2) => (),
            e => panic!("proof should not fail with {e:?}"),
        }
    }

    #[test]
    fn passing_p256() {
        passing_test::<generic_ec::curves::Secp256r1, sha2::Sha256>()
    }

    #[test]
    fn passing_million() {
        passing_test::<crate::curve::C, sha2::Sha256>()
    }
    #[test]
    fn failing_check_1_p256() {
        failing_check_lambda_::<generic_ec::curves::Secp256r1, sha2::Sha256>()
    }

    #[test]
    fn failing_check_1_million() {
        failing_check_lambda_::<crate::curve::C, sha2::Sha256>()
    }
    #[test]
    fn failing_check_2_p256() {
        failing_check_y_::<generic_ec::curves::Secp256r1, sha2::Sha256>()
    }

    #[test]
    fn failing_check_2_million() {
        failing_check_y_::<crate::curve::C, sha2::Sha256>()
    }
}

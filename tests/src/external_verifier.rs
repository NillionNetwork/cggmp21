use cggmp21::signing::Signature;
use generic_ec::{Curve, Point};

/// Verifies signature produced by cggmp21 implemenation using external library
pub trait ExternalVerifier<E: Curve> {
    fn verify(
        public_key: &Point<E>,
        signature: &Signature<E>,
        message: &[u8],
    ) -> anyhow::Result<()>;
}

/// Doesn't do any external verification
pub struct Noop;

impl<E: Curve> ExternalVerifier<E> for Noop {
    fn verify(
        _public_key: &Point<E>,
        _signature: &Signature<E>,
        _message: &[u8],
    ) -> anyhow::Result<()> {
        Ok(())
    }
}

pub mod blockchains {
    use anyhow::Context;
    use cggmp21::supported_curves::Secp256k1;

    use crate::external_verifier::ExternalVerifier;

    /// Verifies ECDSA signature using the same library as used in Bitcoin
    pub struct Bitcoin;

    impl ExternalVerifier<Secp256k1> for Bitcoin {
        fn verify(
            public_key: &generic_ec::Point<Secp256k1>,
            signature: &cggmp21::signing::Signature<Secp256k1>,
            message: &[u8],
        ) -> anyhow::Result<()> {
            let public_key = secp256k1::PublicKey::from_slice(&public_key.to_bytes(true))
                .context("public key is not valid")?;
            let message =
                secp256k1::Message::from_hashed_data::<secp256k1::hashes::sha256::Hash>(message);

            let mut signature_bytes = [0u8; 64];
            signature.write_to_slice(&mut signature_bytes);
            let signature = secp256k1::ecdsa::Signature::from_compact(&signature_bytes)
                .context("malformed signature")?;

            signature
                .verify(&message, &public_key)
                .context("invalid siganture")
        }
    }
}

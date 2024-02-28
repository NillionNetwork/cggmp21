use digest::Digest;
use futures::SinkExt;
use generic_ec::{Curve, NonZero, Point, Scalar, SecretScalar};
use generic_ec_zkp::schnorr_pok;
use rand_core::{CryptoRng, RngCore};
use round_based::{
    rounds_router::simple_store::RoundInput, rounds_router::RoundsRouter, Delivery, Mpc, MpcParty,
    Outgoing, ProtocolMessage,
};
use serde::{Deserialize, Serialize};

use crate::progress::Tracer;
use crate::{
    errors::IoError,
    key_share::{CoreKeyShare, DirtyCoreKeyShare, DirtyKeyInfo, Validate},
    security_level::SecurityLevel,
    utils, ExecutionId,
};

use super::{Bug, KeygenAborted, KeygenError};

/// Message of key generation protocol
#[derive(ProtocolMessage, Clone, Serialize, Deserialize)]
#[serde(bound = "")]
pub enum Msg<E: Curve, L: SecurityLevel, D: Digest> {
    /// Round 1 message
    Round1(MsgRound1<D>),
    /// Reliability check message (optional additional round)
    ReliabilityCheck(MsgReliabilityCheck<D>),
    /// Round 2 message
    Round2(MsgRound2<E, L>),
    /// Round 3 message
    Round3(MsgRound3<E>),
}

/// Message from round 1
#[derive(Clone, Serialize, Deserialize, udigest::Digestable)]
#[serde(bound = "")]
#[udigest(bound = "")]
#[udigest(tag = "dfns.cggmp21.keygen.non_threshold.round1")]
pub struct MsgRound1<D: Digest> {
    /// $V_i$
    #[udigest(as_bytes)]
    pub commitment: digest::Output<D>,
}
/// Message from round 2
#[serde_with::serde_as]
#[derive(Clone, Serialize, Deserialize, udigest::Digestable)]
#[serde(bound = "")]
#[udigest(bound = "")]
#[udigest(tag = "dfns.cggmp21.keygen.non_threshold.round2")]
pub struct MsgRound2<E: Curve, L: SecurityLevel> {
    /// `rid_i`
    #[serde_as(as = "utils::HexOrBin")]
    #[udigest(as_bytes)]
    pub rid: L::Rid,
    /// $X_i$
    pub X: NonZero<Point<E>>,
    /// $A_i$
    pub sch_commit: schnorr_pok::Commit<E>,
    /// Party contribution to chain code
    #[cfg(feature = "hd-wallets")]
    #[serde_as(as = "Option<utils::HexOrBin>")]
    #[udigest(with = utils::encoding::maybe_bytes)]
    pub chain_code: Option<slip_10::ChainCode>,
    /// $u_i$
    #[serde(with = "hex::serde")]
    #[udigest(as_bytes)]
    pub decommit: L::Rid,
}
/// Message from round 3
#[derive(Clone, Serialize, Deserialize)]
#[serde(bound = "")]
pub struct MsgRound3<E: Curve> {
    /// $\psi_i$
    pub sch_proof: schnorr_pok::Proof<E>,
}
/// Message parties exchange to ensure reliability of broadcast channel
#[derive(Clone, Serialize, Deserialize)]
#[serde(bound = "")]
pub struct MsgReliabilityCheck<D: Digest>(pub digest::Output<D>);

#[derive(udigest::Digestable)]
#[udigest(tag = "dfns.cggmp21.keygen.non_threshold.tag")]
enum Tag<'a> {
    /// Tag that includes the prover index
    Indexed {
        party_index: u16,
        #[udigest(as_bytes)]
        sid: &'a [u8],
    },
    /// Tag w/o party index
    Unindexed {
        #[udigest(as_bytes)]
        sid: &'a [u8],
    },
}

pub async fn run_keygen<E, R, M, L, D>(
    mut tracer: Option<&mut dyn Tracer>,
    i: u16,
    n: u16,
    reliable_broadcast_enforced: bool,
    execution_id: ExecutionId<'_>,
    rng: &mut R,
    party: M,
    #[cfg(feature = "hd-wallets")] hd_enabled: bool,
) -> Result<CoreKeyShare<E>, KeygenError>
where
    E: Curve,
    L: SecurityLevel,
    D: Digest + Clone + 'static,
    R: RngCore + CryptoRng,
    M: Mpc<ProtocolMessage = Msg<E, L, D>>,
{
    tracer.protocol_begins();

    tracer.stage("Setup networking");
    let MpcParty { delivery, .. } = party.into_party();
    let (incomings, mut outgoings) = delivery.split();

    let mut rounds = RoundsRouter::<Msg<E, L, D>>::builder();
    let round1 = rounds.add_round(RoundInput::<MsgRound1<D>>::broadcast(i, n));
    let round1_sync = rounds.add_round(RoundInput::<MsgReliabilityCheck<D>>::broadcast(i, n));
    let round2 = rounds.add_round(RoundInput::<MsgRound2<E, L>>::broadcast(i, n));
    let round3 = rounds.add_round(RoundInput::<MsgRound3<E>>::broadcast(i, n));
    let mut rounds = rounds.listen(incomings);

    // Round 1
    tracer.round_begins();

    tracer.stage("Compute execution id");
    let sid = execution_id.as_bytes();
    let tag = |j| {
        udigest::Tag::<D>::new_structured(Tag::Indexed {
            party_index: j,
            sid,
        })
    };
    let tag_i = tag(i);

    tracer.stage("Sample x_i, rid_i, chain_code");
    let x_i = NonZero::<SecretScalar<E>>::random(rng);
    let X_i = Point::generator() * &x_i;

    let mut rid = L::Rid::default();
    rng.fill_bytes(rid.as_mut());

    #[cfg(feature = "hd-wallets")]
    let chain_code_local = if hd_enabled {
        let mut chain_code = slip_10::ChainCode::default();
        rng.fill_bytes(&mut chain_code);
        Some(chain_code)
    } else {
        None
    };

    tracer.stage("Sample schnorr commitment");
    let (sch_secret, sch_commit) = schnorr_pok::prover_commits_ephemeral_secret::<E, _>(rng);

    tracer.stage("Commit to public data");
    let my_decommitment = MsgRound2 {
        rid,
        X: X_i,
        sch_commit,
        #[cfg(feature = "hd-wallets")]
        chain_code: chain_code_local,
        decommit: {
            let mut nonce = L::Rid::default();
            rng.fill_bytes(nonce.as_mut());
            nonce
        },
    };
    let hash_commit = tag_i.clone().digest(&my_decommitment);
    let my_commitment = MsgRound1 {
        commitment: hash_commit,
    };

    tracer.send_msg();
    outgoings
        .send(Outgoing::broadcast(Msg::Round1(my_commitment.clone())))
        .await
        .map_err(IoError::send_message)?;
    tracer.msg_sent();

    // Round 2
    tracer.round_begins();

    tracer.receive_msgs();
    let commitments = rounds
        .complete(round1)
        .await
        .map_err(IoError::receive_message)?;
    tracer.msgs_received();

    // Optional reliability check
    if reliable_broadcast_enforced {
        tracer.stage("Hash received msgs (reliability check)");
        let h_i = udigest::Tag::<D>::new_structured(Tag::Unindexed { sid })
            .digest_iter(commitments.iter_including_me(&my_commitment));

        tracer.send_msg();
        outgoings
            .send(Outgoing::broadcast(Msg::ReliabilityCheck(
                MsgReliabilityCheck(h_i.clone()),
            )))
            .await
            .map_err(IoError::send_message)?;
        tracer.msg_sent();

        tracer.round_begins();

        tracer.receive_msgs();
        let round1_hashes = rounds
            .complete(round1_sync)
            .await
            .map_err(IoError::receive_message)?;
        tracer.msgs_received();

        tracer.stage("Assert other parties hashed messages (reliability check)");
        let parties_have_different_hashes = round1_hashes
            .into_iter_indexed()
            .filter(|(_j, _msg_id, hash_j)| hash_j.0 != h_i)
            .map(|(j, msg_id, _)| (j, msg_id))
            .collect::<Vec<_>>();
        if !parties_have_different_hashes.is_empty() {
            return Err(KeygenAborted::Round1NotReliable(parties_have_different_hashes).into());
        }
    }

    tracer.send_msg();
    outgoings
        .send(Outgoing::broadcast(Msg::Round2(my_decommitment.clone())))
        .await
        .map_err(IoError::send_message)?;
    tracer.msg_sent();

    // Round 3
    tracer.round_begins();

    tracer.receive_msgs();
    let decommitments = rounds
        .complete(round2)
        .await
        .map_err(IoError::receive_message)?;
    tracer.msgs_received();

    tracer.stage("Validate decommitments");
    let blame = utils::collect_blame(&commitments, &decommitments, |j, com, decom| {
        let com_expected = tag(j).digest(decom);
        com.commitment != com_expected
    });
    if !blame.is_empty() {
        return Err(KeygenAborted::InvalidDecommitment(blame).into());
    }

    #[cfg(feature = "hd-wallets")]
    let chain_code = if hd_enabled {
        tracer.stage("Calculate chain_code");
        let blame = utils::collect_simple_blame(&decommitments, |decom| decom.chain_code.is_none());
        if !blame.is_empty() {
            return Err(KeygenAborted::MissingChainCode(blame).into());
        }
        Some(decommitments.iter_including_me(&my_decommitment).try_fold(
            slip_10::ChainCode::default(),
            |acc, decom| {
                Ok::<_, Bug>(utils::xor_array(
                    acc,
                    decom.chain_code.ok_or(Bug::NoChainCode)?,
                ))
            },
        )?)
    } else {
        None
    };

    tracer.stage("Calculate challege rid");
    let rid = decommitments
        .iter_including_me(&my_decommitment)
        .map(|d| &d.rid)
        .fold(L::Rid::default(), utils::xor_array);
    let challenge = {
        let hash = |d: D| {
            d.chain_update(sid)
                .chain_update(i.to_be_bytes())
                .chain_update(rid.as_ref())
                .finalize()
        };
        let mut rng = crate::rng::HashRng::new(hash);
        Scalar::random(&mut rng)
    };
    let challenge = schnorr_pok::Challenge { nonce: challenge };

    tracer.stage("Prove knowledge of `x_i`");
    let sch_proof = schnorr_pok::prove(&sch_secret, &challenge, &x_i);

    tracer.send_msg();
    let my_sch_proof = MsgRound3 { sch_proof };
    outgoings
        .send(Outgoing::broadcast(Msg::Round3(my_sch_proof.clone())))
        .await
        .map_err(IoError::send_message)?;
    tracer.msg_sent();

    // Round 4
    tracer.round_begins();

    tracer.receive_msgs();
    let sch_proofs = rounds
        .complete(round3)
        .await
        .map_err(IoError::receive_message)?;
    tracer.msgs_received();

    tracer.stage("Validate schnorr proofs");
    let blame = utils::collect_blame(&decommitments, &sch_proofs, |j, decom, sch_proof| {
        let challenge = {
            let hash = |d: D| {
                d.chain_update(sid)
                    .chain_update(j.to_be_bytes())
                    .chain_update(rid.as_ref())
                    .finalize()
            };
            let mut rng = crate::rng::HashRng::new(hash);
            Scalar::random(&mut rng)
        };
        let challenge = schnorr_pok::Challenge { nonce: challenge };
        sch_proof
            .sch_proof
            .verify(&decom.sch_commit, &challenge, &decom.X)
            .is_err()
    });
    if !blame.is_empty() {
        return Err(KeygenAborted::InvalidSchnorrProof(blame).into());
    }

    tracer.protocol_ends();

    Ok(DirtyCoreKeyShare {
        i,
        key_info: DirtyKeyInfo {
            curve: Default::default(),
            shared_public_key: NonZero::from_point(
                decommitments
                    .iter_including_me(&my_decommitment)
                    .map(|d| d.X)
                    .sum(),
            )
            .ok_or(Bug::ZeroPk)?,
            public_shares: decommitments
                .iter_including_me(&my_decommitment)
                .map(|d| d.X)
                .collect(),
            vss_setup: None,
            #[cfg(feature = "hd-wallets")]
            chain_code,
        },
        x: x_i,
    }
    .validate()
    .map_err(|e| Bug::InvalidKeyShare(e.into_error()))?)
}

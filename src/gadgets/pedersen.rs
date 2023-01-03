use bls12_381::{Scalar};
use jubjub::Fr;
use crate::{ConstraintSystem, SynthesisError, Circuit};

use super::boolean::{self, Boolean};
use super::ecc;

use super::constants::{
    VALUE_COMMITMENT_RANDOMNESS_GENERATOR, VALUE_COMMITMENT_VALUE_GENERATOR,
    _VALUE_COMMITMENT_RANDOMNESS_GENERATOR, _VALUE_COMMITMENT_VALUE_GENERATOR,
};
pub struct PedersenCommit {
    pub msg: Fr,
    pub rand: Fr,
}

impl PedersenCommit {
    pub fn commit(&self) -> jubjub::SubgroupPoint {
        (_VALUE_COMMITMENT_VALUE_GENERATOR * self.msg) + 
            (_VALUE_COMMITMENT_RANDOMNESS_GENERATOR * self.rand)
    }
}

pub struct PedersenCommitCircuit {
    pub pedersen_commitment: Option<PedersenCommit>,
}

impl Circuit<Scalar> for PedersenCommitCircuit {
    fn synthesize<CS: ConstraintSystem<Scalar>>(self, cs: &mut CS) -> Result<(), SynthesisError> {
        
        let _value_bits = expose_pedersen_commitment(cs.namespace(|| "Pedersen Commitment"), self.pedersen_commitment)?;
        
        Ok(())
    }
}

fn expose_pedersen_commitment<CS>(
    mut cs: CS,
    pedersen_commitment: Option<PedersenCommit>,
) -> Result<Vec<Boolean>, SynthesisError>
where CS: ConstraintSystem<Scalar>,
{
    let msg_bits = boolean::field_into_boolean_vec_le(
        cs.namespace(||"msg"),
        pedersen_commitment.as_ref().map(|c| c.msg),
    )?;
    
    let msg = ecc::fixed_base_multiplication(
        cs.namespace(|| "compute the msg in the exponent"),
    &VALUE_COMMITMENT_VALUE_GENERATOR,
        &msg_bits
    )?;

    let r = boolean::field_into_boolean_vec_le(
        cs.namespace(|| "r"),
        pedersen_commitment.as_ref().map(|c| c.rand),
    )?;

    let r = ecc::fixed_base_multiplication(
        cs.namespace(|| "r"),
        &VALUE_COMMITMENT_RANDOMNESS_GENERATOR,
    &r,
    )?;

    let res = msg.add(cs.namespace(|| "computation of result"), &r)?;

    res.inputize(cs.namespace(|| "commitment point"))?;

    Ok(msg_bits)
}
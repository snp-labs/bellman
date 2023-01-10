use std::ops::Mul;

use bls12_381::{Bls12, Scalar};
use ff::PrimeField;
use group::Group;
use jubjub::{AffinePoint, ExtendedPoint, Fr, SubgroupPoint};
use crate::{ConstraintSystem, SynthesisError, Circuit};

use super::boolean::{self, Boolean};
use super::ecc::{EdwardsPoint, MontgomeryPoint, fixed_base_multiplication};
use super::mimc7::{self, get_mimc_constants, mimc7, mimc7_cs};
use super::num::{AllocatedNum, Num};

use super::constants::{
    VALUE_COMMITMENT_RANDOMNESS_GENERATOR, VALUE_COMMITMENT_VALUE_GENERATOR,
    _VALUE_COMMITMENT_RANDOMNESS_GENERATOR, _VALUE_COMMITMENT_VALUE_GENERATOR,
};
pub struct Encryption {
    pub msg : Vec<Scalar>,
    pub user_key : SubgroupPoint,
    pub auditor_key : SubgroupPoint,
    pub rand_r: Fr,
    pub rand_k: Fr,
}

pub struct PCT {
    pub keys : Vec<jubjub::SubgroupPoint>,
    pub ct : Vec<Scalar>
}

impl Encryption {
    pub fn encrypt(&self) -> PCT {
        let K = (_VALUE_COMMITMENT_VALUE_GENERATOR * self.rand_k);
        let enc_pk = self.user_key * self.rand_r;
        let enc_apk = self.auditor_key * self.rand_r;

        let affine_K = AffinePoint::from(ExtendedPoint::from(K));

        let mut keys : Vec<jubjub::SubgroupPoint> = Vec::new();

        keys.push(_VALUE_COMMITMENT_VALUE_GENERATOR * self.rand_r);
        keys.push(K + enc_pk);
        keys.push(K + enc_apk);

        let mut enc_msg : Vec<Scalar> = Vec::new();
        let constants = get_mimc_constants().to_vec();
        for (idx, e) in self.msg.iter().enumerate() {
            let bias : Scalar = PrimeField::from_str_vartime(idx.to_string().as_str()).unwrap();
            let input = affine_K.get_u() + bias;
            let res = mimc7(input, input, &constants);
            enc_msg.push(e + res);
        } 

        PCT {
            keys : keys,
            ct : enc_msg
        }
    }
}

pub struct EncryptionCircuit {
    pub encryption: Option<Encryption>,
}

impl Circuit<Scalar> for EncryptionCircuit {
    fn synthesize<CS: ConstraintSystem<Scalar>>(self, cs: &mut CS) -> Result<(), SynthesisError> {
        
        let _value_bits = expose_encryption(cs.namespace(|| "Encryption"), self.encryption)?;
        
        Ok(())
    }
}

fn expose_encryption<CS>(
    mut cs: CS,
    encryption: Option<Encryption>,
    pct: Option<PCT>
) -> Result<(), SynthesisError>
where CS: ConstraintSystem<Scalar>,
{
    let rand_k_bits = boolean::field_into_boolean_vec_le(
        cs.namespace(||"Random k"),
        encryption.as_ref().map(|c| c.rand_k),
    )?;

    let curve_k = fixed_base_multiplication(
        cs.namespace(|| "compute the k in the exponent"),
    &VALUE_COMMITMENT_VALUE_GENERATOR,
        &rand_k_bits
    )?;

    let rand_r_bits = boolean::field_into_boolean_vec_le(
        cs.namespace(|| "Random r"),
        encryption.as_ref().map(|c| c.rand_r),
    )?;

    let enc = encryption.unwrap();
    let affine_pk = AffinePoint::from(ExtendedPoint::from(enc.user_key));
    let raw_pk_u = affine_pk.get_u();
    let raw_pk_v = affine_pk.get_v();

    let alloc_pk_u = AllocatedNum::alloc(cs.namespace(|| "pk u"), || Ok(raw_pk_u)).unwrap();
    let alloc_pk_v = AllocatedNum::alloc(cs.namespace(|| "pk v"), || Ok(raw_pk_v)).unwrap();

    let mon_pk = MontgomeryPoint::interpret_unchecked(alloc_pk_u.into(), alloc_pk_v.into());
    let ed_pk = mon_pk.into_edwards(&mut cs).unwrap();

    ed_pk.mul(cs.namespace(||"compute enc pk"), &rand_r_bits);

    let affine_apk = AffinePoint::from(ExtendedPoint::from(enc.auditor_key));
    let raw_apk_u = affine_apk.get_u();
    let raw_apk_v = affine_apk.get_v();

    let alloc_apk_u = AllocatedNum::alloc(cs.namespace(|| "apk u"), || Ok(raw_apk_u)).unwrap();
    let alloc_apk_v = AllocatedNum::alloc(cs.namespace(|| "apk v"), || Ok(raw_apk_v)).unwrap();

    let mon_apk = MontgomeryPoint::interpret_unchecked(alloc_apk_u.into(), alloc_apk_v.into());
    let ed_apk = mon_apk.into_edwards(&mut cs).unwrap();

    ed_apk.mul(cs.namespace(||"compute enc apk"), &rand_r_bits);

    
    // let constants = get_mimc_constants();
    // let constants = get_mimc_constants().to_vec();
    // for (idx, e) in self.msg.iter().enumerate() {
    //     let bias : Scalar = PrimeField::from_str_vartime(idx.to_string().as_str()).unwrap();
    //     bias.alloc(annotation, f)
    //     let input = curve_k.get_u() + bias;
    //     let res = mimc7_cs(
    //         &mut cs,
    //         input,
    //         input,
    //         &constants
    //     );
    //     enc_msg.push(e + res);
    // } 

    

    Ok(())
}
use super::constants::ROUND_CONSTANTS;
use bls12_381::Scalar;
use crate::{ConstraintSystem, SynthesisError, Variable};
use ff::{PrimeField, Field};

pub const MIMC7_ROUNDS: usize = 91;

pub fn mimc7<S: PrimeField>(xl: S, xr: S, constants: &[S]) -> S{
    assert_eq!(constants.len(), MIMC7_ROUNDS);
    let mut res = mimc7_round(xl, xr, &Field::zero());
    for i in 1..MIMC7_ROUNDS {
        res = mimc7_round(res, xr, &constants[i]);
    }
    res.add_assign(xr);
    res.add_assign(xl);
    res.add_assign(xr);

    res
}

pub fn mimc7_round<S: PrimeField>(mut msg: S, key: S, constant: &S) -> S{
    msg.add_assign(key);
    msg.add_assign(constant);
    let mut tmp = msg;
    let mut res = msg;
    tmp = tmp.square(); 
    let tmp2 = tmp.square();
    res.mul_assign(&tmp); 
    res.mul_assign(&tmp2);

    res
}

pub(crate) fn get_mimc_constants() -> Vec<Scalar> {
    let constants = (0..91)
        .map(|idx| Scalar::from_bytes(&ROUND_CONSTANTS[idx]).unwrap())
        .collect::<Vec<_>>();

    constants
}

pub(crate) fn mimc7_cs<S: PrimeField, CS: ConstraintSystem<S>>(
    cs: &mut CS,
    mut xl_value: Option<S>,
    mut xr_value: Option<S>,
    round_constants: &[S],
) -> Option<S> {
    assert_eq!(round_constants.len(), MIMC7_ROUNDS);

    // let mut xl_value = self.xl;
    let mut xl = cs
        .alloc(
            || "preimage xl",
            || xl_value.ok_or(SynthesisError::AssignmentMissing),
        )
        .unwrap();

    let copy_xl_value = xl_value.clone();
    let copy_xl = cs
        .alloc(
            || "copy xl",
            || xl_value.ok_or(SynthesisError::AssignmentMissing),
        )
        .unwrap();

    cs.enforce(
        || "copy xl",
        |lc| lc + xl,
        |lc| lc + CS::one(),
        |lc| lc + copy_xl,
    );

    // let xr_value = self.xr;
    let xr = cs
        .alloc(
            || "preimage xr",
            || xr_value.ok_or(SynthesisError::AssignmentMissing),
        )
        .unwrap();

    // MiMC7 ROUND
    for i in 0..MIMC7_ROUNDS {
        let x1_value;
        if i == 0 {
            x1_value = xl_value.map(|mut e| {
                e.add_assign(&xr_value.unwrap());
                e
            });
        } else {
            x1_value = xl_value.map(|mut e| {
                e.add_assign(&xr_value.unwrap());
                e.add_assign(&round_constants[i]);
                e
            });
        };

        let x1 = cs
            .alloc(|| "x", || x1_value.ok_or(SynthesisError::AssignmentMissing))
            .unwrap();

        // x1 = (xl + xr) * 1
        if i == 0 {
            cs.enforce(
                || "x1 = (xl + xr) * 1",
                |lc| lc + xl + xr,
                |lc| lc + CS::one(),
                |lc| lc + x1,
            );
        } else {
            cs.enforce(
                || "x1 = (xl + xr) * 1",
                |lc| lc + xl + xr + (round_constants[i], CS::one()),
                |lc| lc + CS::one(),
                |lc| lc + x1,
            );
        }

        // MiMC7 ROUND: x^2
        let x2_value = x1_value.map(|e| e.square());

        let x2 = cs
            .alloc(
                || "x2 = (xl + xr) * (xl + xr)",
                || x2_value.ok_or(SynthesisError::AssignmentMissing),
            )
            .unwrap();

        cs.enforce(
            || "x2 = (xl + xr) * (xl + xr)",
            |lc| lc + x1,
            |lc| lc + x1,
            |lc| lc + x2,
        );

        // MiMC7 ROUND: x^4
        let x4_value = x2_value.map(|e| e.square());

        let x4 = cs
            .alloc(
                || "x^4 = x2^2 * x2^2",
                || x4_value.ok_or(SynthesisError::AssignmentMissing),
            )
            .unwrap();

        cs.enforce(
            || "x^4 = x2^2 * x2^2",
            |lc| lc + x2,
            |lc| lc + x2,
            |lc| lc + x4,
        );

        let x3_value = x1_value.map(|mut e| {
            e.mul_assign(&x2_value.unwrap());
            e
        });

        let x3 = cs
            .alloc(
                || "x3 = x * x^2",
                || x3_value.ok_or(SynthesisError::AssignmentMissing),
            )
            .unwrap();

        cs.enforce(|| "x3 = x * x^2", |lc| lc + x1, |lc| lc + x2, |lc| lc + x3);

        let x7_value = x3_value.map(|mut e| {
            e.mul_assign(&x4_value.unwrap());
            e
        });

        let x7 = cs
            .alloc(
                || "x7 = x3 * x4",
                || x7_value.ok_or(SynthesisError::AssignmentMissing),
            )
            .unwrap();

        cs.enforce(|| "x7 = x3 * x4", |lc| lc + x3, |lc| lc + x4, |lc| lc + x7);

        xl = x7;
        xl_value = x7_value;
    }

    let output_value = xl_value.map(|mut e| {
        e.add_assign(&xr_value.unwrap());
        e.add_assign(&copy_xl_value.unwrap());
        e.add_assign(&xr_value.unwrap());
        e
    });

    let output = cs
        .alloc(
            || "image",
            || output_value.ok_or(SynthesisError::AssignmentMissing),
        )
        .unwrap();

    cs.enforce(
        || "image = xl + ((xl + xr + c)^7 + xr) + xr",
        |lc| lc + xl + xr + xr + copy_xl,
        |lc| lc + CS::one(),
        |lc| lc + output,
    );

    output_value
}

use std::ops::MulAssign;

use ff::{PrimeField, Field};

use bellman::{Circuit, ConstraintSystem, SynthesisError};

pub const MIMC_ROUNDS: usize = 322;

pub const MIMC7_ROUNDS: usize = 91;
/// This is an implementation of MiMC, specifically a
/// variant named `LongsightF322p3` for BLS12-381.
/// See http://eprint.iacr.org/2016/492 for more
/// information about this construction.
///
/// ```
/// function LongsightF322p3(xL ⦂ Fp, xR ⦂ Fp) {
///     for i from 0 up to 321 {
///         xL, xR := xR + (xL + Ci)^3, xL
///     }
///     return xL
/// }
/// ```
pub fn mimc<S: PrimeField>(mut xl: S, mut xr: S, constants: &[S]) -> S {
    assert_eq!(constants.len(), MIMC_ROUNDS);

    for c in constants {
        let mut tmp1 = xl;
        tmp1.add_assign(c);
        let mut tmp2 = tmp1.square();
        tmp2.mul_assign(&tmp1);
        tmp2.add_assign(&xr);
        xr = xl;
        xl = tmp2;
    }

    xl
}

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

#[allow(clippy::upper_case_acronyms)]
pub struct MiMC7Demo<'a, S: PrimeField> {
    pub xl: Option<S>,
    pub xr: Option<S>,
    pub constants: &'a [S],
}

/// This is our demo circuit for proving knowledge of the
/// preimage of a MiMC hash invocation.
#[allow(clippy::upper_case_acronyms)]
pub struct MiMCDemo<'a, S: PrimeField> {
    pub xl: Option<S>,
    pub xr: Option<S>,
    pub constants: &'a [S],
}



impl<'a, S: PrimeField> Circuit<S> for MiMC7Demo<'a, S> {
    fn synthesize<CS: ConstraintSystem<S>>(self, cs: &mut CS) -> Result<(), SynthesisError> {
        assert_eq!(self.constants.len(), MIMC7_ROUNDS);

        let mut xl_value = self.xl;
        let mut xl = cs.alloc(
            || "preimage xl",
            || xl_value.ok_or(SynthesisError::AssignmentMissing),
        )?;

        let copy_xl_value =self.xl;
        let copy_xl = cs.alloc(
            || "copy xl",
            || copy_xl_value.ok_or(SynthesisError::AssignmentMissing),
        )?;

        cs.enforce(
            ||"copy xl",
            |lc| lc + xl,
            |lc| lc + CS::one(),
            |lc| lc + copy_xl,
        );

        let xr_value = self.xr;
        let xr = cs.alloc(
            || "preimage xr",
            || xr_value.ok_or(SynthesisError::AssignmentMissing),
        )?;

        // MiMC7 ROUND
        for i in 0..MIMC7_ROUNDS {

            let x1_value;
            if i == 0 {
                x1_value = xl_value.map(|mut e| {
                    e.add_assign(&xr_value.unwrap());
                    e
                });
            }
            else {
                x1_value = xl_value.map(|mut e| {
                    e.add_assign(&xr_value.unwrap());
                    e.add_assign(&self.constants[i]);
                    e
                });
            };
            
            let x1 = cs.alloc(
                || "x",
                || x1_value.ok_or(SynthesisError::AssignmentMissing),
            )?;

            // x1 = (xl + xr) * 1
            if i ==0 {
                cs.enforce(
                    || "x1 = (xl + xr) * 1",
                    |lc| lc + xl + xr,
                    |lc| lc + CS::one(),
                    |lc |lc + x1,
                );    
            }
            else {
                cs.enforce(
                    || "x1 = (xl + xr) * 1",
                    |lc| lc + xl + xr + (self.constants[i], CS::one()),
                    |lc| lc + CS::one(),
                    |lc |lc + x1,
                );
            }
            
            // MiMC7 ROUND: x^2
            let x2_value = x1_value.map(|e|{
                e.square()
            });

            let x2 = cs.alloc(
                || "x2 = (xl + xr) * (xl + xr)",
                || x2_value.ok_or(SynthesisError::AssignmentMissing),
            )?;

            cs.enforce(
                || "x2 = (xl + xr) * (xl + xr)",
                |lc| lc + x1,
                |lc| lc + x1,
                |lc| lc + x2,
            );

            // MiMC7 ROUND: x^4
            let x4_value = x2_value.map(|e|{
                e.square()
            });

            let x4 = cs.alloc(
                || "x^4 = x2^2 * x2^2",
                || x4_value.ok_or(SynthesisError::AssignmentMissing),
            )?;

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

            let x3 = cs.alloc(
                || "x3 = x * x^2",
            || x3_value.ok_or(SynthesisError::AssignmentMissing),
            )?;

            cs.enforce(
                ||"x3 = x * x^2",
                |lc| lc + x1,
                |lc| lc + x2,
                |lc| lc + x3,
            );

            let x7_value = x3_value.map(|mut e| {
                e.mul_assign(&x4_value.unwrap());
                e
            });

            let x7 = cs.alloc(
                ||"x7 = x3 * x4",
                ||x7_value.ok_or(SynthesisError::AssignmentMissing),
            )?;

            cs.enforce(
                ||"x7 = x3 * x4",
                |lc| lc + x3,
                |lc| lc + x4,
                |lc|lc + x7,
            );

            xl = x7;
            xl_value = x7_value;
        }

        let output_value = xl_value.map(|mut e|{
            e.add_assign(&xr_value.unwrap());
            e.add_assign(&copy_xl_value.unwrap());
            e.add_assign(&xr_value.unwrap());
            e
        });

        let output = cs.alloc_input(
            || "image",
        ||output_value.ok_or(SynthesisError::AssignmentMissing),
        )?;

        cs.enforce(
            || "image = xl + ((xl + xr + c)^7 + xr) + xr",
            |lc| lc + xl + xr + xr + copy_xl,
            |lc| lc + CS::one(),
            |lc| lc + output,
        );

        Ok(())
    }
}


/// Our demo circuit implements this `Circuit` trait which
/// is used during paramgen and proving in order to
/// synthesize the constraint system.
impl<'a, S: PrimeField> Circuit<S> for MiMCDemo<'a, S> {
    fn synthesize<CS: ConstraintSystem<S>>(self, cs: &mut CS) -> Result<(), SynthesisError> {
        assert_eq!(self.constants.len(), MIMC_ROUNDS);

        // Allocate the first component of the preimage.
        let mut xl_value = self.xl;
        let mut xl = cs.alloc(
            || "preimage xl",
            || xl_value.ok_or(SynthesisError::AssignmentMissing),
        )?;

        // Allocate the second component of the preimage.
        let mut xr_value = self.xr;
        let mut xr = cs.alloc(
            || "preimage xr",
            || xr_value.ok_or(SynthesisError::AssignmentMissing),
        )?;

        for i in 0..MIMC_ROUNDS {
            // xL, xR := xR + (xL + Ci)^3, xL
            let cs = &mut cs.namespace(|| format!("round {}", i));

            // tmp = (xL + Ci)^2
            let tmp_value = xl_value.map(|mut e| {
                e.add_assign(&self.constants[i]);
                e.square()
            });
            let tmp = cs.alloc(
                || "tmp",
                || tmp_value.ok_or(SynthesisError::AssignmentMissing),
            )?;

            cs.enforce(
                || "tmp = (xL + Ci)^2",
                |lc| lc + xl + (self.constants[i], CS::one()),
                |lc| lc + xl + (self.constants[i], CS::one()),
                |lc| lc + tmp,
            );

            // new_xL = xR + (xL + Ci)^3
            // new_xL = xR + tmp * (xL + Ci)
            // new_xL - xR = tmp * (xL + Ci)
            let new_xl_value = xl_value.map(|mut e| {
                e.add_assign(&self.constants[i]);
                e.mul_assign(&tmp_value.unwrap());
                e.add_assign(&xr_value.unwrap());
                e
            });

            let new_xl = if i == (MIMC_ROUNDS - 1) {
                // This is the last round, xL is our image and so
                // we allocate a public input.
                cs.alloc_input(
                    || "image",
                    || new_xl_value.ok_or(SynthesisError::AssignmentMissing),
                )?
            } else {
                cs.alloc(
                    || "new_xl",
                    || new_xl_value.ok_or(SynthesisError::AssignmentMissing),
                )?
            };

            cs.enforce(
                || "new_xL = xR + (xL + Ci)^3",
                |lc| lc + tmp,
                |lc| lc + xl + (self.constants[i], CS::one()),
                |lc| lc + new_xl - xr,
            );

            // xR = xL
            xr = xl;
            xr_value = xl_value;

            // xL = new_xL
            xl = new_xl;
            xl_value = new_xl_value;
        }

        Ok(())
    }
}

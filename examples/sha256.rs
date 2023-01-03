use bellman::ConstraintSystem;
// Bring in some tools for using finite fiels
use ff::{PrimeField};

// We're going to use the BLS12-381 pairing-friendly elliptic curve.
use bls12_381::{Bls12};

use bellman::{
    Circuit, SynthesisError, gadgets::{
        boolean::{AllocatedBit, Boolean}, sha256::sha256, multipack
    },
    groth16
};

use rand::rngs::OsRng;
use sha2::{Digest, Sha256};
use std::{time::SystemTime};

fn sha256t<Scalar: PrimeField, CS: ConstraintSystem<Scalar>>(
        mut cs: CS,
        data: &[Boolean],
    ) -> Result<Vec<Boolean>, SynthesisError> {
        // Flip endianness of each input byte
        let input: Vec<_> = data
            .chunks(8)
            .map(|c| c.iter().rev())
            .flatten()
            .cloned()
            .collect();
        let res = sha256(cs.namespace(|| "SHA-256(input)"), &input)?;
        // Flip endianness of each output byte
        Ok(res
            .chunks(8)
            .map(|c| c.iter().rev())
            .flatten()
            .cloned()
            .collect())
}

struct Sha256Circuit {
    preimage: Option<[u8; 80]>,
}

impl<Scalar: PrimeField> Circuit<Scalar> for Sha256Circuit {
    fn synthesize<CS: ConstraintSystem<Scalar>>(self, cs: &mut CS) ->Result<(), SynthesisError>{
        let bit_values = if let Some(preimage) = self.preimage {
                         preimage
                             .into_iter()
                             .map(|byte| (0..8).map(move |i| (byte >> i) & 1u8 == 1u8))
                             .flatten()
                             .map(|b| Some(b))
                             .collect()
                     } else {
                         vec![None; 80 * 8]
                     };
                     assert_eq!(bit_values.len(), 80 * 8);
            
                     // Witness the bits of the preimage.
                     let preimage_bits = bit_values
                         .into_iter()
                         .enumerate()
                         // Allocate each bit.
                         .map(|(i, b)| {
                             AllocatedBit::alloc(cs.namespace(|| format!("preimage bit {}", i)), b)
                         })
                         // Convert the AllocatedBits into Booleans (required for the sha256 gadget).
                         .map(|b| b.map(Boolean::from))
                         .collect::<Result<Vec<_>, _>>()?;
            
                     // Compute hash = SHA-256d(preimage).
                     let hash = sha256t(cs.namespace(|| "SHA-256d(preimage)"), &preimage_bits)?;
            
                     // Expose the vector of 32 boolean variables as compact public inputs.
                     multipack::pack_into_inputs(cs.namespace(|| "pack hash"), &hash)

    }
}

fn main() {
     let params = {
     let c = Sha256Circuit { preimage: None };
     groth16::generate_random_parameters::<Bls12, _, _>(c, &mut OsRng).unwrap()
 };

    // Prepare the verification key (for proof verification).
    let pvk = groth16::prepare_verifying_key(&params.vk);

    // Pick a preimage and compute its hash.
    let preimage = [42; 80];
    let hash = Sha256::digest(&preimage);
    // Create an instance of our circuit (with the preimage as a witness).
    
    
    let now = SystemTime::now();
    
    let c = Sha256Circuit {
        preimage: Some(preimage),
    };

    // Create a Groth16 proof with our parameters.
    let proof = groth16::create_random_proof(c, &params, &mut OsRng).unwrap();

    match now.elapsed() {
        Ok(elapsed) => {
            println!("Time elapsed, {}", elapsed.as_millis());
        }
        Err(err) => {
            println!("Error calculating time taken, err: {}", err);
        }
    }

    // Pack the hash as inputs for proof verification.
    let hash_bits = multipack::bytes_to_bits_le(&hash);
    let inputs = multipack::compute_multipacking(&hash_bits);

    // Check the proof
    assert!(groth16::verify_proof(&pvk, &proof, &inputs).is_ok());
}
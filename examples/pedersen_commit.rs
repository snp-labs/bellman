
use std::time::SystemTime;

// Bring in some tools for using finite fiels
use ff::{Field};
use group::Curve;
use rand_xorshift::XorShiftRng;
use rand_core::{RngCore, SeedableRng};
// We're going to use the BLS12-381 pairing-friendly elliptic curve.
use bls12_381::{Bls12};
use jubjub::Fr;

use bellman::{
    Circuit, SynthesisError, gadgets::{ boolean::{Boolean}},
    ConstraintSystem, gadgets::pedersen::*,
    groth16
};

use rand::rngs::OsRng;

fn main() {
    let mut rng = XorShiftRng::from_seed([
        0x58, 0x62, 0xbe, 0x3d, 0x76, 0x3d, 0x31, 0x8d, 0x17, 0xdb, 0x37, 0x32, 0x54, 0x06, 0xbc,
        0xe5,
    ]);

    let pedersen_commitment = PedersenCommit {
        msg: Fr::random(&mut rng), 
        rand: Fr::random(&mut rng),
    };

    println!(
        "msg: {:?} \nrandomness: {:?}",
        pedersen_commitment.msg,
        pedersen_commitment.rand,
    );

    let expected_value_commitment =
        jubjub::ExtendedPoint::from(pedersen_commitment.commit()).to_affine();

    let params = {
        let c = PedersenCommitCircuit { 
            pedersen_commitment : None,
        };
        groth16::generate_random_parameters::<Bls12, _, _>(c, &mut OsRng).unwrap()
    };

    let pvk = groth16::prepare_verifying_key(&params.vk);

    let now = SystemTime::now();

    let c = PedersenCommitCircuit { 
        pedersen_commitment : Some(pedersen_commitment),
    };

    let proof = groth16::create_random_proof(c, &params, &mut OsRng).unwrap();
    
    match now.elapsed() {
        Ok(elapsed) => {
            println!("Time elapsed, {}", elapsed.as_millis());
        }
        Err(err) => {
            println!("Error calculating time taken, err: {}", err);
        }
    }

    let input = [
        expected_value_commitment.get_u(),
        expected_value_commitment.get_v(),
    ];

    assert!(groth16::verify_proof(&pvk, &proof, &input).is_ok());
}
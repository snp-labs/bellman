use lazy_static::lazy_static;
use group::{Curve};
use jubjub::{self, SubgroupPoint};
use bls12_381::Scalar;

/// Reference to a circuit version of a generator for fixed-base salar multiplication.
pub type FixedGenerator = &'static [Vec<(Scalar, Scalar)>];

/// Circuit version of a generator for fixed-base salar multiplication.
pub type FixedGeneratorOwned = Vec<Vec<(Scalar, Scalar)>>;

/// The number of chunks needed to represent a full scalar during fixed-base
/// exponentiation.
const FIXED_BASE_CHUNKS_PER_GENERATOR: usize = 84;

/// The `d` constant of the twisted Edwards curve.
pub(crate) const EDWARDS_D: Scalar = Scalar::from_raw([
    0x0106_5fd6_d634_3eb1,
    0x292d_7f6d_3757_9d26,
    0xf5fd_9207_e6bd_7fd4,
    0x2a93_18e7_4bfa_2b48,
]);

/// The scaling factor used for conversion to and from the Montgomery form.
pub(crate) const MONTGOMERY_SCALE: Scalar = Scalar::from_raw([
    0x8f45_35f7_cf82_b8d9,
    0xce40_6970_3da8_8abd,
    0x31de_341e_77d7_64e5,
    0x2762_de61_e862_645e,
]);

/// The `A` constant of the birationally equivalent Montgomery curve.
pub(crate) const MONTGOMERY_A: Scalar = Scalar::from_raw([
    0x0000_0000_0000_a002,
    0x0000_0000_0000_0000,
    0x0000_0000_0000_0000,
    0x0000_0000_0000_0000,
]);

/// Creates the 3-bit window table `[0, 1, ..., 8]` for different magnitudes of a fixed
/// generator.
pub fn generate_circuit_generator(mut gen: jubjub::SubgroupPoint) -> FixedGeneratorOwned {
    let mut windows = vec![];

    for _ in 0..FIXED_BASE_CHUNKS_PER_GENERATOR {
        let mut coeffs = vec![(Scalar::zero(), Scalar::one())];
        let mut g = gen;
        for _ in 0..7 {
            let g_affine = jubjub::ExtendedPoint::from(g).to_affine();
            coeffs.push((g_affine.get_u(), g_affine.get_v()));
            g += gen;
        }
        windows.push(coeffs);

        // gen = gen * 8
        gen = g;
    }

    windows
}

/// The value commitment is used to check balance between inputs and outputs. The value is
/// placed over this generator.
pub const _VALUE_COMMITMENT_VALUE_GENERATOR: SubgroupPoint = SubgroupPoint::from_raw_unchecked(
    Scalar::from_raw([
        0x3618_3b2c_b4d7_ef51,
        0x9472_c89a_c043_042d,
        0xd861_8ed1_d15f_ef4e,
        0x273f_910d_9ecc_1615,
    ]),
    Scalar::from_raw([
        0xa77a_81f5_0667_c8d7,
        0xbc33_32d0_fa1c_cd18,
        0xd322_94fd_8977_4ad6,
        0x466a_7e3a_82f6_7ab1,
    ]),
);

/// The value commitment is randomized over this generator, for privacy.
pub const _VALUE_COMMITMENT_RANDOMNESS_GENERATOR: SubgroupPoint = SubgroupPoint::from_raw_unchecked(
    Scalar::from_raw([
        0x3bce_3b77_9366_4337,
        0xd1d8_da41_af03_744e,
        0x7ff6_826a_d580_04b4,
        0x6800_f4fa_0f00_1cfc,
    ]),
    Scalar::from_raw([
        0x3cae_fab9_380b_6a8b,
        0xad46_f1b0_473b_803b,
        0xe6fb_2a6e_1e22_ab50,
        0x6d81_d3a9_cb45_dedb,
    ]),
);

lazy_static! {
    pub static ref VALUE_COMMITMENT_VALUE_GENERATOR: FixedGeneratorOwned =
        generate_circuit_generator(_VALUE_COMMITMENT_VALUE_GENERATOR);
    pub static ref VALUE_COMMITMENT_RANDOMNESS_GENERATOR: FixedGeneratorOwned =
        generate_circuit_generator(_VALUE_COMMITMENT_RANDOMNESS_GENERATOR);
}

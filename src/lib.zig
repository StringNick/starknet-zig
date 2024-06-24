pub const curve = struct {
    pub usingnamespace @import("math/curve/short_weierstrass/affine.zig");
    pub usingnamespace @import("math/curve/short_weierstrass/projective.zig");
    pub usingnamespace @import("math/fields/elliptic_curve.zig");
    pub usingnamespace @import("math/curve/short_weierstrass/projective_jacobian.zig");
    pub usingnamespace @import("math/curve/curve_params.zig");
};

pub const fields = struct {
    pub usingnamespace @import("math/fields/fields.zig");
    pub usingnamespace @import("math/fields/montgomery.zig");
    pub usingnamespace @import("math/fields/starknet.zig");
    pub usingnamespace @import("math/fields/arithmetic.zig");
    pub usingnamespace @import("math/fields/biginteger.zig");
    pub usingnamespace @import("math/fields/prime.zig");
    pub usingnamespace @import("math/fields/const_choice.zig");
};

pub const prime = struct {
    pub usingnamespace @import("math/numprime/prime.zig");
};

pub const bench = struct {
    pub usingnamespace @import("bench/bench.zig");
    // pub usingnamespace @import("bench/bench_field.zig");
};

pub const crypto = struct {
    pub usingnamespace @import("crypto/ecdsa.zig");
    pub usingnamespace @import("crypto/pedersen_hash.zig");
    pub usingnamespace @import("crypto/poseidon_hash.zig");
};

pub const core = struct {
    pub usingnamespace @import("core/types/eth_address.zig");
    pub usingnamespace @import("core/types/hash256.zig");
    pub usingnamespace @import("core/types/message.zig");
};

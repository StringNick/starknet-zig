const std = @import("std");

const bigInteger = @import("biginteger.zig").bigInt;

/// Compute CIOS multiplication of `a` * `b`
/// `q` is the modulus
/// `mu` is the inverse of -q modulo 2^{64}
/// Notice CIOS stands for Coarsely Integrated Operand Scanning
/// For more information see section 2.3.2 of Tolga Acar's thesis
/// https://www.microsoft.com/en-us/research/wp-content/uploads/1998/06/97Acar.pdf
/// ported from rust lambdaworks-math
pub fn cios(
    comptime N: usize,
    a: *const bigInteger(N),
    b: *const bigInteger(N),
    q: *const bigInteger(N),
    mu: u64,
) bigInteger(N) {
    var t = [_]u64{0} ** N;
    var t_extra = [_]u64{0} ** 2;

    inline for (0..N) |i| {
        // C := 0
        var c: u128 = 0;

        // for j=N-1 to 0
        //    (C,t[j]) := t[j] + a[j]*b[i] + C
        var cs: u128 = undefined;
        inline for (0..N) |j| {
            cs = t[j] + (@as(u128, a.limbs[j])) * (@as(u128, b.limbs[i])) + c;
            c = cs >> 64;
            t[j] = @truncate(cs);
        }

        // (t_extra[0],t_extra[1]) := t_extra[1] + C
        cs = @as(u128, t_extra[1]) + c;
        t_extra[0] = @truncate(cs >> 64);
        t_extra[1] = @truncate(cs);

        // m := t[N-1]*q'[N-1] mod D
        const m = ((@as(u128, t[0]) * @as(u128, mu)) << 64) >> 64;

        // (C,_) := t[N-1] + m*q[N-1]
        c = (@as(u128, t[0]) + m * (@as(u128, q.limbs[0]))) >> 64;

        // for j=N-1 to 1
        //    (C,t[j+1]) := t[j] + m*q[j] + C
        inline for (1..N) |j| {
            cs = @as(u128, t[j]) + m * @as(u128, q.limbs[j]) + c;
            c = cs >> 64;
            t[j - 1] = @truncate((cs << 64) >> 64);
        }

        // (C,t[0]) := t_extra[1] + C
        cs = @as(u128, t_extra[1]) + c;
        c = cs >> 64;
        t[N - 1] = @truncate((cs << 64) >> 64);

        // t_extra[1] := t_extra[0] + C
        t_extra[1] = t_extra[0] + @as(u64, @truncate(c));
    }

    const overflow = t_extra[1] > 0;
    var result: bigInteger(N) = .{ .limbs = t };

    if (overflow or q.cmp(&result).compare(.lte)) {
        _ = result.subWithBorrowAssign(q);
    }

    return result;
}

/// Compute CIOS multiplication of `a` * `b`
/// This is the Algorithm 2 described in the paper
/// "EdMSM: Multi-Scalar-Multiplication for SNARKs and Faster Montgomery multiplication"
/// https://eprint.iacr.org/2022/1400.pdf.
/// It is only suited for moduli with `q[0]` smaller than `2^63 - 1`.
/// `q` is the modulus
/// `mu` is the inverse of -q modulo 2^{64}
pub fn ciosOptimizedForModuliWithOneSpareBit(
    comptime N: usize,
    a: *const bigInteger(N),
    b: *const bigInteger(N),
    q: *const bigInteger(N),
    mu: u64,
) bigInteger(N) {
    var t = [_]u64{0} ** N;
    var t_extra: u64 = undefined;

    inline for (0..N) |i| {
        // C := 0
        var c: u128 = 0;

        // for j=N-1 to 0
        //    (C,t[j]) := t[j] + a[j]*b[i] + C
        var cs: u128 = undefined;
        inline for (0..N) |j| {
            cs = @as(u128, t[j]) + @as(u128, a.limbs[j]) * @as(u128, b.limbs[i]) + c;
            c = cs >> 64;
            t[j] = @truncate(cs);
        }

        t_extra = @truncate(c);

        // m := t[N-1]*q'[N-1] mod D
        const m = ((@as(u128, t[0]) * @as(u128, mu)) << 64) >> 64;

        // (C,_) := t[0] + m*q[0]
        c = (@as(u128, t[0]) + m * @as(u128, q.limbs[0])) >> 64;

        // for j=N-1 to 1
        //    (C,t[j+1]) := t[j] + m*q[j] + C
        inline for (1..N) |j| {
            cs = @as(u128, t[j]) + m * @as(u128, q.limbs[j]) + c;
            c = cs >> 64;
            t[j - 1] = @truncate((cs << 64) >> 64);
        }

        // (C,t[0]) := t_extra + C
        cs = @as(u128, t_extra) + c;
        t[N - 1] = @truncate((cs << 64) >> 64);
    }

    // std.mem.reverse(u64, &t);
    var result = bigInteger(N).init(t);

    if (q.cmp(&result).compare(.lte)) {
        _ = result.subWithBorrowAssign(q);
    }

    return result;
}

test "Montgomery: from and to montgomery equal" {
    const Felt252 = @import("starknet.zig").Felt252;

    const expected = bigInteger(4).fromInt(u256, 255);

    const montg = cios(4, &expected, &Felt252.R2, &Felt252.Modulus, Felt252.Inv);

    const from_montg = cios(4, &montg, &bigInteger(4).fromInt(u8, 1), &Felt252.Modulus, Felt252.Inv);

    try std.testing.expectEqual(expected, from_montg);
}

test "cios vs cios optimized" {
    const a = [6]u64{ 12432, 1241, 34343434, 3434343, 3434343, 34343434 };
    const b = [6]u64{ 12432, 1241, 34343434, 3434343, 3434343, 34343434 };
    const U384 = bigInteger(6);

    const x = U384.init(a);
    const y = U384.init(b);
    const m = U384.fromInt(u384, 0xcdb061954fdd36e5176f50dbdcfd349570a29ce1); // this is prime

    const mu: u64 = 16085280245840369887; // negative of the inverse of `m` modulo 2^{64}
    try std.testing.expectEqual(cios(6, &x, &y, &m, mu), ciosOptimizedForModuliWithOneSpareBit(6, &x, &y, &m, mu));
}

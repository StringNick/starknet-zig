const std = @import("std");

const bigInt = @import("biginteger.zig").bigInt;
const arithmetic = @import("arithmetic.zig");

/// Represents a finite field with runtime modulus.
///
/// This finite field struct encapsulates operations and properties related to arithmetic operations modulo a given modulus.
/// It provides methods for arithmetic operations on field elements within the finite field.
///
/// # Parameters
/// - `n_limbs`: The number of limbs used to represent a field element.
///
/// # Returns
/// A finite field struct with the specified modulus.
pub fn Field(comptime N: usize) type {
    return struct {
        const bigInteger = bigInt(N);
        const Self = @This();

        modulus: bigInteger,
        r2: bigInteger,
        fe: bigInteger,
        inv: u64,

        pub fn toBigInt(self: Self) bigInteger {
            // Initialize an array to store the limbs of the resulting value
            var r = self.fe.limbs;

            // Iterate over the limbs of the field element
            inline for (0..N) |i| {
                // Compute the Montgomery factor k
                const k: u64 = r[i] *% self.inv;
                var carry: u64 = 0;

                // Multiply the current limb with k and the modulus, adding the carry
                _ = arithmetic.macWithCarry(r[i], k, self.modulus.limbs[0], &carry);

                // Iterate over the remaining limbs and perform multiplication with carry
                inline for (1..N) |j|
                    r[(i + j) % N] = arithmetic.macWithCarry(r[(i + j) % N], k, self.modulus.limbs[j], &carry);

                // Store the final carry
                r[i % N] = carry;
            }

            // Return the resulting `bigInt` value
            return .{ .limbs = r };
        }

        pub fn modulusHasSpareBit(modulus: bigInteger) bool {
            // Check if the highest bit of the modulus is zero
            return modulus.limbs[N - 1] >> 63 == 0;
        }

        pub fn mulAssign(self: *Self, rhs: *const Self) void {
            // Dereference the pointer to obtain the actual field element

            if (modulusHasSpareBit(self.modulus)) {
                self.fe = ciosOptimizedForModuliWithOneSpareBit(N, self.fe, rhs.fe, self.modulus, self.inv);
            } else {
                self.fe = cios(N, self.fe, rhs.fe, self.modulus, self.inv);
            }
        }

        pub fn squareAssign(self: *Self) void {
            self.fe = sosSquare(N, self.fe, self.modulus, self.inv);
        }

        pub fn subModulusAssign(self: *Self) void {
            if (self.fe.gte(&self.modulus)) {
                _ = self.fe.subWithBorrowAssign(&self.modulus);
            }
        }

        pub fn subWithBorrowAssign(self: *Self, rhs: *const Self) bool {
            // Initialize a variable to hold the borrow
            var borrow: u8 = 0;

            // Iterate through each limb and perform subtraction with borrow
            inline for (0..N) |i|
                // Perform subtraction with borrow for the current limb
                borrow = arithmetic.sbb(u8, &self.limbs[i], rhs.limbs[i], borrow);

            // Return a flag indicating whether there was a borrow during subtraction
            return borrow != 0;
        }

        pub fn powToInt(self: *const Self, exponent: anytype) Self {
            if (@typeInfo(@TypeOf(exponent)) != .Int) @compileError("expected type int");
            // Initialize result as the identity element
            var res = fromIntWithData(u64, 1, self.r2, self.modulus, self.inv);

            // Copy the exponent for manipulation
            var exp = exponent;

            // Copy the base field element for manipulation
            var base = self.*;

            // Perform binary exponentiation algorithm
            while (exp > 0) : (exp /= 2) {
                // If current bit of exponent is 1, multiply result by base
                if (exp & 1 == 1) res.mulAssign(&base);
                // Square the base for the next iteration
                base.mulAssign(&base);
            }
            // Return the computed result
            return res;
        }

        pub fn isZero(self: *const Self) bool {
            return self.fe.eql(.{});
        }

        pub fn zero(base: *const Self) Self {
            return fromIntWithData(u64, 0, base.r2, base.modulus, base.inv);
        }

        pub fn negAssign(self: *Self) void {
            // Check if the field element is non-zero
            if (!self.isZero()) {
                // Create a temporary big integer representing the modulus
                var tmp = self.modulus;
                // Subtract the field element from the modulus with borrow and assign the result to tmp
                _ = tmp.subWithBorrowAssign(&self.fe);
                // Update the original field element with the negated value
                self.*.fe = tmp;
            }
        }

        pub fn inverse(self: *const Self) ?Self {
            // Check if the input is zero
            if (self.isZero()) return null;

            // Constant representing the value 1 in the field
            const o = bigInt.one();

            var u = self.fe;
            var v = self.modulus;
            var b: Self = .{ .fe = self.r2, .r2 = self.r2, .modulus = self.modulus, .inv = self.inv };
            var c = zero(self);

            const hasSpareBit = modulusHasSpareBit(self.fe);

            // Iterate while both u and v are not one
            while (!u.eql(o) and !v.eql(o)) {
                // Perform operations while u is even
                while (u.isEven()) {
                    u.div2Assign();

                    if (b.fe.isEven()) {
                        b.fe.div2Assign();
                    } else {
                        const carry = b.fe.addWithCarryAssign(&self.modulus);
                        b.fe.div2Assign();
                        // Handle overflow if necessary
                        if (!hasSpareBit and carry) {
                            b.fe.limbs[N - 1] |= 1 << 63;
                        }
                    }
                }

                // Perform operations while v is even
                while (v.isEven()) {
                    v.div2Assign();
                    if (c.fe.isEven()) {
                        c.fe.div2Assign();
                    } else {
                        const carry = c.fe.addWithCarryAssign(&self.modulus);
                        c.fe.div2Assign();
                        // Handle overflow if necessary
                        if (!hasSpareBit and carry) {
                            c.fe.limbs[N - 1] |= 1 << 63;
                        }
                    }
                }

                // Update based on u vs v values
                if (v.cmp(&u).compare(.lte)) {
                    _ = u.subWithBorrowAssign(&v);
                    b.subAssign(&c);
                } else {
                    _ = v.subWithBorrowAssign(&u);
                    c.subAssign(&b);
                }
            }

            // Return the inverse based on the final values of u and v
            return if (u.eql(o)) b else c;
        }

        pub fn subAssign(self: *Self, rhs: *const Self) void {
            // If b > a, add the modulus to `a` first.
            if (rhs.fe.cmp(&self.fe) == .gt)
                _ = self.fe.addWithCarryAssign(&self.modulus);

            // Perform the subtraction operation
            _ = self.fe.subWithBorrowAssign(&rhs.fe);
        }

        pub fn fromIntWithData(comptime T: type, value: T, R2: bigInteger, modulus: bigInteger, inv: u64) Self {
            return .{
                .modulus = modulus,
                .inv = inv,
                .r2 = R2,
                .fe = if (modulusHasSpareBit(modulus))
                    ciosOptimizedForModuliWithOneSpareBit(N, bigInteger.fromInt(T, value), R2, modulus, inv)
                else
                    cios(N, bigInteger.fromInt(T, value), R2, modulus, inv),
            };
        }

        pub fn toBigInt2(self: Self) bigInt {
            return cios(N, self.fe, comptime bigInt.one(), self.modulus, self.inv);
        }

        pub fn fromInt(comptime T: type, value: T, modulus: T) Self {
            const modulus_big = bigInteger.fromInt(T, modulus);
            const R2 = computeR2Montgomery(N, modulus_big);
            const Inv = computeInvMontgomery(N, modulus_big);

            return fromIntWithData(T, value, R2, modulus_big, Inv);
        }
    };
}
/// Compute CIOS multiplication of `a` * `b`
/// `q` is the modulus
/// `mu` is the inverse of -q modulo 2^{64}
/// Notice CIOS stands for Coarsely Integrated Operand Scanning
/// For more information see section 2.3.2 of Tolga Acar's thesis
/// https://www.microsoft.com/en-us/research/wp-content/uploads/1998/06/97Acar.pdf
/// ported from rust lambdaworks-math
pub fn cios(
    comptime N: usize,
    a: bigInt(N),
    b: bigInt(N),
    q: bigInt(N),
    mu: u64,
) bigInt(N) {
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
    var result: bigInt(N) = .{ .limbs = t };

    if (overflow or q.cmp(&result).compare(.lte)) {
        _ = result.subWithBorrowAssign(&q);
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
    a: bigInt(N),
    b: bigInt(N),
    q: bigInt(N),
    mu: u64,
) bigInt(N) {
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
    var result = bigInt(N).init(t);

    if (q.cmp(&result).compare(.lte)) {
        _ = result.subWithBorrowAssign(&q);
    }

    return result;
}

/// Computes the value of -M^{-1} mod 2^64.
///
/// This function calculates the modular inverse of `-MODULUS` modulo 2^64.
/// The computation involves exponentiating by the multiplicative group order, which is Euler's totient function (φ(2^64)) - 1.
///
/// # Returns:
/// The modular inverse of `-MODULUS` modulo 2^64.
///
/// Remarks:
/// - This function is used in Montgomery exponentiation to compute the R value.
/// - The computation involves exponentiating by the multiplicative group order, which is Euler's totient function (φ(2^64)) - 1.
pub fn computeInvMontgomery(comptime N: usize, modulus: bigInt(N)) u64 {
    // Initialize the modular inverse.
    var inv: u64 = 1;

    // Iterate through 63 times.
    for (0..63) |_| {
        // Square the modular inverse.
        inv *%= inv;
        // Multiply the modular inverse by the least significant limb of the modulus.
        inv *%= modulus.limbs[0];
    }
    // Negate the computed modular inverse.
    return -%inv;
}

/// Computes R^2 in Montgomery domain.
///
/// Montgomery multiplication requires precomputing R^2 mod modulus, where R is a power of 2
/// such that R > modulus. R^2 is used to convert a product back to the Montgomery domain.
///
/// Returns:
///     A big integer representing R^2 in the Montgomery domain.
pub fn computeR2Montgomery(comptime N: usize, modulus: bigInt(N)) bigInt(N) {
    // Initialize the loop counter
    var l: u32 = 0;

    // Define `c` as the largest power of 2 smaller than `modulus`
    while (l < N * @bitSizeOf(usize)) {
        // Check if modulus shifted right by `l` bits is not equal to zero
        if (modulus.shr(l).ne(.{})) break;
        l += 1;
    }
    var c = bigInt(N).one().shl(l);

    // Double `c` and reduce modulo `modulus` until getting
    // `2^{2 * number_limbs * word_size}` mod `modulus`
    var i: usize = 1;
    while (i <= 2 * N * @bitSizeOf(usize) - l) {
        // Double `c`
        const double_c = c.addWithCarry(&c);

        // Update `c` using modular reduction
        c = if (modulus.lte(&double_c[0]) or double_c[1])
            double_c[0].subWithBorrow(&modulus)[0]
        else
            double_c[0];

        i += 1;
    }

    return c;
}

// Separated Operand Scanning Method (2.3.1)
pub fn sosSquare(
    comptime N: usize,
    a: bigInt(N),
    q: bigInt(N),
    mu: u64,
) bigInt(N) {
    // NOTE: we use explicit `while` loops in this function because profiling pointed
    // at iterators of the form `(<x>..<y>).rev()` as the main performance bottleneck.

    // Step 1: Compute `(hi, lo) = a * a`
    var hi, var lo = a.square();

    // Step 2: Add terms to `(hi, lo)` until multiple it
    // is a multiple of both `2^{NUM_LIMBS * 64}` and
    // `q`.
    var c: u128 = 0;

    inline for (0..N) |i| {
        c = 0;
        const m: u64 = @truncate(@as(u128, lo.limbs[i]) * @as(u128, mu));
        inline for (0..N) |j| {
            if (i + j <= N - 1) {
                const index = i + j;
                const cs = @as(u128, lo.limbs[index]) + @as(u128, m) * @as(u128, q.limbs[j]) + c;
                c = cs >> 64;
                lo.limbs[index] = @truncate(cs);
            } else {
                const index = i + j - N + 1;
                const cs = @as(u128, hi.limbs[index]) + @as(u128, m) * @as(u128, q.limbs[j]) + c;
                c = cs >> 64;
                hi.limbs[index] = @truncate(cs);
            }
        }

        // Carry propagation to `hi`
        var t: usize = 0;
        while (c > 0 and i >= t) {
            const cs = @as(u128, hi.limbs[i - t]) + c;
            c = cs >> 64;
            hi.limbs[i - t] = @truncate(cs);
            t += 1;
        }
    }

    // Step 3: At this point `overflow * 2^{2 * NUM_LIMBS * 64} + (hi, lo)` is a multiple
    // of `2^{NUM_LIMBS * 64}` and the result is obtained by dividing it by `2^{NUM_LIMBS * 64}`.
    // In other words, `lo` is zero and the result is
    // `overflow * 2^{NUM_LIMBS * 64} + hi`.
    // That number is always strictly smaller than `2 * q`. To normalize it we substract
    // `q` whenever it is larger than `q`.
    // The easy case is when `overflow` is zero. We just use the `sub` function.
    // If `overflow` is 1, then `hi` is smaller than `q`. The function `sub(hi, q)` wraps
    // around `2^{NUM_LIMBS * 64}`. This is the result we need.
    const overflow = c > 0;
    if (overflow or q.lte(&hi)) {
        _ = hi.subWithBorrowAssign(&q);
    }
    return hi;
}

test "Montgomery: from and to montgomery equal" {
    const Felt252 = @import("starknet.zig").Felt252;

    const expected = bigInt(4).fromInt(u256, 255);

    const modulus = bigInt(4).fromInt(u256, 0x800000000000011000000000000000000000000000000000000000000000001);
    const R2 = computeR2Montgomery(4, modulus);
    const inv = computeInvMontgomery(4, modulus);

    try std.testing.expectEqualSlices(u64, modulus.limbs[0..], Felt252.Modulus.limbs[0..]);
    try std.testing.expectEqualSlices(u64, R2.limbs[0..], Felt252.R2.limbs[0..]);
    try std.testing.expectEqual(inv, Felt252.Inv);

    const montg = cios(4, expected, R2, modulus, Felt252.Inv);

    const from_montg = cios(4, montg, bigInt(4).fromInt(u8, 1), Felt252.Modulus, Felt252.Inv);

    try std.testing.expectEqual(expected, from_montg);
}

test "cios vs cios optimized" {
    const a = [6]u64{ 12432, 1241, 34343434, 3434343, 3434343, 34343434 };
    const b = [6]u64{ 12432, 1241, 34343434, 3434343, 3434343, 34343434 };
    const U384 = bigInt(6);

    const x = U384.init(a);
    const y = U384.init(b);
    const m = U384.fromInt(u384, 0xcdb061954fdd36e5176f50dbdcfd349570a29ce1); // this is prime

    const mu: u64 = 16085280245840369887; // negative of the inverse of `m` modulo 2^{64}
    try std.testing.expectEqual(cios(6, x, y, m, mu), ciosOptimizedForModuliWithOneSpareBit(6, x, y, m, mu));
}

test "sosSquare" {}

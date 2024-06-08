const std = @import("std");
const table = @import("table.zig");
const RndGen = std.rand.DefaultPrng;
const field = @import("../fields/fields.zig").Field;

pub const PrimalityTestConfig = struct {
    // TODO: add option to divides small primes in the table
    //       and this option should be enabled if the probabilistic test is used for strict config
    /// Number of strong probable prime test, starting from base 2
    sprp_trials: usize = 2,

    /// Number of strong probable prime test with random bases
    sprp_random_trials: usize = 3,

    /// Whether perform strong lucas probable prime test (with automatically selected parameters)
    slprp_test: bool = false,

    /// Whether perform extra strong lucas probable prime test (with automatically selected parameters)
    eslprp_test: bool = false,
};

pub fn mulm(comptime T: type, lhs: T, rhs: T, m: T) T {
    const dT = std.meta.Int(.unsigned, @typeInfo(T).Int.bits * 2);
    const result, const over = @mulWithOverflow(lhs, rhs);
    if (over == 0) {
        return result;
    }

    const val = @as(dT, lhs) * @as(dT, rhs) % @as(dT, m);
    return @truncate(val);
}

pub fn sqm(comptime T: type, self: T, m: T) T {
    return mulm(T, self, self, m);
}

/// NaiveBuffer implements a very simple Sieve of Eratosthenes
pub const NaiveBuffer = struct {
    list: std.ArrayList(u64), // list of found prime numbers
    next: u64, // all primes smaller than this value has to be in the prime list, should be an odd number

    pub fn init(allocator: std.mem.Allocator) !NaiveBuffer {
        var list = try std.ArrayList(u64).initCapacity(allocator, table.SMALL_PRIMES.len);
        errdefer list.deinit();

        for (table.SMALL_PRIMES) |p| {
            try list.append(@intCast(p));
        }

        return .{
            .list = list,
            .next = table.SMALL_PRIMES_NEXT,
        };
    }
};

// fn is_sprp(&self, base: Self) -> bool {
//     self.test_sprp(base).either(|v| v, |_| false)
// }
pub fn Either(comptime leftT: type, comptime rightT: type) type {
    return union(enum) {
        left: leftT,
        right: rightT,
    };
}

/// Modular multiplication
pub fn mulmod(comptime T: type, a: T, b: T, m: T) !T {
    if (T == comptime_int) return @mod(a * b, m);
    if (@typeInfo(T) != .Int) {
        @compileError("mulmod not implemented for " ++ @typeName(T));
    }

    const modulus = m;
    if (modulus == 0) return error.DivisionByZero;
    if (modulus < 0) return error.NegativeModulus;

    // On overflow, falling back on the multiplication property first
    if (std.math.mul(T, a, b) catch std.math.mul(T, @mod(a, m), @mod(b, m))) |product| {
        return @mod(product, modulus);
    } else |_| {
        const WideInt = std.meta.Int(@typeInfo(T).Int.signedness, @typeInfo(T).Int.bits * 2);
        return @intCast(@mod(@as(WideInt, a) * @as(WideInt, b), @as(WideInt, m)));
    }
}

/// Modular exponentiation
///
/// wikipedia.org/wiki/Modular_exponentiation#Right-to-left_binary_method
pub fn powmod(comptime T: type, base: T, exponent: T, modulus: T) !T {
    if (@typeInfo(T) != .Int) {
        @compileError("powmod not implemented for " ++ @typeName(T));
    }

    // zig fmt: off
    if (modulus == 0) return error.DivisionByZero;
    if (modulus  < 0) return error.NegativeModulus;
    // TODO Perform by finding the modular multiplicative inverse
    if (exponent < 0) return error.NegativeExponent;

    if (modulus  == 1) return 0;
    if (exponent == 2) return mulmod(T, base, base, modulus) catch unreachable;
    // zig fmt: on

    const Field = std.crypto.ff.Modulus(@bitSizeOf(T));

    var ff = try Field.fromPrimitive(T, modulus);

    const b = try Field.Fe.fromPrimitive(T, ff, base);

    const e = try Field.Fe.fromPrimitive(T, ff, exponent);

    const res = try ff.pow(
        b,
        e,
    );

    return try res.toPrimitive(T);
}

fn isSprp(comptime T: type, self: T, base: T) bool {
    return switch (testSprp(T, self, base)) {
        .left => |v| v,
        .right => |_| false,
    };
}

fn testSprp(comptime T: type, self: T, base: T) Either(bool, T) {
    if (self < 1) {
        return .{ .left = false };
    }

    // find 2^shift*u + 1 = n
    const tm1 = self - 1;
    const shift = @ctz(tm1);
    const u = std.math.shr(T, tm1, shift);

    // prevent reduction if the input is in montgomery form
    const m1 = 1 % self;
    const mm1 = (val: {
        const x = m1 % self;
        if (x == 0) {
            break :val 0;
        } else {
            break :val self - x;
        }
    });

    // var x: T = @intCast(powModulus(std.meta.Int(.unsigned, @typeInfo(T).Int.bits * 2), base, u, self));
    var x: T = powmod(T, base, u, self) catch unreachable;

    if (x == m1 or x == mm1) {
        return .{ .left = true };
    }

    for (0..shift) |_| {
        // const y: T = @intCast(powModulus(std.meta.Int(.unsigned, @typeInfo(T).Int.bits * 2), x, 2, self));
        const y: T = mulmod(T, x, x, self) catch unreachable;

        if (y == m1) {
            return .{ .right = std.math.gcd(self, x - 1) };
        }

        if (y == mm1) {
            return .{ .left = true };
        }

        x = y;
    }

    return .{ .left = x == m1 };
}

pub fn isPrime(comptime T: type, comptime cfg: ?PrimalityTestConfig, number: T) !bool {
    // // shortcuts
    //        if target.is_even() {
    //            return if target == &T::from_u8(2u8).unwrap() {
    //                Primality::Yes
    //            } else {
    //                Primality::No
    //            };
    //        }

    // // do deterministic test if target is under 2^64
    // if let Some(x) = target.to_u64() {
    //     return match is_prime64(x) {
    //         true => Primality::Yes,
    //         false => Primality::No,
    //     };
    // }

    // TODO: config as argument
    const config = cfg orelse PrimalityTestConfig{};
    // var probability: f32 = 1.0;

    // miller-rabin test
    var witness_list = [_]u64{0} ** (config.sprp_trials + config.sprp_random_trials);

    if (config.sprp_trials > 0) {
        witness_list[0..config.sprp_trials].* = table.SMALL_PRIMES[0..config.sprp_trials].*;
        // probability *= (1.0 - std.math.pow(f32, 0.25, @floatFromInt(config.sprp_trials)));
    }

    var rnd = RndGen.init(0);
    if (config.sprp_random_trials > 0) {
        for (0..config.sprp_random_trials) |i| {
            // we have ensured target is larger than 2^64

            var w: u64 = rnd.random().int(u64);

            while (std.mem.indexOfScalar(u64, &witness_list, w) != null) {
                w = rnd.random().int(u64);
            }

            witness_list[i + config.sprp_trials] = w;
        }

        // probability *= (1.0 - std.math.pow(f32, 0.25, @floatFromInt(config.sprp_random_trials)));
    }

    // we dont need another algos to check probability, only sprp

    for (witness_list) |x| {
        if (isSprp(T, number, x) == false) return false;
    }

    return true;
}

test "is_prime" {
    try std.testing.expectEqual(true, try isPrime(u64, null, 7));
    try std.testing.expectEqual(true, try isPrime(u64, null, 18446744069414584321));
    try std.testing.expectEqual(true, try isPrime(u256, null, 1489313108020924784844819367773615431304754137524579622245743070945963));
    try std.testing.expectEqual(false, try isPrime(u256, null, 1489313108020924784844819367773615431304754137524579622245743070945961));
}

test "sprp_test" {
    const spsp = [5]u16{ 2047, 3277, 4033, 4681, 8321 };
    inline for (spsp) |psp|
        try std.testing.expect(isSprp(u16, psp, 2));

    try std.testing.expectEqual(Either(bool, u16){ .right = 31 }, testSprp(u16, 341, 2));
}

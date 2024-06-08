const std = @import("std");
const table = @import("table.zig");
const RndGen = std.rand.DefaultPrng;

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

fn ipow(comptime T: type, b: T, e: T, modulus: T) T {
    var result: T = 1;
    var base = b;

    var exp = e;

    while (true) {
        if (exp & 1 != 0)
            result = mulm(T, result, base, modulus);
        exp >>= 1;

        if (exp != 0)
            break;

        base = mulm(T, base, base, modulus);
    }

    return result;
}

fn powModulus(comptime T: type, b: T, e: T, modulus: T) T {
    var base: T = b;
    var exponent: T = e;

    if (modulus == 1) return 0;

    base = base % modulus;

    var result: T = 1;

    while (exponent > 0) {
        if ((exponent & 1) == 1) {
            result = @rem(result * base, modulus);
        }

        base = @rem(base * base, modulus);
        exponent = exponent >> 1;
    }

    return result;
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
    var x: T = ipow(T, base, u, self);

    if (x == m1 or x == mm1) {
        return .{ .left = true };
    }

    for (0..shift) |_| {
        // const y: T = @intCast(powModulus(std.meta.Int(.unsigned, @typeInfo(T).Int.bits * 2), x, 2, self));
        const y: T = sqm(T, x, self);

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

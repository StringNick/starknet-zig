const std = @import("std");
const table = @import("table.zig");

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
    if (exponent == 2) return mulmod(T, base, base, modulus) ;
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

fn binpow(comptime T: type, _a: T, _b: T, m: T) T {
    var b = _b;
    var a = _a % m;
    var res: T = 1;

    while (b > 0) {
        if (b & 1 != 0)
            res = mulmod(T, res, a, m) catch unreachable;

        a = mulmod(T, a, a, m) catch unreachable;
        b >>= 1;
    }

    return res;
}

fn checkComposite(comptime T: type, n: T, a: T, d: T, s: usize) bool {
    var x = powmod(T, a, d, n) catch unreachable;

    if (x == 1 or x == n - 1)
        return false;

    for (1..s) |_| {
        x = mulmod(T, x, x, n) catch unreachable;
        if (x == n - 1)
            return false;
    }

    return true;
}

pub fn MillerRabin(comptime T: type, n: T) bool { // returns true if n is prime, else returns false.
    if (n < 2)
        return false;

    var r: usize = 0;
    var d = n - 1;

    while ((d & 1) == 0) {
        d >>= 1;
        r += 1;
    }

    inline for (table.SMALL_PRIMES[0..15]) |a| {
        if (n == a)
            return true;
        if (checkComposite(T, n, a, d, r))
            return false;
    }

    return true;
}

test "check primality" {
    std.log.err("check primality {any}", .{MillerRabin(u256, 1489313108020924784844819367773615431304754137524579622245743070945963)});
}

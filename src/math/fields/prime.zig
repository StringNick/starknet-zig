const std = @import("std");

const prime = @cImport(@cInclude("prime.h"));

pub fn isPrimeStdBigInt(t: std.math.big.int.Managed) bool {
    return prime.is_prime(@ptrCast(&t.limbs[0]), t.len());
}

pub fn isPrime(comptime T: type, t: T) bool {
    if (@typeInfo(T).Int.signedness == .signed) @compileError("type is signed, need unsigned");

    if (comptime @typeInfo(T).Int.bits % @bitSizeOf(u64) != 0) @compileError("type is not div exact on u64");

    return prime.is_prime(@ptrCast(@constCast(&t)), comptime @typeInfo(T).Int.bits / @bitSizeOf(u64));
}

test "test_is_prime" {
    try std.testing.expectEqual(true, isPrime(u256, 20988936657440586486151264256610222593863921));
}

test "test_prime" {
    var t = try std.math.big.int.Managed.initSet(std.testing.allocator, 179000439371151441295299346079990763017);
    defer t.deinit();

    try std.testing.expectEqual(true, isPrimeStdBigInt(t));
}

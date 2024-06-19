const zul = @import("zul");
const std = @import("std");
const Felt252 = @import("math/fields/starknet.zig").Felt252;

const a = Felt252.fromInt(
    u256,
    0x6606d7dccf23a0f61182da8d1149497f01b909036384bedb3e4c3284e2f2c1e1,
);
const b = Felt252.fromInt(
    u256,
    0x4cd366c0feadabcd6c61a395f6d9f91484ac4e51c3f8aede6c0ab49e2a55446a,
);

pub fn benchMul(_: std.mem.Allocator, _: *std.time.Timer) !void {
    _ = a.mul(&b);
}

pub fn benchMul2(_: std.mem.Allocator, _: *std.time.Timer) !void {
    _ = a.mul2(&b);
}

pub fn benchMulAssign(_: std.mem.Allocator, t: *std.time.Timer) !void {
    var c = Felt252.fromInt(u256, 12345);
    t.reset();
    c.mulAssign(&b);
}
pub fn benchMulAssign2(_: std.mem.Allocator, t: *std.time.Timer) !void {
    var c = Felt252.fromInt(u256, 12345);
    t.reset();
    c.mulAssign2(&b);
}

pub fn from2(context: Context, _: std.mem.Allocator, t: *std.time.Timer) !void {
    const idx = context.rand.random().intRangeLessThan(usize, 0, context.randomNumbers.len);

    const val = context.randomNumbers[idx];
    t.reset();
    _ = Felt252.fromInt2(u256, val);
}

pub fn from(context: Context, _: std.mem.Allocator, t: *std.time.Timer) !void {
    const idx = context.rand.random().intRangeLessThan(usize, 0, context.randomNumbers.len);

    const val = context.randomNumbers[idx];
    t.reset();

    _ = Felt252.fromInt(u256, val);
}

const Context = struct {
    randomNumbers: [1000]u256,
    rand: *std.Random.Xoshiro256,
};

pub fn main() !void {
    var rand = std.Random.DefaultPrng.init(13414);

    const ctxtes: []const []const Context = &.{
        ([_]Context{.{ .randomNumbers = undefined, .rand = undefined }} ** 3)[0..],
    };

    _ = ctxtes; // autofix
    var ctx = Context{
        .randomNumbers = undefined,
        .rand = &rand,
    };

    rand.fill(std.mem.asBytes(ctx.randomNumbers[0..]));

    (try zul.benchmark.runC(ctx, from, .{})).print("from");
    (try zul.benchmark.runC(ctx, from2, .{})).print("from2");

    (try zul.benchmark.run(benchMul, .{})).print("mul");
    (try zul.benchmark.run(benchMul2, .{})).print("mul2");
    (try zul.benchmark.run(benchMulAssign, .{})).print("mulAssign");
    (try zul.benchmark.run(benchMulAssign2, .{})).print("mulAssign2");
}

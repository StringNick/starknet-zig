const zul = @import("zul");
const std = @import("std");
const Felt252 = @import("math/fields/starknet.zig").Felt252;
const poseidon = @import("crypto/poseidon_hash.zig");
const pedersen = @import("crypto/pedersen_hash.zig");

const isPrime = @import("math/fields/prime.zig").isPrime;
const isPrime2 = @import("math/numprime/prime.zig").isPrime;

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

pub fn powToInt(context: Context, _: std.mem.Allocator, t: *std.time.Timer) !void {
    const idx = context.rand.random().intRangeLessThan(usize, 0, context.randomNumbers.len);

    const val = context.randomNumbers[idx];
    const d = Felt252.fromInt(u256, val);
    t.reset();

    _ = d.powToInt(3);
}

pub fn powToIntConst(context: Context, _: std.mem.Allocator, t: *std.time.Timer) !void {
    const idx = context.rand.random().intRangeLessThan(usize, 0, context.randomNumbers.len);

    const val = context.randomNumbers[idx];
    const d = Felt252.fromInt(u256, val);
    t.reset();

    _ = d.powToIntConst(3);
}

pub fn poseidonBench(context: Context, _: std.mem.Allocator, t: *std.time.Timer) !void {
    const idx = context.rand.random().intRangeLessThan(usize, 0, context.randomNumbers.len - 3);

    const val = context.randomNumbers[idx .. idx + 3];

    var hasher = poseidon.PoseidonHasher{ .state = .{
        Felt252.fromInt(u256, val[0]),
        Felt252.fromInt(u256, val[1]),
        Felt252.fromInt(u256, val[2]),
    } };
    t.reset();

    hasher.permuteComp();
}

pub fn pedersenBench(context: Context, _: std.mem.Allocator, t: *std.time.Timer) !void {
    const idx = context.rand.random().intRangeLessThan(usize, 0, context.randomNumbers.len - 2);

    const val1 = Felt252.fromInt(u256, context.randomNumbers[idx]);
    const val2 = Felt252.fromInt(u256, context.randomNumbers[idx + 1]);

    t.reset();

    _ = pedersen.pedersenHash(val1, val2);
}

pub fn isPrimeBench(context: Context, _: std.mem.Allocator, t: *std.time.Timer) !void {
    const idx = context.rand.random().intRangeLessThan(usize, 0, context.randomNumbers.len - 1);
    _ = idx; // autofix

    // const val1 = context.randomNumbers[idx];

    t.reset();
    _ = isPrime(u256, 18446744069414584321);
}

fn isPrime2Bench(context: Context, _: std.mem.Allocator, t: *std.time.Timer) !void {
    const idx = context.rand.random().intRangeLessThan(usize, 0, context.randomNumbers.len - 1);

    const val1 = context.randomNumbers[idx];
    _ = val1; // autofix

    t.reset();
    // big prime 1489313108020924784844819367773615431304754137524579622245743070945963
    _ = try isPrime2(u64, .{}, 18446744069414584321);
}

const Context = struct {
    randomNumbers: [1000]u256,
    rand: *std.Random.Xoshiro256,
};

pub fn main() !void {
    var rand = std.Random.DefaultPrng.init(13414);

    var ctx = Context{
        .randomNumbers = undefined,
        .rand = &rand,
    };

    rand.fill(std.mem.asBytes(ctx.randomNumbers[0..]));

    (try zul.benchmark.runC(ctx, isPrimeBench, .{})).print("isPrime");
    (try zul.benchmark.runC(ctx, isPrime2Bench, .{})).print("isPrime2");

    (try zul.benchmark.run(benchMulAssign, .{})).print("mulAssign");
    (try zul.benchmark.run(benchMulAssign2, .{})).print("mulAssign2");

    (try zul.benchmark.run(benchMul, .{})).print("mul");
    (try zul.benchmark.run(benchMul2, .{})).print("mul2");

    (try zul.benchmark.runC(ctx, from, .{})).print("from");
    (try zul.benchmark.runC(ctx, poseidonBench, .{})).print("poseidonPermuteComp");
    // (try zul.benchmark.runC(ctx, pedersenBench, .{})).print("pedersenHash");
    (try zul.benchmark.runC(ctx, from2, .{})).print("from2");
    (try zul.benchmark.runC(ctx, powToInt, .{})).print("powToInt");
    (try zul.benchmark.runC(ctx, powToIntConst, .{})).print("powToIntConst");
}

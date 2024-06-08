const benchmark = @import("bench.zig").benchmark;
const bigInt = @import("../math/fields/biginteger.zig").bigInt;
const Felt252 = @import("../math/fields/starknet.zig").Felt252;
const montgomery = @import("../math/fields/montgomery.zig");
const field = @import("../math/fields/fields.zig");
const std = @import("std");
const isPrime = @import("../math/fields/prime.zig").isPrime;
const isPrime2 = @import("../math/numprime/prime.zig").isPrime;

test "benchmark field multiplication" {
    try benchmark(struct {
        // How many iterations to run each benchmark.
        // If not present then a default will be used.
        pub const iterations = 1000;

        const a = Felt252.fromInt(
            u256,
            0x6606d7dccf23a0f61182da8d1149497f01b909036384bedb3e4c3284e2f2c1e1,
        );
        const b = Felt252.fromInt(
            u256,
            0x4cd366c0feadabcd6c61a395f6d9f91484ac4e51c3f8aede6c0ab49e2a55446a,
        );

        const v = std.mem.toBytes(@as(u256, 0x4cd366c0feadabcd6c61a395f6d9f91484ac4e51c3f8aede6c0ab49e2a55446a));

        pub fn benchIsPrime() void {
            _ = isPrime(u256, 1489313108020924784844819367773615431304754137524579622245743070945963);
        }

        pub fn benchIsPrime2() void {
            _ = isPrime2(u256, null, 1489313108020924784844819367773615431304754137524579622245743070945963) catch unreachable;
        }

        pub fn modFloor2() void {
            _ = a.modFloor2(b);
        }

        pub fn modFloor() void {
            _ = a.modFloor(b);
        }

        pub fn fromBytesLe() void {
            _ = bigInt(4).fromBytesLe(v);
        }
        pub fn fromBytesLe2() void {
            _ = bigInt(4).fromBytesLe2(v);
        }

        pub fn toDigitsLe() void {
            _ = a.toLeDigits();
        }
        pub fn toDigitsBe() void {
            _ = a.toBeDigits();
        }

        pub fn bench_div_rem() void {
            _, _ = a.divRem(b);
        }
        pub fn bench_div_rem2() void {
            _, _ = a.divRem2(b);
        }

        pub fn bench_from2() void {
            for (0..100) |_| {
                _ = Felt252.fromInt2(u256, 0x4cd366c0feadabcd6c61a395f6d9f91484ac4e51c3f8aede6c0ab49e2a55446a);
            }
        }

        pub fn bench_from() void {
            for (0..100) |_| {
                _ = Felt252.fromInt(u256, 0x4cd366c0feadabcd6c61a395f6d9f91484ac4e51c3f8aede6c0ab49e2a55446a);
            }
        }
    });
}

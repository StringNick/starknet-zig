const std = @import("std");
const table = @import("table.zig");
const RndGen = std.Random.DefaultPrng;
const field = @import("../fields/fields.zig").Field;

const bigInteger = @import("../fields/biginteger.zig").bigInt;
const montgomery = @import("../fields/montgomery.zig");
const arithmetic = @import("../fields/arithmetic.zig");

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

pub fn Either(comptime leftT: type, comptime rightT: type) type {
    return union(enum) {
        left: leftT,
        right: rightT,
    };
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

    const limb_count =
        comptime @typeInfo(T).Int.bits / 64;

    const F = montgomery.Field(limb_count);

    var x = F.fromInt(T, base, self);

    x = x.powToInt(u);

    const m1 = F.fromIntWithData(u64, 1, x.r2, x.modulus, x.inv);

    var mm1 = m1;
    mm1.negAssign();

    if (x.fe.eql(m1.fe) or x.fe.eql(mm1.fe)) {
        return .{ .left = true };
    }

    for (0..shift) |_| {
        var y = x;
        y.mulAssign(&y);

        if (y.fe.eql(m1.fe)) {
            return .{ .right = 0 };
        }

        if (y.fe.eql(mm1.fe)) {
            return .{ .left = true };
        }

        x = y;
    }

    return .{ .left = x.fe.eql(m1.fe) };
}

/// Very fast primality test on a u64 integer is a prime number. It's based on
/// deterministic Miller-rabin tests with hashing. if target is larger than 2^64 or more controlled
/// primality tests are desired, please use [is_prime()]
pub fn isPrime64(target: u64) bool {
    // shortcuts
    if (target < 2) {
        return false;
    }
    if (target & 1 == 0) {
        return target == 2;
    }

    const S = struct {
        fn order_u64(context: u64, rhs: u64) std.math.Order {
            return std.math.order(context, rhs);
        }
    };

    // remove small factors
    if (target < table.SMALL_PRIMES_NEXT) {
        // find in the prime list if the target is small enough
        return if (std.sort.binarySearch(u64, table.SMALL_PRIMES[0..], target, S.order_u64) != null) true else false;
    } else {
        // check remainder against the wheel table
        // this step eliminates any number that is not coprime to WHEEL_SIZE
        const pos: usize = @intCast(target % table.WHEEL_SIZE);
        if (pos == 0 or table.WHEEL_NEXT[pos] < table.WHEEL_NEXT[pos - 1]) {
            return false;
        }
    }

    // Then do a deterministic Miller-rabin test
    return isPrime64Miller(target);
}

fn isPrime32Miller(target: u32) bool {
    var h = target;
    h = ((h >> 16) ^ h) *% 0x45d9f3b;
    h = ((h >> 16) ^ h) *% 0x45d9f3b;
    h = ((h >> 16) ^ h) & 255;
    return isSprp(u64, target, table.MILLER_RABIN_BASE32[h]);
}

fn isPrime64Miller(target: u64) bool {
    if (target <= std.math.maxInt(u32)) {
        return isPrime32Miller(@truncate(target));
    }

    if (!isSprp(u64, target, 2)) {
        return false;
    }

    var h = target;
    h = ((h >> 32) ^ h) *% 0x45d9f3b3335b369;
    h = ((h >> 32) ^ h) *% 0x3335b36945d9f3b;
    h = ((h >> 32) ^ h) & 16383;
    const b = table.MILLER_RABIN_BASE64[h];

    return isSprp(u64, target, b & 4095) and isSprp(u64, target, b >> 12);
}

pub fn isPrime(comptime T: type, comptime cfg: ?PrimalityTestConfig, number: T) !bool {
    // shortcuts
    if ((number & 1) == 0)
        return if (number == 2)
            true
        else
            false;

    // do deterministic test if target is under 2^64
    if (comptime @typeInfo(T).Int.bits <= 64) {
        return isPrime64(number);
    }

    if (number <= std.math.maxInt(u64)) {
        return isPrime64(@truncate(number));
    }
    // TODO: config as argument
    const config = cfg orelse PrimalityTestConfig{};
    // var probability: f32 = 1.0;

    // miller-rabin test
    var witness_list = [_]u64{0} ** (config.sprp_trials + config.sprp_random_trials);

    if (config.sprp_trials > 0) {
        witness_list[0..config.sprp_trials].* = table.SMALL_PRIMES[0..config.sprp_trials].*;
        // probability *= (1.0 - std.math.pow(f32, 0.25, @floatFromInt(config.sprp_trials)));
    }

    var rnd = RndGen.init(@bitCast(std.time.timestamp()));
    if (config.sprp_random_trials > 0) {
        for (0..config.sprp_random_trials) |i| {
            // we have ensured target is larger than 2^64

            var w: u64 = rnd.random().int(u64);

            while (std.mem.indexOfScalar(u64, &witness_list, w) != null) {
                w = rnd.random().int(u64);
            }

            witness_list[i + config.sprp_trials] = w;
        }
    }

    // we dont need another algos to check probability, only sprp
    for (witness_list) |x| {
        if (isSprp(T, number, x) == false) return false;
    }

    return true;
}

test "is_prime" {
    try std.testing.expectEqual(true, try isPrime(u192, null, 3160627764927513289620384882792821105838947206441849759813));

    try std.testing.expectEqual(true, try isPrime(u192, null, 5764832813333353515948744950633024025186264784801315604993));

    try std.testing.expectEqual(true, try isPrime(u128, null, 3459023292465768231769));

    // less than equal 64 bit prime check
    try std.testing.expectEqual(false, try isPrime(u64, null, 4));
    try std.testing.expectEqual(true, try isPrime(u64, null, 7));
    try std.testing.expectEqual(true, try isPrime(u64, null, 193));
    try std.testing.expectEqual(true, try isPrime(u64, null, 63521));
    try std.testing.expectEqual(true, try isPrime(u64, null, 2876412481));
    try std.testing.expectEqual(true, try isPrime(u64, null, 758105553599));
    try std.testing.expectEqual(true, try isPrime(u64, null, 144260754288829));
    try std.testing.expectEqual(true, try isPrime(u64, null, 48156831561473449));
    try std.testing.expectEqual(true, try isPrime(u64, null, 10241082060017620583));
    try std.testing.expectEqual(true, try isPrime(u64, null, 12654157));
    try std.testing.expectEqual(true, try isPrime(u64, null, 18446744069414584321));
    // end of it

    try std.testing.expectEqual(false, try isPrime(u128, null, 174333704537955174084498843097741352677));
    try std.testing.expectEqual(false, try isPrime(u192, null, 13632814354644356817627273389719377694617167));

    try std.testing.expectEqual(true, try isPrime(u256, null, 1489313108020924784844819367773615431304754137524579622245743070945963));
    try std.testing.expectEqual(true, try isPrime(u256, null, 58439895430190581190666953963550055989211516870257943404603512697700968118731));
    try std.testing.expectEqual(true, try isPrime(u256, null, 94767248925316862866242448706799956186377412712579780364762034688283235353763));

    try std.testing.expectEqual(false, try isPrime(u256, null, 77371252455336267181195165));
    try std.testing.expectEqual(false, try isPrime(u256, null, 1489313108020924784844819367773615431304754137524579622245743070945961));
}

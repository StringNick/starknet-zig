// Core imports.
const std = @import("std");
const builtin = @import("builtin");

// Local imports.
const customlogFn = @import("utils/log.zig").logFn;
const Felt252 = @import("math/fields/starknet.zig").Felt252;
const isPrime2 = @import("math/numprime/prime.zig").isPrime;

// *****************************************************************************
// *                     GLOBAL CONFIGURATION                                  *
// *****************************************************************************

/// Standard library options.
/// log_level and log_scope_levels make it configurable.
pub const std_options = .{
    .logFn = customlogFn,
    .log_level = .debug,
    .log_scope_levels = &[_]std.log.ScopeLevel{},
};

pub fn main() !void {
    for (0..10000) |_| {
        _ = try isPrime2(u256, null, 1489313108020924784844819367773615431304754137524579622245743070945963);
    }
    std.log.debug("done", .{});
    // std.log.debug("starknet-zig\n", .{});
    // std.debug.print(
    //     \\Let's add two field elements together.
    //     \\We will use the Starknet prime field 0x800000000000011000000000000000000000000000000000000000000000001.
    //     \\We will add 0x800000000000011000000000000000000000000000000000000000000000000 and 0x4.
    //     \\The result should be 3.
    // , .{});
    // const a = Felt252.fromInt(u256, 0x800000000000011000000000000000000000000000000000000000000000000);
    // const b = Felt252.fromInt(u256, 0x4);
    // const c = a.add(&b);
    // std.debug.print("\nResult: {}\n", .{c.toU256()});
}

pub const TEST_ITERATIONS = 1;

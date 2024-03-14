/// Generating constants for poseidon hashing function
/// All memory allocation is on arena allocator
///
const std = @import("std");
const Felt252 = @import("./math/fields/starknet.zig").Felt252;
const Allocator = std.mem.Allocator;

const round_constants_block = "/// based on https://github.com/starkware-industries/poseidon/blob/5403dff9ff4eadb07deb5c0a43e88bedb011deb8/poseidon3.txt\n" ++
    "///\n///\n/// This file is autogenerated by poseidon_consts_gen.zig\n///\n///\n" ++
    "const Felt252 = @import(\"../../math/fields/starknet.zig\").Felt252;\n\npub const POSEIDON_COMPRESSED_ROUND_CONSTS: [{d}]Felt252 = .{{\n {s}\n }};" ++
    "\n\npub const POSEIDON_FULL_ROUNDS: usize = {d};\n\npub const POSEIDON_PARTIAL_ROUNDS: usize = {d};\n\n";

const round_constants_block_item = ".{{ .fe = .{{.limbs = .{{ {}, {}, {}, {}, }}, }}, }},\n";

const ConfigJSON = struct {
    full_rounds: usize,
    partial_rounds: usize,
    round_keys: [][3]u256,
};

// generateRoundConstantBlock - injecting compressed round constants and config into template
// result slice owner is caller, so it should be deinit by caller
fn generateRoundConstantBlock(allocator: Allocator, config: ConfigJSON, round_keys: []Felt252) ![]const u8 {
    var array_tpl = std.ArrayList(u8).init(allocator);
    defer array_tpl.deinit();

    for (round_keys) |round_key| {
        const value = round_key.fe;
        // writing array felt item
        try std.fmt.format(
            array_tpl.writer(),
            round_constants_block_item,
            .{
                value.limbs[0],
                value.limbs[1],
                value.limbs[2],
                value.limbs[3],
            },
        );
    }

    var result = std.ArrayList(u8).init(allocator);

    try std.fmt.format(
        result.writer(),
        round_constants_block,
        .{
            round_keys.len,
            try array_tpl.toOwnedSlice(),
            config.full_rounds,
            config.partial_rounds,
        },
    );

    return try result.toOwnedSlice();
}

// parseConfig - parsing config from json, allocator should be arena allocator
fn parseConfig(allocator: Allocator, json_spec: []const u8) !ConfigJSON {
    return try std.json.parseFromSliceLeaky(
        ConfigJSON,
        allocator,
        json_spec,
        .{ .allocate = std.json.AllocWhen.alloc_always },
    );
}

// compressRoundConstants - compressing round constants
// caller is owner of result slice and should deinit it
fn compressRoundConstants(allocator: Allocator, config: ConfigJSON, round_constants: [][3]Felt252) ![]Felt252 {
    var result = std.ArrayList(Felt252).init(allocator);

    for (round_constants[0 .. config.full_rounds / 2]) |rk| {
        inline for (0..3) |i| try result.append(rk[i]);
    }

    var idx = config.full_rounds / 2;

    var state = [_]Felt252{.{}} ** 3;

    // Add keys for partial rounds
    for (0..config.partial_rounds) |_| {
        inline for (0..3) |i|
            state[i].addAssign(&round_constants[idx][i]);

        // Add last state
        try result.append(state[2]);

        // Reset last state
        state[2] = .{};

        const st = state[0].add(&state[1]).add(&state[2]);

        // MixLayer
        state[0] = st.add(&Felt252.two().mul(&state[0]));
        state[1] = st.sub(&Felt252.two().mul(&state[1]));
        state[2] = st.sub(&Felt252.two().mul(&state[2]));

        idx += 1;
    }

    // Add keys for first of the last full rounds
    inline for (0..3) |i| {
        state[i].addAssign(&round_constants[idx][i]);
        try result.append(state[i]);
    }

    for (round_constants[config.full_rounds / 2 + config.partial_rounds + 1 ..]) |rk| {
        try result.appendSlice(&rk);
    }

    return try result.toOwnedSlice();
}

// parseNumbersToFieldElement - parsing numbers to field element
// caller is owner of result slice and should deinit it
fn parseNumbersToFieldElement(allocator: Allocator, keys: [][3]u256) ![][3]Felt252 {
    var result = try allocator.alloc([3]Felt252, keys.len);

    for (keys, 0..) |key, idx| {
        for (key, 0..) |k, idy| {
            result[idx][idy] = Felt252.fromInt(u256, k);
        }
    }

    return result;
}

pub fn main() !void {
    var arena = std.heap.ArenaAllocator.init(std.heap.page_allocator);
    defer arena.deinit();
    const allocator = arena.allocator();

    // writing constants for poseidon
    const config = try parseConfig(allocator, @embedFile("./crypto/poseidon/config.json"));

    const round_consts = try parseNumbersToFieldElement(allocator, config.round_keys);

    const compressed_round_consts = try compressRoundConstants(allocator, config, round_consts);

    const result = try generateRoundConstantBlock(allocator, config, compressed_round_consts);

    var file = try std.fs.cwd().openFile("./src/crypto/poseidon/constants.zig", .{ .mode = .write_only });

    try file.writer().writeAll(result);
}

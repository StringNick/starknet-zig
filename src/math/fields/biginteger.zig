const std = @import("std");
const arithmetic = @import("./arithmetic.zig");
const TEST_ITERATIONS = @import("../../main.zig").TEST_ITERATIONS;

const expectEqual = std.testing.expectEqual;
const expect = std.testing.expect;

/// This function generates a new type representing an arbitrary-precision integer with a fixed number of limbs.
/// Each limb is a 64-bit unsigned integer, allowing the representation of integers larger than the native machine word size.
///
/// Parameters:
///   - N: The number of limbs for the big integer.
///
/// Returns:
///   - A struct representing a big integer with the specified number of limbs.
pub fn bigInt(comptime N: usize) type {
    return struct {
        const Self = @This();

        /// The `limbs` field is an array of u64 integers that store the individual limbs of the big integer.
        ///
        /// Each limb contributes to the overall magnitude of the integer, allowing representation of large numbers.
        limbs: [N]u64 = [_]u64{0} ** N,

        /// Initializes a big integer with the specified limbs.
        ///
        /// This function initializes a big integer with the given array of limbs.
        /// It creates a new instance of the `bigInt` struct and assigns the provided limbs to the `limbs` field.
        ///
        /// Parameters:
        ///   - limbs: An array of u64 integers representing the limbs of the big integer.
        ///
        /// Returns:
        ///   - A new instance of the `bigInt` struct with the specified limbs.
        pub fn init(limbs: [N]u64) Self {
            return .{ .limbs = limbs };
        }

        /// Returns a big integer representing the value one.
        ///
        /// This function generates a big integer with a value of one. It creates a new instance of the `bigInt` struct
        /// where all limbs are initialized to zero, except for the first limb which is set to one, representing the value one.
        ///
        /// Returns:
        ///   - A new instance of the `bigInt` struct representing the value one.
        pub fn one() Self {
            // Compile-time computation to generate a big integer with a value of one
            comptime {
                // Initialize an array of limbs with all zeros
                var o = [_]u64{0} ** N;
                // Set the first limb to one
                o[0] = 1;
                // Return a new instance of the `bigInt` struct with the generated limbs
                return .{ .limbs = o };
            }
        }

        /// Generates a random big integer with the specified limbs.
        ///
        /// This function generates a random big integer with the specified number of limbs using a provided random number generator.
        /// It creates a new instance of the `bigInt` struct and assigns randomly generated values to each limb within the specified range.
        ///
        /// Parameters:
        ///   - r: An instance of the `std.Random` random number generator used to generate random limb values.
        ///
        /// Returns:
        ///   - A new instance of the `bigInt` struct representing a randomly generated big integer.
        pub fn rand(r: std.Random) Self {
            // Create a new instance of the `bigInt` struct
            var big_int: Self = .{};
            // Iterate over each limb
            for (0..N) |i|
                // Generate a random u64 value and assign it to the current limb
                big_int.limbs[i] = r.int(u64);

            // Return the randomly generated big integer
            return big_int;
        }

        /// Checks if two big integers are equal.
        ///
        /// This function compares two big integers for equality by performing a byte-wise memory comparison of their limbs.
        /// It returns true if the big integers have identical limb values, indicating equality, and false otherwise.
        ///
        /// Parameters:
        ///   - self: The first big integer to compare.
        ///   - rhs: The second big integer to compare.
        ///
        /// Returns:
        ///   - true if the big integers are equal, false otherwise.
        pub fn eql(self: Self, rhs: Self) bool {
            return std.mem.eql(u64, &self.limbs, &rhs.limbs);
        }

        /// Checks if a big integer is odd.
        ///
        /// This function determines whether a given big integer is odd by examining the least significant bit of its first limb.
        /// It returns true if the integer is odd (i.e., the LSB is set), and false otherwise.
        ///
        /// Parameters:
        ///   - self: A pointer to the big integer to be checked.
        ///
        /// Returns:
        ///   - true if the big integer is odd, false otherwise.
        pub fn isOdd(self: *const Self) bool {
            return self.limbs[0] & 1 == 1;
        }

        /// Checks if a big integer is even.
        ///
        /// This function determines whether a given big integer is even by negating the result of the `isOdd` function.
        /// It returns true if the integer is even (i.e., not odd), and false otherwise.
        ///
        /// Parameters:
        ///   - self: A pointer to the big integer to be checked.
        ///
        /// Returns:
        ///   - true if the big integer is even, false otherwise.
        pub fn isEven(self: *const Self) bool {
            return !self.isOdd();
        }

        /// Divides a big integer by two in place.
        ///
        /// This function divides a big integer by two in place, effectively performing a right shift operation on each limb.
        /// It iterates through each limb of the big integer, starting from the most significant limb, and performs the division operation.
        /// During the iteration, it propagates the carry bit from the higher-order bits to the lower-order bits to maintain precision.
        ///
        /// Parameters:
        ///   - self: A pointer to the big integer to be divided by two.
        ///
        /// Returns:
        ///   - void
        ///
        /// Notes:
        ///   - This function modifies the big integer in place, effectively dividing it by two.
        ///   - The carry bit from the higher-order bits is propagated to the lower-order bits to maintain precision during the division operation.
        ///   - Inline loops are used for performance optimization.
        ///   - The operation effectively performs a right shift of each limb by one bit, equivalent to division by two.
        ///   - The division is performed in place, and the result is updated in the original big integer.
        pub inline fn div2Assign(self: *Self) void {
            // Initialize a variable to hold the carry
            var t: u64 = 0;

            // Compile-time variable to iterate through limbs
            comptime var i = N;

            // Iterate through limbs starting from the most significant
            inline while (i > 0) {
                // Decrement the iterator
                i -= 1;
                // Get a pointer to the current limb
                const a = &self.limbs[i];
                // Shift the carry bit to the left
                const t2 = a.* << 63;
                // Divide the limb by 2 (right shift)
                a.* >>= 1;
                // Add the carry bit to the current limb
                a.* |= t;
                // Update the carry with the previous carry bit
                t = t2;
            }
        }

        /// Adds a big integer to another big integer with carry and returns the result along with a carry flag.
        ///
        /// This function adds a big integer `rhs` to another big integer `self` with carry and returns the result along with a carry flag.
        /// It dereferences the pointers to obtain the actual big integers, performs the addition operation using the `addWithCarryAssign` function,
        /// and then returns a tuple containing the result of the addition operation (`a`) and a boolean flag indicating whether there was a carry.
        ///
        /// Parameters:
        ///   - self: A pointer to the first big integer to which the second big integer will be added.
        ///   - rhs: A pointer to the second big integer that will be added to the first big integer.
        ///
        /// Returns:
        ///   - A tuple containing the result of the addition operation and a boolean flag indicating whether there was a carry.
        ///
        /// Notes:
        ///   - This function does not modify the original big integers; it only performs the addition operation.
        ///   - Addition with carry is performed using the `addWithCarryAssign` function.
        ///   - The result of the addition operation and the carry flag are returned as a tuple.
        ///   - The carry flag can be used to detect overflow during addition operations.
        pub fn addWithCarry(self: *const Self, rhs: *const Self) std.meta.Tuple(&.{ Self, bool }) {
            // Dereference the pointer to obtain the actual big integer
            var a = self.*;
            // Perform addition with carry and obtain the carry flag
            const carry = a.addWithCarryAssign(rhs);
            // Return a tuple containing the result of the addition operation and the carry flag
            return .{ a, carry };
        }

        /// Adds a big integer to another big integer with carry and assigns the result to the first big integer.
        ///
        /// This function performs addition with carry between two big integers and assigns the result to the first big integer.
        /// It iterates through each limb of the big integers, adding the corresponding limbs along with the carry bit.
        /// At the end of the addition operation, it returns a flag indicating whether there was a carry during the addition process.
        ///
        /// Parameters:
        ///   - self: A pointer to the first big integer to which the addition result will be assigned.
        ///   - rhs: A pointer to the second big integer that will be added to the first big integer.
        ///
        /// Returns:
        ///   - A boolean flag indicating whether there was a carry during the addition operation.
        ///
        /// Notes:
        ///   - This function modifies the first big integer in place, assigning the result of the addition operation to it.
        ///   - Addition with carry is performed for each corresponding limb of the big integers, starting from the least significant limb.
        ///   - The carry bit from each addition operation is propagated to the next higher-order limb to ensure correct arithmetic.
        ///   - Inline loops are used for performance optimization.
        ///   - At the end of the addition operation, the function returns true if there was a carry, and false otherwise.
        ///   - The carry flag can be used to detect overflow during addition operations.
        ///
        pub fn addWithCarryAssign(self: *Self, rhs: *const Self) bool {
            // Initialize a variable to hold the carry
            var carry: u8 = 0;

            // Iterate through each limb and perform addition with carry
            inline for (0..N) |i|
                // Perform addition with carry for the current limb
                carry = arithmetic.adc(u8, &self.limbs[i], rhs.limbs[i], carry);

            // Return a flag indicating whether there was a carry during addition
            return carry != 0;
        }

        /// Subtracts a big integer from another big integer with borrow and returns the result along with a borrow flag.
        ///
        /// This function subtracts a big integer `rhs` from another big integer `self` with borrow and returns the result along with a borrow flag.
        /// It dereferences the pointers to obtain the actual big integers, performs the subtraction operation using the `subWithBorrowAssign` function,
        /// and then returns a tuple containing the result of the subtraction operation (`a`) and a boolean flag indicating whether there was a borrow.
        ///
        /// Parameters:
        ///   - self: A pointer to the first big integer from which the second big integer will be subtracted.
        ///   - rhs: A pointer to the second big integer that will be subtracted from the first big integer.
        ///
        /// Returns:
        ///   - A tuple containing the result of the subtraction operation and a boolean flag indicating whether there was a borrow.
        ///
        /// Notes:
        ///   - This function does not modify the original big integers; it only performs the subtraction operation.
        ///   - Subtraction with borrow is performed using the `subWithBorrowAssign` function.
        ///   - The result of the subtraction operation and the borrow flag are returned as a tuple.
        ///   - The borrow flag can be used to detect underflow during subtraction operations.
        pub fn subWithBorrow(self: *const Self, rhs: *const Self) std.meta.Tuple(&.{ Self, bool }) {
            // Dereference the pointer to obtain the actual big integer
            var a = self.*;
            // Perform subtraction with borrow and obtain the borrow flag
            const borrow = a.subWithBorrowAssign(rhs);
            // Return a tuple containing the result of the subtraction operation and the borrow flag
            return .{ a, borrow };
        }

        /// Subtracts a big integer from another big integer with borrow and assigns the result to the first big integer.
        ///
        /// This function performs subtraction with borrow between two big integers and assigns the result to the first big integer.
        /// It iterates through each limb of the big integers, subtracting the corresponding limbs along with the borrow bit.
        /// At the end of the subtraction operation, it returns a flag indicating whether there was a borrow during the subtraction process.
        ///
        /// Parameters:
        ///   - self: A pointer to the first big integer from which the second big integer will be subtracted, and the result will be assigned.
        ///   - rhs: A pointer to the second big integer that will be subtracted from the first big integer.
        ///
        /// Returns:
        ///   - A boolean flag indicating whether there was a borrow during the subtraction operation.
        ///
        /// Notes:
        ///   - This function modifies the first big integer in place, assigning the result of the subtraction operation to it.
        ///   - Subtraction with borrow is performed for each corresponding limb of the big integers, starting from the least significant limb.
        ///   - The borrow bit from each subtraction operation is propagated to the next higher-order limb to ensure correct arithmetic.
        ///   - Inline loops are used for performance optimization.
        ///   - At the end of the subtraction operation, the function returns true if there was a borrow, and false otherwise.
        ///   - The borrow flag can be used to detect underflow during subtraction operations.
        pub fn subWithBorrowAssign(self: *Self, rhs: *const Self) bool {
            // Initialize a variable to hold the borrow
            var borrow: u8 = 0;

            // Iterate through each limb and perform subtraction with borrow
            inline for (0..N) |i|
                // Perform subtraction with borrow for the current limb
                borrow = arithmetic.sbb(u8, &self.limbs[i], rhs.limbs[i], borrow);

            // Return a flag indicating whether there was a borrow during subtraction
            return borrow != 0;
        }

        /// Compares two big integers and returns their relative order.
        ///
        /// This function compares two big integers `self` and `rhs` and returns their relative order.
        /// It first reverses the order of limbs from little-endian to big-endian to ensure correct comparison.
        /// Then, it compares the big integers using byte-wise comparison to determine their relative order.
        ///
        /// Parameters:
        ///   - self: A pointer to the first big integer to be compared.
        ///   - rhs: A pointer to the second big integer to be compared.
        ///
        /// Returns:
        ///   - An enum value representing the relative order of the two big integers:
        ///     - `std.math.Order.lt` if `self` is less than `rhs`.
        ///     - `std.math.Order.eq` if `self` is equal to `rhs`.
        ///     - `std.math.Order.gt` if `self` is greater than `rhs`.
        ///
        /// Notes:
        ///   - This function does not modify the original big integers; it only performs the comparison operation.
        ///   - The big integers are compared in a byte-wise manner after converting their limbs to big-endian order.
        ///   - The comparison result is returned as an enum value indicating the relative order.
        ///   - This function can be used to determine the relative order of big integers for sorting or comparison purposes.
        pub fn cmp(self: *const Self, rhs: *const Self) std.math.Order {
            // Obtain pointers to the limbs of the big integers
            var a = self.limbs;
            var b = rhs.limbs;

            // Reverse the order of limbs from little-endian to big-endian
            _ = std.mem.reverse(u64, a[0..]);
            _ = std.mem.reverse(u64, b[0..]);

            // Compare the big integers using byte-wise comparison
            return std.mem.order(u64, &a, &b);
        }

        /// Converts a big integer to a little-endian bit representation.
        ///
        /// This function converts a big integer to its little-endian bit representation.
        /// It iterates through each limb of the big integer, starting from the least significant limb,
        /// and generates a bit representation where each bit corresponds to a single limb.
        ///
        /// Parameters:
        ///   - self: A pointer to the big integer to be converted to a little-endian bit representation.
        ///
        /// Returns:
        ///   - An array of boolean values representing the little-endian bit representation of the big integer.
        ///
        /// Notes:
        ///   - The function generates a bit representation where the least significant bit of the big integer corresponds to the first bit of the array.
        ///   - Each limb of the big integer contributes 64 bits to the overall bit representation.
        ///   - The resulting bit representation is little-endian, with the least significant bits appearing first.
        ///   - Inline loops are used for performance optimization.
        ///   - The function returns an array of boolean values representing the bit representation of the big integer.
        pub fn toBitsLe(self: *const Self) [N * 64]bool {
            // Initialize an array to hold the bit representation
            var bits = [_]bool{false} ** (N * 64);

            // Iterate through each limb of the big integer
            inline for (0..N) |idx_limb| {
                // Calculate the starting index for the current limb
                const i = idx_limb * 64;
                // Iterate through each bit of the current limb
                inline for (0..64) |ind_bit|
                    // Extract the bit value and assign it to the corresponding position in the bit representation array
                    bits[i + ind_bit] = (self.limbs[idx_limb] >> ind_bit) & 1 == 1;
            }

            // Return the little-endian bit representation of the big integer
            return bits;
        }

        /// Converts a big integer to a big-endian bit representation.
        ///
        /// This function converts a big integer to its big-endian bit representation.
        /// It iterates through each limb of the big integer, starting from the most significant limb,
        /// and generates a bit representation where each bit corresponds to a single limb.
        ///
        /// Parameters:
        ///   - self: A pointer to the big integer to be converted to a big-endian bit representation.
        ///
        /// Returns:
        ///   - An array of boolean values representing the big-endian bit representation of the big integer.
        ///
        /// Notes:
        ///   - The function generates a bit representation where the most significant bit of the big integer corresponds to the first bit of the array.
        ///   - Each limb of the big integer contributes 64 bits to the overall bit representation.
        ///   - The resulting bit representation is big-endian, with the most significant bits appearing first.
        ///   - Inline loops are used for performance optimization.
        ///   - The function returns an array of boolean values representing the bit representation of the big integer.
        pub fn toBitsBe(self: *const Self) [N * 64]bool {
            // Initialize an array to hold the bit representation
            var bits = [_]bool{false} ** (N * 64);

            // Iterate through each limb of the big integer
            inline for (0..N) |idx_limb| {
                // Calculate the starting index for the current limb in the big-endian bit representation
                const pre_index = (N - idx_limb - 1) * 64 + 63;
                // Iterate through each bit of the current limb
                inline for (0..64) |ind_bit|
                    // Calculate the index in the bit representation array and extract the bit value
                    bits[pre_index - ind_bit] = (self.limbs[idx_limb] >> ind_bit) & 1 == 1;
            }

            // Return the big-endian bit representation of the big integer
            return bits;
        }
    };
}

test "bigInt: div2 function should divide by 2" {
    // Test case: Divide a big integer by 2
    // Initialize a big integer with a value
    var a = bigInt(4).init(.{ 8446744, 0, 0, 0 });
    // Perform division by 2 in place
    a.div2Assign();
    // Assert that the result matches the expected value
    try expectEqual([_]u64{ 4223372, 0, 0, 0 }, a.limbs);

    // Test case: Divide zero by 2
    // Initialize a big integer with zero value
    var b = bigInt(4).init(.{ 0, 0, 0, 0 });
    // Perform division by 2 in place
    b.div2Assign();
    // Assert that the result is still zero
    try expectEqual([_]u64{ 0, 0, 0, 0 }, b.limbs);
}

test "bigInt: constants" {
    // Test case: Verify initialization of zero constant
    // Assert that the zero constant matches the expected value
    try expectEqual(bigInt(4).init(.{ 0, 0, 0, 0 }), comptime bigInt(4){});

    // Test case: Verify initialization of one constant
    // Assert that the one constant matches the expected value
    try expectEqual(bigInt(4).init(.{ 1, 0, 0, 0 }), comptime bigInt(4).one());
}

test "bigInt: fuzzing test for equality" {
    // Test case: Fuzzing test for equality
    // Initialize a pseudo-random number generator (PRNG) with a seed of 0.
    var prng = std.Random.DefaultPrng.init(0);
    // Generate a random number using the PRNG.
    const random = prng.random();

    // Iterate over the test iterations.
    for (0..TEST_ITERATIONS) |_| {
        // Test case: Verify equality of randomly generated big integers
        // Generate random unsigned integers of different sizes.
        var a = bigInt(4).rand(random);
        var b = bigInt(4).rand(random);
        var c = bigInt(4).rand(random);
        // Obtain constant big integers for comparison
        var one = comptime bigInt(4).one();
        var zero = comptime bigInt(4){};

        // Assert that each big integer is equal to itself
        try expect(a.eql(a));
        try expect(b.eql(b));
        try expect(c.eql(c));

        // Assert that constant big integers are equal to themselves
        try expect(one.eql(one));
        try expect(zero.eql(zero));
    }
}

test "bigInt: fuzzing test for add and sub operations" {
    // Test case: Fuzzing test for addition and subtraction operations
    // Initialize a pseudo-random number generator (PRNG) with a seed of 0.
    var prng = std.Random.DefaultPrng.init(0);
    // Generate a random number using the PRNG.
    const random = prng.random();

    // Iterate over the test iterations.
    for (0..TEST_ITERATIONS) |_| {
        // Test case: Verify addition and subtraction operations with random big integers
        // Generate random unsigned integers of different sizes.
        const a = bigInt(4).rand(random);
        const b = bigInt(4).rand(random);
        const c = bigInt(4).rand(random);

        // Constant zero big integer
        const zero = comptime bigInt(4){};

        // Test addition with zero
        // Assert that adding zero to a big integer results in the same big integer
        try expectEqual(
            @as(std.meta.Tuple(&.{ bigInt(4), bool }), .{ a, false }),
            a.addWithCarry(&zero),
        );
        try expectEqual(
            @as(std.meta.Tuple(&.{ bigInt(4), bool }), .{ b, false }),
            b.addWithCarry(&zero),
        );
        try expectEqual(
            @as(std.meta.Tuple(&.{ bigInt(4), bool }), .{ c, false }),
            c.addWithCarry(&zero),
        );

        // Test subtraction with zero
        // Assert that subtracting zero from a big integer results in the same big integer
        try expectEqual(
            @as(std.meta.Tuple(&.{ bigInt(4), bool }), .{ a, false }),
            a.subWithBorrow(&zero),
        );
        try expectEqual(
            @as(std.meta.Tuple(&.{ bigInt(4), bool }), .{ b, false }),
            b.subWithBorrow(&zero),
        );
        try expectEqual(
            @as(std.meta.Tuple(&.{ bigInt(4), bool }), .{ c, false }),
            c.subWithBorrow(&zero),
        );

        // Test subtraction resulting in zero
        // Assert that subtracting a big integer from itself results in zero
        try expectEqual(
            @as(std.meta.Tuple(&.{ bigInt(4), bool }), .{ zero, false }),
            a.subWithBorrow(&a),
        );
        try expectEqual(
            @as(std.meta.Tuple(&.{ bigInt(4), bool }), .{ zero, false }),
            b.subWithBorrow(&b),
        );
        try expectEqual(
            @as(std.meta.Tuple(&.{ bigInt(4), bool }), .{ zero, false }),
            c.subWithBorrow(&c),
        );
    }
}

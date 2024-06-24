const std = @import("std");

const bigInt = @import("biginteger.zig").bigInt;
const arithmetic = @import("arithmetic.zig");

/// Compute CIOS multiplication of `a` * `b`
/// `q` is the modulus
/// `mu` is the inverse of -q modulo 2^{64}
/// Notice CIOS stands for Coarsely Integrated Operand Scanning
/// For more information see section 2.3.2 of Tolga Acar's thesis
/// https://www.microsoft.com/en-us/research/wp-content/uploads/1998/06/97Acar.pdf
/// ported from rust lambdaworks-math
pub fn cios(
    comptime N: usize,
    a: bigInt(N),
    b: bigInt(N),
    q: bigInt(N),
    mu: u64,
) bigInt(N) {
    var t = [_]u64{0} ** N;
    var t_extra = [_]u64{0} ** 2;

    inline for (0..N) |i| {
        // C := 0
        var c: u128 = 0;

        // for j=N-1 to 0
        //    (C,t[j]) := t[j] + a[j]*b[i] + C
        var cs: u128 = undefined;
        inline for (0..N) |j| {
            cs = t[j] + (@as(u128, a.limbs[j])) * (@as(u128, b.limbs[i])) + c;
            c = cs >> 64;
            t[j] = @truncate(cs);
        }

        // (t_extra[0],t_extra[1]) := t_extra[1] + C
        cs = @as(u128, t_extra[1]) + c;
        t_extra[0] = @truncate(cs >> 64);
        t_extra[1] = @truncate(cs);

        // m := t[N-1]*q'[N-1] mod D
        const m = ((@as(u128, t[0]) * @as(u128, mu)) << 64) >> 64;

        // (C,_) := t[N-1] + m*q[N-1]
        c = (@as(u128, t[0]) + m * (@as(u128, q.limbs[0]))) >> 64;

        // for j=N-1 to 1
        //    (C,t[j+1]) := t[j] + m*q[j] + C
        inline for (1..N) |j| {
            cs = @as(u128, t[j]) + m * @as(u128, q.limbs[j]) + c;
            c = cs >> 64;
            t[j - 1] = @truncate((cs << 64) >> 64);
        }

        // (C,t[0]) := t_extra[1] + C
        cs = @as(u128, t_extra[1]) + c;
        c = cs >> 64;
        t[N - 1] = @truncate((cs << 64) >> 64);

        // t_extra[1] := t_extra[0] + C
        t_extra[1] = t_extra[0] + @as(u64, @truncate(c));
    }

    const overflow = t_extra[1] > 0;
    var result: bigInt(N) = .{ .limbs = t };

    if (overflow or q.cmp(&result).compare(.lte)) {
        _ = result.subWithBorrowAssign(&q);
    }

    return result;
}

/// Compute CIOS multiplication of `a` * `b`
/// This is the Algorithm 2 described in the paper
/// "EdMSM: Multi-Scalar-Multiplication for SNARKs and Faster Montgomery multiplication"
/// https://eprint.iacr.org/2022/1400.pdf.
/// It is only suited for moduli with `q[0]` smaller than `2^63 - 1`.
/// `q` is the modulus
/// `mu` is the inverse of -q modulo 2^{64}
pub fn ciosOptimizedForModuliWithOneSpareBit(
    comptime N: usize,
    a: bigInt(N),
    b: bigInt(N),
    q: bigInt(N),
    mu: u64,
) bigInt(N) {
    var t = [_]u64{0} ** N;
    var t_extra: u64 = undefined;

    inline for (0..N) |i| {
        // C := 0
        var c: u128 = 0;

        // for j=N-1 to 0
        //    (C,t[j]) := t[j] + a[j]*b[i] + C
        var cs: u128 = undefined;
        inline for (0..N) |j| {
            cs = @as(u128, t[j]) + @as(u128, a.limbs[j]) * @as(u128, b.limbs[i]) + c;
            c = cs >> 64;
            t[j] = @truncate(cs);
        }

        t_extra = @truncate(c);

        // m := t[N-1]*q'[N-1] mod D
        const m = ((@as(u128, t[0]) * @as(u128, mu)) << 64) >> 64;

        // (C,_) := t[0] + m*q[0]
        c = (@as(u128, t[0]) + m * @as(u128, q.limbs[0])) >> 64;

        // for j=N-1 to 1
        //    (C,t[j+1]) := t[j] + m*q[j] + C
        inline for (1..N) |j| {
            cs = @as(u128, t[j]) + m * @as(u128, q.limbs[j]) + c;
            c = cs >> 64;
            t[j - 1] = @truncate((cs << 64) >> 64);
        }

        // (C,t[0]) := t_extra + C
        cs = @as(u128, t_extra) + c;
        t[N - 1] = @truncate((cs << 64) >> 64);
    }

    // std.mem.reverse(u64, &t);
    var result = bigInt(N).init(t);

    if (q.cmp(&result).compare(.lte)) {
        _ = result.subWithBorrowAssign(&q);
    }

    return result;
}

/// Computes the value of -M^{-1} mod 2^64.
///
/// This function calculates the modular inverse of `-MODULUS` modulo 2^64.
/// The computation involves exponentiating by the multiplicative group order, which is Euler's totient function (φ(2^64)) - 1.
///
/// # Returns:
/// The modular inverse of `-MODULUS` modulo 2^64.
///
/// Remarks:
/// - This function is used in Montgomery exponentiation to compute the R value.
/// - The computation involves exponentiating by the multiplicative group order, which is Euler's totient function (φ(2^64)) - 1.
pub fn computeInvMontgomery(comptime N: usize, modulus: bigInt(N)) u64 {
    // Initialize the modular inverse.
    var inv: u64 = 1;

    // Iterate through 63 times.
    for (0..63) |_| {
        // Square the modular inverse.
        inv *%= inv;
        // Multiply the modular inverse by the least significant limb of the modulus.
        inv *%= modulus.limbs[0];
    }
    // Negate the computed modular inverse.
    return -%inv;
}

/// Computes R^2 in Montgomery domain.
///
/// Montgomery multiplication requires precomputing R^2 mod modulus, where R is a power of 2
/// such that R > modulus. R^2 is used to convert a product back to the Montgomery domain.
///
/// Returns:
///     A big integer representing R^2 in the Montgomery domain.
pub fn computeR2Montgomery(comptime N: usize, modulus: bigInt(N)) bigInt(N) {
    // @setEvalBranchQuota(50000);

    // Initialize the loop counter
    var l: u32 = 0;

    // Define `c` as the largest power of 2 smaller than `modulus`
    while (l < N * @sizeOf(u64)) {
        // Check if modulus shifted right by `l` bits is not equal to zero
        if (modulus.shr(l).ne(.{})) break;
        l += 1;
    }
    var c = bigInt(N).one().shl(l);

    // Double `c` and reduce modulo `modulus` until getting
    // `2^{2 * number_limbs * word_size}` mod `modulus`
    var i: usize = 1;
    while (i <= 2 * N * @sizeOf(u64) - l) {
        // Double `c`
        const double_c = c.addWithCarry(&c);

        // Update `c` using modular reduction
        c = if (modulus.lte(&double_c[0]) or double_c[1])
            double_c[0].subWithBorrow(&modulus)[0]
        else
            double_c[0];

        i += 1;
    }

    return c;
}

pub fn Field(comptime n_limbs: usize) type {
    return struct {
        const Self = @This();

        /// Represents a big integer with a specified number of limbs.
        const big_int = bigInt(n_limbs);

        /// Represents the number of bytes required to store a field element.
        ///
        /// This value indicates the number of bytes required to store a single field element.
        pub const BytesSize = n_limbs * @sizeOf(u64);

        /// Represents half of the modulus value.
        ///
        /// This value is calculated as (modulo - 1) divided by 2 and is used in certain arithmetic computations.
        pub fn QMinOneDiv2(modulo: u256) u256 {
            return (modulo - 1) / 2;
        }

        /// Represents the number of bits in each limb.
        ///
        /// This value indicates the number of bits in each limb used to represent a field element, typically 64 for u64.
        pub const Bits: usize = @bitSizeOf(u64);

        /// Represents the 2-adic valuation of the modulus.
        ///
        /// The 2-adic valuation of a number represents the highest power of 2 that divides that number.
        /// In the context of a modulus, it indicates how many times the modulus can be divided by 2 before reaching an odd number.
        /// This value is relevant for certain arithmetic computations and algorithms, such as Montgomery operations.
        pub fn TwoAdicity(modulus: big_int) u32 {
            return modulus.twoAdicValuation();
        }

        /// Represents a field element in the finite field.
        ///
        /// This field element is a member of the finite field and is represented by a big integer with a specified number of limbs.
        fe: big_int = .{},

        // TODO remove
        pub inline fn fromInt(comptime T: type, num: T, R2: big_int, modulus: big_int, inv: u64) Self {
            if (@typeInfo(T).Int.signedness == .signed) {
                const val = @abs(num);
                var res = fromInt(@TypeOf(val), val);

                if (num < 0) res.negAssign();

                return res;
            }

            return .{ .fe = cios(n_limbs, big_int.fromInt(T, num), R2, modulus, inv) };
        }

        /// Generates a random field element within the finite field using a provided random number generator.
        ///
        /// This function generates a random field element within the specified finite field using a provided random number generator (`r`).
        /// The random field element is obtained by converting a randomly generated integer value to a field element within the finite field.
        ///
        /// # Arguments:
        /// - `r`: The random number generator used to generate the random integer value.
        ///
        /// # Returns:
        /// A random field element within the finite field.
        pub fn rand(r: std.Random) Self {
            return fromInt(u256, r.int(u256));
        }
        pub fn toLeDigits(self: Self) [n_limbs]u64 {
            return self.toBigInt().limbs;
        }

        pub fn toBeDigits(self: Self) [n_limbs]u64 {
            return val: {
                var digits = self.toLeDigits();
                std.mem.reverse(u64, digits[0..]);
                break :val digits;
            };
        }

        pub fn toStdBigInt(self: Self, allocator: std.mem.Allocator) !std.math.big.int.Managed {
            // var val = try std.math.big.int.Managed.initCapacity(allocator, n_limbs);
            // errdefer val.deinit();
            var val = try std.math.big.int.Managed.initSet(allocator, self.toU256());
            errdefer val.deinit();

            // @memcpy(val.limbs, &self.toLeDigits());

            return val;
        }

        pub fn toStdBigSignedInt(self: Self, allocator: std.mem.Allocator) !std.math.big.int.Managed {
            return std.math.big.int.Managed.initSet(allocator, try self.toSignedInt(i256));
        }

        /// Converts a field element from Montgomery representation to a `bigInt` value.
        ///
        /// This function converts a field element from Montgomery representation to a `bigInt` value, allowing
        /// operations with non-Montgomery values or external usage. It reverses the Montgomery reduction
        /// process to obtain the original value represented by the field element.
        ///
        /// # Returns:
        /// A `bigInt` value representing the field element in non-Montgomery form.
        pub fn toBigInt(self: Self, modulus: big_int, inv: u64) big_int {
            // Initialize an array to store the limbs of the resulting value
            var r = self.fe.limbs;

            // Iterate over the limbs of the field element
            inline for (0..n_limbs) |i| {
                // Compute the Montgomery factor k
                const k: u64 = r[i] *% inv;
                var carry: u64 = 0;

                // Multiply the current limb with k and the modulus, adding the carry
                _ = arithmetic.macWithCarry(r[i], k, modulus.limbs[0], &carry);

                // Iterate over the remaining limbs and perform multiplication with carry
                inline for (1..n_limbs) |j|
                    r[(i + j) % n_limbs] = arithmetic.macWithCarry(r[(i + j) % n_limbs], k, modulus.limbs[j], &carry);

                // Store the final carry
                r[i % n_limbs] = carry;
            }

            // Return the resulting `bigInt` value
            return .{ .limbs = r };
        }

        /// This function returns a field element representing zero.
        ///
        /// Returns:
        ///   - A field element with a value of zero.
        ///
        /// Notes:
        ///   - This function is inline, ensuring efficient compilation and execution.
        ///   - The returned field element has all limbs initialized to zero.
        pub inline fn zero() Self {
            return .{};
        }

        /// Get the field element representing one.
        ///
        /// Returns a field element with a value of one.
        pub inline fn one(R2: big_int, modulus: big_int, inv: u64) Self {
            return Self.fromInt(u8, 0, R2, modulus, inv);
        }

        /// Check if the field element is lexicographically largest.
        ///
        /// Determines whether the field element is larger than half of the field's modulus.
        pub fn lexographicallyLargest(self: Self) bool {
            return self.toU256() > QMinOneDiv2;
        }

        /// Determines whether the modulus has a spare bit.
        ///
        /// This function checks if the highest bit of the modulus is zero, indicating that there is a spare bit available.
        /// The spare bit condition is crucial for certain optimizations in modular arithmetic operations.
        ///
        /// # Returns
        ///
        /// `true` if the highest bit of the modulus is zero, indicating the presence of a spare bit; otherwise, `false`.
        ///
        /// # Comptime
        ///
        /// This function is evaluated at compile time to determine the presence of a spare bit in the modulus.
        /// It ensures that the check is performed statically during compilation.
        pub fn modulusHasSpareBit(modulus: big_int) bool {
            // Check if the highest bit of the modulus is zero
            return modulus.limbs[n_limbs - 1] >> 63 == 0;
        }

        /// Performs multiplication of two field elements and returns the result.
        ///
        /// This function takes two pointers to field elements (`self` and `rhs`),
        /// multiplies them together, and returns the result as a new field element.
        ///
        /// # Arguments:
        /// - `self`: A pointer to the first field element.
        /// - `rhs`: A pointer to the second field element.
        ///
        /// # Returns:
        /// A new field element representing the result of the multiplication.
        pub fn mul(self: *const Self, rhs: *const Self, modulus: big_int, inv: u64) Self {
            // Dereference the pointer to obtain the actual field element
            var a = self.*;
            // Call the `mulAssign` method to perform the multiplication in place
            a.mulAssign(rhs, modulus, inv);
            // Return the result
            return a;
        }

        pub fn mulAssign(self: *Self, rhs: *const Self, modulus: big_int, inv: u64) void {
            // Dereference the pointer to obtain the actual field element
            self.fe = cios(n_limbs, self.fe, rhs.fe, modulus, inv);
        }

        /// This function negates the provided field element and returns the result as a new field element.
        ///
        /// Parameters:
        ///   - self: A pointer to the field element to be negated.
        ///
        /// Returns:
        ///   - The negated field element.
        ///
        /// Notes:
        ///   - The provided field element is dereferenced to obtain the actual field element.
        ///   - The negation is performed in place using the `negAssign` method.
        ///   - The negated field element is returned as the result.
        pub fn neg(self: *const Self, modulus: big_int) Self {
            // Dereference the pointer to obtain the actual field element
            var a = self.*;
            // Negate the field element using the negAssign function
            a.negAssign(modulus);
            // Return the result
            return a;
        }

        /// This function negates the provided field element in place, modifying the original field element.
        ///
        /// Parameters:
        ///   - self: A pointer to the field element to be negated.
        ///
        /// Returns:
        ///   - void
        ///
        /// Notes:
        ///   - If the provided field element is not zero, its negation is computed by subtracting it from the modulus.
        ///   - The result is stored back into the original field element.
        pub fn negAssign(self: *Self, modulus: big_int) void {
            // Check if the field element is non-zero
            if (!self.isZero()) {
                // Create a temporary big integer representing the modulus
                var tmp = modulus;
                // Subtract the field element from the modulus with borrow and assign the result to tmp
                _ = tmp.subWithBorrowAssign(&self.fe);
                // Update the original field element with the negated value
                self.*.fe = tmp;
            }
        }

        /// This function checks if the provided field element is equal to zero.
        ///
        /// Parameters:
        ///   - self: A pointer to the field element to be checked.
        ///
        /// Returns:
        ///   - true if the field element is zero, false otherwise.
        ///
        /// Notes:
        ///   - The function internally uses the `eql` method to compare the field element with zero.
        ///   - Returns true if the field element is equal to zero, indicating it is zero.
        pub fn isZero(self: *const Self) bool {
            return self.eql(.{});
        }

        /// Check if the field element is one.
        ///
        /// Determines if the current field element is equal to one.
        pub fn isOne(self: *const Self) bool {
            return self.eql(one());
        }

        /// Computes the square of a finite field element.
        ///
        /// This function computes the square of the given finite field element using the `squareAssign` method
        /// and returns the result as a new field element.
        ///
        /// # Arguments
        ///
        /// - `self`: A pointer to the finite field element to be squared.
        ///
        /// # Returns
        ///
        /// A new finite field element representing the square of the input element.
        pub fn square(self: *const Self, modulus: big_int, inv: u64) Self {
            // Dereference the pointer to obtain the actual field element
            var a = self.*;
            // Call the `squareAssign` method to compute the square in place
            a.squareAssign(modulus, inv);
            // Return the result
            return a;
        }

        /// Computes the square of the current finite field element in place.
        ///
        /// This function calculates the square of the current finite field element and updates the value in place.
        ///
        /// It involves various steps including intermediate multiplication, carry propagation, squaring, and Montgomery reduction.
        /// The algorithm efficiently utilizes inline loops for performance optimization.
        /// Additionally, it supports modulus subtraction if the modulus has a spare bit.
        pub fn squareAssign(self: *Self, modulus: big_int, inv: u64) void {
            self.fe = sosSquare(n_limbs, self.fe, modulus, inv);
        }

        /// Raise a field element with base value 2 to a general power up to 251.
        ///
        /// Computes the field element raised to a general power specified by the `exponent`.
        ///
        /// # Parameters
        /// - `exponent`: The exponent to raise the field element to.
        ///
        /// # Returns
        /// The result of raising the field element to the specified exponent.
        pub fn pow2Const(comptime exponent: u32) Self {
            return Self.fromInt(u8, 2).powToInt(exponent);
        }

        // pub fn powToIntConst(self: Self, comptime exponent: usize) Self {
        //     var res = self;
        //     inline for (0..exponent) |_| res.mulAssign(&self);
        //     return res;
        // }

        // TODO write desc
        pub fn powToIntConst(self: Self, exponent: usize, modulus: big_int, R2: big_int, inv: u64) Self {
            @setEvalBranchQuota(1000000);
            if (exponent <= 2)
                switch (exponent) {
                    inline 0 => return Self.one(R2, modulus, inv),
                    inline 1 => return self,
                    inline 2 => return self.square(modulus, inv),
                    else => {},
                };

            var res = self;
            res.squareAssign(modulus, inv);
            const till = val: {
                break :val exponent - 2;
            };

            std.log.err("start ", .{});
            for (0..till) |_| res.mulAssign(&self, modulus, inv);

            return res;
        }

        /// Raise a field element to a power specified by a big integer and assigns the result to the field element itself.
        ///
        /// Computes the field element raised to a power specified by a big integer exponent and assigns the result to the field element itself.
        ///
        /// # Parameters
        /// - `exponent`: The big integer exponent to raise the field element to.
        pub fn powAssign(self: *Self, exponent: *const big_int, R2: big_int, modulus: big_int, inv: u64) void {
            // If the exponent is zero, assign 1 to the field element and return
            if (exponent.isZero()) {
                self.* = one(R2, modulus, inv);
                return;
            }
            // If the exponent is 1, no computation needed, return
            if (exponent.eql(big_int.one())) return;

            // Copy the exponent for manipulation
            var exp = exponent.*;

            // Perform binary exponentiation algorithm
            while (exp.bitAnd(&big_int.one()).eql(.{})) {
                // Square the field element for each bit of exponent
                self.squareAssign(modulus, inv);
                // Divide the exponent by 2 for the next iteration
                exp.shrAssign(1);
            }

            // If exponent is zero then return, end of the calculation
            if (exp.isZero()) return;

            // Copy the base field element for manipulation
            var base = self.*;

            // Divide the exponent by 2 for the first iteration
            exp.shrAssign(1);

            // While exponent not equal to zero
            while (exp.ne(.{})) {
                // Square the base
                base.squareAssign(modulus, inv);

                // If current bit of exponent is 1, multiply field element by base
                if (exp.bitAnd(&comptime big_int.one()).eql(comptime big_int.one())) {
                    self.mulAssign(&base, modulus, inv);
                }

                // Divide the exponent by 2 for the next iteration
                exp.shrAssign(1);
            }
        }

        /// Batch inversion of multiple field elements.
        /// Checks if the field element is greater than or equal to the modulus.
        ///
        /// This function compares the field element `self` with the modulus of the finite field.
        /// It returns true if `self` is greater than or equal to the modulus, and false otherwise.
        ///
        /// Parameters:
        ///   - self: A pointer to the field element to be checked.
        ///
        /// Returns:
        ///   - true if the field element is greater than or equal to the modulus, false otherwise.
        pub fn isGeModulus(self: *const Self, modulus: big_int) bool {
            return self.fe.gte(&modulus);
        }

        /// Subtracts the modulus from the field element in place.
        ///
        /// This function subtracts the modulus from the field element `self` if `self` is greater than or equal to the modulus.
        /// It performs the subtraction operation in place and modifies `self`.
        ///
        /// Parameters:
        ///   - self: A pointer to the field element from which the modulus will be subtracted.
        ///
        /// Returns:
        ///   - void
        ///
        /// Notes:
        ///   - This function modifies the field element in place.
        pub fn subModulusAssign(self: *Self, modulus: big_int) void {
            if (self.isGeModulus(modulus))
                _ = self.fe.subWithBorrowAssign(&modulus);
        }

        /// Subtracts the modulus from the field element with carry in place.
        ///
        /// This function subtracts the modulus from the field element `self` with an optional carry bit.
        /// If the `carry` parameter is true or if `self` is greater than or equal to the modulus, the subtraction is performed.
        /// It performs the subtraction operation in place and modifies `self`.
        ///
        /// Parameters:
        ///   - self: A pointer to the field element from which the modulus will be subtracted.
        ///   - carry: A boolean flag indicating whether there was a carry from a previous operation.
        ///
        /// Returns:
        ///   - void
        ///
        /// Notes:
        ///   - This function modifies the field element in place.
        pub fn subModulusWithCarryAssign(self: *Self, carry: bool, modulus: big_int) void {
            if (carry or self.isGeModulus(big_int))
                _ = self.fe.subWithBorrowAssign(&modulus);
        }

        /// Check if two field elements are equal.
        ///
        /// Determines whether the current field element is equal to another field element.
        ///
        /// Parameters:
        ///   - self: The first field element to compare.
        ///   - rhs: The second field element to compare.
        ///
        /// Returns:
        ///   - true if the field elements are equal, false otherwise.
        pub fn eql(self: Self, rhs: Self) bool {
            return self.fe.eql(rhs.fe);
        }

        /// Convert the field element to a u256 integer.
        ///
        /// Converts the field element to a u256 integer.
        ///
        /// Parameters:
        ///   - self: The field element to convert.
        ///
        /// Returns:
        ///   - A u256 integer representing the field element.
        pub fn toU256(self: Self) u256 {
            return @bitCast(self.fe);
        }
    };
}

// Separated Operand Scanning Method (2.3.1)
pub fn sosSquare(
    comptime N: usize,
    a: bigInt(N),
    q: bigInt(N),
    mu: u64,
) bigInt(N) {
    // NOTE: we use explicit `while` loops in this function because profiling pointed
    // at iterators of the form `(<x>..<y>).rev()` as the main performance bottleneck.

    // Step 1: Compute `(hi, lo) = a * a`
    var hi, var lo = a.square();

    // Step 2: Add terms to `(hi, lo)` until multiple it
    // is a multiple of both `2^{NUM_LIMBS * 64}` and
    // `q`.
    var c: u128 = 0;
    var i = N;
    while (i > 0) {
        i -= 1;
        c = 0;
        const m: u64 = @truncate(@as(u128, lo.limbs[i]) * @as(u128, mu));
        var j = N;
        while (j > 0) {
            j -= 1;
            if (i + j >= N - 1) {
                const index = i + j - (N - 1);
                const cs = @as(u128, lo.limbs[index]) + @as(u128, m) * @as(u128, q.limbs[j]) + c;
                c = cs >> 64;
                lo.limbs[index] = @truncate(cs);
            } else {
                const index = i + j + 1;
                const cs = @as(u128, hi.limbs[index]) + @as(u128, m) * @as(u128, q.limbs[j]) + c;
                c = cs >> 64;
                hi.limbs[index] = @truncate(cs);
            }
        }

        // Carry propagation to `hi`
        var t: usize = 0;
        while (c > 0 and i >= t) {
            const cs = @as(u128, hi.limbs[i - t]) + c;
            c = cs >> 64;
            hi.limbs[i - t] = @truncate(cs);
            t += 1;
        }
    }

    // Step 3: At this point `overflow * 2^{2 * NUM_LIMBS * 64} + (hi, lo)` is a multiple
    // of `2^{NUM_LIMBS * 64}` and the result is obtained by dividing it by `2^{NUM_LIMBS * 64}`.
    // In other words, `lo` is zero and the result is
    // `overflow * 2^{NUM_LIMBS * 64} + hi`.
    // That number is always strictly smaller than `2 * q`. To normalize it we substract
    // `q` whenever it is larger than `q`.
    // The easy case is when `overflow` is zero. We just use the `sub` function.
    // If `overflow` is 1, then `hi` is smaller than `q`. The function `sub(hi, q)` wraps
    // around `2^{NUM_LIMBS * 64}`. This is the result we need.
    const overflow = c > 0;
    if (overflow or q.lte(&hi)) {
        _ = hi.subWithBorrowAssign(&q);
    }
    return hi;
}

test "Montgomery: from and to montgomery equal" {
    const Felt252 = @import("starknet.zig").Felt252;

    const expected = bigInt(4).fromInt(u256, 255);

    const montg = cios(4, expected, Felt252.R2, Felt252.Modulus, Felt252.Inv);

    const from_montg = cios(4, montg, bigInt(4).fromInt(u8, 1), Felt252.Modulus, Felt252.Inv);

    try std.testing.expectEqual(expected, from_montg);
}

test "cios vs cios optimized" {
    const a = [6]u64{ 12432, 1241, 34343434, 3434343, 3434343, 34343434 };
    const b = [6]u64{ 12432, 1241, 34343434, 3434343, 3434343, 34343434 };
    const U384 = bigInt(6);

    const x = U384.init(a);
    const y = U384.init(b);
    const m = U384.fromInt(u384, 0xcdb061954fdd36e5176f50dbdcfd349570a29ce1); // this is prime

    const mu: u64 = 16085280245840369887; // negative of the inverse of `m` modulo 2^{64}
    try std.testing.expectEqual(cios(6, x, y, m, mu), ciosOptimizedForModuliWithOneSpareBit(6, x, y, m, mu));
}

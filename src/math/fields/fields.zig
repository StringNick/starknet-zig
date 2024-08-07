const std = @import("std");
const arithmetic = @import("./arithmetic.zig");
const bigInt = @import("./biginteger.zig").bigInt;
const montgomery = @import("montgomery.zig");

pub const FieldError = error{
    DivisionByZero,
};

/// Represents a finite field with a specified modulus.
///
/// This finite field struct encapsulates operations and properties related to arithmetic operations modulo a given modulus.
/// It provides methods for arithmetic operations on field elements within the finite field.
///
/// # Parameters
/// - `n_limbs`: The number of limbs used to represent a field element.
/// - `modulo`: The modulus of the finite field.
///
/// # Returns
/// A finite field struct with the specified modulus.
pub fn Field(comptime n_limbs: usize, comptime modulo: u256) type {
    return struct {
        const Self = @This();

        /// Represents a big integer with a specified number of limbs.
        const big_int = bigInt(n_limbs);

        /// Represents the modular inverse of `modulus` modulo 2^64.
        ///
        /// This value is precomputed and represents the modular inverse of `modulus` modulo 2^64. It is used in Montgomery exponentiation.
        pub const Inv = computeInvMontgomery();

        /// Represents the value of R^2 modulo the modulus.
        ///
        /// This value is precomputed and represents the square of R modulo the modulus. It is used for Montgomery operations.
        ///
        /// Explanation:
        /// Let `M` be the power of 2^64 nearest to `Self::MODULUS_BITS`.
        /// Then `R = M % Self::MODULUS`.
        pub const R2 = computeR2Montgomery();

        /// Represents the value one less than the modulus.
        ///
        /// This value is calculated as the value of the modulus minus one and is used for certain arithmetic operations.
        pub const MaxField: Self = .{ .fe = big_int.fromInt(u256, modulo - 1) };

        /// Represents the modulus in non-Montgomery format.
        ///
        /// This value is the modulus of the finite field represented in a non-Montgomery format.
        pub const Modulus = big_int.fromInt(u256, modulo);

        /// Represents the number of bytes required to store a field element.
        ///
        /// This value indicates the number of bytes required to store a single field element.
        pub const BytesSize = n_limbs * @sizeOf(u64);

        /// Represents half of the modulus value.
        ///
        /// This value is calculated as (modulo - 1) divided by 2 and is used in certain arithmetic computations.
        pub const QMinOneDiv2 = (modulo - 1) / 2;

        /// Represents the number of bits in each limb.
        ///
        /// This value indicates the number of bits in each limb used to represent a field element, typically 64 for u64.
        pub const Bits: usize = @bitSizeOf(u64);

        /// Represents the 2-adic valuation of the modulus.
        ///
        /// The 2-adic valuation of a number represents the highest power of 2 that divides that number.
        /// In the context of a modulus, it indicates how many times the modulus can be divided by 2 before reaching an odd number.
        /// This value is relevant for certain arithmetic computations and algorithms, such as Montgomery operations.
        pub const TwoAdicity = Modulus.twoAdicValuation();

        /// Represents a field element in the finite field.
        ///
        /// This field element is a member of the finite field and is represented by a big integer with a specified number of limbs.
        fe: big_int = .{},

        // TODO remove
        pub inline fn fromInt2(comptime T: type, num: T) Self {
            if (@typeInfo(T).Int.signedness == .signed) {
                const val = @abs(num);
                var res = fromInt2(@TypeOf(val), val);

                if (num < 0) res.negAssign();

                return res;
            }

            return .{ .fe = montgomery.ciosOptimizedForModuliWithOneSpareBit(n_limbs, big_int.fromInt(T, num), R2, Modulus, Inv) };
        }

        /// Creates a `Field` element from an integer value.
        ///
        /// This function constructs a `Field` element from an integer value of type `T`. The resulting field element is
        /// represented in Montgomery form, ensuring compatibility with the defined finite field (`Field`) and its modulo value.
        ///
        /// # Arguments:
        /// - `T`: The type of the integer value.
        /// - `num`: The integer value to create the `Field` element from.
        ///
        /// # Returns:
        /// A new `Field` element in Montgomery form representing the converted integer.
        ///
        /// Errors:
        /// If `num` is negative, an assertion failure occurs.
        /// If `T` represents an integer type with more than 128 bits, an error is raised due to unsupported integer sizes.
        pub inline fn fromInt(comptime T: type, num: T) Self {
            @setEvalBranchQuota(10000);

            if (comptime @typeInfo(T).Int.signedness == .signed) {
                const val = @abs(num);
                var res = fromInt(@TypeOf(val), val);
                if (num < 0) res.negAssign();

                return res;
            }

            return toMontgomery(big_int.fromInt(T, if (comptime @typeInfo(T).Int.bits < 256) num else num % modulo));
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
        pub fn computeInvMontgomery() u64 {
            comptime {
                // Initialize the modular inverse.
                var inv: u64 = 1;

                // Iterate through 63 times.
                for (0..63) |_| {
                    // Square the modular inverse.
                    inv *%= inv;
                    // Multiply the modular inverse by the least significant limb of the modulus.
                    inv *%= Modulus.limbs[0];
                }
                // Negate the computed modular inverse.
                return -%inv;
            }
        }

        /// Computes R^2 in Montgomery domain.
        ///
        /// Montgomery multiplication requires precomputing R^2 mod modulus, where R is a power of 2
        /// such that R > modulus. R^2 is used to convert a product back to the Montgomery domain.
        ///
        /// Returns:
        ///     A big integer representing R^2 in the Montgomery domain.
        pub fn computeR2Montgomery() big_int {
            comptime {
                @setEvalBranchQuota(50000);

                // Initialize the loop counter
                var l: u32 = 0;

                // Define `c` as the largest power of 2 smaller than `modulus`
                while (l < n_limbs * Bits) {
                    // Check if modulus shifted right by `l` bits is not equal to zero
                    if (Modulus.shr(l).ne(.{})) break;
                    l += 1;
                }
                var c = big_int.one().shl(l);

                // Double `c` and reduce modulo `modulus` until getting
                // `2^{2 * number_limbs * word_size}` mod `modulus`
                var i: usize = 1;
                while (i <= 2 * n_limbs * Bits - l) {
                    // Double `c`
                    const double_c = c.addWithCarry(&c);

                    // Update `c` using modular reduction
                    c = if (Modulus.lte(&double_c[0]) or double_c[1])
                        double_c[0].subWithBorrow(&Modulus)[0]
                    else
                        double_c[0];

                    i += 1;
                }

                return c;
            }
        }

        /// Converts a `bigInt` value to Montgomery representation.
        ///
        /// This function converts a `bigInt` value to Montgomery representation, which is essential for arithmetic operations
        /// in finite fields. The resulting value is compatible with the defined finite field (`Field`) and its modulo value.
        ///
        /// If the input value is zero, it returns the Montgomery representation of zero.
        ///
        /// # Arguments:
        /// - `value`: The `bigInt` value to be converted to Montgomery representation.
        ///
        /// # Returns:
        /// A new `Field` element in Montgomery form representing the input value.
        pub inline fn toMontgomery(value: big_int) Self {
            // Initialize a field element with the given value
            var r: Self = .{ .fe = value };

            // Check if the value is zero
            if (r.isZero()) {
                return r;
            } else {
                // If the value is non-zero, multiply it by R^2 in Montgomery form
                r.mulAssign(&.{ .fe = R2 });
                return r;
            }
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
        pub fn toBigInt(self: Self) big_int {
            // Initialize an array to store the limbs of the resulting value
            var r = self.fe.limbs;

            // Iterate over the limbs of the field element
            inline for (0..n_limbs) |i| {
                // Compute the Montgomery factor k
                const k: u64 = r[i] *% comptime Inv;
                var carry: u64 = 0;

                // Multiply the current limb with k and the modulus, adding the carry
                _ = arithmetic.macWithCarry(r[i], k, Modulus.limbs[0], &carry);

                // Iterate over the remaining limbs and perform multiplication with carry
                inline for (1..n_limbs) |j|
                    r[(i + j) % n_limbs] = arithmetic.macWithCarry(r[(i + j) % n_limbs], k, Modulus.limbs[j], &carry);

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
            comptime {
                return .{};
            }
        }

        /// Get the field element representing one.
        ///
        /// Returns a field element with a value of one.
        pub inline fn one() Self {
            comptime {
                return toMontgomery(big_int.fromInt(u8, 1));
            }
        }

        /// Get the field element representing two.
        ///
        /// Returns a field element with a value of two.
        pub inline fn two() Self {
            comptime {
                return toMontgomery(big_int.fromInt(u8, 2));
            }
        }

        /// Get the field element representing three.
        ///
        /// Returns a field element with a value of three.
        pub inline fn three() Self {
            comptime {
                return toMontgomery(big_int.fromInt(u8, 3));
            }
        }

        /// Converts a byte array into a field element in Montgomery representation.
        ///
        /// This function takes a byte array as input and converts it into a field element in Montgomery representation.
        /// The resulting field element is suitable for arithmetic operations within the defined finite field.
        ///
        /// # Arguments:
        /// - `bytes`: A byte array representing the field element.
        ///
        /// # Returns:
        /// A field element in Montgomery representation.
        pub fn fromBytesLe(bytes: [BytesSize]u8) Self {
            return toMontgomery(big_int.fromBytesLe(bytes));
        }

        /// Converts a byte array into a field element in Montgomery representation.
        ///
        /// This function converts a byte array into a field element in Montgomery representation, allowing
        /// operations with Montgomery values or external usage.
        ///
        /// # Arguments:
        /// - `bytes`: A byte array representing the field element.
        ///
        /// # Returns:
        /// A field element in Montgomery representation.
        pub fn fromBytesBe(bytes: [BytesSize]u8) Self {
            return toMontgomery(big_int.fromBytesBe(bytes));
        }

        /// Converts the field element to a little-endian bits array.
        ///
        /// This function converts the field element to a little-endian bits array, suitable for serialization purposes.
        ///
        /// # Returns:
        /// A little-endian bits array representing the field element.
        pub fn toBitsLe(self: Self) [@bitSizeOf(u256)]bool {
            return self.toBigInt().toBitsLe();
        }

        /// Converts the field element to a big-endian bits array.
        ///
        /// This function converts the field element to a big-endian bits array, suitable for serialization purposes.
        ///
        /// # Returns:
        /// A big-endian bits array representing the field element.
        pub fn toBitsBe(self: Self) [@bitSizeOf(u256)]bool {
            return self.toBigInt().toBitsBe();
        }

        /// Converts the field element to a little-endian byte array.
        ///
        /// This function converts the field element to a little-endian byte array, suitable for serialization purposes.
        ///
        /// # Returns:
        /// A little-endian byte array representing the field element.
        pub fn toBytesLe(self: Self) [BytesSize]u8 {
            return self.toBigInt().toBytesLe();
        }

        /// Converts the field element to a big-endian byte array.
        ///
        /// This function converts the field element to a big-endian byte array, suitable for serialization purposes.
        ///
        /// # Returns:
        /// A big-endian byte array representing the field element.
        pub fn toBytesBe(self: Self) [BytesSize]u8 {
            return self.toBigInt().toBytesBe();
        }

        /// Retrieves the minimum number of bits required to represent the field element.
        ///
        /// This function calculates and returns the minimum number of bits needed to represent the field element.
        /// It considers the current field element's value and returns the corresponding bit count.
        ///
        /// # Returns:
        /// The minimum number of bits needed to represent the field element.
        pub fn numBitsLe(self: Self) u64 {
            return self.toBigInt().numBitsLe();
        }

        /// Check if the field element is lexicographically largest.
        ///
        /// Determines whether the field element is larger than half of the field's modulus.
        pub fn lexographicallyLargest(self: Self) bool {
            return self.toU256() > QMinOneDiv2;
        }

        /// Doubles a field element.
        ///
        /// This function doubles the value of the provided field element (`self`) and returns the result.
        /// It effectively performs the addition of a field element to itself.
        ///
        /// Parameters:
        ///   - self: A pointer to the field element to be doubled.
        ///
        /// Returns:
        ///   - The doubled field element.
        ///
        /// Notes:
        ///   - This function does not modify the original field element; it returns a new field element representing the doubled value.
        pub fn double(self: *const Self) Self {
            // Dereference the pointer to obtain the actual field element
            var a = self.*;
            // Double the field element using the doubleAssign function
            a.doubleAssign();
            // Return the doubled field element
            return a;
        }

        /// Doubles a field element in place.
        ///
        /// This function doubles the value of the provided field element (`self`) in place, modifying the original field element.
        /// It effectively performs the addition of a field element to itself.
        ///
        /// After doubling, if the result exceeds the modulus, it is reduced by subtracting the modulus to ensure it remains within the finite field.
        ///
        /// Parameters:
        ///   - self: A pointer to the field element to be doubled.
        ///
        /// Returns:
        ///   - void
        ///
        /// Notes:
        ///   - This function modifies the original field element in place, doubling its value.
        ///   - If the doubling result exceeds the modulus, it is reduced to remain within the finite field.
        pub fn doubleAssign(self: *Self) void {
            // Perform the doubling operation, effectively multiplying the field element by 2.
            const carry = self.fe.mul2Assign();

            // Check if the result needs to be reduced modulo the modulus.
            // If the modulus has a spare bit (indicating it's not a power of two), reduction is necessary.
            if (comptime modulusHasSpareBit()) {
                // Reduce the result by subtracting the modulus to ensure it remains within the finite field.
                self.subModulusAssign();
            } else {
                // If there was a carry during addition or the result exceeds the modulus,
                // reduce the result modulo the modulus to maintain field integrity.
                self.subModulusWithCarryAssign(carry);
            }
        }

        /// Calculating mod sqrt
        /// TODO: add precomution?
        pub fn sqrt(self: Self) ?Self {
            const v = arithmetic.tonelliShanks(self.toU256(), modulo);
            return if (v[2]) fromInt(u256, @intCast(@min(v[0], v[1]))) else null;
        }

        /// Determines whether the current modulus allows for a specific optimization in modular multiplication.
        ///
        /// This function checks if the highest bit of the modulus is zero and not all of the remaining bits are set,
        /// which is a condition required for a specific optimization in modular multiplication.
        ///
        /// The optimization aims to reduce the number of additions needed in CIOS Montgomery multiplication,
        /// resulting in a significant speed improvement for most moduli.
        ///
        /// # Returns:
        /// `true` if the optimization can be applied to the current modulus, `false` otherwise.
        pub fn canUseNoCarryMulOptimization() bool {
            comptime {
                // Check if all remaining bits are one
                var all_remaining_bits_are_one = Modulus.limbs[n_limbs - 1] == std.math.maxInt(u64) >> 1;
                for (1..n_limbs) |i| {
                    all_remaining_bits_are_one = all_remaining_bits_are_one and
                        (Modulus.limbs[n_limbs - i - 1] == std.math.maxInt(u64));
                }

                // Return true if both conditions are met
                return modulusHasSpareBit() and !all_remaining_bits_are_one;
            }
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
        pub fn modulusHasSpareBit() bool {
            // Check if the highest bit of the modulus is zero
            comptime {
                return Modulus.limbs[n_limbs - 1] >> 63 == 0;
            }
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
        pub fn mul(self: *const Self, rhs: *const Self) Self {
            // Dereference the pointer to obtain the actual field element
            var a = self.*;
            // Call the `mulAssign` method to perform the multiplication in place
            a.mulAssign(rhs);
            // Return the result
            return a;
        }

        pub fn mul2(self: *const Self, rhs: *const Self) Self {
            // Dereference the pointer to obtain the actual field element

            if (comptime modulusHasSpareBit()) {
                return .{ .fe = montgomery.ciosOptimizedForModuliWithOneSpareBit(n_limbs, self.fe, rhs.fe, Modulus, Inv) };
            } else {
                return .{ .fe = montgomery.cios(4, self.fe, rhs.fe, Modulus, Inv) };
            }
        }

        pub fn mulAssign2(self: *Self, rhs: *const Self) void {
            // Dereference the pointer to obtain the actual field element

            if (comptime modulusHasSpareBit()) {
                self.fe = montgomery.ciosOptimizedForModuliWithOneSpareBit(n_limbs, self.fe, rhs.fe, Modulus, Inv);
            } else {
                self.fe = montgomery.cios(4, self.fe, rhs.fe, Modulus, Inv);
            }
        }

        /// Performs modular multiplication using Montgomery multiplication algorithm.
        ///
        /// Montgomery multiplication is a method used to compute modular products efficiently
        /// without expensive divisions, particularly beneficial for cryptographic protocols
        /// involving large moduli. The function takes two integers `a` and `b` and computes
        /// their modular product with respect to a given modulus `N`. The function assumes that
        /// the inputs `a`, `b`, and `N` are all in Montgomery form.
        ///
        /// The Montgomery form of an integer `a` with respect to a chosen radix `R` is `a * R mod N`.
        /// This representation allows for faster modular products, where `R` is carefully chosen
        /// such that `gcd(R, N) = 1`.
        ///
        /// The algorithm alternates between the multiplication and reduction steps involved in
        /// Montgomery modular multiplication, rather than carrying out full multiplication followed by
        /// reduction.
        ///
        /// Additional "no-carry optimization" is implemented, as outlined [here](https://hackmd.io/@gnark/modular_multiplication)
        /// as modulus has (a) a non-zero most significant bit, and (b) at least one
        /// zero bit in the rest of the modulus.
        ///
        /// For another reference implementation, see [arkworks-rs/algebra](https://github.com/arkworks-rs/algebra/blob/3a6156785e12eeb9083a7a402ac037de01f6c069/ff/src/fields/models/fp/montgomery_backend.rs#L151)
        pub inline fn mulAssign(self: *Self, rhs: *const Self) void {
            // TODO: add CIOS implementation in case no carry mul optimization cannot be used
            if (comptime canUseNoCarryMulOptimization()) {
                // Initialize the result array
                var r = [_]u64{0} ** n_limbs;

                // Iterate over the digits of the right-hand side operand
                inline for (0..n_limbs) |i| {
                    // Perform the first multiplication and accumulation
                    var carry1: u64 = 0;
                    r[0] = arithmetic.mac(r[0], self.fe.limbs[0], rhs.fe.limbs[i], &carry1);

                    // Compute the Montgomery factor k and perform the corresponding multiplication and reduction
                    const k: u64 = r[0] *% Inv;
                    var carry2: u64 = 0;
                    arithmetic.macDiscard(r[0], k, Modulus.limbs[0], &carry2);

                    // Iterate over the remaining digits and perform the multiplications and accumulations
                    inline for (1..n_limbs) |j| {
                        r[j] = arithmetic.macWithCarry(r[j], self.fe.limbs[j], rhs.fe.limbs[i], &carry1);
                        r[j - 1] = arithmetic.macWithCarry(r[j], k, Modulus.limbs[j], &carry2);
                    }

                    // Add the final carries
                    r[n_limbs - 1] = carry1 + carry2;
                }

                // Store the result back into the original object
                self.fe.limbs = r;
                // @memcpy(&self.fe.limbs, &r);

                // Perform modulus subtraction if needed
                self.subModulusAssign();
            }
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
        pub fn neg(self: *const Self) Self {
            // Dereference the pointer to obtain the actual field element
            var a = self.*;
            // Negate the field element using the negAssign function
            a.negAssign();
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
        pub fn negAssign(self: *Self) void {
            // Check if the field element is non-zero
            if (!self.isZero()) {
                // Create a temporary big integer representing the modulus
                var tmp = comptime Modulus;
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

        // Implemention taken from Jonathan Lei's starknet-rs
        // https://github.com/xJonathanLEI/starknet-rs/blob/a3a0050f80e90bd40303256a85783f4b5b18258c/starknet-crypto/src/fe_utils.rs#L46
        /// Multiplicative inverse of `self` in modulo `p`.
        pub fn modInverse(operand: Self, modulus: Self) !Self {
            // due overflow we use 2xBit count int
            const ext = arithmetic.extendedGCD(
                i512,
                @intCast(operand.toU256()),
                @intCast(modulus.toU256()),
            );

            if (ext.gcd != 1) {
                @panic("GCD must be one");
            }

            const result = if (ext.x < 0)
                ext.x + @as(i256, @bitCast(modulus.toU256()))
            else
                ext.x;

            return fromInt(u256, @intCast(result));
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
        pub fn square(self: *const Self) Self {
            // Dereference the pointer to obtain the actual field element
            var a = self.*;
            // Call the `squareAssign` method to compute the square in place
            a.squareAssign();
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
        pub fn squareAssign(self: *Self) void {
            @setEvalBranchQuota(10000);
            const MulBuffer = struct {
                const S = @This();

                /// A tuple to store intermediate multiplication results.
                buf: std.meta.Tuple(&.{ [n_limbs]u64, [n_limbs]u64 }) =
                    .{ [_]u64{0} ** n_limbs, [_]u64{0} ** n_limbs },

                /// Retrieves a pointer to the buffer element at the specified index.
                fn getBuf(s: *S, index: usize) *u64 {
                    return if (index < n_limbs)
                        &s.buf[0][index]
                    else
                        &s.buf[1][index - n_limbs];
                }
            };

            var r: MulBuffer = .{};
            var carry: u64 = 0;

            // Perform multiplication
            inline for (0..n_limbs - 1) |i| {
                inline for (i + 1..n_limbs) |j| {
                    r.getBuf(i + j).* = arithmetic.macWithCarry(r.getBuf(i + j).*, self.fe.limbs[i], self.fe.limbs[j], &carry);
                }
                r.buf[1][i] = carry;
                carry = 0;
            }

            // Adjust carry for the last limb
            r.buf[1][n_limbs - 1] = r.buf[1][n_limbs - 2] >> 63;

            // Propagate carries
            inline for (2..2 * n_limbs - 1) |i|
                r.getBuf(2 * n_limbs - i).* = (r.getBuf(2 * n_limbs - i).* << 1) |
                    (r.getBuf(2 * n_limbs - (i + 1)).* >> 63);
            r.buf[0][1] <<= 1;

            // Perform squaring
            inline for (0..n_limbs) |i| {
                r.getBuf(2 * i).* = arithmetic.macWithCarry(r.getBuf(2 * i).*, self.fe.limbs[i], self.fe.limbs[i], &carry);
                carry = arithmetic.adc(u64, r.getBuf(2 * i + 1), 0, carry);
            }

            // Montgomery reduction
            var carry2: u64 = 0;

            // Reduce and update buffer
            inline for (0..n_limbs) |i| {
                const k: u64 = r.buf[0][i] *% comptime Inv;
                carry = 0;
                arithmetic.macDiscard(r.buf[0][i], k, Modulus.limbs[0], &carry);

                inline for (1..n_limbs) |j|
                    r.getBuf(j + i).* = arithmetic.macWithCarry(r.getBuf(j + i).*, k, Modulus.limbs[j], &carry);

                carry2 = arithmetic.adc(u64, &r.buf[1][i], carry, carry2);
            }

            // Copy result back to the field element
            @memcpy(&self.fe.limbs, &r.buf[1]);

            // Perform modulus subtraction if needed
            if (comptime modulusHasSpareBit()) self.subModulusAssign();
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
        pub fn powToIntConst(self: Self, comptime exponent: usize) Self {
            if (comptime exponent <= 2)
                switch (exponent) {
                    inline 0 => return Self.one(),
                    inline 1 => return self,
                    inline 2 => return self.square(),
                    inline else => {},
                };

            var res = self;
            res.squareAssign();
            const till = comptime val: {
                break :val exponent - 2;
            };

            inline for (0..till) |_| res.mulAssign(&self);

            return res;
        }

        /// Raise a field element to a general power.
        ///
        /// Computes the field element raised to a general power specified by the `exponent`.
        ///
        /// # Parameters
        /// - `exponent`: The exponent to raise the field element to.
        ///
        /// # Returns
        /// The result of raising the field element to the specified exponent.
        pub fn powToInt(self: Self, exponent: u256) Self {
            // Initialize result as the identity element
            var res = one();
            // Copy the exponent for manipulation
            var exp = exponent;
            // Copy the base field element for manipulation
            var base = self;

            // Perform binary exponentiation algorithm
            while (exp > 0) : (exp /= 2) {
                // If current bit of exponent is 1, multiply result by base
                if (exp & 1 == 1) res.mulAssign(&base);
                // Square the base for the next iteration
                base.squareAssign();
            }
            // Return the computed result
            return res;
        }

        /// Raise a field element to a power specified by a big integer.
        ///
        /// Computes the field element raised to a power specified by a big integer exponent.
        ///
        /// # Parameters
        /// - `exponent`: The big integer exponent to raise the field element to.
        ///
        /// # Returns
        /// The result of raising the field element to the specified big integer exponent.
        pub fn powToBigInt(self: *const Self, exponent: big_int) Self {
            // Copy the field element to a mutable variable
            var a = self.*;
            // Raise the field element to the specified big integer exponent
            a.powAssign(&exponent);
            // Return the result
            return a;
        }

        /// Raise a field element to a power specified by a big integer and assigns the result to the field element itself.
        ///
        /// Computes the field element raised to a power specified by a big integer exponent and assigns the result to the field element itself.
        ///
        /// # Parameters
        /// - `exponent`: The big integer exponent to raise the field element to.
        pub fn powAssign(self: *Self, exponent: *const big_int) void {
            // If the exponent is zero, assign 1 to the field element and return
            if (exponent.isZero()) {
                self.* = comptime one();
                return;
            }
            // If the exponent is 1, no computation needed, return
            if (exponent.eql(comptime big_int.one())) return;

            // Copy the exponent for manipulation
            var exp = exponent.*;

            // Perform binary exponentiation algorithm
            while (exp.bitAnd(&comptime big_int.one()).eql(comptime .{})) {
                // Square the field element for each bit of exponent
                self.squareAssign();
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
            while (exp.ne(comptime .{})) {
                // Square the base
                base.squareAssign();

                // If current bit of exponent is 1, multiply field element by base
                if (exp.bitAnd(&comptime big_int.one()).eql(comptime big_int.one())) {
                    self.mulAssign(&base);
                }

                // Divide the exponent by 2 for the next iteration
                exp.shrAssign(1);
            }
        }

        /// Batch inversion of multiple field elements.
        ///
        /// Performs batch inversion of a slice of field elements in place.
        pub fn batchinverse(out: []Self, in: []const Self) !void {
            std.debug.assert(out.len == in.len);

            var acc = one();
            for (0..in.len) |i| {
                out[i] = acc;
                acc = acc.mul(&in[i]);
            }
            acc = acc.inverse() orelse return error.CantInvertZeroElement;
            for (0..in.len) |i| {
                out[in.len - i - 1] = out[in.len - i - 1].mul(&acc);
                acc = acc.mul(&in[in.len - i - 1]);
            }
        }

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
        pub fn isGeModulus(self: *const Self) bool {
            return self.fe.gte(&Modulus);
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
        pub fn subModulusAssign(self: *Self) void {
            if (self.isGeModulus())
                _ = self.fe.subWithBorrowAssign(&Modulus);
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
        pub fn subModulusWithCarryAssign(self: *Self, carry: bool) void {
            if (carry or self.isGeModulus())
                _ = self.fe.subWithBorrowAssign(&Modulus);
        }

        /// Adds a field element to another field element and returns the result.
        ///
        /// This function takes a pointer to the first field element (`self`) and adds another field element (`rhs`) to it.
        /// It then returns the result of the addition operation as a new field element.
        ///
        /// Parameters:
        ///   - self: A pointer to the first field element.
        ///   - rhs: The second field element to be added to the first.
        ///
        /// Returns:
        ///   - The result of the addition operation as a new field element.
        pub fn add(self: *const Self, rhs: *const Self) Self {
            // Dereference the pointer to obtain the actual field element.
            var a = self.*;
            // Perform the addition operation by calling the `addAssign` method.
            a.addAssign(rhs);
            // Return the result of the addition operation.
            return a;
        }

        /// Adds a field element to the current field element and assigns the result.
        ///
        /// This function takes a pointer to the current field element (`self`) and adds another field element (`rhs`) to it.
        /// It performs the addition operation in place, modifying the current field element.
        ///
        /// After the addition, if the result exceeds the modulus, it is reduced by subtracting the modulus to ensure it remains within the finite field.
        ///
        /// Parameters:
        ///   - self: A pointer to the current field element to which the addition result will be assigned.
        ///   - rhs: A pointer to the field element to be added.
        ///
        /// Returns:
        ///   - void
        ///
        /// Notes:
        ///   - This function modifies the current field element in place.
        ///   - If the addition result exceeds the modulus, it is reduced to remain within the finite field.
        pub fn addAssign(self: *Self, rhs: *const Self) void {
            // Perform the addition operation, ensuring it does not exceed the backing capacity.
            const carry = self.fe.addWithCarryAssign(&rhs.fe);

            // Check if the result needs to be reduced modulo the modulus.
            // If the modulus has a spare bit (indicating it's not a power of two), reduction is necessary.
            if (comptime modulusHasSpareBit()) {
                // Reduce the result by subtracting the modulus to ensure it remains within the finite field.
                self.subModulusAssign();
            } else {
                // If there was a carry during addition or the result exceeds the modulus,
                // reduce the result modulo the modulus to maintain field integrity.
                self.subModulusWithCarryAssign(carry);
            }
        }

        /// Subtracts a field element from another field element.
        ///
        /// This function subtracts a field element `rhs` from another field element `self`.
        /// It dereferences the pointers `self` and `rhs` to obtain the actual field elements,
        /// performs the subtraction operation using the `subAssign` function, and returns the result.
        ///
        /// Parameters:
        ///   - self: A pointer to the first field element from which the second field element will be subtracted.
        ///   - rhs: A pointer to the second field element that will be subtracted from the first field element.
        ///
        /// Returns:
        ///   - The result of the subtraction operation as a new field element.
        ///
        /// Notes:
        ///   - This function does not modify the original field elements; it only performs the subtraction operation.
        ///   - The subtraction operation is performed using the `subAssign` function, which modifies the first field element in place.
        pub fn sub(self: *const Self, rhs: *const Self) Self {
            // Dereference the pointer to obtain the actual field element
            var a = self.*;
            // Perform the subtraction operation in place
            a.subAssign(rhs);
            // Return the result
            return a;
        }

        /// Subtracts a bigint from another bigint and assigns the result to the first bigint.
        ///
        /// This function subtracts a bigint `b` from another bigint `a` and assigns the result to `a`.
        /// If `b` is greater than `a`, it adds the modulus to `a` first to ensure correct subtraction in a finite field.
        ///
        /// Parameters:
        ///   - a: A pointer to the bigint from which `b` will be subtracted, and the result will be assigned.
        ///   - b: A pointer to the bigint that will be subtracted from `a`.
        ///
        /// Returns:
        ///   - void
        ///
        /// Notes:
        ///   - If `b` is greater than `a`, the modulus of the finite field is added to `a` before subtraction to ensure correct arithmetic in a finite field context.
        ///   - The subtraction operation is performed in place, and the result is assigned to `a`.
        pub fn subAssign(self: *Self, rhs: *const Self) void {
            // If b > a, add the modulus to `a` first.
            if (rhs.fe.cmp(&self.fe) == .gt)
                _ = self.fe.addWithCarryAssign(&Modulus);

            // Perform the subtraction operation
            _ = self.fe.subWithBorrowAssign(&rhs.fe);
        }

        pub fn inverse2(self: *const Self) ?Self {
            if (self.isZero()) return null;

            const _one = comptime big_int.fromInt(u64, 1);
            const modulus = comptime Modulus;

            var u = self.fe;
            var v = Modulus;
            var b = R2; // Avoids unnecessary reduction step.
            var c = Self.zero().fe;

            while (!u.eql(_one) and !v.eql(_one)) {
                while (u.limbs[0] & 1 == 0) {
                    u.shrAssign(1);

                    if (b.limbs[0] & 1 == 0) {
                        b.shrAssign(1);
                    } else {
                        const carry = b.addWithCarryAssign(&modulus);
                        b.shrAssign(1);

                        if (comptime !modulusHasSpareBit() and carry) {
                            b.limbs[n_limbs - 1] |= 1 << 63;
                        }
                    }
                }

                while (v.limbs[0] & 1 == 0) {
                    v.shrAssign(1);

                    if (c.limbs[0] & 1 == 0) {
                        c.shrAssign(1);
                    } else {
                        const carry = c.addWithCarryAssign(&modulus);
                        c.shrAssign(1);
                        if (comptime !modulusHasSpareBit() and carry) {
                            c.limbs[n_limbs - 1] |= 1 << 63;
                        }
                    }
                }

                if (v.cmp(&u).compare(.lte)) {
                    _ = u.subWithBorrowAssign(&v);
                    if (b.cmp(&c).compare(.lt)) {
                        _ = b.addWithCarryAssign(&modulus);
                    }

                    _ = b.subWithBorrowAssign(&c);
                } else {
                    _ = v.subWithBorrowAssign(&u);
                    if (c.cmp(&b).compare(.lt)) {
                        _ = c.addWithCarryAssign(&modulus);
                    }
                    _ = c.subWithBorrowAssign(&b);
                }
            }

            if (u.eql(_one)) return .{ .fe = b };

            return .{ .fe = c };
        }

        /// Computes the multiplicative inverse of a given element in a finite field using the binary Extended Euclidean Algorithm (BEA).
        ///
        /// Reference: Efficient Software-Implementation of Finite Fields with Applications to Cryptography
        /// DOI: DOI: 10.1007/s10440-006-9046-1
        ///
        /// This function implements the binary Extended Euclidean Algorithm (BEA) to compute the multiplicative inverse of a given element in a finite field.
        /// It follows the steps outlined in the BEA, including successive division and modular arithmetic operations, to determine the inverse.
        ///
        ///  BEA does not require integer divisions but only simple operations such as shifts and additions
        ///
        /// Parameters:
        ///   - self: A pointer to the element for which the inverse is to be computed.
        ///
        /// Returns:
        ///   - On success: A structure containing the computed inverse.
        ///   - On failure (if the input is zero): `null`.
        ///
        /// Notes:
        ///   - The binary Extended Euclidean Algorithm (BEA) is a general and efficient method for computing multiplicative inverses in finite fields.
        ///   - Montgomery parameters are used to optimize the computation, improving performance on digital computers.
        ///   - Overflow handling is performed to ensure correct arithmetic operations during the inversion process.
        pub fn inverse(self: *const Self) ?Self {
            // Check if the input is zero
            if (self.isZero()) return null;

            // Constant representing the value 1 in the field
            const o = comptime big_int.one();

            var u = self.fe;
            var v = Modulus;
            var b: Self = .{ .fe = R2 };
            var c = comptime zero();

            // Iterate while both u and v are not one
            while (!u.eql(o) and !v.eql(o)) {
                // Perform operations while u is even
                while (u.isEven()) {
                    u.div2Assign();

                    if (b.fe.isEven()) {
                        b.fe.div2Assign();
                    } else {
                        const carry = b.fe.addWithCarryAssign(&Modulus);
                        b.fe.div2Assign();
                        // Handle overflow if necessary
                        if (comptime !modulusHasSpareBit() and carry) {
                            b.fe.limbs[n_limbs - 1] |= 1 << 63;
                        }
                    }
                }

                // Perform operations while v is even
                while (v.isEven()) {
                    v.div2Assign();
                    if (c.fe.isEven()) {
                        c.fe.div2Assign();
                    } else {
                        const carry = c.fe.addWithCarryAssign(&Modulus);
                        c.fe.div2Assign();
                        // Handle overflow if necessary
                        if (comptime !modulusHasSpareBit() and carry) {
                            c.fe.limbs[n_limbs - 1] |= 1 << 63;
                        }
                    }
                }

                // Update based on u vs v values
                if (v.cmp(&u).compare(.lte)) {
                    _ = u.subWithBorrowAssign(&v);
                    b.subAssign(&c);
                } else {
                    _ = v.subWithBorrowAssign(&u);
                    c.subAssign(&b);
                }
            }

            // Return the inverse based on the final values of u and v
            return if (u.eql(o)) b else c;
        }

        /// Divides one field element by another field element.
        ///
        /// This function divides the current field element (`self`) by another field element (`den`).
        /// If the denominator (`den`) is zero, it returns an error indicating a division by zero.
        ///
        /// # Arguments:
        /// - `den`: The field element by which to divide the current field element.
        ///
        /// # Returns:
        /// - Result of the division operation as a new field element.
        ///
        /// Errors:
        /// - If the denominator is zero, a division by zero error is returned.
        pub fn div(self: Self, den: Self) FieldError!Self {
            return self.mul(&(den.inverse() orelse return FieldError.DivisionByZero));
        }

        pub fn modFloor2(self: Self, den: Self) Self {
            _, const rem = self.divRem(den);
            return rem;
        }

        pub fn modFloor(self: Self, den: Self) Self {
            return Self.toMontgomery(big_int.fromInt(u256, @mod(self.toU256(), den.toU256())));
        }

        // TODO: make desc
        pub fn divRem2(self: Self, den: Self) struct { Self, Self } {
            const q, const r = self.toBigInt().divRem(&den.toBigInt());
            return .{
                Self.toMontgomery(q),
                Self.toMontgomery(r),
            };
        }

        pub fn divRem(self: Self, den: Self) struct { Self, Self } {
            const s = self.toU256();
            const d = den.toU256();
            const q, const r = .{
                @divFloor(s, d),
                @mod(s, d),
            };
            return .{
                Self.toMontgomery(big_int.fromInt(u256, q)),
                Self.toMontgomery(big_int.fromInt(u256, r)),
            };
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
            return std.mem.readInt(
                u256,
                &self.toBigInt().toBytesLe(),
                .little,
            );
        }

        /// Try to convert the field element to a u64 if its value is small enough.
        ///
        /// Attempts to convert the field element to a u64 if its value is within the representable range.
        ///
        /// Parameters:
        ///   - self: The field element to convert.
        ///   - T: The target type for conversion (must be u64 or smaller).
        ///
        /// Returns:
        ///   - A u64 representation of the field element if conversion succeeds.
        ///   - Error(ValueTooLarge) if the value exceeds the representable range of the target type.
        pub fn toInt(self: Self, comptime T: type) !T {
            // TODO maybe exist a better way to convert
            const asU256 = self.toU256();
            // Check if the value is small enough to fit into a type T integer
            if (asU256 > std.math.maxInt(T)) return error.ValueTooLarge;

            // Otherwise, it's safe to cast
            return @intCast(asU256);
        }

        pub fn toSignedInt(self: Self, comptime T: type) !T {
            // TODO: improve this code
            const asU256 = self.toU256();

            if (asU256 > comptime value: {
                break :value modulo >> @as(u32, 1);
            }) return -@as(T, @intCast(modulo - asU256));

            return @intCast(asU256);
        }

        /// Calculate the Legendre symbol of a field element.
        ///
        /// Computes the Legendre symbol of the field element using Euler's criterion.
        /// The Legendre symbol is a mathematical function commonly used in number theory
        /// to determine whether a given integer is a quadratic residue modulo a prime number.
        ///
        /// # Arguments:
        /// - `a`: The field element for which the Legendre symbol is calculated.
        ///
        /// # Returns:
        /// - `1` if `a` has a square root modulo the modulus (`p`),
        /// - `-1` if `a` does not have a square root modulo `p`,
        /// - `0` if `a` is zero modulo `p`.
        pub fn legendre(a: Self) i2 {
            // Compute the Legendre symbol a|p using
            // Euler's criterion. p is a prime, a is
            // relatively prime to p (if p divides
            // a, then a|p = 0)
            // Returns 1 if a has a square root modulo
            // p, -1 otherwise.

            // Calculate a^(p-1)/2 modulo p
            const ls = a.powToBigInt(comptime big_int.fromInt(u256, QMinOneDiv2));

            // Check if a^(p-1)/2 is equivalent to -1 modulo p
            if (ls.toU256() == comptime modulo - 1) return -1;

            // Check if a^(p-1)/2 is equivalent to 0 modulo p
            if (ls.isZero()) return 0;

            // Otherwise, a^(p-1)/2 is equivalent to 1 modulo p
            return 1;
        }

        /// Compare two field elements and return the ordering result.
        ///
        /// # Parameters
        /// - `self` - The first field element to compare.
        /// - `rhs` - The second field element to compare.
        ///
        /// # Returns
        /// A `std.math.Order` enum indicating the ordering relationship.
        pub fn cmp(self: *const Self, rhs: *const Self) std.math.Order {
            return self.toBigInt().cmp(&rhs.toBigInt());
        }

        /// Check if this field element is less than the rhs.
        ///
        /// # Parameters
        /// - `self` - The field element to check.
        /// - `rhs` - The field element to compare against.
        ///
        /// # Returns
        /// `true` if `self` is less than `rhs`, `false` otherwise.
        pub fn lt(self: *const Self, rhs: *const Self) bool {
            return self.cmp(rhs) == .lt;
        }

        /// Check if this field element is less than or equal to the rhs.
        ///
        /// # Parameters
        /// - `self` - The field element to check.
        /// - `rhs` - The field element to compare against.
        ///
        /// # Returns
        /// `true` if `self` is less than or equal to `rhs`, `false` otherwise.
        pub fn lte(self: *const Self, rhs: *const Self) bool {
            return self.cmp(rhs).compare(.lte);
        }

        /// Check if this field element is greater than the rhs.
        ///
        /// # Parameters
        /// - `self` - The field element to check.
        /// - `rhs` - The field element to compare against.
        ///
        /// # Returns
        /// `true` if `self` is greater than `rhs`, `false` otherwise.
        pub fn gt(self: *const Self, rhs: *const Self) bool {
            return self.cmp(rhs) == .gt;
        }

        /// Check if this field element is greater than or equal to the rhs.
        ///
        /// # Parameters
        /// - `self` - The field element to check.
        /// - `rhs` - The field element to compare against.
        ///
        /// # Returns
        /// `true` if `self` is greater than or equal to `rhs`, `false` otherwise.
        pub fn gte(self: *const Self, rhs: *const Self) bool {
            return self.cmp(rhs).compare(.gte);
        }
    };
}

// use num_prime;

extern crate libc;

// C representation of a bit array: a raw pointer to a mutable unsigned 8 bits integer.
type Bytes = *mut u8;

use num_bigint::BigUint;

#[no_mangle]
extern "C" fn is_prime(number: Bytes, num_limbs: usize) -> bool {
    let n = unsafe {
        let slice: &mut [u8] = std::slice::from_raw_parts_mut(number, num_limbs * 8);
        BigUint::from_bytes_le(slice)
    };

    num_prime::nt_funcs::is_prime(&n, None).probably()
}

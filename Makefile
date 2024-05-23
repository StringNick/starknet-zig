build_prime_library:
	@rm -f src/math/fields/prime/libstarknet_crypto.a
	@cd src/math/fields/prime/prime && cargo build --release
	@mv src/math/fields/prime/prime/target/release/libprime.a src/math/fields/prime

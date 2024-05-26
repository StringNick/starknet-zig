build_prime_library:
	@rm -f src/math/fields/prime/libprime.a
	@cd src/math/fields/prime/prime && cargo build --release
	@mv src/math/fields/prime/prime/target/release/libprime.a src/math/fields/prime

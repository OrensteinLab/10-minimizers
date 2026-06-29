
#define _CRT_SECURE_NO_WARNINGS

#include "sampling_runtimes.h"
#include <ctime>
#include <malloc.h> // For _alloca on Windows or alloca on Linux

/**
 * PCG RXS-M-XS 64-bit Mixer
 * Highly robust pseudorandom mapping for single 64-bit integers.
 * Source: https://en.wikipedia.org/wiki/Permuted_congruential_generator
 */
inline uint64_t pcg_mix_64(uint64_t x) {
	uint64_t state = x * 6364136223846793005ULL + 1442695040888963407ULL;
	uint64_t word = ((state >> ((state >> 59u) + 5u)) ^ state) * 12605985483714917081ULL;
	return (word >> 43u) ^ word;
}

inline uint32_t pcg_mix_32(uint32_t x) {
	uint32_t state = x * 747796405U + 2891336453U;
	uint32_t word = ((state >> ((state >> 28u) + 4u)) ^ state) * 277803737U;
	return (word >> 22u) ^ word;
}



std::vector<bool> get_bits(uint64_t n_bits) {
	std::vector<bool> bits(n_bits);

	std::random_device rd;
	std::mt19937 gen(rd());
	std::uniform_int_distribution<int> dist(0, 1);

	size_t extra_bits = 1000; // extra bits for algorithms that need it

	for (size_t i = 0; i < (n_bits + extra_bits); ++i) {
		bits[i] = dist(gen); // Either 0 or 1
	}

	return bits;

}


double test_runtime_gm(uint32_t w, uint32_t k, uint64_t N, const std::vector<bool>& bits) {

	uint64_t k_mask = (1ULL << k) - 1;

	// 1. Create a vector "order" of size 2^k
	std::vector<std::uint64_t> order(1ULL << k);

	// 2. Set up random number generation
	std::random_device rd;              // Non-deterministic seed
	std::mt19937_64 gen(rd());          // Mersenne Twister 64-bit
	std::uniform_int_distribution<std::uint64_t> dist(
		0, std::numeric_limits<std::uint64_t>::max()
	);

	// 3. Fill the vector with random data
	for (auto& elem : order) {
		elem = dist(gen);
	}




	//std::vector<uint64_t> order = load_gm_order(w, k);

	// multiply all ranks in order by 2^k
	for (uint64_t i = 0; i < order.size(); i++) {
		order[i] = order[i] << k;
	}


	uint64_t sum = 0;




	// start timing 
	auto start = std::chrono::high_resolution_clock::now();

	uint64_t odd_kmer = 0;
	uint64_t even_kmer = 0;

	for (uint64_t i = 0; i < N; i += 2) {
		even_kmer = (even_kmer << 1) & k_mask;
		even_kmer |= bits[i];
		odd_kmer = (odd_kmer << 1) & k_mask;
		odd_kmer |= bits[i + 1];
		uint64_t kmer_rank = order[odd_kmer] + even_kmer;
		sum += kmer_rank;
	}


	auto end = std::chrono::high_resolution_clock::now();

	std::cout << "sum:\t" << sum << std::endl;

	// Measure elapsed time in nanoseconds
	auto elapsed_ns = std::chrono::duration_cast<std::chrono::nanoseconds>(end - start).count();
	// Or microseconds
	auto elapsed_us = std::chrono::duration_cast<std::chrono::microseconds>(end - start).count();

	std::cout << "Time taken gm (ms):  " << elapsed_us / 1000.0 << "\n";
	return elapsed_us / 1000.0;
}





double test_runtime_baseline_xor(uint32_t w, uint32_t k, uint64_t N, const std::vector<bool>& bits) {
	uint64_t dna_k_mask = (1 << (2 * k)) - 1;
	const uint64_t random_number = 0x9E3779B97F4A7C15ULL;

	
	uint64_t sum = 0;

	// start timing
	auto start = std::chrono::high_resolution_clock::now();

	uint64_t current_kmer = 0;
	for (uint64_t i = 0; i < N; i += 2) {
		current_kmer = (current_kmer << 1) & dna_k_mask;
		current_kmer |= (bits[i] << 1) | bits[i + 1];
		sum += current_kmer ^ random_number;
	}



	auto end = std::chrono::high_resolution_clock::now();

	std::cout << "sum:\t" << sum << std::endl;

	// Measure elapsed time in nanoseconds
	auto elapsed_ns = std::chrono::duration_cast<std::chrono::nanoseconds>(end - start).count();
	// Or microseconds
	auto elapsed_us = std::chrono::duration_cast<std::chrono::microseconds>(end - start).count();

	std::cout << "Time taken baseline xor (ms):  " << elapsed_us / 1000.0 << "\n";

	return elapsed_us / 1000.0;

}

double test_runtime_one_zero_most_optimized_lazy(uint32_t w, uint32_t k, uint64_t N, const std::vector<bool>& bits) {
	uint64_t dna_k_mask = (1ull << (2 * k)) - 1;
	uint64_t proj_k_mask = (1ull << k) - 1;
	uint64_t mask_after_prefix = (1ull << (k - 1)) - 1;

	const uint64_t is_10_match = 2ull << (k - 2);
	const uint64_t is_10_mask = 3ull << (k - 2);

	uint64_t sum1 = 0;
	uint64_t sum2 = 0;

	// Precompute Lemma 1 projections for the current k
	std::unordered_map<uint64_t, int> lemma1_ranks;
	std::vector<uint64_t> tau;
	tau.push_back((1ull << k) - 2);
	tau.push_back((1ull << (k - 1)) - 2);
	for (uint64_t i = 2; i < k - 1; ++i) tau.push_back(2 * ((1ull << (k - 1 - i)) - 1));
	tau.push_back(1);
	tau.push_back((1ull << k) - 1);
	tau.push_back(0);

	uint64_t class_bits = (1ull << (62));
	uint64_t half_shift = (1ull << (32));

	for (int i = 0; i < (int)tau.size(); ++i) lemma1_ranks[tau[i]] = i;


	// Circular Buffer for Lazy Ranking
	struct PendingKmer { uint64_t kmer; uint64_t projection; };
	constexpr uint64_t BUF_SIZE = 128; // Fixed power of 2
	constexpr uint64_t BUF_MASK = BUF_SIZE - 1;
	PendingKmer buffer[BUF_SIZE]; // Allocated on the stack

	uint64_t buf_start = 0;
	uint64_t buf_count = 0;
	uint64_t class1_constant = (1 * class_bits);
	uint64_t class2_constant = (2 * class_bits);

	auto start = std::chrono::high_resolution_clock::now();

	uint64_t current_kmer = 0;
	uint64_t projection = 0;
	uint64_t sixty_four_minus_k = 64 - k;


	for (uint64_t i = 0; i < N; i += 2) {
		current_kmer = (current_kmer << 2) & dna_k_mask;
		uint8_t incoming_nuc = (bits[i] << 1) | bits[i + 1];
		current_kmer |= incoming_nuc;

		projection = (projection << 1) & proj_k_mask;
		projection |= (incoming_nuc >> 1) & incoming_nuc;

		// Use the mask logic to check for Tier 0 (10-kmer)
		if ((projection & is_10_mask) == is_10_match) [[unlikely]] {

			// Filter for 10 patterns
			uint64_t next_10 = (projection & ~(projection << 1)) ^ (1ull << (k - 1));


			uint64_t tail_len = std::countl_zero(next_10) - sixty_four_minus_k;// technically counts 32-k more zeros, but it should still fit in 
			uint64_t tail_val = projection & ((1ull << tail_len) - 1);

			// assume k<=32
			sum1 += (tail_len << 32ULL) - tail_val;
			sum2 += current_kmer;



			// RESET: Clear the buffer index because a Tier 0 k-mer makes previous Tiers irrelevant
			buf_count = 0;
			buf_start = 0;
			}
		else [[likely]] {
			// Window management for non-10-kmers
			if (buf_count == (w - 1)) [[unlikely]] {
				// Rank the oldest element before overwriting
				const PendingKmer& oldest = buffer[buf_start];
				uint64_t fallback_rank = 0;

				auto it = lemma1_ranks.find(oldest.projection);
				if (it != lemma1_ranks.end()) [[unlikely]] { // there are ~k such kmers and an expo amount of other kmers
					sum1 += (class1_constant)+(it->second);
					sum2 += oldest.kmer; // current_kmer
					}
				else [[likely]] {
					sum1 += class2_constant;
					sum2 -= oldest.kmer; // current_kmer
					}

				buffer[buf_start] = { current_kmer, projection };
				buf_start = (buf_start + 1) & BUF_MASK; // since buffer size is 128
				}
			else [[likely]] {
				uint64_t next_pos = (buf_start + buf_count) & BUF_MASK;
				buffer[next_pos] = { current_kmer, projection };
				buf_count++;
				sum1 += 1; // infty
				sum2 += 1; // simulate not ranking
				}
			}
	}

	auto end = std::chrono::high_resolution_clock::now();
	std::cout << "sum:\t" << sum1 << sum2 << std::endl;
	auto elapsed_us = std::chrono::duration_cast<std::chrono::microseconds>(end - start).count();
	std::cout << "Time taken 10-Min Lazy (ms): " << elapsed_us / 1000.0 << "\n";
	return elapsed_us / 1000.0;
}


double test_runtime_baseline_lex(uint32_t w, uint32_t k, uint64_t N, const std::vector<bool>& bits) {
	uint64_t dna_k_mask = (1 << (2 * k)) - 1;
	uint64_t sum = 0;

	// start timing
	auto start = std::chrono::high_resolution_clock::now();

	uint64_t current_kmer = 0;
	for (uint64_t i = 0; i < N; i += 2) {
		current_kmer = (current_kmer << 1) & dna_k_mask;
		current_kmer |= (bits[i] << 1) | bits[i + 1];
		sum += current_kmer;
	}



	auto end = std::chrono::high_resolution_clock::now();

	std::cout << "sum:\t" << sum << std::endl;

	// Measure elapsed time in nanoseconds
	auto elapsed_ns = std::chrono::duration_cast<std::chrono::nanoseconds>(end - start).count();
	// Or microseconds
	auto elapsed_us = std::chrono::duration_cast<std::chrono::microseconds>(end - start).count();

	std::cout << "Time taken baseline lex (ms):  " << elapsed_us / 1000.0 << "\n";

	return elapsed_us / 1000.0;

}

double test_runtime_abb_plus_lex(uint32_t w, uint32_t k, uint64_t N, const std::vector<bool>& bits) {
	uint64_t dna_k_mask = (1ull << (2 * k)) - 1;
	const uint64_t a_mask_limit = (1ull << (k - 1)) - 1; // We only care about k-1 A-bits
	uint64_t sum1 = 0;
	uint64_t sum2 = 0;
	uint64_t L_shift = 2 * k;

	auto start = std::chrono::high_resolution_clock::now();

	uint64_t current_kmer = 0;
	uint64_t a_bits_window = 0; // Sliding window of (nuc == 0) bits
	for (uint64_t i = 0; i < N; i += 2) {
		// 1. Extract the new nucleotide from the bitstream
		uint8_t nuc = (bits[i] << 1) | bits[i + 1];

		// 2. Update the DNA k-mer
		current_kmer = ((current_kmer << 2) | nuc) & dna_k_mask;

		// 3. Update the "A-bits" window (1 if nuc is 'A' (0), else 0)
		// We shift in the new bit and mask to size k
		a_bits_window = ((a_bits_window << 1) | (nuc == 0)) & ((1ull << k) - 1);

		// 4. Construct the Key on the fly
		// Part A: Prefix (The first nucleotide of the k-mer)
		uint64_t nuc0 = (current_kmer >> (2 * (k - 1))) & 0b11;

		// Part B: The A-bits for positions 1 to k-1
		// These are exactly the lower k-1 bits of our sliding a_bits_window
		uint64_t key = (nuc0 << (k - 1)) | (a_bits_window & a_mask_limit);

		// 5. Calculate final rank
		sum1 += key;
		sum2 += current_kmer;
	}

	auto end = std::chrono::high_resolution_clock::now();
	std::cout << "sum:\t" << sum1 << sum2 << std::endl;
	auto elapsed_us = std::chrono::duration_cast<std::chrono::microseconds>(end - start).count();
	std::cout << "Time taken ABB+ lex (ms): " << elapsed_us / 1000.0 << "\n";
	return elapsed_us / 1000.0;
}

double test_runtime_one_zero_states_lex(uint32_t w, uint32_t k, uint64_t N, const std::vector<bool>& bits) {
	// Primary Key Constants
	const uint64_t X1 = 1ULL << 62; // p ends in 10
	const uint64_t X2 = 2ULL << 62; // p = 0...01
	const uint64_t X3 = 3ULL << 62; // p = 1...1
	const uint64_t X4 = 4ULL << 62; // p = 0...0
	const uint64_t INF = 5ULL << 62; // Others



	const uint64_t dna_k_mask = (1ull << (2 * k)) - 1;
	const uint64_t proj_k_mask = (1ull << (k-1)) - 1; // only remember last k-1 bits, and not k - to discard leading 1
	const uint64_t k_prime = k - 1;
	const uint64_t sixty_four_minus_k = 64ULL - k;

	uint64_t sum1 = 0, sum2 = 0;
	uint64_t current_kmer = 0, projection = 0;

	auto start = std::chrono::high_resolution_clock::now();

	// --- Manual Initialization of first two nucleotides to be 10
	
	current_kmer = (3ULL << 2);
	projection = (1ULL << 2);

	uint64_t i = 0;
	// Read the remaining k-3 nucleotides to fill the initial k-mer window aside from last bit
	uint64_t startup_nucs = (k > 4) ? k - 4 : 0;
	for (uint64_t j = 0; j < startup_nucs && i < N; ++j, i += 2) {
		uint8_t nuc = (bits[i] << 1) | bits[i + 1];
		current_kmer = ((current_kmer << 2) | nuc) & dna_k_mask;
		projection = ((projection << 1) | ((nuc >> 1) & nuc)) & proj_k_mask;
	}

	// For simplicity, we jump to STATE_10 and let it handle the first logic gate
	goto STATE_10;

STATE_10:
	while (i < N) {
		uint8_t nuc = (bits[i] << 1) | bits[i + 1];
		current_kmer = ((current_kmer << 2) | nuc) & dna_k_mask;
		projection = ((projection << 1) | ((nuc >> 1) & nuc)) & proj_k_mask;
		i += 2;

		uint64_t next_10 = (projection & ~(projection << 1)) ; // no need to do ^ (1ull << k_prime), since we save only k-1 bits (the leading 1 is already discarded)
		uint64_t t = 64ULL - std::countl_zero(next_10);
		uint64_t tail_val = projection & ((1ull << t) - 1);

		sum1 += (t << 32) - tail_val;
		sum2 += current_kmer;

		if (t > 1) {
			// Trivial iterations: next k'-t kmers are INF
			//trivial_limit = k_prime - t;
			uint64_t trivial_limit = i + 2*(k_prime - t);
			for (;i < trivial_limit; i += 2)
			{
				uint8_t tnuc = (bits[i] << 1) | bits[i + 1];
				current_kmer = ((current_kmer << 2) | tnuc) & dna_k_mask;
				projection = ((projection << 1) | ((tnuc >> 1) & tnuc)) & proj_k_mask;
				sum1 += INF;
				//sum2 += current_kmer;
			}
			continue; // Remain in STATE_10
		}
		else if (t == 0) {
			goto STATE_0;
		}
		else {
			// t == 1. Compute i_count by counting leading zeros relative to k_prime
			// projection looks like 0^{k'-i} 1^i
			uint64_t i_count = k_prime - (std::countl_zero(projection ^ (1ull << k_prime)) - sixty_four_minus_k);
			// Move to state1 with i_count
			goto STATE_1_ENTRY;
		}
	}
	goto END;

STATE_0:
	while (i < N) {
		uint8_t nuc = (bits[i] << 1) | bits[i + 1];
		current_kmer = ((current_kmer << 2) | nuc) & dna_k_mask;
		projection = ((projection << 1) | ((nuc >> 1) & nuc)) & proj_k_mask;
		i += 2;

		if ((projection & 1) == 0) { // proj(bp) == 0
			sum1 += X4;
			sum2 += current_kmer;
			// Remain in STATE_0
		}
		else {
			sum1 += X2;
			//sum2 += current_kmer;
			// Move to state1(1)
			uint64_t i_count = 1;
			goto STATE_1_ENTRY;
		}
	}
	goto END;

STATE_1_ENTRY:
	{
		uint64_t i_count = 1;
		while (i < N) {
			uint8_t nuc = (bits[i] << 1) | bits[i + 1];
			current_kmer = ((current_kmer << 2) | nuc) & dna_k_mask;
			projection = ((projection << 1) | ((nuc >> 1) & nuc)) & proj_k_mask;
			i += 2;

			if ((projection & 1) == 0) { // proj(bp) == 0
				sum1 += X1;
				//sum2 += current_kmer;
				// Trivial iterations
				for (uint64_t j = 0; j < k_prime - 2 && i < N; ++j, i += 2) {
					uint8_t tnuc = (bits[i] << 1) | bits[i + 1];
					current_kmer = ((current_kmer << 2) | tnuc) & dna_k_mask;
					projection = ((projection << 1) | ((tnuc >> 1) & tnuc)) & proj_k_mask;
					sum1 += INF;
					//sum2 += current_kmer;
				}
				goto STATE_10;
			}
			else {
				if (i_count < k_prime) {
					sum1 += INF;
					//sum2 += current_kmer;
					i_count++;
				}
				else {
					sum1 += X3;
					//sum2 += current_kmer;
				}
			}
		}
	}

END:
	auto end = std::chrono::high_resolution_clock::now();
	std::cout << "sum:\t" << sum1 << sum2 << std::endl;
	auto elapsed_us = std::chrono::duration_cast<std::chrono::microseconds>(end - start).count();
	std::cout << "Time taken 10-min lex (ms): " << elapsed_us / 1000.0 << "\n";
	return elapsed_us / 1000.0;
}

double test_runtime_one_zero_states_random_lex(uint32_t w, uint32_t k, uint64_t N, const std::vector<bool>& bits) {
	// Primary Key Constants
	const uint64_t X1 = 1ULL << 62; // p ends in 10
	const uint64_t X2 = 2ULL << 62; // p = 0...01
	const uint64_t X3 = 3ULL << 62; // p = 1...1
	const uint64_t X4 = 4ULL << 62; // p = 0...0
	const uint64_t INF = 5ULL << 62; // Others


	


	const uint64_t dna_k_mask = (1ull << (2 * k)) - 1;
	const uint64_t proj_k_mask = (1ull << (k - 1)) - 1; // only remember last k-1 bits, and not k - to discard leading 1
	const uint64_t k_prime = k - 1;
	const uint64_t sixty_four_minus_k = 64 - k;

	uint64_t sum1 = 0, sum2 = 0;
	uint64_t current_kmer = 0;
	uint64_t projection = 0;

	auto start = std::chrono::high_resolution_clock::now();

	// --- Manual Initialization of first two nucleotides to be 10

	current_kmer = (3ULL << 2);
	projection = (1ULL << 2);

	uint64_t i = 0;
	// Read the remaining k-3 nucleotides to fill the initial k-mer window aside from last bit
	uint64_t startup_nucs = (k > 4) ? k - 4 : 0;
	for (uint64_t j = 0; j < startup_nucs && i < N; ++j, i += 2) {
		uint8_t nuc = (bits[i] << 1) | bits[i + 1];
		current_kmer = ((current_kmer << 2) | nuc) & dna_k_mask;
		projection = ((projection << 1) | (nuc >> 1)) & proj_k_mask;
	}

	// For simplicity, we jump to STATE_10 and let it handle the first logic gate
	goto STATE_10;

STATE_10:
	while (i < N) {
		uint8_t nuc = (bits[i] << 1) | bits[i + 1];
		current_kmer = ((current_kmer << 2) | nuc) & dna_k_mask;
		projection = ((projection << 1) | (nuc >> 1)) & proj_k_mask;
		i += 2;

		uint64_t next_10 = (projection & ~(projection << 1));
		uint64_t t = 64 - std::countl_zero(next_10);

		sum1 += projection;
		sum2 += current_kmer;

		if (t > 1) {
			// Trivial iterations: next k'-t kmers are INF
			//uint32_t trivial_limit = k_prime - t;
			uint64_t trivial_limit = i + 2*(k_prime - t);
			//for (uint32_t j = 0; j < trivial_limit && i < N; ++j, i += 2) {
			for (; i < trivial_limit; i += 2)
			{
				uint8_t tnuc = (bits[i] << 1) | bits[i + 1];
				current_kmer = ((current_kmer << 2) | tnuc) & dna_k_mask;
				projection = ((projection << 1) | ((tnuc >> 1) & tnuc)) & proj_k_mask;
				sum1 += INF;
				//sum2 += current_kmer;
			}
			continue; // Remain in STATE_10
		}
		else if (t != 0) {
			// t == 1. Compute i_count by counting leading zeros relative to k_prime
			// projection looks like 0^{k'-i} 1^i
			uint64_t i_count = k_prime - (std::countl_zero(projection ^ (1ull << k_prime)) - sixty_four_minus_k);
			// Move to state1 with i_count
			goto STATE_1_ENTRY;
		}
		else{
			goto STATE_0;
		}

	}
	goto END;

STATE_0:
	while (i < N) {
		uint8_t nuc = (bits[i] << 1) | bits[i + 1];
		current_kmer = ((current_kmer << 2) | nuc) & dna_k_mask;
		projection = ((projection << 1) | (nuc >> 1)) & proj_k_mask;
		i += 2;

		if ((projection & 1) == 0) { // proj(bp) == 0
			sum1 += X4;
			sum2 += current_kmer;
			// Remain in STATE_0
		}
		else {
			sum1 += X2;
			//sum2 += current_kmer;
			// Move to state1(1)
			uint64_t i_count = 1;
			goto STATE_1_ENTRY;
		}
	}
	goto END;

STATE_1_ENTRY:
	{
		uint64_t i_count = 1;
		while (i < N) {
			uint8_t nuc = (bits[i] << 1) | bits[i + 1];
			current_kmer = ((current_kmer << 2) | nuc) & dna_k_mask;
			projection = ((projection << 1) | (nuc >> 1)) & proj_k_mask;
			i += 2;

			if ((projection & 1) == 0) { // proj(bp) == 0
				sum1 += X1;
				//sum2 += current_kmer;
				// Trivial iterations
				for (uint64_t j = 0; j < k_prime - 2 && i < N; ++j, i += 2) {
					uint8_t tnuc = (bits[i] << 1) | bits[i + 1];
					current_kmer = ((current_kmer << 2) | tnuc) & dna_k_mask;
					projection = ((projection << 1) | (tnuc >> 1)) & proj_k_mask;
					sum1 += INF;
					//sum2 += current_kmer;
				}
				goto STATE_10;
			}
			else {
				if (i_count < k_prime) {
					sum1 += INF;
					//sum2 += current_kmer;
					i_count++;
				}
				else {
					sum1 += X3;
					//sum2 += current_kmer;
				}
			}
		}
	}

END:
	auto end = std::chrono::high_resolution_clock::now();
	std::cout << "sum:\t" << sum1 << sum2 << std::endl;
	auto elapsed_us = std::chrono::duration_cast<std::chrono::microseconds>(end - start).count();
	std::cout << "Time taken 10-min random lex (ms): " << elapsed_us / 1000.0 << "\n";
	return elapsed_us / 1000.0;
}

double test_runtime_one_zero_states_random_hybrid_lex(uint32_t w, uint32_t k, uint64_t N, const std::vector<bool>& bits) {
	// Primary Key Constants
	const uint64_t X1 = 1ULL << 62; // p ends in 10
	const uint64_t X2 = 2ULL << 62; // p = 0...01
	const uint64_t X3 = 3ULL << 62; // p = 1...1
	const uint64_t X4 = 4ULL << 62; // p = 0...0
	const uint64_t INF = 5ULL << 62; // Others





	const uint64_t dna_k_mask = (1ull << (2 * k)) - 1;
	const uint64_t proj_k_mask = (1ull << (k - 1)) - 1; // only remember last k-1 bits, and not k - to discard leading 1
	const uint64_t k_prime = k - 1;
	const uint64_t sixty_four_minus_k = 64 - k;

	uint64_t sum1 = 0, sum2 = 0;
	uint64_t current_kmer = 0;
	uint64_t projection = 0;

	auto start = std::chrono::high_resolution_clock::now();

	// --- Manual Initialization of first two nucleotides to be 10

	current_kmer = (3ULL << 2);
	projection = (1ULL << 2);

	uint64_t i = 0;
	// Read the remaining k-3 nucleotides to fill the initial k-mer window aside from last bit
	uint64_t startup_nucs = (k > 4) ? k - 4 : 0;
	for (uint64_t j = 0; j < startup_nucs && i < N; ++j, i += 2) {
		uint8_t nuc = (bits[i] << 1) | bits[i + 1];
		current_kmer = ((current_kmer << 2) | nuc) & dna_k_mask;
		projection = ((projection << 1) | (nuc >> 1)) & proj_k_mask;
	}

	// For simplicity, we jump to STATE_10 and let it handle the first logic gate
	goto STATE_10;

STATE_10:
	while (i < N) {
		uint8_t nuc = (bits[i] << 1) | bits[i + 1];
		current_kmer = ((current_kmer << 2) | nuc) & dna_k_mask;
		projection = ((projection << 1) | (nuc >> 1)) & proj_k_mask;
		i += 2;

		uint64_t next_10 = (projection & ~(projection << 1));
		uint64_t t = 64 - std::countl_zero(next_10);

		uint64_t tail_val = projection & ((1ull << t) - 1);

		sum1 += (t << 32) - tail_val;
		sum2 += current_kmer;


		if (t > 1) {
			// Trivial iterations: next k'-t kmers are INF
			//uint32_t trivial_limit = k_prime - t;
			uint64_t trivial_limit = i + 2 * (k_prime - t);
			//for (uint32_t j = 0; j < trivial_limit && i < N; ++j, i += 2) {
			for (; i < trivial_limit; i += 2)
			{
				uint8_t tnuc = (bits[i] << 1) | bits[i + 1];
				current_kmer = ((current_kmer << 2) | tnuc) & dna_k_mask;
				projection = ((projection << 1) | ((tnuc >> 1) & tnuc)) & proj_k_mask;
				sum1 += INF;
				//sum2 += current_kmer;
			}
			continue; // Remain in STATE_10
		}
		else if (t != 0) {
			// t == 1. Compute i_count by counting leading zeros relative to k_prime
			// projection looks like 0^{k'-i} 1^i
			uint64_t i_count = k_prime - (std::countl_zero(projection ^ (1ull << k_prime)) - sixty_four_minus_k);
			// Move to state1 with i_count
			goto STATE_1_ENTRY;
		}
		else {
			goto STATE_0;
		}

	}
	goto END;

STATE_0:
	while (i < N) {
		uint8_t nuc = (bits[i] << 1) | bits[i + 1];
		current_kmer = ((current_kmer << 2) | nuc) & dna_k_mask;
		projection = ((projection << 1) | (nuc >> 1)) & proj_k_mask;
		i += 2;

		if ((projection & 1) == 0) { // proj(bp) == 0
			sum1 += X4;
			sum2 += current_kmer;
			// Remain in STATE_0
		}
		else {
			sum1 += X2;
			//sum2 += current_kmer;
			// Move to state1(1)
			uint64_t i_count = 1;
			goto STATE_1_ENTRY;
		}
	}
	goto END;

STATE_1_ENTRY:
	{
		uint64_t i_count = 1;
		while (i < N) {
			uint8_t nuc = (bits[i] << 1) | bits[i + 1];
			current_kmer = ((current_kmer << 2) | nuc) & dna_k_mask;
			projection = ((projection << 1) | (nuc >> 1)) & proj_k_mask;
			i += 2;

			if ((projection & 1) == 0) { // proj(bp) == 0
				sum1 += X1;
				//sum2 += current_kmer;
				// Trivial iterations
				for (uint64_t j = 0; j < k_prime - 2 && i < N; ++j, i += 2) {
					uint8_t tnuc = (bits[i] << 1) | bits[i + 1];
					current_kmer = ((current_kmer << 2) | tnuc) & dna_k_mask;
					projection = ((projection << 1) | (tnuc >> 1)) & proj_k_mask;
					sum1 += INF;
					//sum2 += current_kmer;
				}
				goto STATE_10;
			}
			else {
				if (i_count < k_prime) {
					sum1 += INF;
					//sum2 += current_kmer;
					i_count++;
				}
				else {
					sum1 += X3;
					//sum2 += current_kmer;
				}
			}
		}
	}

END:
	auto end = std::chrono::high_resolution_clock::now();
	std::cout << "sum:\t" << sum1 << sum2 << std::endl;
	auto elapsed_us = std::chrono::duration_cast<std::chrono::microseconds>(end - start).count();
	std::cout << "Time taken 10-min random lex (ms): " << elapsed_us / 1000.0 << "\n";
	return elapsed_us / 1000.0;
}

double test_runtime_one_zero_states_random_hybrid_lex_2(uint32_t w, uint32_t k, uint64_t N, const std::vector<bool>& bits) {
	// Primary Key Constants
	const uint64_t X1 = 1ULL << 62; // p ends in 10
	const uint64_t X2 = 2ULL << 62; // p = 0...01
	const uint64_t X3 = 3ULL << 62; // p = 1...1
	const uint64_t X4 = 4ULL << 62; // p = 0...0
	const uint64_t INF = 5ULL << 62; // Others





	const uint64_t dna_k_mask = (1ull << (2 * k)) - 1;
	const uint64_t proj_k_mask = (1ull << (k - 1)) - 1; // only remember last k-1 bits, and not k - to discard leading 1
	const uint64_t k_prime = k - 1;
	const uint64_t sixty_four_minus_k = 64 - k;

	uint64_t sum1 = 0, sum2 = 0;
	uint64_t current_kmer = 0;
	uint64_t projection = 0;

	auto start = std::chrono::high_resolution_clock::now();

	// --- Manual Initialization of first two nucleotides to be 10

	current_kmer = (3ULL << 2);
	projection = (1ULL << 2);

	uint64_t i = 0;
	// Read the remaining k-3 nucleotides to fill the initial k-mer window aside from last bit
	uint64_t startup_nucs = (k > 4) ? k - 4 : 0;
	for (uint64_t j = 0; j < startup_nucs && i < N; ++j, i += 2) {
		uint8_t nuc = (bits[i] << 1) | bits[i + 1];
		current_kmer = ((current_kmer << 2) | nuc) & dna_k_mask;
		projection = ((projection << 1) | (nuc >> 1)) & proj_k_mask;
	}

	// For simplicity, we jump to STATE_10 and let it handle the first logic gate
	goto STATE_10;

STATE_10:
	while (i < N) {
		uint8_t nuc = (bits[i] << 1) | bits[i + 1];
		current_kmer = ((current_kmer << 2) | nuc) & dna_k_mask;
		projection = ((projection << 1) | (nuc >> 1)) & proj_k_mask;
		i += 2;

		uint64_t next_10 = (projection & ~(projection << 1));
		uint64_t t = 64 - std::countl_zero(next_10);

		uint64_t tail_val = projection & ((1ull << t) - 1);

		sum1 += (t << 32) - tail_val;
		sum2 += current_kmer;


		if (t > 1) {
			// Trivial iterations: next k'-t kmers are INF
			//uint32_t trivial_limit = k_prime - t;
			uint64_t trivial_limit = i + 2 * (k_prime - t);
			//for (uint32_t j = 0; j < trivial_limit && i < N; ++j, i += 2) {
			for (; i < trivial_limit; i += 2)
			{
				uint8_t tnuc = (bits[i] << 1) | bits[i + 1];
				current_kmer = ((current_kmer << 2) | tnuc) & dna_k_mask;
				projection = ((projection << 1) | ((tnuc >> 1) & tnuc)) & proj_k_mask;
				sum1 += INF;
				//sum2 += current_kmer;
			}
			continue; // Remain in STATE_10
		}
		else if (t == 0) {
			goto STATE_0;
		}
		else {
			// t == 1. Compute i_count by counting leading zeros relative to k_prime
			// projection looks like 0^{k'-i} 1^i
			uint64_t i_count = k_prime - (std::countl_zero(projection ^ (1ull << k_prime)) - sixty_four_minus_k);
			// Move to state1 with i_count
			goto STATE_1_ENTRY;
		}

	}
	goto END;

STATE_0:
	while (i < N) {
		uint8_t nuc = (bits[i] << 1) | bits[i + 1];
		current_kmer = ((current_kmer << 2) | nuc) & dna_k_mask;
		projection = ((projection << 1) | (nuc >> 1)) & proj_k_mask;
		i += 2;

		if ((projection & 1) == 0) { // proj(bp) == 0
			sum1 += X4;
			sum2 += current_kmer;
			// Remain in STATE_0
		}
		else {
			sum1 += X2;
			//sum2 += current_kmer;
			// Move to state1(1)
			uint64_t i_count = 1;
			goto STATE_1_ENTRY;
		}
	}
	goto END;

STATE_1_ENTRY:
	{
		uint64_t i_count = 1;
		while (i < N) {
			uint8_t nuc = (bits[i] << 1) | bits[i + 1];
			current_kmer = ((current_kmer << 2) | nuc) & dna_k_mask;
			projection = ((projection << 1) | (nuc >> 1)) & proj_k_mask;
			i += 2;

			if ((projection & 1) == 0) { // proj(bp) == 0
				sum1 += X1;
				//sum2 += current_kmer;
				// Trivial iterations
				for (uint64_t j = 0; j < k_prime - 2 && i < N; ++j, i += 2) {
					uint8_t tnuc = (bits[i] << 1) | bits[i + 1];
					current_kmer = ((current_kmer << 2) | tnuc) & dna_k_mask;
					projection = ((projection << 1) | (tnuc >> 1)) & proj_k_mask;
					sum1 += INF;
					//sum2 += current_kmer;
				}
				goto STATE_10;
			}
			else {
				if (i_count < k_prime) {
					sum1 += INF;
					//sum2 += current_kmer;
					i_count++;
				}
				else {
					sum1 += X3;
					//sum2 += current_kmer;
				}
			}
		}
	}

END:
	auto end = std::chrono::high_resolution_clock::now();
	std::cout << "sum:\t" << sum1 << sum2 << std::endl;
	auto elapsed_us = std::chrono::duration_cast<std::chrono::microseconds>(end - start).count();
	std::cout << "Time taken 10-min random lex (ms): " << elapsed_us / 1000.0 << "\n";
	return elapsed_us / 1000.0;
}

double test_runtime_one_zero_states_random_hybrid_lex_3(uint32_t w, uint32_t k, uint64_t N, const std::vector<bool>& bits) {
	// Primary Key Constants
	const uint64_t X1 = 1ULL << 62; // p ends in 10
	const uint64_t X2 = 2ULL << 62; // p = 0...01
	const uint64_t X3 = 3ULL << 62; // p = 1...1
	const uint64_t X4 = 4ULL << 62; // p = 0...0
	const uint64_t INF = 5ULL << 62; // Others





	const uint64_t dna_k_mask = (1ull << (2 * k)) - 1;
	const uint64_t proj_k_mask = (1ull << (k - 1)) - 1; // only remember last k-1 bits, and not k - to discard leading 1
	const uint64_t k_prime = k - 1;
	const uint64_t sixty_four_minus_k = 64 - k;

	uint64_t sum1 = 0, sum2 = 0;
	uint64_t current_kmer = 0;
	uint64_t projection = 0;

	auto start = std::chrono::high_resolution_clock::now();

	// --- Manual Initialization of first two nucleotides to be 10

	current_kmer = (3ULL << 2);
	projection = (1ULL << 2);

	uint64_t i = 0;
	// Read the remaining k-3 nucleotides to fill the initial k-mer window aside from last bit
	uint64_t startup_nucs = (k > 4) ? k - 4 : 0;
	for (uint64_t j = 0; j < startup_nucs && i < N; ++j, i += 2) {
		uint8_t nuc = (bits[i] << 1) | bits[i + 1];
		current_kmer = ((current_kmer << 2) | nuc) & dna_k_mask;
		projection = ((projection << 1) | ((nuc >> 1) & nuc)) & proj_k_mask;
	}

	// For simplicity, we jump to STATE_10 and let it handle the first logic gate
	goto STATE_10;

STATE_10:
	while (i < N) {
		uint8_t nuc = (bits[i] << 1) | bits[i + 1];
		current_kmer = ((current_kmer << 2) | nuc) & dna_k_mask;
		projection = ((projection << 1) | ((nuc >> 1) & nuc)) & proj_k_mask;
		i += 2;

		uint64_t next_10 = (projection & ~(projection << 1));
		uint64_t t = 64 - std::countl_zero(next_10);

		uint64_t tail_val = projection & ((1ull << t) - 1);

		sum1 += (t << 32) - tail_val;
		sum2 += current_kmer;


		if (t > 1) {
			// Trivial iterations: next k'-t kmers are INF
			//uint32_t trivial_limit = k_prime - t;
			uint64_t trivial_limit = i + 2 * (k_prime - t);
			//for (uint32_t j = 0; j < trivial_limit && i < N; ++j, i += 2) {
			for (; i < trivial_limit; i += 2)
			{
				uint8_t tnuc = (bits[i] << 1) | bits[i + 1];
				current_kmer = ((current_kmer << 2) | tnuc) & dna_k_mask;
				projection = ((projection << 1) | ((tnuc >> 1) & tnuc)) & proj_k_mask;
				sum1 += INF;
				//sum2 += current_kmer;
			}
			continue; // Remain in STATE_10
		}
		else if (t != 0) {
			// t == 1. Compute i_count by counting leading zeros relative to k_prime
			// projection looks like 0^{k'-i} 1^i
			uint64_t i_count = k_prime - (std::countl_zero(projection ^ (1ull << k_prime)) - sixty_four_minus_k);
			// Move to state1 with i_count
			goto STATE_1_ENTRY;
		}
		else {
			goto STATE_0;
		}

	}
	goto END;

STATE_0:
	while (i < N) {
		uint8_t nuc = (bits[i] << 1) | bits[i + 1];
		current_kmer = ((current_kmer << 2) | nuc) & dna_k_mask;
		projection = ((projection << 1) | ((nuc >> 1) & nuc)) & proj_k_mask;
		i += 2;

		if ((projection & 1) == 0) { // proj(bp) == 0
			sum1 += X4;
			sum2 += current_kmer;
			// Remain in STATE_0
		}
		else {
			sum1 += X2;
			//sum2 += current_kmer;
			// Move to state1(1)
			uint64_t i_count = 1;
			goto STATE_1_ENTRY;
		}
	}
	goto END;

STATE_1_ENTRY:
	{
		uint64_t i_count = 1;
		while (i < N) {
			uint8_t nuc = (bits[i] << 1) | bits[i + 1];
			current_kmer = ((current_kmer << 2) | nuc) & dna_k_mask;
			projection = ((projection << 1) | ((nuc >> 1) & nuc)) & proj_k_mask;
			i += 2;

			if ((projection & 1) == 0) { // proj(bp) == 0
				sum1 += X1;
				//sum2 += current_kmer;
				// Trivial iterations
				for (uint64_t j = 0; j < k_prime - 2 && i < N; ++j, i += 2) {
					uint8_t tnuc = (bits[i] << 1) | bits[i + 1];
					current_kmer = ((current_kmer << 2) | tnuc) & dna_k_mask;
					projection = ((projection << 1) | ((tnuc >> 1) & tnuc)) & proj_k_mask;
					sum1 += INF;
					//sum2 += current_kmer;
				}
				goto STATE_10;
			}
			else {
				if (i_count < k_prime) {
					sum1 += INF;
					//sum2 += current_kmer;
					i_count++;
				}
				else {
					sum1 += X3;
					//sum2 += current_kmer;
				}
			}
		}
	}

END:
	auto end = std::chrono::high_resolution_clock::now();
	std::cout << "sum:\t" << sum1 << sum2 << std::endl;
	auto elapsed_us = std::chrono::duration_cast<std::chrono::microseconds>(end - start).count();
	std::cout << "Time taken 10-min random lex (ms): " << elapsed_us / 1000.0 << "\n";
	return elapsed_us / 1000.0;
}

double test_runtime_double_decycling_optimized_lex(uint32_t w, uint32_t k, uint64_t N, const std::vector<bool>& bits) {
	uint64_t dna_k_mask = (k == 32) ? ~0ULL : (1ULL << (2 * k)) - 1;

	double M_PI = 3.14159265358979323846;
	const double EPS = 1e-8;

	uint64_t sum1 = 0;
	uint64_t sum2 = 0;

	// 1. Precompute weights for O(1) sliding window
	// weights[nucleotide][position_in_kmer]
	double u = 2.0 * M_PI / (double)k;
	std::vector<std::vector<double>> weights(4, std::vector<double>(k));
	for (int j = 0; j < 4; j++) {
		for (int i = 0; i < k; i++) {
			weights[j][i] = std::sin(i * u) * j;
		}
	}

	std::vector<uint8_t> kmer_window(k, 0);
	uint64_t window_ptr = 0; // Circular buffer pointer
	double current_weight_sum = 0.0;
	uint64_t current_kmer = 0;

	auto start = std::chrono::high_resolution_clock::now();

	for (uint64_t n = 0; n < N; n += 2) {
		uint8_t incoming_nuc = (bits[n] << 1) | bits[n + 1];

		current_kmer = ((current_kmer << 2) | incoming_nuc) & dna_k_mask;
		kmer_window[window_ptr] = incoming_nuc;
		window_ptr = (window_ptr + 1 == k) ? 0 : window_ptr + 1;

		// --- Core Membership Logic ---
		// We calculate 'sum' and 'shiftsum'. 
		// 'shiftsum' is just the sum of the k-mer minus the oldest nuc.
		double sum = 0.0;
		double shiftsum = 0.0;

		for (uint64_t i = 0; i < k - 1; ++i) {
			uint64_t idx = (window_ptr + i) % k;
			sum += weights[kmer_window[idx]][i];
		}

		for (uint64_t i = 0; i < k - 1; ++i) {
			uint64_t idx = (window_ptr + i + 1) % k;
			shiftsum += weights[kmer_window[idx]][i];
		}

		unsigned setMembership = 0;
		if (sum > EPS && shiftsum < -EPS) {
			sum1+= 1; // assign some set membership
		}
		else if (sum < -EPS && shiftsum > EPS) {
			sum1 += 2; // assign some set membership
		}
		else if (std::abs(sum) <= EPS && std::abs(shiftsum) <= EPS) {
			// can handle this case, for simplicity just assign some set membersip and ingore claculation
			sum1 += 3;  // assign some set membership
		}
		else {
			sum1 += 4;  // assign some set membership
		}

		sum2 += current_kmer;
	}

	auto end = std::chrono::high_resolution_clock::now();
	double elapsed_ms = std::chrono::duration_cast<std::chrono::microseconds>(end - start).count() / 1000.0;

	std::cout << "Sum : " << sum1 << sum2 << std::endl;
	std::cout << "Time taken DD lex (ms): " << elapsed_ms  << "\n";
	return elapsed_ms;
}

double test_runtime_oc_minimizer_lex(uint32_t w, uint32_t k, uint64_t N, const std::vector<bool>& bits) {
	const uint64_t s = 4;
	const uint64_t v = (k - s) / 2;

	const uint64_t w_prime = k - s + 1;

	uint64_t dna_k_mask = (k == 32) ? ~0ULL : (1ULL << (2 * k)) - 1;
	uint64_t s_mask = (1ULL << (2 * s)) - 1;


	uint64_t sum1 = 0, sum2 = 0;

	// Circular buffers
	 uint64_t smer_hashes[128] = { 0 };
	const uint64_t BUF_MASK = 127;
	uint64_t ptr = 0;

	// Track the current minimum s-mer rank and its absolute index
	uint64_t min_s_val = 0xFFFFFFFFFFFFFFFFULL;
	uint64_t min_s_idx = 0; // Absolute position in the stream

	uint64_t current_kmer = 0;
	uint64_t current_smer = 0;



	auto start = std::chrono::high_resolution_clock::now();


	for (uint64_t n = 0; n < N; n += 2) {
		uint64_t current_idx = n / 2;
		uint8_t nuc = (bits[n] << 1) | bits[n + 1];
		current_kmer = ((current_kmer << 2) | nuc) & dna_k_mask;
		current_smer = ((current_smer << 2) | nuc) & s_mask;

		// 1. Update the circular buffer with the new s-mer hash
		smer_hashes[ptr] = current_smer;

		// 2. Optimized Sliding Window Minimum
		// Check if the old minimum has expired (slid out of the k-mer window)
		if (min_s_idx <= (current_idx >= w_prime ? current_idx - w_prime : 0)) {
			// Expired: Must re-scan the entire w' window
			min_s_val = 0xFFFFFFFFFFFFFFFFULL;
			for (uint64_t j = 0; j < w_prime; ++j) {
				uint64_t lookup_idx = (ptr - (w_prime - 1 - j)) & BUF_MASK;
				if (smer_hashes[lookup_idx] < min_s_val) {
					min_s_val = smer_hashes[lookup_idx];
					min_s_idx = current_idx - (w_prime - 1 - j);
				}
			}
		}
		else {
			// Not expired: Just compare with the incoming s-mer
			if (current_smer < min_s_val) {
				min_s_val = current_smer;
				min_s_idx = current_idx;
			}
		}

		// 3. Classify as Open/Closed Syncmer based on relative position
		// The relative position in the k-mer is: (min_s_idx) - (current_idx - w_prime + 1)
		uint64_t relative_pos = min_s_idx - (current_idx - w_prime + 1);
		if (relative_pos == v) {
			sum1 += 1; // Open Syncmer
		}
		else if (relative_pos == 0 || relative_pos == w_prime - 1) {
			sum1 += 2; // Closed Syncmer
		}
		else {
			sum1 += 3;
		}

		sum2 += current_kmer;

		ptr = (ptr + 1) & BUF_MASK;
	}

	auto end = std::chrono::high_resolution_clock::now();
	double elapsed_ms = std::chrono::duration_cast<std::chrono::microseconds>(end - start).count() / 1000.0;

	std::cout << "Sums: " << sum1 << " " << sum2 << std::endl;
	std::cout << "Time taken OC lex (ms): " << elapsed_ms << "\n";

	return elapsed_ms;
}

double test_runtime_miniception_lex(uint32_t w, uint32_t k, uint64_t N, const std::vector<bool>& bits) {
	const uint64_t s = 4;
	const uint64_t w_prime = k - s + 1;

	uint64_t dna_k_mask = (k == 32) ? ~0ULL : (1ULL << (2 * k)) - 1;
	uint64_t s_mask = (1ULL << (2 * s)) - 1;


	uint64_t sum1 = 0, sum2 = 0;

	// Circular buffers
	uint64_t smer_hashes[128] = { 0 };
	const uint64_t BUF_MASK = 127;
	uint64_t ptr = 0;

	// Track the current minimum s-mer rank and its absolute index
	uint64_t min_s_val = 0xFFFFFFFFFFFFFFFFULL;
	uint64_t min_s_idx = 0; // Absolute position in the stream

	uint64_t current_kmer = 0;
	uint64_t current_smer = 0;



	auto start = std::chrono::high_resolution_clock::now();


	for (uint64_t n = 0; n < N; n += 2) {
		uint64_t current_idx = n / 2;
		uint8_t nuc = (bits[n] << 1) | bits[n + 1];
		current_kmer = ((current_kmer << 2) | nuc) & dna_k_mask;
		current_smer = ((current_smer << 2) | nuc) & s_mask;

		// 1. Update the circular buffer with the new s-mer hash
		smer_hashes[ptr] = current_smer;

		// 2. Optimized Sliding Window Minimum
		// Check if the old minimum has expired (slid out of the k-mer window)
		if (min_s_idx <= (current_idx >= w_prime ? current_idx - w_prime : 0)) {
			// Expired: Must re-scan the entire w' window
			min_s_val = 0xFFFFFFFFFFFFFFFFULL;
			for (uint64_t j = 0; j < w_prime; ++j) {
				uint64_t lookup_idx = (ptr - (w_prime - 1 - j)) & BUF_MASK;
				if (smer_hashes[lookup_idx] < min_s_val) {
					min_s_val = smer_hashes[lookup_idx];
					min_s_idx = current_idx - (w_prime - 1 - j);
				}
			}
		}
		else {
			// Not expired: Just compare with the incoming s-mer
			if (current_smer < min_s_val) {
				min_s_val = current_smer;
				min_s_idx = current_idx;
			}
		}

		// 3. Classify as Closed Syncmer based on relative position
		uint64_t relative_pos = min_s_idx - (current_idx - w_prime + 1);
		if (relative_pos == 0 || relative_pos == w_prime - 1) {
			sum1 +=  1; // Closed Syncmer
		}
		else {
			sum1 += 2;
		}

		sum2 += current_kmer;

		ptr = (ptr + 1) & BUF_MASK;
	}

	auto end = std::chrono::high_resolution_clock::now();
	double elapsed_ms = std::chrono::duration_cast<std::chrono::microseconds>(end - start).count() / 1000.0;

	std::cout << "sums : " << sum1 << " " << sum2 << std::endl;
	std::cout << "Time taken MC lex (ms): " << elapsed_ms << "\n";

	return elapsed_ms;
}



double test_runtime_baseline_std(uint32_t w, uint32_t k, uint64_t N, const std::vector<bool>& bits) {
	uint64_t dna_k_mask = (1ull << (2 * k)) - 1;
	static constexpr uint64_t random_number = 0x5B8A1C4D9F2E7360ULL;
	std::hash<uint64_t> hash_fn; // Create a hash function instance
	uint64_t sum = 0;

	// Start timing
	auto start = std::chrono::high_resolution_clock::now();

	uint64_t current_kmer = 0;
	for (uint64_t i = 0; i < N; i += 2) {
		current_kmer = (current_kmer << 1) & dna_k_mask;
		current_kmer |= (bits[i] << 1) | bits[i + 1];
		sum += hash_fn(current_kmer ^ random_number);
	}

	auto end = std::chrono::high_resolution_clock::now();

	std::cout << "sum:\t" << sum << std::endl;

	// Measure elapsed time in nanoseconds
	auto elapsed_ns = std::chrono::duration_cast<std::chrono::nanoseconds>(end - start).count();
	// Or microseconds
	auto elapsed_us = std::chrono::duration_cast<std::chrono::microseconds>(end - start).count();

	std::cout << "Time taken baseline std (ms):  " << elapsed_us / 1000.0 << "\n";

	return elapsed_us / 1000.0;
}

double test_runtime_one_zero_states_random_std(uint32_t w, uint32_t k, uint64_t N, const std::vector<bool>& bits) {
	// Primary Key Constants
	const uint64_t X1 = 1ULL << 62; // p ends in 10
	const uint64_t X2 = 2ULL << 62; // p = 0...01
	const uint64_t X3 = 3ULL << 62; // p = 1...1
	const uint64_t X4 = 4ULL << 62; // p = 0...0
	const uint64_t INF = 5ULL << 62; // Others

	//std::hash<uint32_t> hash_fn_small;
	static constexpr uint64_t random_number = 0x5B8A1C4D9F2E7360ULL;
	std::hash<uint64_t> hash_fn;





	const uint64_t dna_k_mask = (1ull << (2 * k)) - 1;
	const uint64_t proj_k_mask = (1ull << (k - 1)) - 1; // only remember last k-1 bits, and not k - to discard leading 1
	const uint64_t k_prime = k - 1;
	const uint64_t sixty_four_minus_k = 64 - k;

	uint64_t sum1 = 0, sum2 = 0;
	uint64_t current_kmer = 0;
	uint64_t projection = 0;

	auto start = std::chrono::high_resolution_clock::now();

	// --- Manual Initialization of first two nucleotides to be 10

	current_kmer = (3ULL << 2);
	projection = (1ULL << 2);

	uint64_t i = 0;
	// Read the remaining k-3 nucleotides to fill the initial k-mer window aside from last bit
	uint64_t startup_nucs = (k > 4) ? k - 4 : 0;
	for (uint64_t j = 0; j < startup_nucs && i < N; ++j, i += 2) {
		uint8_t nuc = (bits[i] << 1) | bits[i + 1];
		current_kmer = ((current_kmer << 2) | nuc) & dna_k_mask;
		projection = ((projection << 1) | (nuc >> 1)) & proj_k_mask;
	}

	// For simplicity, we jump to STATE_10 and let it handle the first logic gate
	goto STATE_10;

STATE_10:
	while (i < N) {
		uint8_t nuc = (bits[i] << 1) | bits[i + 1];
		current_kmer = ((current_kmer << 2) | nuc) & dna_k_mask;
		projection = ((projection << 1) | (nuc >> 1)) & proj_k_mask;
		i += 2;

		uint64_t next_10 = (projection & ~(projection << 1));
		uint64_t t = 64ULL - std::countl_zero(next_10);

		sum1 += hash_fn(projection ^ random_number); //// In real implementation use std::hash 32 bits, on 32 bit projection so first 32 bits are 0
		sum2 += current_kmer;

		if (t > 1) {
			// Trivial iterations: next k'-t kmers are INF
			//uint32_t trivial_limit = k_prime - t;
			uint64_t trivial_limit = i + 2*(k_prime - t);
			//for (uint32_t j = 0; j < trivial_limit && i < N; ++j, i += 2) {
			for (; i < trivial_limit; i += 2)
			{
				uint8_t tnuc = (bits[i] << 1) | bits[i + 1];
				current_kmer = ((current_kmer << 2) | tnuc) & dna_k_mask;
				projection = ((projection << 1) | ((tnuc >> 1) & tnuc)) & proj_k_mask;
				sum1 += INF;
				//sum2 += current_kmer;
			}
			continue; // Remain in STATE_10
		}
		else if (t != 0) {
			// t == 1. Compute i_count by counting leading zeros relative to k_prime
			// projection looks like 0^{k'-i} 1^i
			uint64_t i_count = k_prime - (std::countl_zero(projection ^ (1ull << k_prime)) - sixty_four_minus_k);
			// Move to state1 with i_count
			goto STATE_1_ENTRY;
		}
		else {
			goto STATE_0;
		}

	}
	goto END;

STATE_0:
	while (i < N) {
		uint8_t nuc = (bits[i] << 1) | bits[i + 1];
		current_kmer = ((current_kmer << 2) | nuc) & dna_k_mask;
		projection = ((projection << 1) | (nuc >> 1)) & proj_k_mask;
		i += 2;

		if ((projection & 1) == 0) { // proj(bp) == 0
			sum1 += X4;
			sum2 += current_kmer;
			// Remain in STATE_0
		}
		else {
			sum1 += X2;
			//sum2 += current_kmer;
			// Move to state1(1)
			uint64_t i_count = 1;
			goto STATE_1_ENTRY;
		}
	}
	goto END;

STATE_1_ENTRY:
	{
		uint64_t i_count = 1;
		while (i < N) {
			uint8_t nuc = (bits[i] << 1) | bits[i + 1];
			current_kmer = ((current_kmer << 2) | nuc) & dna_k_mask;
			projection = ((projection << 1) | (nuc >> 1)) & proj_k_mask;
			i += 2;

			if ((projection & 1) == 0) { // proj(bp) == 0
				sum1 += X1;
				//sum2 += current_kmer;
				// Trivial iterations
				for (uint64_t j = 0; j < k_prime - 2 && i < N; ++j, i += 2) {
					uint8_t tnuc = (bits[i] << 1) | bits[i + 1];
					current_kmer = ((current_kmer << 2) | tnuc) & dna_k_mask;
					projection = ((projection << 1) | (tnuc >> 1)) & proj_k_mask;
					sum1 += INF;
					//sum2 += current_kmer;
				}
				goto STATE_10;
			}
			else {
				if (i_count < k_prime) {
					sum1 += INF;
					//sum2 += current_kmer;
					i_count++;
				}
				else {
					sum1 += X3;
					//sum2 += current_kmer;
				}
			}
		}
	}

END:
	auto end = std::chrono::high_resolution_clock::now();
	std::cout << "sum:\t" << sum1 << sum2 << std::endl;
	auto elapsed_us = std::chrono::duration_cast<std::chrono::microseconds>(end - start).count();
	std::cout << "Time taken 10-min random std (ms): " << elapsed_us / 1000.0 << "\n";
	return elapsed_us / 1000.0;
}

double test_runtime_double_decycling_optimized_std(uint32_t w, uint32_t k, uint64_t N, const std::vector<bool>& bits) {
	uint64_t dna_k_mask = (k == 32) ? ~0ULL : (1ULL << (2 * k)) - 1;
	static constexpr uint64_t random_number = 0x5B8A1C4D9F2E7360ULL;

	double M_PI = 3.14159265358979323846;
	const double EPS = 1e-8;
	std::hash<uint64_t> hash_fn;
	uint64_t sum1 = 0;
	uint64_t sum2 = 0;

	// 1. Precompute weights for O(1) sliding window
	// weights[nucleotide][position_in_kmer]
	double u = 2.0 * M_PI / (double)k;
	std::vector<std::vector<double>> weights(4, std::vector<double>(k));
	for (int j = 0; j < 4; j++) {
		for (int i = 0; i < k; i++) {
			weights[j][i] = std::sin(i * u) * j;
		}
	}

	std::vector<uint8_t> kmer_window(k, 0);
	uint64_t window_ptr = 0; // Circular buffer pointer
	double current_weight_sum = 0.0;
	uint64_t current_kmer = 0;

	auto start = std::chrono::high_resolution_clock::now();

	for (uint64_t n = 0; n < N; n += 2) {
		uint8_t incoming_nuc = (bits[n] << 1) | bits[n + 1];

		current_kmer = ((current_kmer << 2) | incoming_nuc) & dna_k_mask;
		kmer_window[window_ptr] = incoming_nuc;
		window_ptr = (window_ptr + 1 == k) ? 0 : window_ptr + 1;

		// --- Core Membership Logic ---
		// We calculate 'sum' and 'shiftsum'. 
		// 'shiftsum' is just the sum of the k-mer minus the oldest nuc.
		double sum = 0.0;
		double shiftsum = 0.0;

		for (uint64_t i = 0; i < k - 1; ++i) {
			uint64_t idx = (window_ptr + i) % k;
			sum += weights[kmer_window[idx]][i];
		}

		for (uint64_t i = 0; i < k - 1; ++i) {
			uint64_t idx = (window_ptr + i + 1) % k;
			shiftsum += weights[kmer_window[idx]][i];
		}

		unsigned setMembership = 0;
		if (sum > EPS && shiftsum < -EPS) {
			sum1 += 1; // assign some set membership
		}
		else if (sum < -EPS && shiftsum > EPS) {
			sum1 += 2; // assign some set membership
		}
		else if (std::abs(sum) <= EPS && std::abs(shiftsum) <= EPS) {
			// can handle this case, for simplicity just assign some set membersip and ingore claculation
			sum1 += 3;  // assign some set membership
		}
		else {
			sum1 += 4;  // assign some set membership
		}

		sum2 += hash_fn(current_kmer ^ random_number);
	}

	auto end = std::chrono::high_resolution_clock::now();
	double elapsed_ms = std::chrono::duration_cast<std::chrono::microseconds>(end - start).count() / 1000.0;

	std::cout << "Sum : " << sum1 << sum2 << std::endl;
	std::cout << "Time taken DD std: " << elapsed_ms << "\n";
	return elapsed_ms;
}

double test_runtime_oc_minimizer_std(uint32_t w, uint32_t k, uint64_t N, const std::vector<bool>& bits) {
	const uint64_t s = 4;
	const uint64_t v = (k - s) / 2;

	const uint64_t w_prime = k - s + 1;

	uint64_t dna_k_mask = (k == 32) ? ~0ULL : (1ULL << (2 * k)) - 1;
	uint64_t s_mask = (1ULL << (2 * s)) - 1;


	std::hash<uint64_t> hash_fn;
	static constexpr uint64_t random_number = 0x5B8A1C4D9F2E7360ULL;
	uint64_t sum1 = 0, sum2 = 0;

	// Circular buffers
	uint64_t smer_hashes[128] = { 0 };
	const uint64_t BUF_MASK = 127;
	uint64_t ptr = 0;

	// Track the current minimum s-mer rank and its absolute index
	uint64_t min_s_val = 0xFFFFFFFFFFFFFFFFULL;
	uint64_t min_s_idx = 0; // Absolute position in the stream

	uint64_t current_kmer = 0;
	uint64_t current_smer = 0;



	auto start = std::chrono::high_resolution_clock::now();


	for (uint64_t n = 0; n < N; n += 2) {
		uint64_t current_idx = n / 2;
		uint8_t nuc = (bits[n] << 1) | bits[n + 1];
		current_kmer = ((current_kmer << 2) | nuc) & dna_k_mask;
		current_smer = ((current_smer << 2) | nuc) & s_mask;

		// 1. Update the circular buffer with the new s-mer hash
		uint64_t current_h = hash_fn(current_smer ^ random_number);
		smer_hashes[ptr] = current_h;

		// 2. Optimized Sliding Window Minimum
		// Check if the old minimum has expired (slid out of the k-mer window)
		if (min_s_idx <= (current_idx >= w_prime ? current_idx - w_prime : 0)) {
			// Expired: Must re-scan the entire w' window
			min_s_val = 0xFFFFFFFFFFFFFFFFULL;
			for (uint64_t j = 0; j < w_prime; ++j) {
				uint64_t lookup_idx = (ptr - (w_prime - 1 - j)) & BUF_MASK;
				if (smer_hashes[lookup_idx] < min_s_val) {
					min_s_val = smer_hashes[lookup_idx];
					min_s_idx = current_idx - (w_prime - 1 - j);
				}
			}
		}
		else {
			// Not expired: Just compare with the incoming s-mer
			if (current_h < min_s_val) {
				min_s_val = current_h;
				min_s_idx = current_idx;
			}
		}

		// 3. Classify as Open/Closed Syncmer based on relative position
		// The relative position in the k-mer is: (min_s_idx) - (current_idx - w_prime + 1)
		uint64_t relative_pos = min_s_idx - (current_idx - w_prime + 1);
		if (relative_pos == v) {
			sum1 += 1; // Open Syncmer
		}
		else if (relative_pos == 0 || relative_pos == w_prime - 1) {
			sum1 += 2; // Closed Syncmer
		}
		else {
			sum1 += 3;
		}

		sum2 += hash_fn(current_kmer ^ random_number);

		ptr = (ptr + 1) & BUF_MASK;
	}

	auto end = std::chrono::high_resolution_clock::now();
	double elapsed_ms = std::chrono::duration_cast<std::chrono::microseconds>(end - start).count() / 1000.0;

	std::cout << "Sums: " << sum1 << " " << sum2 << std::endl;
	std::cout << "Time taken OC std (ms): " << elapsed_ms << "\n";

	return elapsed_ms;
}

double test_runtime_miniception_std(uint32_t w, uint32_t k, uint64_t N, const std::vector<bool>& bits) {
	const uint64_t s = 4;
	const uint64_t w_prime = k - s + 1;

	uint64_t dna_k_mask = (k == 32) ? ~0ULL : (1ULL << (2 * k)) - 1;
	uint64_t s_mask = (1ULL << (2 * s)) - 1;

	std::hash<uint64_t> hash_fn;
	static constexpr uint64_t random_number = 0x5B8A1C4D9F2E7360ULL;

	uint64_t sum1 = 0, sum2 = 0;

	// Circular buffers
	uint64_t smer_hashes[128] = { 0 };
	const uint64_t BUF_MASK = 127;
	uint64_t ptr = 0;

	// Track the current minimum s-mer rank and its absolute index
	uint64_t min_s_val = 0xFFFFFFFFFFFFFFFFULL;
	uint64_t min_s_idx = 0; // Absolute position in the stream

	uint64_t current_kmer = 0;
	uint64_t current_smer = 0;



	auto start = std::chrono::high_resolution_clock::now();


	for (uint64_t n = 0; n < N; n += 2) {
		uint64_t current_idx = n / 2;
		uint8_t nuc = (bits[n] << 1) | bits[n + 1];
		current_kmer = ((current_kmer << 2) | nuc) & dna_k_mask;
		current_smer = ((current_smer << 2) | nuc) & s_mask;

		// 1. Update the circular buffer with the new s-mer hash
		uint64_t current_h = hash_fn(current_smer ^ random_number);
		smer_hashes[ptr] = current_h;

		// 2. Optimized Sliding Window Minimum
		// Check if the old minimum has expired (slid out of the k-mer window)
		if (min_s_idx <= (current_idx >= w_prime ? current_idx - w_prime : 0)) {
			// Expired: Must re-scan the entire w' window
			min_s_val = 0xFFFFFFFFFFFFFFFFULL;
			for (uint64_t j = 0; j < w_prime; ++j) {
				uint64_t lookup_idx = (ptr - (w_prime - 1 - j)) & BUF_MASK;
				if (smer_hashes[lookup_idx] < min_s_val) {
					min_s_val = smer_hashes[lookup_idx];
					min_s_idx = current_idx - (w_prime - 1 - j);
				}
			}
		}
		else {
			// Not expired: Just compare with the incoming s-mer
			if (current_h < min_s_val) {
				min_s_val = current_h;
				min_s_idx = current_idx;
			}
		}

		// 3. Classify as Closed Syncmer based on relative position
		uint64_t relative_pos = min_s_idx - (current_idx - w_prime + 1);
		if (relative_pos == 0 || relative_pos == w_prime - 1) {
			sum1 += 1; // Closed Syncmer
		}
		else {
			sum1 += 2;
		}

		sum2 += hash_fn(current_kmer ^ random_number);

		ptr = (ptr + 1) & BUF_MASK;
	}

	auto end = std::chrono::high_resolution_clock::now();
	double elapsed_ms = std::chrono::duration_cast<std::chrono::microseconds>(end - start).count() / 1000.0;

	std::cout << "sums : " << sum1 << " " << sum2 << std::endl;
	std::cout << "Time taken MC std (ms): " << elapsed_ms << "\n";

	return elapsed_ms;
}



double test_runtime_baseline_pcg(uint32_t w, uint32_t k, uint64_t N, const std::vector<bool>& bits) {
	uint64_t dna_k_mask = (1ull << (2 * k)) - 1;
	static constexpr uint64_t random_number = 0x5B8A1C4D9F2E7360ULL;

	uint64_t sum = 0;

	// Start timing
	auto start = std::chrono::high_resolution_clock::now();

	uint64_t current_kmer = 0;
	for (uint64_t i = 0; i < N; i += 2) {
		current_kmer = (current_kmer << 1) & dna_k_mask;
		current_kmer |= (bits[i] << 1) | bits[i + 1];
		sum += pcg_mix_64(current_kmer ^ random_number);
	}

	auto end = std::chrono::high_resolution_clock::now();

	std::cout << "sum:\t" << sum << std::endl;

	// Measure elapsed time in nanoseconds
	auto elapsed_ns = std::chrono::duration_cast<std::chrono::nanoseconds>(end - start).count();
	// Or microseconds
	auto elapsed_us = std::chrono::duration_cast<std::chrono::microseconds>(end - start).count();

	std::cout << "Time taken baseline pcg (ms):  " << elapsed_us / 1000.0 << "\n";

	return elapsed_us / 1000.0;
}

double test_runtime_one_zero_states_random_pcg(uint32_t w, uint32_t k, uint64_t N, const std::vector<bool>& bits) {
	// Primary Key Constants
	const uint64_t X1 = 1ULL << 62; // p ends in 10
	const uint64_t X2 = 2ULL << 62; // p = 0...01
	const uint64_t X3 = 3ULL << 62; // p = 1...1
	const uint64_t X4 = 4ULL << 62; // p = 0...0
	const uint64_t INF = 5ULL << 62; // Others


	static constexpr uint64_t random_number = 0x5B8A1C4D9F2E7360ULL;



	const uint64_t dna_k_mask = (1ull << (2 * k)) - 1;
	const uint64_t proj_k_mask = (1ull << (k - 1)) - 1; // only remember last k-1 bits, and not k - to discard leading 1
	const uint64_t k_prime = k - 1;
	const uint64_t sixty_four_minus_k = 64ULL - k;

	uint64_t sum1 = 0, sum2 = 0;
	uint64_t current_kmer = 0;
	uint64_t projection = 0;

	auto start = std::chrono::high_resolution_clock::now();

	// --- Manual Initialization of first two nucleotides to be 10

	current_kmer = (3ULL << 2);
	projection = (1ULL << 2);

	uint64_t i = 0;
	// Read the remaining k-3 nucleotides to fill the initial k-mer window aside from last bit
	uint64_t startup_nucs = (k > 4) ? k - 4 : 0;
	for (uint64_t j = 0; j < startup_nucs && i < N; ++j, i += 2) {
		uint8_t nuc = (bits[i] << 1) | bits[i + 1];
		current_kmer = ((current_kmer << 2) | nuc) & dna_k_mask;
		projection = ((projection << 1) | (nuc >> 1)) & proj_k_mask;
	}

	// For simplicity, we jump to STATE_10 and let it handle the first logic gate
	goto STATE_10;

STATE_10:
	while (i < N) {
		uint8_t nuc = (bits[i] << 1) | bits[i + 1];
		current_kmer = ((current_kmer << 2) | nuc) & dna_k_mask;
		projection = ((projection << 1) | (nuc >> 1)) & proj_k_mask;
		i += 2;

		uint64_t next_10 = (projection & ~(projection << 1));
		uint64_t t = 64ULL - std::countl_zero(next_10);

		sum1 += pcg_mix_64(projection ^ random_number); // In real implementation use pcg_mix_32, on 32 bit projection so first 32 bits are 0
		sum2 += current_kmer;

		if (t > 1) {
			// Trivial iterations: next k'-t kmers are INF
			uint64_t trivial_limit = i + 2*(k_prime - t);
			//for (uint32_t j = 0; j < trivial_limit && i < N; ++j, i += 2) {
			for (; i < trivial_limit; i += 2)
			{
				uint8_t tnuc = (bits[i] << 1) | bits[i + 1];
				current_kmer = ((current_kmer << 2) | tnuc) & dna_k_mask;
				projection = ((projection << 1) | ((tnuc >> 1) & tnuc)) & proj_k_mask;
				sum1 += INF;
				//sum2 += current_kmer;
			}
			continue; // Remain in STATE_10
		}
		else if (t != 0) {
			// t == 1. Compute i_count by counting leading zeros relative to k_prime
			// projection looks like 0^{k'-i} 1^i
			uint64_t i_count = k_prime - (std::countl_zero(projection ^ (1ull << k_prime)) - sixty_four_minus_k);
			// Move to state1 with i_count
			goto STATE_1_ENTRY;
		}
		else {
			goto STATE_0;
		}

	}
	goto END;

STATE_0:
	while (i < N) {
		uint8_t nuc = (bits[i] << 1) | bits[i + 1];
		current_kmer = ((current_kmer << 2) | nuc) & dna_k_mask;
		projection = ((projection << 1) | (nuc >> 1)) & proj_k_mask;
		i += 2;

		if ((projection & 1) == 0) { // proj(bp) == 0
			sum1 += X4;
			sum2 += current_kmer;
			// Remain in STATE_0
		}
		else {
			sum1 += X2;
			//sum2 += current_kmer;
			// Move to state1(1)
			uint64_t i_count = 1;
			goto STATE_1_ENTRY;
		}
	}
	goto END;

STATE_1_ENTRY:
	{
		uint64_t i_count = 1;
		while (i < N) {
			uint8_t nuc = (bits[i] << 1) | bits[i + 1];
			current_kmer = ((current_kmer << 2) | nuc) & dna_k_mask;
			projection = ((projection << 1) | (nuc >> 1)) & proj_k_mask;
			i += 2;

			if ((projection & 1) == 0) { // proj(bp) == 0
				sum1 += X1;
				//sum2 += current_kmer;
				// Trivial iterations
				for (uint64_t j = 0; j < k_prime - 2 && i < N; ++j, i += 2) {
					uint8_t tnuc = (bits[i] << 1) | bits[i + 1];
					current_kmer = ((current_kmer << 2) | tnuc) & dna_k_mask;
					projection = ((projection << 1) | (tnuc >> 1)) & proj_k_mask;
					sum1 += INF;
					//sum2 += current_kmer;
				}
				goto STATE_10;
			}
			else {
				if (i_count < k_prime) {
					sum1 += INF;
					//sum2 += current_kmer;
					i_count++;
				}
				else {
					sum1 += X3;
					//sum2 += current_kmer;
				}
			}
		}
	}

END:
	auto end = std::chrono::high_resolution_clock::now();
	std::cout << "sum:\t" << sum1 << sum2 << std::endl;
	auto elapsed_us = std::chrono::duration_cast<std::chrono::microseconds>(end - start).count();
	std::cout << "Time taken 10-min random pcg (ms): " << elapsed_us / 1000.0 << "\n";
	return elapsed_us / 1000.0;
}

double test_runtime_double_decycling_optimized_pcg(uint32_t w, uint32_t k, uint64_t N, const std::vector<bool>& bits) {
	uint64_t dna_k_mask = (k == 32) ? ~0ULL : (1ULL << (2 * k)) - 1;

	static constexpr uint64_t random_number = 0x5B8A1C4D9F2E7360ULL;

	double M_PI = 3.14159265358979323846;
	const double EPS = 1e-8;
	uint64_t sum1 = 0;
	uint64_t sum2 = 0;

	// 1. Precompute weights for O(1) sliding window
	// weights[nucleotide][position_in_kmer]
	double u = 2.0 * M_PI / (double)k;
	std::vector<std::vector<double>> weights(4, std::vector<double>(k));
	for (int j = 0; j < 4; j++) {
		for (int i = 0; i < k; i++) {
			weights[j][i] = std::sin(i * u) * j;
		}
	}

	std::vector<uint8_t> kmer_window(k, 0);
	uint64_t window_ptr = 0; // Circular buffer pointer
	double current_weight_sum = 0.0;
	uint64_t current_kmer = 0;

	auto start = std::chrono::high_resolution_clock::now();

	for (uint64_t n = 0; n < N; n += 2) {
		uint8_t incoming_nuc = (bits[n] << 1) | bits[n + 1];

		current_kmer = ((current_kmer << 2) | incoming_nuc) & dna_k_mask;
		kmer_window[window_ptr] = incoming_nuc;
		window_ptr = (window_ptr + 1 == k) ? 0 : window_ptr + 1;

		// --- Core Membership Logic ---
		// We calculate 'sum' and 'shiftsum'. 
		// 'shiftsum' is just the sum of the k-mer minus the oldest nuc.
		double sum = 0.0;
		double shiftsum = 0.0;

		for (uint64_t i = 0; i < k - 1; ++i) {
			uint64_t idx = (window_ptr + i) % k;
			sum += weights[kmer_window[idx]][i];
		}

		for (uint64_t i = 0; i < k - 1; ++i) {
			uint64_t idx = (window_ptr + i + 1) % k;
			shiftsum += weights[kmer_window[idx]][i];
		}

		unsigned setMembership = 0;
		if (sum > EPS && shiftsum < -EPS) {
			sum1 += 1; // assign some set membership
		}
		else if (sum < -EPS && shiftsum > EPS) {
			sum1 += 2; // assign some set membership
		}
		else if (std::abs(sum) <= EPS && std::abs(shiftsum) <= EPS) {
			// can handle this case, for simplicity just assign some set membersip and ingore claculation
			sum1 += 3;  // assign some set membership
		}
		else {
			sum1 += 4;  // assign some set membership
		}

		sum2 += pcg_mix_64(current_kmer ^ random_number);
	}

	auto end = std::chrono::high_resolution_clock::now();
	double elapsed_ms = std::chrono::duration_cast<std::chrono::microseconds>(end - start).count() / 1000.0;

	std::cout << "Sum : " << sum1 << sum2 << std::endl;
	std::cout << "Time taken DD pcg (ms): " << elapsed_ms << "\n";
	return elapsed_ms;
}

double test_runtime_oc_minimizer_pcg(uint32_t w, uint32_t k, uint64_t N, const std::vector<bool>& bits) {
	const uint64_t s = 4;
	const uint64_t v = (k - s) / 2;

	const uint64_t w_prime = k - s + 1;

	uint64_t dna_k_mask = (k == 32) ? ~0ULL : (1ULL << (2 * k)) - 1;
	uint64_t s_mask = (1ULL << (2 * s)) - 1;

	static constexpr uint64_t random_number = 0x5B8A1C4D9F2E7360ULL;

	uint64_t sum1 = 0, sum2 = 0;

	// Circular buffers
	 uint64_t smer_hashes[128] = { 0 };
	const uint64_t BUF_MASK = 127;
	uint64_t ptr = 0;

	// Track the current minimum s-mer rank and its absolute index
	uint64_t min_s_val = 0xFFFFFFFFFFFFFFFFULL;
	uint64_t min_s_idx = 0; // Absolute position in the stream

	uint64_t current_kmer = 0;
	uint64_t current_smer = 0;



	auto start = std::chrono::high_resolution_clock::now();


	for (uint64_t n = 0; n < N; n += 2) {
		uint64_t current_idx = n / 2;
		uint8_t nuc = (bits[n] << 1) | bits[n + 1];
		current_kmer = ((current_kmer << 2) | nuc) & dna_k_mask;
		current_smer = ((current_smer << 2) | nuc) & s_mask;

		// 1. Update the circular buffer with the new s-mer hash
		uint64_t current_h = pcg_mix_64(current_smer ^ random_number);
		smer_hashes[ptr] = current_h;

		// 2. Optimized Sliding Window Minimum
		// Check if the old minimum has expired (slid out of the k-mer window)
		if (min_s_idx <= (current_idx >= w_prime ? current_idx - w_prime : 0)) {
			// Expired: Must re-scan the entire w' window
			min_s_val = 0xFFFFFFFFFFFFFFFFULL;
			for (uint64_t j = 0; j < w_prime; ++j) {
				uint64_t lookup_idx = (ptr - (w_prime - 1 - j)) & BUF_MASK;
				if (smer_hashes[lookup_idx] < min_s_val) {
					min_s_val = smer_hashes[lookup_idx];
					min_s_idx = current_idx - (w_prime - 1 - j);
				}
			}
		}
		else {
			// Not expired: Just compare with the incoming s-mer
			if (current_h < min_s_val) {
				min_s_val = current_h;
				min_s_idx = current_idx;
			}
		}

		// 3. Classify as Open/Closed Syncmer based on relative position
		// The relative position in the k-mer is: (min_s_idx) - (current_idx - w_prime + 1)
		uint64_t relative_pos = min_s_idx - (current_idx - w_prime + 1);
		if (relative_pos == v) {
			sum1 += 1; // Open Syncmer
		}
		else if (relative_pos == 0 || relative_pos == w_prime - 1) {
			sum1 += 2; // Closed Syncmer
		}
		else {
			sum1 += 3;
		}

		sum2 += pcg_mix_64(current_kmer ^ random_number);

		ptr = (ptr + 1) & BUF_MASK;
	}

	auto end = std::chrono::high_resolution_clock::now();
	double elapsed_ms = std::chrono::duration_cast<std::chrono::microseconds>(end - start).count() / 1000.0;

	std::cout << "Sums: " << sum1 << " " << sum2 << std::endl;
	std::cout << "Time taken OC pcg (ms): " << elapsed_ms << "\n";

	return elapsed_ms;
}

double test_runtime_miniception_pcg(uint32_t w, uint32_t k, uint64_t N, const std::vector<bool>& bits) {
	const uint64_t s = 4;
	const uint64_t w_prime = k - s + 1;

	uint64_t dna_k_mask = (k == 32) ? ~0ULL : (1ULL << (2 * k)) - 1;
	uint64_t s_mask = (1ULL << (2 * s)) - 1;

	static constexpr uint64_t random_number = 0x5B8A1C4D9F2E7360ULL;

	uint64_t sum1 = 0, sum2 = 0;

	// Circular buffers
	 uint64_t smer_hashes[128] = { 0 };
	const uint64_t BUF_MASK = 127;
	uint64_t ptr = 0;

	// Track the current minimum s-mer rank and its absolute index
	uint64_t min_s_val = 0xFFFFFFFFFFFFFFFFULL;
	uint64_t min_s_idx = 0; // Absolute position in the stream

	uint64_t current_kmer = 0;
	uint64_t current_smer = 0;



	auto start = std::chrono::high_resolution_clock::now();


	for (uint64_t n = 0; n < N; n += 2) {
		uint64_t current_idx = n / 2;
		uint8_t nuc = (bits[n] << 1) | bits[n + 1];
		current_kmer = ((current_kmer << 2) | nuc) & dna_k_mask;
		current_smer = ((current_smer << 2) | nuc) & s_mask;

		// 1. Update the circular buffer with the new s-mer hash
		uint64_t current_h = pcg_mix_64(current_smer ^ random_number);
		smer_hashes[ptr] = current_h;

		// 2. Optimized Sliding Window Minimum
		// Check if the old minimum has expired (slid out of the k-mer window)
		if (min_s_idx <= (current_idx >= w_prime ? current_idx - w_prime : 0)) {
			// Expired: Must re-scan the entire w' window
			min_s_val = 0xFFFFFFFFFFFFFFFFULL;
			for (uint64_t j = 0; j < w_prime; ++j) {
				uint64_t lookup_idx = (ptr - (w_prime - 1 - j)) & BUF_MASK;
				if (smer_hashes[lookup_idx] < min_s_val) {
					min_s_val = smer_hashes[lookup_idx];
					min_s_idx = current_idx - (w_prime - 1 - j);
				}
			}
		}
		else {
			// Not expired: Just compare with the incoming s-mer
			if (current_h < min_s_val) {
				min_s_val = current_h;
				min_s_idx = current_idx;
			}
		}

		// 3. Classify as Closed Syncmer based on relative position
		uint64_t relative_pos = min_s_idx - (current_idx - w_prime + 1);
		if (relative_pos == 0 || relative_pos == w_prime - 1) {
			sum1 += 1; // Closed Syncmer
		}
		else {
			sum1 += 2;
		}

		sum2 += pcg_mix_64(current_kmer ^ random_number);

		ptr = (ptr + 1) & BUF_MASK;
	}

	auto end = std::chrono::high_resolution_clock::now();
	double elapsed_ms = std::chrono::duration_cast<std::chrono::microseconds>(end - start).count() / 1000.0;

	std::cout << "sums : " << sum1 << " " << sum2 << std::endl;
	std::cout << "Time taken MC pcg (ms): " << elapsed_ms << "\n";

	return elapsed_ms;
}











// Simple struct to hold cumulative results
struct Stats {
	std::vector<double> results;
	double average() const {
		if (results.empty()) return 0.0;
		return std::accumulate(results.begin(), results.end(), 0.0) / results.size();
	}
};

// 1. Define the Function Signature type for readability
using SchemeFunc = std::function<double(uint32_t, uint32_t, uint64_t, const std::vector<bool>&)>;

// 2. Structure to hold everything about a specific sampling method
struct BenchmarkScheme {
	std::string name;
	SchemeFunc function;
	bool requires_w_limit; // If true, only runs if w >= k - 2
};

// 3. Unit of work for the shuffled pool
struct KTask { uint32_t k; const BenchmarkScheme* scheme; };
struct WTask { uint32_t w; const BenchmarkScheme* scheme; };

void run_k_scaling_benchmark(uint32_t num_runs, uint64_t N, uint32_t w_fixed, const std::vector<bool>& bits) {
	// --- CONSTANTS ---
	const uint32_t NUM_CORES = 6;
	const uint32_t K_START = 5;
	const uint32_t K_END = 25; // we put 32 for most stuff, but 25 if using all since DD-based ones were VERY slow for 32

	// --- STEP 1: REGISTRY ---
	std::vector<BenchmarkScheme> schemes = {
		{"baseline_lex", test_runtime_baseline_lex, false},
		{"abb_lex", test_runtime_abb_plus_lex, false},
		{"dd_lex", test_runtime_double_decycling_optimized_lex, false},
		{"oc_lex", test_runtime_oc_minimizer_lex, false},
		{"miniception_lex", test_runtime_miniception_lex, false},
		{"states_10_lex", test_runtime_one_zero_states_lex, true},
		{"states_10_random_lex", test_runtime_one_zero_states_random_lex, true},
		{"baseline_std", test_runtime_baseline_std, false},
		{"dd_std", test_runtime_double_decycling_optimized_std, false},
		{"oc_std", test_runtime_oc_minimizer_std, false},
		{"miniception_std", test_runtime_miniception_std, false},
		{"states_10_random_std", test_runtime_one_zero_states_random_std, true},
		{"baseline_pcg", test_runtime_baseline_pcg, false},
		{"dd_pcg", test_runtime_double_decycling_optimized_pcg, false},
		{"oc_pcg", test_runtime_oc_minimizer_pcg, false},
		{"miniception_pcg", test_runtime_miniception_pcg, false},
		{"states_10_random_pcg", test_runtime_one_zero_states_random_pcg, true},
		{"baseline_xor", test_runtime_baseline_xor, false},
		{"prev_states_10_lex", test_runtime_one_zero_most_optimized_lazy, false},
		{"GreedyMini", test_runtime_gm, false}
	};



	// --- STEP 2: TASK POOL GENERATION ---
	std::vector<KTask> task_pool;
	for (uint32_t k = K_START; k <= K_END; ++k) {
		for (const auto& scheme : schemes) {
			for (uint32_t i = 0; i < num_runs; ++i) {
				task_pool.push_back({ k, &scheme });
			}
		}
	}

	std::random_device rd;
	std::mt19937 g(rd());
	std::shuffle(task_pool.begin(), task_pool.end(), g);

	// --- STEP 3: MULTI-THREADED EXECUTION ---
	// We store results in thread-local maps to avoid locking overhead
	using ResultMap = std::map<std::string, std::map<uint32_t, Stats>>;
	std::vector<ResultMap> thread_local_results(NUM_CORES);

	std::atomic<size_t> next_task_idx(0);
	std::vector<std::thread> workers;

	std::cout << "Launching benchmark on " << NUM_CORES << " cores..." << std::endl;
	std::cout << "Total tasks: " << task_pool.size() << std::endl;

	for (uint32_t t = 0; t < NUM_CORES; ++t) {
		workers.emplace_back([&, t]() {
			while (true) {
				size_t idx = next_task_idx.fetch_add(1);
				if (idx >= task_pool.size()) break;

				const KTask& task = task_pool[idx];
				const BenchmarkScheme& s = *(task.scheme);

				// Periodic progress (only one thread handles the print)
				if (t == 0 && idx % 100 == 0) {
					std::cout << "Progress: " << idx << "/" << task_pool.size() << "\r" << std::flush;
				}

				if (s.requires_w_limit && w_fixed < (task.k - 2)) {
					thread_local_results[t][s.name][task.k].results.push_back(-1.0);
				}
				else {
					double time = s.function(w_fixed, task.k, N, bits);
					thread_local_results[t][s.name][task.k].results.push_back(time);
				}
			}
			});
	}

	// Wait for all cores to finish
	for (auto& worker : workers) worker.join();

	// --- STEP 4: MERGE RESULTS ---
	std::map<std::string, std::map<uint32_t, Stats>> all_results;
	for (const auto& local_map : thread_local_results) {
		for (const auto& scheme_entry : local_map) {
			const std::string& scheme_name = scheme_entry.first;
			for (const auto& k_entry : scheme_entry.second) {
				uint32_t k = k_entry.first;
				const auto& stats = k_entry.second;

				// Append thread-local vectors to the global vector
				all_results[scheme_name][k].results.insert(
					all_results[scheme_name][k].results.end(),
					stats.results.begin(),
					stats.results.end()
				);
			}
		}
	}

	// --- STEP 5: CSV EXPORT ---
	std::ofstream ofs("sampling_runtime_k_scaling.csv");
	ofs << "k";
	for (const auto& s : schemes) ofs << "," << s.name;
	ofs << "\n";

	for (uint32_t k = K_START; k <= K_END; ++k) {
		ofs << k;
		for (const auto& s : schemes) {
			ofs << "," << all_results[s.name][k].average();
		}
		ofs << "\n";
	}

	std::cout << "\nMulti-threaded benchmark complete! Results saved." << std::endl;
}

void run_w_scaling_benchmark(uint32_t num_runs, uint64_t N, uint32_t k_fixed, const std::vector<bool>& bits) {
	// 1. Define the registry for W-scaling
	// You can add more schemes here (e.g., OC, MC) to compare them across different window sizes
	std::vector<BenchmarkScheme> schemes = {
		{"states_10_lex", test_runtime_one_zero_states_lex, true}
	};

	uint32_t w_start = (k_fixed > 2) ? (k_fixed - 2) : 3;
	uint32_t w_end = 100;

	// 2. Create the Task Pool (Randomizing both W and the Scheme)
	std::vector<WTask> task_pool;
	for (uint32_t w = w_start; w <= w_end; ++w) {
		for (const auto& scheme : schemes) {
			// Respect the w >= k-2 limit during task generation
			if (scheme.requires_w_limit && w < (k_fixed - 2)) continue;

			for (uint32_t i = 0; i < num_runs; ++i) {
				task_pool.push_back({ w, &scheme });
			}
		}
	}

	std::random_device rd;
	std::mt19937 g(rd());
	std::shuffle(task_pool.begin(), task_pool.end(), g);

	// 3. Storage
	std::map<std::string, std::map<uint32_t, Stats>> all_results;

	std::cout << "Starting randomized w-scaling benchmark (" << task_pool.size() << " tasks)..." << std::endl;

	// 4. Execution
	for (size_t i = 0; i < task_pool.size(); ++i) {
		const WTask& task = task_pool[i];
		const BenchmarkScheme& s = *(task.scheme);

		if (i % 50 == 0) std::cout << "Progress: " << i << "/" << task_pool.size() << "\r" << std::flush;

		double time = s.function(task.w, k_fixed, N, bits);
		all_results[s.name][task.w].results.push_back(time);
	}

	// 5. Automatic CSV Export
	std::ofstream ofs("one_zero_runtime_w_scaling.csv");

	// Header
	ofs << "w";
	for (const auto& s : schemes) ofs << "," << s.name;
	ofs << "\n";

	// Rows
	for (uint32_t w = w_start; w <= w_end; ++w) {
		ofs << w;
		for (const auto& s : schemes) {
			// Using a default of 0.0 if a specific w/scheme combo wasn't run
			double avg = all_results[s.name].count(w) ? all_results[s.name][w].average() : 0.0;
			ofs << "," << avg;
		}
		ofs << "\n";
	}

	std::cout << "\nBenchmark complete. Results saved to one_zero_runtime_w_scaling.csv" << std::endl;
}

namespace fs = std::filesystem;



// Helper to calculate Wilcoxon Rank-Sum p-values (Normal Approximation)
// returns {p_value_oc_gt_mc, p_value_mc_gt_oc}
std::pair<double, double> calculate_wilcoxon(const std::vector<double>& oc, const std::vector<double>& mc) {
	size_t n1 = oc.size();
	size_t n2 = mc.size();
	if (n1 == 0 || n2 == 0) return { 1.0, 1.0 };

	struct Observed { double val; int group; }; // 0 for OC, 1 for MC
	std::vector<Observed> combined;
	for (double v : oc) combined.push_back({ v, 0 });
	for (double v : mc) combined.push_back({ v, 1 });

	std::sort(combined.begin(), combined.end(), [](const Observed& a, const Observed& b) {
		return a.val < b.val;
		});

	std::vector<double> ranks(combined.size());
	for (size_t i = 0; i < combined.size(); ) {
		size_t j = i;
		while (j < combined.size() && combined[j].val == combined[i].val) j++;
		double rank = 1.0 + (i + (j - 1)) / 2.0;
		for (size_t k = i; k < j; k++) ranks[k] = rank;
		i = j;
	}

	double rank_sum_oc = 0;
	for (size_t i = 0; i < combined.size(); i++) {
		if (combined[i].group == 0) rank_sum_oc += ranks[i];
	}

	double U_oc = rank_sum_oc - (n1 * (n1 + 1.0)) / 2.0;
	double mu_U = (n1 * n2) / 2.0;
	double sigma_U = std::sqrt((n1 * n2 * (n1 + n2 + 1.0)) / 12.0);

	// Continuity correction
	double z_oc_gt_mc = (U_oc - 0.5 - mu_U) / sigma_U; // Testing if OC ranks are higher (smaller U)
	double z_mc_gt_oc = (U_oc + 0.5 - mu_U) / sigma_U;

	auto normal_cdf = [](double z) { return 0.5 * std::erfc(-z * std::sqrt(0.5)); };

	// We want the probability of seeing a U this extreme or more
	// Note: Small U means OC is smaller/faster. Large U means OC is larger/slower.
	double p_oc_slower = 1.0 - normal_cdf((U_oc - 0.5 - mu_U) / sigma_U);
	double p_mc_slower = normal_cdf((U_oc + 0.5 - mu_U) / sigma_U);

	return { p_oc_slower, p_mc_slower };
}


void run_oc_vs_mc_benchmark(uint32_t num_runs, uint64_t N, uint32_t w_fixed, const std::vector<bool>& bits) {
	const uint32_t NUM_CORES = 6;
	const uint32_t K_START = 16;
	const uint32_t K_END = 16;
	const std::string filename = "oc_vs_mc_history.csv";

	std::vector<BenchmarkScheme> schemes = {
		{"oc_std", test_runtime_oc_minimizer_std, false},
		{"miniception_std", test_runtime_miniception_std, false}
	};

	// --- STEP 1: TASK POOL ---
	std::vector<KTask> task_pool;
	for (uint32_t k = K_START; k <= K_END; ++k) {
		for (const auto& scheme : schemes) {
			for (uint32_t i = 0; i < num_runs; ++i) task_pool.push_back({ k, &scheme });
		}
	}
	std::shuffle(task_pool.begin(), task_pool.end(), std::mt19937(std::random_device{}()));

	// --- STEP 2: EXECUTION ---
	using ResultMap = std::map<std::string, std::map<uint32_t, Stats>>;
	std::vector<ResultMap> thread_local_results(NUM_CORES);
	std::atomic<size_t> next_task_idx(0);
	std::vector<std::thread> workers;

	for (uint32_t t = 0; t < NUM_CORES; ++t) {
		workers.emplace_back([&, t]() {
			while (true) {
				size_t idx = next_task_idx.fetch_add(1);
				if (idx >= task_pool.size()) break;
				const KTask& task = task_pool[idx];
				double time = task.scheme->function(w_fixed, task.k, N, bits);
				thread_local_results[t][task.scheme->name][task.k].results.push_back(time);
			}
			});
	}
	for (auto& worker : workers) worker.join();

	// --- STEP 3: MERGE ---
	std::map<std::string, std::map<uint32_t, Stats>> all_results;
	for (const auto& local_map : thread_local_results) {
		for (const auto& s_entry : local_map) {
			for (const auto& k_entry : s_entry.second) {
				auto& v = all_results[s_entry.first][k_entry.first].results;
				v.insert(v.end(), k_entry.second.results.begin(), k_entry.second.results.end());
			}
		}
	}

	// --- STEP 4: APPEND TO CSV ---
	bool exists = fs::exists(filename);
	std::ofstream ofs(filename, std::ios::app);

	if (!exists) {
		ofs << "datetime,k,oc_avg,oc_stdev,mc_avg,mc_stdev,p_oc_gt_mc,p_mc_gt_oc\n";
	}

	// Get current time to the minute
	auto now_time_t = std::chrono::system_clock::to_time_t(std::chrono::system_clock::now());
	std::stringstream ss;
	ss << std::put_time(std::localtime(&now_time_t), "%Y-%m-%d %H:%M");

	for (uint32_t k = K_START; k <= K_END; ++k) {
		const auto& oc_res = all_results["oc_std"][k];
		const auto& mc_res = all_results["miniception_std"][k];

		auto p_values = calculate_wilcoxon(oc_res.results, mc_res.results);

		auto get_stdev = [](const std::vector<double>& v, double avg) {
			if (v.size() < 2) return 0.0;
			double sq_sum = 0;
			for (double x : v) sq_sum += (x - avg) * (x - avg);
			return std::sqrt(sq_sum / (v.size() - 1));
			};

		double oc_avg = oc_res.average();
		double mc_avg = mc_res.average();

		ofs << ss.str() << ","
			<< k << ","
			<< oc_avg << "," << get_stdev(oc_res.results, oc_avg) << ","
			<< mc_avg << "," << get_stdev(mc_res.results, mc_avg) << ","
			<< p_values.first << "," << p_values.second << "\n";
	}

	std::cout << "Results appended to " << filename << std::endl;
}
// Helper to convert k-mer bits to DNA string
std::string kmer_to_dna(uint64_t kmer, uint32_t k) {
	std::string dna = "";
	const char nuc_map[] = { 'A', 'C', 'G', 'T' };
	for (int i = k - 1; i >= 0; --i) {
		uint8_t nuc_bits = (kmer >> (2 * i)) & 0b11;
		dna += nuc_map[nuc_bits];
	}
	return dna;
}

// Helper to convert projection bits to binary string
std::string val_to_bin(uint64_t val, uint32_t k) {
	std::string bin = "";
	for (int i = k - 1; i >= 0; --i) {
		bin += ((val >> i) & 1) ? '1' : '0';
	}
	return bin;
}




void print_state_change(const std::string& new_state, std::string& tracker) {
	if (new_state != tracker) {
		std::cout << "\n==================== " << new_state << " ====================\n";
		tracker = new_state;
	}
}

void simulate_comprehensive_10_states() {
	// --- CONSTANTS ---
	const uint32_t k = 5;
	const uint32_t w = 24;
	const uint64_t N = 350; // Total bits (~175 nucleotides)

	const uint64_t dna_k_mask = (1ull << (2 * k)) - 1;
	const uint64_t proj_k_mask = (1u << (k-1)) - 1;
	const uint64_t k_prime = k - 1;
	const uint64_t sixty_four_minus_k = 64ULL - k;

	std::vector<bool> bits(N);
	std::random_device rd;
	std::mt19937 gen(rd());
	std::uniform_int_distribution<int> dist(0, 1);
	for (size_t i = 0; i < N + 2000; ++i) bits[i] = dist(gen);

	uint64_t current_kmer = 0ULL;
	uint64_t projection = 0ULL;
	uint64_t i = 0ULL;
	std::string current_state_name = "INIT";

	// Lambda to handle the expanded print format
	auto print_step = [&](bool hashed, std::string n10 = "N/A", std::string t_val = "N/A", std::string ic = "N/A") {
		std::cout << "[Step " << std::setw(3) << i / 2 << "] "
			<< "DNA: " << kmer_to_dna(current_kmer, k) << " | "
			<< "PROJ: " << val_to_bin(projection, k) << " | "
			<< "NEXT10: " << std::setw(k) << n10 << " | "
			<< "t: " << std::setw(3) << t_val << " | "
			<< "i_c: " << std::setw(3) << ic << " | "
			<< (hashed ? "\033[1;32m[HASH]\033[0m" : "\033[1;30m[SKIP]\033[0m")
			<< std::endl;
		};

	std::cout << "Starting Detailed Simulation (k=" << k << ", w=" << w << ")\n";

	// --- INITIALIZATION ---
	current_kmer = (3ULL << 2);
	projection = (1U << 2);
	uint64_t startup_nucs = (k > 4) ? k - 4 : 0;
	for (uint64_t j = 0; j < startup_nucs && i < N; ++j, i += 2) {
		uint8_t nuc = (bits[i] << 1) | bits[i + 1];
		current_kmer = ((current_kmer << 2) | nuc) & dna_k_mask;
		projection = ((projection << 1) | (nuc >> 1)) & proj_k_mask;
		print_step(false);
	}

	goto STATE_10;

STATE_10:
	print_state_change("STATE 10", current_state_name);
	while (i < N) {
		uint8_t nuc = (bits[i] << 1) | bits[i + 1];
		current_kmer = ((current_kmer << 2) | nuc) & dna_k_mask;
		projection = ((projection << 1) | (nuc >> 1)) & proj_k_mask;
		i += 2;

		// Calculate metadata for Step 10
		uint64_t next_10 = (projection & ~(projection << 1));
		uint64_t t = 64ULL - std::countl_zero(next_10);

		print_step(true, val_to_bin(next_10, k), std::to_string(t));

		if (t > 1) {
			uint64_t trivial_limit = i + 2*(k_prime - t);
			for (; i < trivial_limit; i += 2) {
				uint8_t tnuc = (bits[i] << 1) | bits[i + 1];
				current_kmer = ((current_kmer << 2) | tnuc) & dna_k_mask;
				projection = ((projection << 1) | (tnuc >> 1)) & proj_k_mask;
				print_step(false); // In trivial loop, next_10/t aren't recalculated
			}
			continue;
		}
		else if (t != 0) goto STATE_1_ENTRY;
		else goto STATE_0;
	}
	goto END;

STATE_0:
	print_state_change("STATE 0", current_state_name);
	while (i < N) {
		uint8_t nuc = (bits[i] << 1) | bits[i + 1];
		current_kmer = ((current_kmer << 2) | nuc) & dna_k_mask;
		projection = ((projection << 1) | (nuc >> 1)) & proj_k_mask;
		i += 2;
		print_step(false);

		if ((projection & 1) != 0) goto STATE_1_ENTRY;
	}
	goto END;

STATE_1_ENTRY:
	print_state_change("STATE 1 ENTRY", current_state_name);
	{
		uint64_t i_count = 1;
		while (i < N) {
			uint8_t nuc = (bits[i] << 1) | bits[i + 1];
			current_kmer = ((current_kmer << 2) | nuc) & dna_k_mask;
			projection = ((projection << 1) | (nuc >> 1)) & proj_k_mask;
			i += 2;

			print_step(false, "N/A", "N/A", std::to_string(i_count));

			if ((projection & 1) == 0) {
				for (uint64_t j = 0; j < k_prime - 2 && i < N; ++j, i += 2) {
					uint8_t tnuc = (bits[i] << 1) | bits[i + 1];
					current_kmer = ((current_kmer << 2) | tnuc) & dna_k_mask;
					projection = ((projection << 1) | (tnuc >> 1)) & proj_k_mask;
					print_step(false);
				}
				goto STATE_10;
			}
			else {
				i_count++;
			}
		}
	}

END:
	std::cout << "\nSimulation Finished.\n";
}







void perform_all_sampling_tests() {
	const uint64_t N = 300'000'000;
	const uint32_t X_runs = 100;  // was 10 or 20
	const uint32_t default_w= 32;


	std::vector<bool> bits = get_bits(N);

	//simulate_comprehensive_10_states();
	run_k_scaling_benchmark(X_runs, N, default_w, bits);



	std::cout << "\nAll tests complete. CSV files generated." << std::endl;
}
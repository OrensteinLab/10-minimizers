#include <iostream>
#include <thread>
#include <mutex>
#include <limits>
#include <algorithm> 
#include <cmath>
#include <random>     
#include <map>
#include <iomanip> 
#include <fstream>
#include <sstream>
#include <string>
#include "config.h"
#include <ctime>    
#include <stdexcept>
#include <unordered_set>
#include <queue>
#include <unordered_map>
#include <vector>
#include <cstdint>
#include <utility> 
#include "sampling_runtimes.h"
#include <functional>
#include <chrono>
#include <bit>
#include <intrin.h>
#include <numeric>


std::vector<bool> get_bits(uint64_t n_bits) {
	std::vector<bool> bits(n_bits);

	std::random_device rd;
	std::mt19937 gen(rd());
	std::uniform_int_distribution<int> dist(0, 1);

	for (size_t i = 0; i < n_bits; ++i) {
		bits[i] = dist(gen); // Either 0 or 1
	}

	return bits;

}






double test_runtime_std_hash(uint32_t w, uint32_t k, uint64_t N, const std::vector<bool>& bits) {
    uint64_t dna_k_mask = (1ull << (2 * k)) - 1;

    std::hash<uint64_t> hash_fn; // Create a hash function instance
    uint64_t sum = 0;

    // Start timing
    auto start = std::chrono::high_resolution_clock::now();

    uint64_t current_kmer = 0;
    for (uint64_t i = 0; i < N; i += 2) {
        current_kmer = (current_kmer << 1) & dna_k_mask;
        current_kmer |= (bits[i] << 1) | bits[i + 1];
        uint64_t kmer_rank = hash_fn(current_kmer); // Use the hash function
        sum += kmer_rank;
    }

    auto end = std::chrono::high_resolution_clock::now();

    std::cout << "sum:\t" << sum << std::endl;

    // Measure elapsed time in nanoseconds
    auto elapsed_ns = std::chrono::duration_cast<std::chrono::nanoseconds>(end - start).count();
    // Or microseconds
    auto elapsed_us = std::chrono::duration_cast<std::chrono::microseconds>(end - start).count();

    std::cout << "Time taken std::hash (ms):  " << elapsed_us / 1000.0 << "\n";

    return elapsed_us / 1000.0;
}

double test_runtime_xor_hash(uint32_t w, uint32_t k, uint64_t N, const std::vector<bool>& bits) {
	uint64_t k_mask = (1 << k) - 1;
	uint64_t dna_k_mask = (1 << (2 * k)) - 1;

	// Create a random number generator and distribution
	std::random_device rd;  // Seed for randomness
	std::mt19937_64 gen(rd()); // Mersenne Twister 64-bit generator
	std::uniform_int_distribution<uint64_t> dist(0, dna_k_mask); // Range based on the dna_k_mask

	uint64_t random_number = dist(gen); // Generate the random number


	uint64_t sum = 0;

	// start timing
	auto start = std::chrono::high_resolution_clock::now();

	uint64_t current_kmer = 0;
	for (uint64_t i = 0; i < N; i += 2) {
		current_kmer = (current_kmer << 1) & dna_k_mask;
		current_kmer |= (bits[i] << 1) | bits[i + 1];
		uint64_t kmer_rank = current_kmer ^ random_number;
		sum += kmer_rank;
	}



	auto end = std::chrono::high_resolution_clock::now();

	std::cout << "sum:\t" << sum << std::endl;

	// Measure elapsed time in nanoseconds
	auto elapsed_ns = std::chrono::duration_cast<std::chrono::nanoseconds>(end - start).count();
	// Or microseconds
	auto elapsed_us = std::chrono::duration_cast<std::chrono::microseconds>(end - start).count();

	std::cout << "Time taken xor hash (ms):  " << elapsed_us / 1000.0 << "\n";

	return elapsed_us / 1000.0;

}

double test_runtime_gm(uint32_t w, uint32_t k, uint64_t N, const std::vector<bool>& bits) {

	uint64_t k_mask = (1 << k) - 1;

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

double test_runtime_gm_no_memory_access(uint32_t w, uint32_t k, uint64_t N, const std::vector<bool>& bits) {

	uint64_t k_mask = (1 << k) - 1;
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
		uint64_t kmer_rank = odd_kmer + even_kmer;
		sum += kmer_rank;
	}


	auto end = std::chrono::high_resolution_clock::now();

	std::cout << "sum:\t" << sum << std::endl;

	// Measure elapsed time in nanoseconds
	auto elapsed_ns = std::chrono::duration_cast<std::chrono::nanoseconds>(end - start).count();
	// Or microseconds
	auto elapsed_us = std::chrono::duration_cast<std::chrono::microseconds>(end - start).count();

	std::cout << "Time taken gm no memory (ms):  " << elapsed_us / 1000.0 << "\n";
	return elapsed_us / 1000.0;
}

double test_runtime_abb_plus(uint32_t w, uint32_t k, uint64_t N, const std::vector<bool>& bits) {
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
		uint64_t kmer_rank = (key << L_shift) | current_kmer;
		sum1 += kmer_rank;
		sum2 += current_kmer;
	}

	auto end = std::chrono::high_resolution_clock::now();
	std::cout << "sum:\t" << sum1 << sum2 << std::endl;
	auto elapsed_us = std::chrono::duration_cast<std::chrono::microseconds>(end - start).count();
	std::cout << "Time taken ABB+ (ms): " << elapsed_us / 1000.0 << "\n";
	return elapsed_us / 1000.0;
}

double test_runtime_double_decycling_64bit(uint32_t k, uint64_t N, const std::vector<bool>& bits) {
	// 1. Setup Constants and Masks
	uint64_t dna_k_mask = (k == 32) ? ~0ULL : (1ULL << (2 * k)) - 1;
	uint64_t nondec_mask = 1ULL << (2 * k + 1); // Logic from your constructor
	uint64_t symdec_mask = 1ULL << (2 * k);
	uint64_t rand_hash = 0x9e3779b97f4a7c15ULL;
	uint64_t checksum = 0;
	const double EPS = 0.00000001;

	// 2. Precompute Weights (Exact replica of getWeights)
	double u = 2.0 * 3.14 / (double)k;
	std::vector<std::vector<double>> weights(4, std::vector<double>(k));
	for (int j = 0; j < 4; j++) {
		for (int i = 0; i < k; i++) {
			// Your logic: sin(i * u) * nucleotide_value
			// We store it such that weights[nuc][i] corresponds to the i-th position
			weights[j][i] = std::sin(i * u) * j;
		}
	}

	// Buffer for I-cycle rotation check
	std::vector<uint8_t> kmerArray(k);

	auto start = std::chrono::high_resolution_clock::now();

	uint64_t current_kmer = 0;
	for (uint64_t n = 0; n < N; n += 2) {
		// A. Slide k-mer
		uint8_t nuc = (bits[n] << 1) | bits[n + 1];
		current_kmer = ((current_kmer << 2) | nuc) & dna_k_mask;

		// B. Logic: inDoubleDec
		double sum = 0.0;
		double shiftsum = 0.0;

		// Compute Current Sum
		uint64_t temp_mer = current_kmer;
		for (int i = 0; i < k - 1; i++) {
			sum += weights[temp_mer & 0b11][i];
			temp_mer >>= 2;
		}

		// Compute Shift Sum (kmer >> 2 logic from your computeShiftSum)
		uint64_t shift_mer = current_kmer >> 2;
		for (int i = 0; i < k - 1; i++) {
			shiftsum += weights[shift_mer & 0b11][i];
			shift_mer >>= 2;
		}

		unsigned setMembership = 0;
		if (sum > EPS && shiftsum < -EPS) {
			setMembership = 2; // First L node
		}
		else if (sum < -EPS && shiftsum > EPS) {
			setMembership = 1; // First R node
		}
		else if (std::abs(sum) <= EPS && std::abs(shiftsum) <= EPS) {
			// I-cycle: Lexicographical Smallest Rotation Check
			uint64_t rot_mer = current_kmer;
			for (int i = k - 1; i >= 0; i--) {
				kmerArray[i] = rot_mer & 0b11;
				rot_mer >>= 2;
			}

			// Your specific I-cycle loop logic
			bool is_smallest = true;
			int i_ptr = 0, j_ptr = 1;
			// First loop pass
			for (; j_ptr < k; j_ptr++) {
				if (kmerArray[j_ptr] < kmerArray[i_ptr]) { is_smallest = false; break; }
				if (kmerArray[j_ptr] > kmerArray[i_ptr]) { i_ptr = 0; }
				else { i_ptr++; }
				if (i_ptr == 0 || i_ptr == k) { is_smallest = true; break; }
			}
			// Second loop pass (if still potentially smallest)
			if (is_smallest && j_ptr == k) {
				for (j_ptr = 0; j_ptr < k; j_ptr++) {
					if (kmerArray[j_ptr] < kmerArray[i_ptr]) { is_smallest = false; break; }
					if (kmerArray[j_ptr] > kmerArray[i_ptr]) { i_ptr = 0; }
					else { i_ptr++; }
					if (i_ptr == 0 || i_ptr == k) { is_smallest = true; break; }
				}
			}
			if (is_smallest) setMembership = 2;
		}

		// C. Apply myMinHash Logic
		uint64_t hashed_mer = current_kmer ^ rand_hash;
		if (setMembership == 0) {
			hashed_mer |= nondec_mask;
		}
		else if (setMembership == 1) {
			hashed_mer |= symdec_mask;
		}

		checksum += hashed_mer;
	}

	auto end = std::chrono::high_resolution_clock::now();
	double elapsed_ms = std::chrono::duration_cast<std::chrono::microseconds>(end - start).count() / 1000.0;

	std::cout << "Double Decycling Checksum: " << checksum << std::endl;
	std::cout << "Time taken (ms): " << elapsed_ms << "\n";
	return elapsed_ms;
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
	for (uint32_t i = 2; i < k - 1; ++i) tau.push_back(2 * ((1ull << (k - 1 - i)) - 1));
	tau.push_back(1);
	tau.push_back((1ull << k) - 1);
	tau.push_back(0);

	uint64_t class_bits = (1ull << (62));
	uint64_t half_shift = (1ull << (32));

	for (int i = 0; i < (int)tau.size(); ++i) lemma1_ranks[tau[i]] = i;


	// Circular Buffer for Lazy Ranking
	struct PendingKmer { uint64_t kmer; uint64_t projection; };
	constexpr uint32_t BUF_SIZE = 128; // Fixed power of 2
	constexpr uint32_t BUF_MASK = BUF_SIZE - 1;
	PendingKmer buffer[BUF_SIZE]; // Allocated on the stack

	uint32_t buf_start = 0;
	uint32_t buf_count = 0;
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
			if (buf_count == (w-1)) [[unlikely]]  {
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
				uint32_t next_pos = (buf_start + buf_count) & BUF_MASK;
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



double test_runtime_double_decycling_optimized(uint32_t w, uint32_t k, uint64_t N, const std::vector<bool>& bits) {
	uint64_t dna_k_mask = (k == 32) ? ~0ULL : (1ULL << (2 * k)) - 1;
	// Masks for hashing flags

	uint64_t rand_hash = 0x9e3779b97f4a7c15ULL;
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
	uint32_t window_ptr = 0; // Circular buffer pointer
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

		for (uint32_t i = 0; i < k - 1; ++i) {
			uint32_t idx = (window_ptr + i) % k;
			sum += weights[kmer_window[idx]][i];
		}

		for (uint32_t i = 0; i < k - 1; ++i) {
			uint32_t idx = (window_ptr + i + 1) % k;
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

		// --- Hashing ---
		uint64_t hashed_mer = current_kmer ^ rand_hash;
		sum2 += hashed_mer;
	}

	auto end = std::chrono::high_resolution_clock::now();
	double elapsed_ms = std::chrono::duration_cast<std::chrono::microseconds>(end - start).count() / 1000.0;

	std::cout << "Double Decycling : " << sum1 << sum2 << std::endl;
	std::cout << "Time taken 10-Min DD (ms): " << elapsed_ms  << "\n";
	return elapsed_ms;
}


double test_runtime_oc_minimizer(uint32_t w, uint32_t k, uint64_t N, const std::vector<bool>& bits) {
	const uint32_t s = 4;
	const uint32_t v = (k - s) / 2;
	const uint32_t w_prime = k - s + 1;

	uint64_t dna_k_mask = (k == 32) ? ~0ULL : (1ULL << (2 * k)) - 1;
	uint64_t s_mask = (1ULL << (2 * s)) - 1;

	const uint64_t s_hash_seed = 0x9e3779b97f4a7c15ULL;
	const uint64_t k_hash_seed = 0xBF58476D1CE4E5B9ULL;

	uint64_t sum1 = 0, sum2 = 0;

	// Circular buffers
	uint64_t smer_hashes[128] = { 0 };
	uint8_t kmer_ranks[128] = { 0 }; // Stores 1 (Open), 2 (Closed), or 3 (None)
	const uint32_t BUF_MASK = 127;
	uint32_t ptr = 0; // Shared pointer for simplicity

	auto start = std::chrono::high_resolution_clock::now();

	uint64_t current_kmer = 0;
	uint64_t current_smer = 0;

	for (uint64_t n = 0; n <  N; n += 2) {
		uint8_t nuc = (bits[n] << 1) | bits[n + 1];
		current_kmer = ((current_kmer << 2) | nuc) & dna_k_mask;
		current_smer = ((current_smer << 2) | nuc) & s_mask;

		// 1. Store the new s-mer hash
		smer_hashes[ptr] = current_smer ^ s_hash_seed;

		// 2. PRE-CALCULATE rank for the NEWEST k-mer only (O(w') instead of O(w*w'))
		uint64_t min_s_val = 0xFFFFFFFFFFFFFFFFULL;
		uint32_t min_s_pos = 0;
		for (uint32_t j = 0; j < w_prime; ++j) {
			uint32_t lookup_idx = (ptr - (w_prime - 1 - j)) & BUF_MASK;
			uint64_t s_h = smer_hashes[lookup_idx];
			if (s_h < min_s_val) {
				min_s_val = s_h;
				min_s_pos = j;
			}
		}

		// Save this rank so we don't re-calculate it for the next w windows
		if (min_s_pos == v) kmer_ranks[ptr] = 1;
		else if (min_s_pos == 0 || min_s_pos == w_prime - 1) kmer_ranks[ptr] = 2;
		else kmer_ranks[ptr] = 3;

		// 3. OC-Minimizer logic
		for (uint32_t i = 0; i < w; ++i) {
			uint32_t k_idx = (ptr - (w - 1 - i)) & BUF_MASK;
			sum1 += kmer_ranks[k_idx]; // accessing the rank directly without recalculating
		}

		sum2 += (current_kmer ^ k_hash_seed); // accessing the hash of the current k-mer
		ptr = (ptr + 1) & BUF_MASK;
	}

	auto end = std::chrono::high_resolution_clock::now();
	double elapsed_ms = std::chrono::duration_cast<std::chrono::microseconds>(end - start).count() / 1000.0;

	std::cout << "OC-Minimizer (Cached) : " << sum1 << " " << sum2 << std::endl;
	std::cout << "Time taken OC-Minimizer (ms): " << elapsed_ms << "\n";

	return elapsed_ms;
}


double test_runtime_miniception(uint32_t w, uint32_t k, uint64_t N, const std::vector<bool>& bits) {
	const uint32_t s = 4;
	const uint32_t w_prime = k - s + 1;

	uint64_t dna_k_mask = (k == 32) ? ~0ULL : (1ULL << (2 * k)) - 1;
	uint64_t s_mask = (1ULL << (2 * s)) - 1;

	const uint64_t s_hash_seed = 0x9e3779b97f4a7c15ULL;
	const uint64_t k_hash_seed = 0xBF58476D1CE4E5B9ULL;

	uint64_t sum1 = 0, sum2 = 0;

	// Circular buffers
	uint64_t smer_hashes[128] = { 0 };
	uint8_t kmer_ranks[128] = { 0 }; // Stores 1 (Closed) or 2 (None)
	const uint32_t BUF_MASK = 127;
	uint32_t ptr = 0;

	auto start = std::chrono::high_resolution_clock::now();

	uint64_t current_kmer = 0;
	uint64_t current_smer = 0;

	for (uint64_t n = 0; n <  N; n += 2) {
		uint8_t nuc = (bits[n] << 1) | bits[n + 1];
		current_kmer = ((current_kmer << 2) | nuc) & dna_k_mask;
		current_smer = ((current_smer << 2) | nuc) & s_mask;

		// 1. Store the new s-mer hash
		smer_hashes[ptr] = current_smer ^ s_hash_seed;

		// 2. PRE-CALCULATE rank for Miniception (Strictly Closed Syncmers)
		uint64_t min_s_val = 0xFFFFFFFFFFFFFFFFULL;
		uint32_t min_s_pos = 0;
		for (uint32_t j = 0; j < w_prime; ++j) {
			uint32_t lookup_idx = (ptr - (w_prime - 1 - j)) & BUF_MASK;
			uint64_t s_h = smer_hashes[lookup_idx];
			if (s_h < min_s_val) {
				min_s_val = s_h;
				min_s_pos = j;
			}
		}

		// Save rank: 1 for Closed Syncmer, 2 otherwise 
		if (min_s_pos == 0 || min_s_pos == w_prime - 1) kmer_ranks[ptr] = 1;
		else kmer_ranks[ptr] = 2;

		// 3. Miniception logic using cached ranks
		for (uint32_t i = 0; i < w; ++i) {
			uint32_t k_idx = (ptr - (w - 1 - i)) & BUF_MASK;
			sum1 += kmer_ranks[k_idx];
		}

		sum2 += (current_kmer ^ k_hash_seed);
		ptr = (ptr + 1) & BUF_MASK;
	}

	auto end = std::chrono::high_resolution_clock::now();
	double elapsed_ms = std::chrono::duration_cast<std::chrono::microseconds>(end - start).count() / 1000.0;

	std::cout << "Miniception (Cached) : " << sum1 << " " << sum2 << std::endl;
	std::cout << "Time taken Miniception (ms): " << elapsed_ms << "\n";

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

void run_k_scaling_benchmark(uint32_t num_runs, uint64_t N, uint32_t w_fixed, const std::vector<bool>& bits) {
	// 1. Create a list of all "tasks" (every k value repeated num_runs times)
	std::vector<uint32_t> task_pool;
	for (uint32_t k = 3; k <= 32; ++k) {
		for (uint32_t i = 0; i < num_runs; ++i) {
			task_pool.push_back(k);
		}
	}

	// 2. Shuffle the tasks
	std::random_device rd;
	std::mt19937 g(rd());
	std::shuffle(task_pool.begin(), task_pool.end(), g);

	// 3. Prepare storage for results
	// We use a map to group results by k since they arrive in random order
	std::map<uint32_t, Stats> rand_stats, gm_no_stats, gm_stats, std_stats, abb_stats, lazy_stats, dd_stats, oc_xor_stats, miniception_xor_stats;

	std::cout << "Starting randomized benchmark (" << task_pool.size() << " total runs)..." << std::endl;

	// 4. Run the shuffled tasks
	for (size_t i = 0; i < task_pool.size(); ++i) {
		uint32_t k = task_pool[i];

		// Progress indicator (optional)
		if (i % 50 == 0) std::cout << "Progress: " << i << "/" << task_pool.size() << "\r" << std::flush;

		// Run each test once for this specific k
		rand_stats[k].results.push_back(test_runtime_xor_hash(w_fixed, k, N, bits));
		//gm_no_stats[k].results.push_back(test_runtime_gm_no_memory_access(w_fixed, k, N, bits));
		//gm_stats[k].results.push_back(test_runtime_gm(w_fixed, k, N, bits));
		std_stats[k].results.push_back(test_runtime_std_hash(w_fixed, k, N, bits));
		abb_stats[k].results.push_back(test_runtime_abb_plus(w_fixed, k, N, bits));
		lazy_stats[k].results.push_back(test_runtime_one_zero_most_optimized_lazy(w_fixed, k, N, bits));
		dd_stats[k].results.push_back(test_runtime_double_decycling_optimized(w_fixed, k, N, bits));
		oc_xor_stats[k].results.push_back(test_runtime_oc_minimizer(w_fixed, k, N, bits));
		miniception_xor_stats[k].results.push_back(test_runtime_miniception(w_fixed, k, N, bits));
	}

	// 5. Write to CSV (sorted by k automatically by the map)
	std::ofstream ofs("sampling_runtime_k_scaling.csv");
	ofs << "k,time_random,time_random_std,time_abb_plus,time_lazy,time_dd,time_oc,time_mc\n";
	//ofs << "k,time_lazy\n";
	//ofs << "k,time_oc,time_mc\n";

	for (uint32_t k = 3; k <= 32; ++k) {
		ofs << k << ","
			<< rand_stats[k].average() << ","
			//<< gm_no_stats[k].average() << ","
			//<< gm_stats[k].average() << ","
			<< std_stats[k].average() << ","
			<< abb_stats[k].average() << ","
			<< lazy_stats[k].average() << ","
			<< dd_stats[k].average() << "," 
			<< oc_xor_stats[k].average() << ","
			<< miniception_xor_stats[k].average() <<
			"\n";
	}
	std::cout << "\nDone!" << std::endl;
}


void run_w_scaling_benchmark(uint32_t num_runs, uint64_t N, uint32_t k_fixed, const std::vector<bool>& bits) {
	// 1. Create a "Task Pool" containing every 'w' value repeated 'num_runs' times
	std::vector<uint32_t> task_pool;
	for (uint32_t w = 3; w <= 100; w++) {
		for (uint32_t i = 0; i < num_runs; ++i) {
			task_pool.push_back(w);
		}
	}

	// 2. Shuffle the execution order
	std::random_device rd;
	std::mt19937 g(rd());
	std::shuffle(task_pool.begin(), task_pool.end(), g);

	// 3. Map to store cumulative times for each 'w'
	// Key: w value, Value: running sum of execution times
	std::map<uint32_t, double> total_times;
	std::cout << "Starting randomized w-scaling benchmark..." << std::endl;

	// 4. Run the shuffled tasks
	for (size_t i = 0; i < task_pool.size(); ++i) {
		uint32_t w = task_pool[i];

		// Run the test exactly once for this specific instance of w
		double run_time = test_runtime_one_zero_most_optimized_lazy(w, k_fixed, N, bits);

		// Accumulate the result
		total_times[w] += run_time;

		if (i % 50 == 0) {
			std::cout << "Progress: " << (i * 100 / task_pool.size()) << "% \r" << std::flush;
		}
	}

	// 5. Calculate averages and write to CSV
	std::ofstream ofs("lazy_runtime_w_scaling.csv");
	ofs << "w,time_lazy_ms\n";

	// Since std::map is ordered, it will write w=2, w=3, etc. automatically
	for (auto const& [w, total_time] : total_times) {
		double avg_time = total_time / static_cast<double>(num_runs);
		ofs << w << "," << avg_time << "\n";
	}

	std::cout << "\nBenchmark complete. Results saved to lazy_runtime_w_scaling.csv" << std::endl;
}


void perform_all_sampling_tests() {
	const uint64_t N = 300'000'000;
	const uint32_t X_runs = 10; 
	const uint32_t default_w = 24;
	const uint32_t default_k = 12;

	std::vector<bool> bits = get_bits(N);


	run_k_scaling_benchmark(X_runs, N, default_w, bits);

	run_w_scaling_benchmark(X_runs, N, default_k, bits);

	std::cout << "\nAll tests complete. CSV files generated." << std::endl;
}
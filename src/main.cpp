/*
 * main.cpp
 *
 * Compile:  make
 * Usage:    ./bin/genome_mapper genome.fasta reads.fasta [OPTIONS]
 *
 * Required:
 *   genome.fasta          genome reference file
 *   reads.fasta           reads to align
 *
 * Options:
 *   --kmer_size   INT     k-mer length used for indexing (default: 21, max: 31)
 *   --table_size  INT     number of hash-table buckets  (default: 4000000)
 *   --heap_size   INT     max candidates per read kept in the BestResultsHeap (default: 15)
 *   --help                print this help message and exit
 *
 * Output columns (TSV, stdout):
 *   #read  chrom  start(0-based)  end(0-based inclusive)  mismatches  matches  candidates
 *   Unmapped reads: chrom=*  start=end=mismatches=matches=-1  candidates=0
 */

#include "BaseLUT.h"
#include "HashTable.h"
#include "GenomicPosition.h"
#include "AlignmentResult.h"
#include "BestResultsHeap.h"

#include <algorithm>
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <stdexcept>
#include <cstdlib>
#include <cstdint>

using namespace dna;

// ─────────────────────────────────────────────────────────────────────────────
// Sanitize a raw sequence string in-place.
// Every byte is run through BASE_LUT; if INVALID, replace with 'X'
// ('X' also maps to INVALID in BASE_LUT — safe sentinel).
// ─────────────────────────────────────────────────────────────────────────────
static void sanitizeSeq(std::string& seq) {
    for (uint32_t i = 0; i < static_cast<uint32_t>(seq.size()); ++i) {
        unsigned char c = static_cast<unsigned char>(seq[i]);
        if (BASE_LUT[c] == INVALID)
            seq[i] = 'X';
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// FASTA reader
// Returns (name, sanitized-sequence) pairs.
// Headers trimmed to first whitespace token. \r stripped for Windows files.
// ─────────────────────────────────────────────────────────────────────────────
typedef std::vector<std::pair<std::string, std::string> > GenomeMap;

static GenomeMap readFasta(const std::string& path) {
    std::ifstream in(path.c_str());
    if (!in.is_open())
        throw std::runtime_error("Cannot open file: " + path);

    GenomeMap records;
    std::string line, name, seq;

    while (std::getline(in, line)) {
        if (!line.empty() && line[line.size() - 1] == '\r')
            line.erase(line.size() - 1);
        if (line.empty()) continue;

        if (line[0] == '>') {
            if (!name.empty()) {
                sanitizeSeq(seq);
                records.push_back(std::make_pair(name, seq));
            }
            name = line.substr(1);
            uint32_t sp = static_cast<uint32_t>(name.find_first_of(" \t"));
            if (sp != static_cast<uint32_t>(std::string::npos))
                name = name.substr(0, sp);
            seq.clear();
        } else {
	    std::transform(line.begin(), line.end(), line.begin(), ::toupper);
            seq += line;
        }
    }
    if (!name.empty()) {
        sanitizeSeq(seq);
        records.push_back(std::make_pair(name, seq));
    }
    return records;
}

// ─────────────────────────────────────────────────────────────────────────────
// Find a chromosome sequence by name (linear scan)
// ─────────────────────────────────────────────────────────────────────────────
static const std::string* findChrom(const GenomeMap& genome,
                                    const std::string& name) {
    for (uint32_t i = 0; i < static_cast<uint32_t>(genome.size()); ++i)
        if (genome[i].first == name) return &genome[i].second;
    return NULL;
}

// ─────────────────────────────────────────────────────────────────────────────
// Score a read against a genomic window.
// BASE_LUT used for every comparison — case differences and 'X' positions
// are handled automatically. Out-of-bounds positions count as mismatches.
// ─────────────────────────────────────────────────────────────────────────────
static void scoreAlignment(const std::string& read,
                            const std::string& chromSeq,
                            long long          genomeStart,
                            int&               mismatches,
                            int&               matches) {
    mismatches = 0;
    matches    = 0;
    const long long glen = static_cast<long long>(chromSeq.size());
    const int       rlen = static_cast<int>(read.size());

    for (int i = 0; i < rlen; ++i) {
        long long gpos = genomeStart + i;
        if (gpos < 0 || gpos >= glen) {
            ++mismatches;
            continue;
        }
        uint8_t r = BASE_LUT[static_cast<unsigned char>(read[i])];
        uint8_t g = BASE_LUT[static_cast<unsigned char>(
                                chromSeq[static_cast<uint32_t>(gpos)])];
        if (r == INVALID || g == INVALID) continue;
        if (r == g) ++matches;
        else        ++mismatches;
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// Build hash-table index from the genome using the 2-bit sliding window.
// Resets on any INVALID base so no kmer spanning an invalid position is stored.
// ─────────────────────────────────────────────────────────────────────────────
static void buildIndex(HashTable&       ht,
                       const GenomeMap& genome,
                       int              k) {
    uint32_t inserted = 0, skipped = 0;
    const uint64_t mask = (k < 32)
                          ? ((uint64_t(1) << (2 * k)) - 1)
                          : ~uint64_t(0);
    
    for (uint32_t ri = 0; ri < static_cast<uint32_t>(genome.size()); ++ri) {
        const std::string& chrom = genome[ri].first;
        const std::string& seq   = genome[ri].second;
        const uint32_t     len   = static_cast<uint32_t>(seq.size());
        if (len < static_cast<uint32_t>(k)) continue;
        
        uint64_t encoded   = 0;
        int      count = 0;

        for (uint32_t i = 0; i < len; ++i) {
            uint8_t base = ht.encodeBase(seq[i]);

            if (base == INVALID) {
                encoded = 0;
                count = 0;
                ++skipped;
                continue;
            }

            encoded = ((encoded << 2) | static_cast<uint64_t>(base)) & mask;
            ++count;

            if (count >= k) {
                uint32_t start = i - static_cast<uint32_t>(k) + 1;
                GenomicPosition pos(chrom, start);

                ht.insert(encoded, pos);   
                ++inserted;
            }            
        }  
    }
    std::cerr << "[index] k=" << k
              << "  inserted=" << inserted
              << "  skipped(non-ACGT)=" << skipped << "\n";
}

// ─────────────────────────────────────────────────────────────────────────────
// Map reads to the genome
// ─────────────────────────────────────────────────────────────────────────────
static void mapReads(const HashTable&  ht,
                     const GenomeMap&  genome,
                     const GenomeMap&  reads,
                     int               k,
                     int               heapSize) {
    uint32_t mapped = 0, unmapped = 0;
    std::cout << "#read\tchrom\tstart\tend\tmismatches\tmatches\tcandidates\n";

    for (uint32_t ri = 0; ri < static_cast<uint32_t>(reads.size()); ++ri) {
        const std::string& readName = reads[ri].first;
        const std::string& seq      = reads[ri].second;
        const int          rlen     = static_cast<int>(seq.size());

        //Declaration du HashTable
        BestResultsHeap heap(heapSize);
        std::vector<std::pair<std::string, long long> > already_scored;

        if (rlen >= k) {
            for (int i = 0; i + k <= rlen; ++i) {
                std::string kmer = seq.substr(static_cast<uint32_t>(i),
                                              static_cast<uint32_t>(k));

                uint64_t modified_base = 0;
                if (!ht.encodeKmer(kmer, modified_base)) {
                    continue; 
                }

                //Rechercher le bucket correspondanr au hash du kmer en lecture (read)
                const LinkedList* bucket = ht.lookup(modified_base);

                //Rechercher 
                const Node* current = bucket->first();
                while (current != nullptr) {

                    //S'assurer qu'on compare seulement avec les nodes du bon seed du genome
                    //ignorant les nodes stockes par collision de hash
                    if (current->m_kmerCode == modified_base) {
                        const GenomicPosition& pos = current->m_position;

                        long long read_start =
                            static_cast<long long>(pos.offset) - static_cast<long long>(i);

                        //Verifier si lu avant
                        bool was_read = false;
                        for (uint32_t j = 0; j < static_cast<uint32_t>(already_scored.size());
                             ++j) {
                            if (already_scored[j].first == pos.chrom &&
                                already_scored[j].second == read_start) {
                                was_read = true;
                                break;
                            }
                        }

                        if (!was_read) {
                            //Sauver dans la liste des reads vus
                            already_scored.push_back(std::make_pair(pos.chrom, read_start));

                            const std::string* chromSeq = findChrom(genome, pos.chrom);
                            if (!chromSeq) { continue; }
                            int mm = 0;
                            int m  = 0;
                            scoreAlignment(seq, *chromSeq, read_start, mm, m);

                            GenomicPosition candidate_pos(pos.chrom,
                                    static_cast<size_t>(read_start < 0 ? 0 : read_start));

                            AlignmentResult result(candidate_pos, mm, m);
                            heap.insert(result);
                            
                        }
                    }
                    current = current->m_next;
                }
            }
        }

        if (!heap.isEmpty()) {
            AlignmentResult best = heap.getBest();
            long long start = static_cast<long long>(best.m_position.offset);
            long long end   = start + static_cast<long long>(rlen) - 1;

            std::cout << readName
                      << "\t" << best.m_position.chrom
                      << "\t" << start
                      << "\t" << end
                      << "\t" << best.m_mismatches
                      << "\t" << best.m_matches
                      << "\t" << heap.getSize()
                      << "\n";
            ++mapped;
        } else {
            std::cout << readName << "\t*\t-1\t-1\t-1\t-1\t0\n";
            ++unmapped;
        }
    }

    std::cerr << "[map] mapped=" << mapped
              << "  unmapped=" << unmapped << "\n";
}

// ─────────────────────────────────────────────────────────────────────────────
// Entry point
// ─────────────────────────────────────────────────────────────────────────────
// ─────────────────────────────────────────────────────────────────────────────
// Print help message to stdout
// ─────────────────────────────────────────────────────────────────────────────
static void printHelp(const char* progName) {
    std::cout
        << "\n"
        << "Usage:\n"
        << "  " << progName << " genome.fasta reads.fasta [OPTIONS]\n"
        << "\n"
        << "Required arguments:\n"
        << "  genome.fasta            FASTA file containing the reference genome\n"
        << "  reads.fasta             FASTA file containing the reads to align\n"
        << "\n"
        << "Options:\n"
        << "  --kmer_size  INT        Length of k-mers used to index the genome\n"
        << "                            default : 21\n"
        << "                            range   : 1 to 31 (must fit in uint64_t)\n"
        << "  --table_size INT        Number of buckets in the hash table\n"
        << "                            default : 4000000\n"
        << "                            tip     : larger = fewer collisions, more memory\n"
        << "  --heap_size  INT        Maximum number of candidate alignments kept\n"
        << "                          per read in the BestResultsHeap\n"
        << "                            default : 15\n"
        << "  --help                  Print this help message and exit\n"
        << "\n"
        << "Output (TSV on stdout):\n"
        << "  #read  chrom  start  end  mismatches  matches  candidates\n"
        << "  Unmapped reads are reported with chrom=* and numeric fields set to -1.\n"
        << "\n"
        << "Examples:\n"
        << "  " << progName << " genome.fasta reads.fasta\n"
        << "  " << progName << " genome.fasta reads.fasta --kmer_size 15\n"
        << "  " << progName << " genome.fasta reads.fasta --kmer_size 15 --table_size 8000000\n"
        << "  " << progName << " genome.fasta reads.fasta --heap_size 32\n"
        << "\n";
}

int main(int argc, char* argv[]) {
    // ── Defaults ──────────────────────────────────────────────────────────────
    std::string genomePath;
    std::string readsPath;
    int      k         = 21;
    uint32_t tableSize = 4000000u;
    int      heapSize  = 15;

    // ── Parse arguments ───────────────────────────────────────────────────────
    // Positional: first two non-flag arguments are genome and reads paths.
    // Named:      --kmer_size, --table_size, --heap_size each consume the next
    //             token as their value.
    for (int i = 1; i < argc; ++i) {
        std::string arg = argv[i];

        if (arg == "--help" || arg == "-h") {
            printHelp(argv[0]);
            return 0;
        } else if (arg == "--kmer_size") {
            if (i + 1 >= argc) {
                std::cerr << "Error: --kmer_size requires a value\n";
                return 1;
            }
            k = std::atoi(argv[++i]);
        } else if (arg == "--table_size") {
            if (i + 1 >= argc) {
                std::cerr << "Error: --table_size requires a value\n";
                return 1;
            }
            tableSize = static_cast<uint32_t>(std::atoi(argv[++i]));
        } else if (arg == "--heap_size") {
            if (i + 1 >= argc) {
                std::cerr << "Error: --heap_size requires a value\n";
                return 1;
            }
            heapSize = std::atoi(argv[++i]);
        } else if (arg.size() > 2 && arg[0] == '-' && arg[1] == '-') {
            std::cerr << "Error: unknown option '" << arg << "'\n";
            std::cerr << "Run with --help for usage.\n";
            return 1;
        } else {
            // Positional argument
            if (genomePath.empty())       genomePath = arg;
            else if (readsPath.empty())   readsPath  = arg;
            else {
                std::cerr << "Error: unexpected argument '" << arg << "'\n";
                std::cerr << "Run with --help for usage.\n";
                return 1;
            }
        }
    }

    // ── Validate ──────────────────────────────────────────────────────────────
    if (genomePath.empty() || readsPath.empty()) {
        std::cerr << "Error: genome.fasta and reads.fasta are required.\n";
        std::cerr << "Run with --help for usage.\n";
        return 1;
    }
    if (k < 1 || k > 31) {
        std::cerr << "Error: --kmer_size must be between 1 and 31\n";
        return 1;
    }
    if (tableSize < 1u) {
        std::cerr << "Error: --table_size must be at least 1\n";
        return 1;
    }
    if (heapSize < 1) {
        std::cerr << "Error: --heap_size must be at least 1\n";
        return 1;
    }

    std::cerr << "[load] genome: " << genomePath << "\n";
    GenomeMap genome = readFasta(genomePath);
    std::cerr << "[load] " << genome.size() << " sequence(s)\n";
    
    HashTable ht(tableSize);   
    buildIndex(ht, genome, k);

    std::cerr << "[load] reads: " << readsPath << "\n";
    GenomeMap reads = readFasta(readsPath);
    std::cerr << "[load] " << reads.size() << " read(s)\n";

    mapReads(ht, genome, reads, k, heapSize);
    return 0;
}
